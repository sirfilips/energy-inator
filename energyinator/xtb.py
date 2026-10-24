"""GFN-xTB Semiempirical Quantum Chemistry Engine Integration.

Provides interface to the xTB software suite (Grimme lab, Bonn) for semiempirical
tight-binding quantum mechanical calculations (GFN2-xTB, GFN1-xTB, GFN-FF).
Enables the hierarchical 3-tier workflow:
    RDKit (ETKDGv3) -> Force Field (MMFF94) -> GFN2-xTB -> Ab Initio DFT (Psi4)
"""

import os
import sys
import shutil
import logging
import tempfile
import subprocess
import re
from typing import Dict, Any, List, Optional, Tuple

from energyinator.models import ConformerResult
from energyinator.geometry import parse_xyz, format_xyz

logger = logging.getLogger("energyinator")

# Mapping of solvent names to xTB ALPB solvent names
XTB_SOLVENT_MAP = {
    "water": "water",
    "chloroform": "chloroform",
    "dmso": "dmso",
    "methanol": "methanol",
    "toluene": "toluene",
    "acetonitrile": "acetonitrile",
    "dcm": "ch2cl2",
    "dichloromethane": "ch2cl2",
    "ether": "diethylether",
    "thf": "thf",
    "acetone": "acetone",
    "benzene": "benzene"
}


def find_xtb_executable() -> Optional[str]:
    """Finds the xtb executable path on the system across PATH and virtual environments."""
    # 1. System PATH
    found = shutil.which("xtb") or shutil.which("xtb.exe")
    if found and os.path.isfile(found):
        return os.path.abspath(found)

    # 2. Local project venv or sys.prefix (common in Conda on Windows)
    potential_dirs = [
        os.path.join(sys.prefix, "Library", "bin"),
        os.path.join(sys.prefix, "bin"),
        os.path.join(sys.prefix, "Scripts"),
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".venv", "Library", "bin"),
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".venv", "bin"),
    ]

    for d in potential_dirs:
        for fname in ["xtb.exe", "xtb"]:
            candidate = os.path.join(d, fname)
            if os.path.isfile(candidate):
                return os.path.abspath(candidate)

    return None


def is_xtb_available() -> bool:
    """Returns True if xTB executable is discovered and operational."""
    return find_xtb_executable() is not None


def get_xtb_version() -> str:
    """Returns xTB version string if available."""
    exe = find_xtb_executable()
    if not exe:
        return "Non installato"
    try:
        res = subprocess.run([exe, "--version"], capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=5)
        for line in res.stdout.splitlines():
            if "xtb version" in line.lower():
                return line.strip().replace("*", "").strip()
        return "xTB (Disponibile)"
    except Exception as e:
        return f"Errore rilevamento: {e}"


def run_xtb_calculation(
    xyz_str: str,
    method: str = "gfn2",
    opt: bool = True,
    charge: int = 0,
    multiplicity: int = 1,
    solvent: Optional[str] = None,
    threads: int = 2
) -> Dict[str, Any]:
    """Executes a GFN-xTB calculation on an XYZ coordinate block.
    
    Args:
        xyz_str: Cartesian coordinates in XYZ format.
        method: xTB method ("gfn2", "gfn1", "gfn0", "gfnff").
        opt: Whether to perform geometry optimization.
        charge: Molecular formal charge.
        multiplicity: Spin multiplicity (2S+1).
        solvent: Optional implicit solvent name.
        threads: OpenMP CPU threads for parallel execution.
        
    Returns:
        Dictionary with energy, gap, dipole, and optimized XYZ.
    """
    exe = find_xtb_executable()
    if not exe:
        raise FileNotFoundError("Eseguibile xtb non trovato nell'ambiente di calcolo.")

    # Calculate number of unpaired electrons (UHF parameter: 2S = mult - 1)
    uhf = max(0, multiplicity - 1)

    cmd = [exe, "mol.xyz"]

    # Parameter flags
    method_lower = method.lower()
    if method_lower == "gfn1":
        cmd.extend(["--gfn", "1"])
    elif method_lower == "gfn0":
        cmd.extend(["--gfn", "0"])
    elif method_lower == "gfnff":
        cmd.append("--gfnff")
    else:
        cmd.extend(["--gfn", "2"])

    if opt:
        cmd.append("--opt")
    else:
        cmd.append("--sp")

    if charge != 0:
        cmd.extend(["--chrg", str(charge)])
    if uhf != 0:
        cmd.extend(["--uhf", str(uhf)])

    if solvent and solvent.lower() != "gas" and solvent.lower() != "vuoto":
        solv_key = solvent.lower().strip()
        xtb_solv = XTB_SOLVENT_MAP.get(solv_key, solv_key)
        cmd.extend(["--alpb", xtb_solv])

    cmd.extend(["--parallel", str(max(1, threads))])

    # Run inside an isolated temporary directory
    with tempfile.TemporaryDirectory() as td:
        mol_path = os.path.join(td, "mol.xyz")
        with open(mol_path, "w", encoding="utf-8") as f:
            f.write(xyz_str.strip() + "\n")

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(threads)
        env["MKL_NUM_THREADS"] = str(threads)

        try:
            res = subprocess.run(
                cmd,
                cwd=td,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=180
            )
        except subprocess.TimeoutExpired:
            raise TimeoutError("Calcolo semiempirico GFN-xTB scaduto per timeout (180 s).")

        stdout = res.stdout

        # Parse Total Energy (Hartree)
        # Matches: | TOTAL ENERGY               -5.070544328472 Eh   |
        energy_hartree: Optional[float] = None
        for line in stdout.splitlines():
            if "TOTAL ENERGY" in line:
                m = re.search(r"TOTAL ENERGY\s+([-\d\.]+)", line)
                if m:
                    energy_hartree = float(m.group(1))

        if energy_hartree is None:
            # Fallback for energy search
            m_alt = re.search(r"total energy\s*:\s*([-\d\.]+)\s*Eh", stdout, re.IGNORECASE)
            if m_alt:
                energy_hartree = float(m_alt.group(1))

        if energy_hartree is None:
            raise RuntimeError(f"xTB non ha restituito l'energia totale. Log:\n{stdout[-1000:]}")

        energy_kcalmol = energy_hartree * 627.509

        # Parse HOMO-LUMO Gap (eV)
        gap_ev: Optional[float] = None
        m_gap = re.search(r"HOMO-LUMO GAP\s+([-\d\.]+)\s+eV", stdout, re.IGNORECASE)
        if m_gap:
            gap_ev = float(m_gap.group(1))

        # Parse Dipole (Debye)
        dipole_debye: Optional[float] = None
        lines = stdout.splitlines()
        for idx, line in enumerate(lines):
            if "molecular dipole:" in line.lower():
                # Following lines contain: 'tot:     X.XXXX'
                for sub_line in lines[idx:idx + 5]:
                    if "tot:" in sub_line.lower():
                        m_dip = re.search(r"tot:\s*([-\d\.]+)", sub_line.lower())
                        if m_dip:
                            dipole_debye = float(m_dip.group(1))
                            break

        # Extract optimized coordinates if opt=True
        opt_xyz_str = xyz_str
        converged = True
        if opt:
            opt_file = os.path.join(td, "xtbopt.xyz")
            if os.path.exists(opt_file):
                with open(opt_file, "r", encoding="utf-8") as f:
                    opt_xyz_str = f.read().strip()
            else:
                converged = False

        return {
            "energy_hartree": energy_hartree,
            "energy_kcalmol": energy_kcalmol,
            "gap_ev": gap_ev,
            "dipole_debye": dipole_debye,
            "xyz_optimized": opt_xyz_str,
            "converged": converged,
            "stdout": stdout
        }


def refine_conformers_with_xtb(
    conformers: List[ConformerResult],
    method: str = "gfn2",
    opt: bool = True,
    charge: int = 0,
    multiplicity: int = 1,
    solvent: Optional[str] = None,
    energy_window_kcalmol: Optional[float] = 6.0,
    top_n: Optional[int] = None,
    threads: int = 2
) -> List[ConformerResult]:
    """Refines and re-ranks conformers using GFN-xTB semiempirical calculations."""
    if not conformers:
        return []

    logger.info("Avvio raffinamento semiempirico GFN-xTB su %d conformeri...", len(conformers))
    refined: List[ConformerResult] = []

    for c in conformers:
        xyz_input = c.xyz_optimized or c.xyz_initial
        if not xyz_input:
            continue

        try:
            res = run_xtb_calculation(
                xyz_str=xyz_input,
                method=method,
                opt=opt,
                charge=charge,
                multiplicity=multiplicity,
                solvent=solvent,
                threads=threads
            )

            c.xtb_energy_hartree = res["energy_hartree"]
            c.xtb_energy_kcalmol = res["energy_kcalmol"]
            c.xtb_gap_ev = res.get("gap_ev")
            c.xtb_dipole_debye = res.get("dipole_debye")

            # Update coordinates to xTB-optimized geometry
            if opt and res.get("converged"):
                c.xyz_optimized = res["xyz_optimized"]

            refined.append(c)

        except Exception as e:
            logger.warning("Conformero %d fallito in GFN-xTB: %s", c.conf_id, e)
            continue

    if not refined:
        logger.warning("Nessun conformero sopravvissuto a GFN-xTB; ripristino conformeri classici.")
        return conformers

    # Sort by xTB energy ascending
    refined.sort(key=lambda x: x.xtb_energy_hartree if x.xtb_energy_hartree is not None else 999999.0)
    min_energy = refined[0].xtb_energy_kcalmol or 0.0

    for c in refined:
        c.delta_xtb_kcalmol = round((c.xtb_energy_kcalmol or min_energy) - min_energy, 2)

    # Filter by energy window if requested
    if energy_window_kcalmol is not None and energy_window_kcalmol > 0:
        filtered = [c for c in refined if c.delta_xtb_kcalmol is not None and c.delta_xtb_kcalmol <= energy_window_kcalmol]
        if filtered:
            refined = filtered

    # Limit to top N
    if top_n is not None and top_n > 0:
        refined = refined[:top_n]

    logger.info("Raffinamento GFN-xTB completato. Conformeri selezionati per DFT: %d", len(refined))
    return refined

