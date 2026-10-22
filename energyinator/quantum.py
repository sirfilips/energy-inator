import logging
import os
import numpy as np
import psi4
from typing import Dict, Any, Optional

logger = logging.getLogger("energyinator")

# NIST / CODATA high-precision physical constants
HARTREE_TO_KCALMOL = 627.5094740631
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KJMOL = 2625.4996395
BOHR_TO_ANGSTROM = 0.529177210903


def setup_psi4(memory: str = "8 GB", threads: int = 8, basis: str = "6-31G(d)", log_file: Optional[str] = None):
    """Configures global Psi4 runtime settings."""
    psi4.set_memory(memory)
    psi4.set_num_threads(threads)
    psi4.set_options({"basis": basis})
    
    if log_file:
        psi4.core.set_output_file(log_file, False)
    else:
        # Default silent Psi4 output if not specifically requested
        psi4.core.be_quiet()


def close_psi4():
    """Flushes and releases Psi4 output files and caches (important on Windows)."""
    try:
        psi4.core.close_outfile()
    except Exception:
        pass
    try:
        psi4.core.clean()
    except Exception:
        pass


def run_psi4_calculation(
    xyz_block: str,
    conf_id: int,
    charge: int = 0,
    multiplicity: int = 1,
    method: str = "b3lyp",
    basis: str = "6-31G(d)",
    sp_only: bool = False,
    smiles_label: str = ""
) -> Optional[Dict[str, Any]]:
    """Executes a single Psi4 calculation (geometry optimization or single point).
    
    Explicitly accounts for formal charge, spin multiplicity, and coordinate units.
    """
    try:
        # Correctly parse standard XYZ block (Line 0: natoms, Line 1: comment/empty, Line 2+: atoms)
        raw_lines = xyz_block.splitlines()
        if len(raw_lines) >= 2 and raw_lines[0].strip().isdigit():
            atom_lines = [l.strip() for l in raw_lines[2:] if l.strip()]
        else:
            atom_lines = [l.strip() for l in raw_lines if l.strip()]

        coords_str = "\n".join(atom_lines)
        geom_input = f"{charge} {multiplicity}\n{coords_str}\nsymmetry c1\nunits angstrom"

        mol = psi4.geometry(geom_input)
        
        # Execute calculation
        effective_method = method.lower().strip()
        wfn = None
        try:
            if sp_only:
                energy_hartree, wfn = psi4.energy(effective_method, molecule=mol, return_wfn=True)
            else:
                energy_hartree, wfn = psi4.optimize(effective_method, molecule=mol, return_wfn=True)
        except Exception as err:
            # Handle missing external dispersion programs (e.g. s-dftd3) gracefully
            if "dftd3" in str(err).lower() or "resourceerror" in str(err).lower():
                logger.warning(
                    f"Metodo '{effective_method}' richiede programma dftd3 esterno non configurato. "
                    "Fallback su functional con dispersione nativa 'wb97x-d'..."
                )
                effective_method = "wb97x-d"
                if sp_only:
                    energy_hartree, wfn = psi4.energy(effective_method, molecule=mol, return_wfn=True)
                else:
                    energy_hartree, wfn = psi4.optimize(effective_method, molecule=mol, return_wfn=True)
            else:
                raise err

        active_mol = psi4.core.get_active_molecule()
        if active_mol is None:
            active_mol = mol

        natoms = active_mol.natom()
        coords = []
        for i in range(natoms):
            sym = active_mol.symbol(i)
            x = active_mol.x(i) * BOHR_TO_ANGSTROM
            y = active_mol.y(i) * BOHR_TO_ANGSTROM
            z = active_mol.z(i) * BOHR_TO_ANGSTROM
            coords.append(f"{sym:<2} {x:14.8f} {y:14.8f} {z:14.8f}")

        # Form formatted XYZ representation
        header = f"{natoms}\n{smiles_label} | Conf {conf_id} | E = {energy_hartree:.8f} Eh"
        xyz_optimized = f"{header}\n" + "\n".join(coords)

        energy_kcalmol = energy_hartree * HARTREE_TO_KCALMOL
        energy_ev = energy_hartree * HARTREE_TO_EV
        energy_kjmol = energy_hartree * HARTREE_TO_KJMOL

        # Extract electronic properties (Dipole, HOMO, LUMO, Gap)
        dipole_debye = None
        homo_ev = None
        lumo_ev = None
        gap_ev = None
        try:
            if psi4.core.has_variable("CURRENT DIPOLE"):
                dip_vec = psi4.core.variable("CURRENT DIPOLE")
                dipole_debye = float(np.linalg.norm(dip_vec) * 2.54174623)
            
            if wfn is not None:
                eps_a = wfn.epsilon_a()
                nalpha = wfn.nalpha()
                if eps_a is not None and eps_a.dim() > nalpha and nalpha > 0:
                    homo_ev = float(eps_a.get(nalpha - 1) * HARTREE_TO_EV)
                    lumo_ev = float(eps_a.get(nalpha) * HARTREE_TO_EV)
                    gap_ev = float(lumo_ev - homo_ev)
        except Exception as prop_err:
            logger.debug(f"Proprietà elettroniche non calcolabili per conf {conf_id}: {prop_err}")

        psi4.core.clean()

        return {
            "conf_id": conf_id,
            "energy_hartree": energy_hartree,
            "energy_kcalmol": energy_kcalmol,
            "energy_ev": energy_ev,
            "energy_kjmol": energy_kjmol,
            "xyz_optimized": xyz_optimized,
            "method_used": effective_method,
            "dipole_debye": dipole_debye,
            "homo_ev": homo_ev,
            "lumo_ev": lumo_ev,
            "gap_ev": gap_ev
        }

    except Exception as e:
        logger.warning(f"Errore durante calcolo Psi4 per conformero {conf_id}: {repr(e)}")
        psi4.core.clean()
        return None

