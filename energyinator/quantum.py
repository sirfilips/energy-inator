import logging
import os
import tempfile
import glob
import shutil
import numpy as np
import psi4
from typing import Dict, Any, Optional

logger = logging.getLogger("energyinator")

# NIST / CODATA high-precision physical constants
HARTREE_TO_KCALMOL = 627.5094740631
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KJMOL = 2625.4996395
BOHR_TO_ANGSTROM = 0.529177210903

# Solvents supported by Psi4 PCM
SOLVENT_MAP = {
    "water": "Water",
    "h2o": "Water",
    "acqua": "Water",
    "chloroform": "Chloroform",
    "cloroformio": "Chloroform",
    "dmso": "Dimethylsulfoxide",
    "methanol": "Methanol",
    "metanolo": "Methanol",
    "toluene": "Toluene",
    "acetonitrile": "Acetonitrile",
    "dcm": "Dichloromethane",
    "dichloromethane": "Dichloromethane",
    "thf": "Tetrahydrofuran",
    "ethanol": "Ethanol"
}


def setup_psi4(memory: str = "8 GB", threads: int = 8, basis: str = "6-31G(d)", log_file: Optional[str] = None):
    """Configures global Psi4 runtime settings."""
    psi4.set_memory(memory)
    psi4.set_num_threads(threads)
    psi4.set_options({"basis": basis})
    
    if log_file:
        psi4.core.set_output_file(log_file, False)
    else:
        # Cross-platform safe null stream (avoids /dev/null error on Windows)
        devnull = "nul" if os.name == "nt" else "/dev/null"
        try:
            psi4.core.set_output_file(devnull, False)
        except Exception:
            try:
                psi4.core.be_quiet()
            except Exception:
                pass


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
    # Reopen null device so subsequent calculations don't write to a closed C++ FILE* pointer
    devnull = "nul" if os.name == "nt" else "/dev/null"
    try:
        psi4.core.set_output_file(devnull, False)
    except Exception:
        pass


def configure_solvent_pcm(solvent: Optional[str]):
    """Configures implicit solvent model (IEFPCM) in Psi4 if requested."""
    if not solvent or solvent.strip().lower() in ["none", "gas", "vacuum", ""]:
        psi4.set_options({"pcm": False})
        return

    solv_clean = solvent.strip().lower()
    solv_name = SOLVENT_MAP.get(solv_clean, solvent.strip().capitalize())

    pcm_spec = f"""Units = Angstrom
Medium {{
  SolverType = IEFPCM
  Solvent = {solv_name}
}}
Cavity {{
  Type = GePol
  Area = 0.3
  Scaling = False
  RadiiSet = UFF
}}"""
    try:
        psi4.set_options({"pcm": True, "pcm_scf_type": "total"})
        psi4.pcm_helper(pcm_spec)
        logger.info(f"Modello di solvatazione implicita IEFPCM configurato con successo per: {solv_name}")
    except Exception as e:
        logger.warning(f"Impossibile inizializzare PCM per solvente '{solvent}': {e}. Continuo in fase gas.")
        psi4.set_options({"pcm": False})


def run_psi4_calculation(
    xyz_block: str,
    conf_id: int,
    charge: int = 0,
    multiplicity: int = 1,
    method: str = "b3lyp",
    basis: str = "6-31G(d)",
    sp_only: bool = False,
    solvent: Optional[str] = None,
    calc_frequencies: bool = False,
    smiles_label: str = ""
) -> Optional[Dict[str, Any]]:
    """Executes a single Psi4 calculation (geometry optimization or single point).
    
    Explicitly accounts for formal charge, spin multiplicity, solvent model, and coordinate units.
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
        
        # Configure solvent if specified
        configure_solvent_pcm(solvent)

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

        # Vibrational Frequencies & Gibbs free energy (optional)
        gibbs_hartree = None
        zpve_hartree = None
        enthalpy_hartree = None
        entropy_cal_mol_k = None

        if calc_frequencies and not sp_only:
            try:
                logger.info(f"Calcolo frequenze vibrazionali e termochimica per conformero {conf_id}...")
                psi4.frequency(effective_method, molecule=active_mol)
                if psi4.core.has_variable("GIBBS FREE ENERGY"):
                    gibbs_hartree = float(psi4.core.variable("GIBBS FREE ENERGY"))
                if psi4.core.has_variable("ZPVE"):
                    zpve_hartree = float(psi4.core.variable("ZPVE"))
                if psi4.core.has_variable("ENTHALPY"):
                    enthalpy_hartree = float(psi4.core.variable("ENTHALPY"))
                if enthalpy_hartree is not None and gibbs_hartree is not None:
                    entropy_cal_mol_k = float((enthalpy_hartree - gibbs_hartree) / 298.15 * 627509.474)
            except Exception as freq_err:
                logger.warning(f"Calcolo frequenze non riuscito per conformero {conf_id}: {freq_err}")

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
            "gibbs_hartree": gibbs_hartree,
            "zpve_hartree": zpve_hartree,
            "enthalpy_hartree": enthalpy_hartree,
            "entropy_cal_mol_k": entropy_cal_mol_k,
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


def generate_orbital_cube(
    xyz_str: str,
    orbital_type: str = "HOMO",
    method: str = "b3lyp",
    basis: str = "6-31G(d)",
    charge: int = 0,
    multiplicity: int = 1,
    solvent: Optional[str] = None,
    grid_spacing: float = 0.35
) -> Optional[Dict[str, Any]]:
    """Genera la griglia volumetrica 3D (Gaussian Cube) per un orbitale molecolare specifico (HOMO, LUMO, ecc.).
    
    Parameters
    ----------
    xyz_str : str
        Coordinate molecolari in formato XYZ.
    orbital_type : str
        Identificativo dell'orbitale: 'HOMO', 'LUMO', 'HOMO-1', 'LUMO+1', 'HOMO-2', 'LUMO+2' o indice 1-based.
    method : str
        Metodo quantistico SCF/DFT (es. 'b3lyp', 'hf', 'wb97x-d').
    basis : str
        Set di basi atomiche (es. '6-31G(d)', 'sto-3g', 'def2-svp').
    charge : int
        Carica formale totale del sistema.
    multiplicity : int
        Molteplicità di spin (2S+1).
    solvent : Optional[str]
        Solvente per il modello PCM implicito.
    grid_spacing : float
        Risoluzione spaziale dei voxel della griglia (default 0.35 Bohr).
        
    Returns
    -------
    Optional[Dict[str, Any]]
        Dizionario con metadati dell'orbitale e dati cubici per il rendering isosuperficie in py3Dmol.
    """
    try:
        devnull = "nul" if os.name == "nt" else "/dev/null"
        try:
            psi4.core.set_output_file(devnull, False)
        except Exception:
            pass

        raw_lines = xyz_str.strip().splitlines()
        if len(raw_lines) >= 2 and raw_lines[0].strip().isdigit():
            atom_lines = [l.strip() for l in raw_lines[2:] if l.strip()]
        else:
            atom_lines = [l.strip() for l in raw_lines if l.strip()]

        coords_str = "\n".join(atom_lines)
        geom_input = f"{charge} {multiplicity}\n{coords_str}\nsymmetry c1\nunits angstrom"

        mol = psi4.geometry(geom_input)
        psi4.set_options({"basis": basis})
        configure_solvent_pcm(solvent)

        effective_method = method.lower().strip()
        try:
            _, wfn = psi4.energy(effective_method, molecule=mol, return_wfn=True)
        except Exception as err:
            if "dftd3" in str(err).lower() or "resourceerror" in str(err).lower():
                logger.warning(f"Fallback su functional con dispersione nativa 'wb97x-d'...")
                effective_method = "wb97x-d"
                _, wfn = psi4.energy(effective_method, molecule=mol, return_wfn=True)
            else:
                raise err

        if wfn is None:
            raise RuntimeError("Psi4 non ha prodotto alcuna wavefunction valida per l'orbitale.")

        nalpha = wfn.nalpha()
        nmo = wfn.nmo()
        orb_label = orbital_type.strip().upper()

        if orb_label == "HOMO":
            orb_idx = nalpha
        elif orb_label == "LUMO":
            orb_idx = nalpha + 1
        elif orb_label in ["HOMO-1", "HOMO - 1"]:
            orb_idx = nalpha - 1
        elif orb_label in ["LUMO+1", "LUMO + 1"]:
            orb_idx = nalpha + 2
        elif orb_label in ["HOMO-2", "HOMO - 2"]:
            orb_idx = nalpha - 2
        elif orb_label in ["LUMO+2", "LUMO + 2"]:
            orb_idx = nalpha + 3
        else:
            try:
                orb_idx = int(orb_label)
            except ValueError:
                orb_idx = nalpha

        if orb_idx < 1 or orb_idx > nmo:
            raise ValueError(f"Indice orbitale {orb_idx} fuori intervallo (totale orbitali disponibili: {nmo}).")

        eps_a = wfn.epsilon_a()
        orbital_energy_ev = None
        if eps_a is not None and eps_a.dim() >= orb_idx:
            orbital_energy_ev = float(eps_a.get(orb_idx - 1) * HARTREE_TO_EV)

        td = tempfile.mkdtemp(prefix="psi4_cube_")
        try:
            psi4.set_options({
                "cubeprop_tasks": ["orbitals"],
                "cubeprop_orbitals": [orb_idx],
                "cubic_grid_spacing": [grid_spacing, grid_spacing, grid_spacing],
                "cubeprop_filepath": td.replace("\\", "/")
            })
            psi4.cubeprop(wfn)

            cube_files = glob.glob(os.path.join(td, "*.cube"))
            if not cube_files:
                raise FileNotFoundError(f"Nessun file .cube prodotto in {td}")

            with open(cube_files[0], "r", encoding="utf-8") as f:
                cube_data = f.read()

        finally:
            shutil.rmtree(td, ignore_errors=True)
            psi4.core.clean()

        return {
            "orbital_type": orbital_type,
            "orbital_idx": orb_idx,
            "orbital_energy_ev": orbital_energy_ev,
            "cube_data": cube_data,
            "method": effective_method,
            "basis": basis
        }

    except Exception as e:
        logger.warning(f"Errore durante generazione orbitale 3D: {repr(e)}")
        psi4.core.clean()
        return None
