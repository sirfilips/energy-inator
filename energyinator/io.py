import os
import logging
import pandas as pd
from typing import List, Tuple
from energyinator.models import ConformerResult, CalculationResult

logger = logging.getLogger("energyinator")


def setup_run_logger(output_dir: str) -> logging.FileHandler:
    """Configures logging to record execution details inside the specific run directory."""
    log_file = os.path.join(output_dir, "run.log")
    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setLevel(logging.INFO)
    file_formatter = logging.Formatter("[%(asctime)s][%(levelname)s] %(message)s", datefmt="%H:%M:%S")
    file_handler.setFormatter(file_formatter)
    
    root_logger = logging.getLogger("energyinator")
    root_logger.setLevel(logging.INFO)
    root_logger.addHandler(file_handler)
    return file_handler


def close_run_logger(file_handler: logging.FileHandler):
    """Closes and removes file handler from logger (vital on Windows to release file locks)."""
    if file_handler:
        root_logger = logging.getLogger("energyinator")
        try:
            file_handler.flush()
            file_handler.close()
            root_logger.removeHandler(file_handler)
        except Exception:
            pass


def save_conformer_xyz(output_dir: str, conf: ConformerResult, timestamp: str) -> str:
    """Saves the optimized 3D geometry of a conformer to an XYZ file."""
    filename = f"geom_ottimizzata_{timestamp}_{conf.conf_id}.xyz"
    filepath = os.path.join(output_dir, filename)
    
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(conf.xyz_optimized.strip() + "\n")
        
    conf.output_path = filepath
    return filepath


def export_results_excel(output_dir: str, conformers: List[ConformerResult], timestamp: str) -> str:
    """Generates a well-formatted Excel spreadsheet summarizing thermodynamic and structural data."""
    rows = []
    for c in conformers:
        xyz_name = os.path.basename(c.output_path) if c.output_path else ""
        rows.append({
            "Conformer_ID": c.conf_id,
            "Degeneracy": c.degeneracy,
            "Energy_Hartree": c.dft_energy_hartree,
            "Energy_kcalmol": c.dft_energy_kcalmol,
            "Energy_eV": c.dft_energy_ev,
            "Delta_E_kcalmol": c.delta_dft_kcalmol,
            "Delta_E_kJmol": c.delta_dft_kjmol,
            "Boltzmann_Pop_%": c.boltzmann_pop,
            "Dipole_Debye": c.dipole_debye,
            "HOMO_eV": c.homo_ev,
            "LUMO_eV": c.lumo_ev,
            "Gap_eV": c.gap_ev,
            "FF_Energy_kcalmol": c.ff_energy,
            "Delta_FF_kcalmol": c.delta_ff_kcalmol,
            "XYZ_File": xyz_name
        })

    df = pd.DataFrame(rows)
    # Sort by Delta E ascending (ground state first)
    if "Delta_E_kcalmol" in df.columns and not df.empty:
        df = df.sort_values(by="Delta_E_kcalmol", ascending=True)

    excel_path = os.path.join(output_dir, f"conformers_energies_{timestamp}.xlsx")
    df.to_excel(excel_path, index=False)
    logger.info(f"Report Excel esportato con successo in: {excel_path}")
    return excel_path

