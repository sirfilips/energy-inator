"""File export and molecular file parsing utilities for GUI with FAIR computational chemistry support."""

import io
import os
import glob
import zipfile
from datetime import datetime
from typing import Optional
from rdkit import Chem

from energyinator.models import CalculationResult
from energyinator.export_chem import (
    generate_orca_input,
    generate_gaussian_input,
    generate_multi_xyz_trajectory,
    generate_run_recipe
)


def mol_file_to_smiles(file_bytes: bytes, filename: str) -> Optional[str]:
    """Converts the uploaded file content (.mol, .sdf, .smi, .txt) into a SMILES string."""
    fname = filename.lower()
    content = file_bytes.decode("utf-8", errors="ignore")
    if fname.endswith((".smi", ".txt")):
        for line in content.splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                return line.split()[0]
    elif fname.endswith((".mol", ".sdf")):
        mol = Chem.MolFromMolBlock(content)
        if mol:
            return Chem.MolToSmiles(mol)
    return None


def create_zip_bundle(res: CalculationResult) -> bytes:
    """Creates a comprehensive in-memory ZIP bundle with:
      - Excel spreadsheet (.xlsx)
      - JSON results and FAIR run recipe (.json)
      - ORCA and Gaussian input templates (.inp, .gjf)
      - Multi-frame trajectory file (.xyz)
      - Individual conformer XYZ files
      - Executive summary text report
    """
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        if res.excel_path and os.path.exists(res.excel_path):
            zf.write(res.excel_path, arcname=os.path.basename(res.excel_path))
        if hasattr(res, "json_path") and res.json_path and os.path.exists(res.json_path):
            zf.write(res.json_path, arcname="result.json")
        elif os.path.isdir(res.output_dir) and os.path.exists(os.path.join(res.output_dir, "result.json")):
            zf.write(os.path.join(res.output_dir, "result.json"), arcname="result.json")

        # Include individual conformers
        for c in res.conformers:
            if c.output_path and os.path.exists(c.output_path):
                zf.write(c.output_path, arcname=os.path.basename(c.output_path))
            elif c.xyz_optimized:
                zf.writestr(f"conformer_{c.conf_id}.xyz", c.xyz_optimized)

        # Include FAIR computational chemistry files if available on disk or generate them
        if res.conformers:
            best_conf = res.conformers[0]
            
            # ORCA input
            orca_file = os.path.join(res.output_dir, "input_orca_minimo.inp")
            if os.path.exists(orca_file):
                zf.write(orca_file, arcname="input_orca_minimo.inp")
            else:
                zf.writestr("input_orca_minimo.inp", generate_orca_input(
                    best_conf, method=res.method, basis=res.basis, solvent=res.solvent,
                    charge=res.formal_charge, multiplicity=res.spin_multiplicity
                ))

            # Gaussian input
            gau_file = os.path.join(res.output_dir, "input_gaussian_minimo.gjf")
            if os.path.exists(gau_file):
                zf.write(gau_file, arcname="input_gaussian_minimo.gjf")
            else:
                zf.writestr("input_gaussian_minimo.gjf", generate_gaussian_input(
                    best_conf, method=res.method, basis=res.basis, solvent=res.solvent,
                    charge=res.formal_charge, multiplicity=res.spin_multiplicity
                ))

            # Multi-XYZ trajectory
            trajs = glob.glob(os.path.join(res.output_dir, "trajectory_ensemble_*.xyz"))
            if trajs and os.path.exists(trajs[0]):
                zf.write(trajs[0], arcname="trajectory_ensemble.xyz")
            else:
                zf.writestr("trajectory_ensemble.xyz", generate_multi_xyz_trajectory(res.conformers, smiles=res.smiles))

            # FAIR Run Recipe
            recipe_file = os.path.join(res.output_dir, "calculation_recipe.json")
            if os.path.exists(recipe_file):
                zf.write(recipe_file, arcname="calculation_recipe.json")
            else:
                import json
                zf.writestr("calculation_recipe.json", json.dumps(generate_run_recipe(res), indent=2))

        # Scientific report
        ens_dip = getattr(res, "ensemble_dipole_debye", None)
        ens_gap = getattr(res, "ensemble_gap_ev", None)
        s_conf = getattr(res, "conformational_entropy_cal_mol_k", None)
        g_conf = getattr(res, "conformational_free_energy_kcalmol", None)

        report_lines = [
            "=" * 65,
            "ENERGY-INATOR - REPORT SCIENTIFICO DI MODELLISTICA CONFORMAZIONALE",
            "=" * 65,
            f"Data di calcolo: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"SMILES bersaglio: {res.smiles}",
            f"Livello di Teoria: {res.method.upper()} / {res.basis}",
            f"Campo di Forze: {res.force_field.upper()}",
            f"Solvente PCM: {res.solvent if res.solvent else 'Fase gassosa (Vuoto)'}",
            f"Conformeri generati: {res.total_generated} | Unici DFT: {res.unique_dft}",
            f"Minimo Globale: {res.min_energy_hartree:.6f} Eh ({res.min_energy_kcalmol:.2f} kcal/mol)",
            f"Tempo di calcolo: {res.elapsed_time:.1f} s",
            "",
            "OSSERVABILI MACROSCOPICI D'INSIEME (BOLTZMANN WEIGHTED):",
            f"  - Momento di Dipolo medio <mu>: {f'{ens_dip:.2f} Debye' if ens_dip is not None else 'N/A'}",
            f"  - Band Gap medio <HOMO-LUMO>:   {f'{ens_gap:.2f} eV' if ens_gap is not None else 'N/A'}",
            f"  - Entropia conformazionale S_conf: {f'{s_conf:.2f} cal/(mol*K)' if s_conf is not None else 'N/A'}",
            f"  - Stabilizzazione libera (-T*S_conf): {f'{g_conf:.3f} kcal/mol' if g_conf is not None else 'N/A'}",
            "",
            "DETTAGLIO CONFORMERI:",
            f"{'Conf':<8} {'Delta E (kcal/mol)':<20} {'Popolazione (%)':<16} {'Degenerazione':<14}"
        ]
        for c in res.conformers:
            report_lines.append(f"Conf {c.conf_id:<3d} {c.delta_dft_kcalmol:<20.2f} {c.boltzmann_pop:<16.1f} {c.degeneracy:<14d}")
        
        report_lines.append("=" * 65)
        zf.writestr("report_scientifico.txt", "\n".join(report_lines))
    
    buf.seek(0)
    return buf.getvalue()
