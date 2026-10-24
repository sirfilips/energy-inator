"""File export and molecular file parsing utilities for GUI."""

import io
import os
import zipfile
from datetime import datetime
from typing import Optional
from rdkit import Chem

from energyinator.models import CalculationResult


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
    """Creates an in-memory ZIP bundle with XYZ geometries, Excel table, JSON, and summary report."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        if res.excel_path and os.path.exists(res.excel_path):
            zf.write(res.excel_path, arcname=os.path.basename(res.excel_path))
        if hasattr(res, "json_path") and res.json_path and os.path.exists(res.json_path):
            zf.write(res.json_path, arcname="result.json")
        elif os.path.isdir(res.output_dir) and os.path.exists(os.path.join(res.output_dir, "result.json")):
            zf.write(os.path.join(res.output_dir, "result.json"), arcname="result.json")

        for c in res.conformers:
            if c.output_path and os.path.exists(c.output_path):
                zf.write(c.output_path, arcname=os.path.basename(c.output_path))
            elif c.xyz_optimized:
                zf.writestr(f"conformer_{c.conf_id}.xyz", c.xyz_optimized)

        report_lines = [
            "=" * 60,
            "ENERGY-INATOR - REPORT DI CALCOLO MOLECOLARE",
            "=" * 60,
            f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"SMILES Target: {res.smiles}",
            f"Metodo / Base: {res.method.upper()} / {res.basis}",
            f"Campo di Forze: {res.force_field.upper()}",
            f"Solvente PCM: {res.solvent if res.solvent else 'Fase gassosa'}",
            f"Conformeri generati: {res.total_generated} | Unici DFT: {res.unique_dft}",
            f"Minimo Globale: {res.min_energy_hartree:.6f} Eh ({res.min_energy_kcalmol:.2f} kcal/mol)",
            f"Tempo di calcolo: {res.elapsed_time:.1f} s",
            "",
            "DETTAGLIO CONFORMERI:",
            f"{'Conf':<8} {'Delta E (kcal/mol)':<20} {'Popolazione (%)':<16} {'Degenerazione':<14}"
        ]
        for c in res.conformers:
            report_lines.append(f"Conf {c.conf_id:<3d} {c.delta_dft_kcalmol:<20.2f} {c.boltzmann_pop:<16.1f} {c.degeneracy:<14d}")
        zf.writestr("report_sintetico.txt", "\n".join(report_lines))
    buf.seek(0)
    return buf.getvalue()

