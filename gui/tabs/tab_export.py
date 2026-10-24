"""Tab 4: Data Table and Export (Excel, XYZ, ZIP Bundle, Output Directory)."""

import os
import subprocess
import platform
from datetime import datetime
import pandas as pd
import streamlit as st

from energyinator.models import CalculationResult
from gui.services.export import create_zip_bundle
from gui.styles import st_dataframe_stretched


def open_folder(path: str):
    """Opens directory in native OS file explorer cross-platform."""
    if not os.path.isdir(path):
        return
    system = platform.system()
    try:
        if system == "Windows":
            subprocess.Popen(["explorer", os.path.normpath(path)])
        elif system == "Darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])
    except Exception as e:
        st.error(f"Impossibile aprire la cartella: {e}")


def render_tab_export(res: CalculationResult):
    """Renders the detailed results table and all download buttons."""
    st.markdown("#### Tabella Risultati Dettagliata")

    has_gibbs = any(c.gibbs_hartree is not None for c in res.conformers)

    table_rows = []
    for c in res.conformers:
        row = {
            "Conformer ID": c.conf_id,
            "Molteplicità (g_i)": c.degeneracy,
            "Energia (Hartree)": f"{c.dft_energy_hartree:.6f}",
            "Energia (kcal/mol)": f"{c.dft_energy_kcalmol:.2f}",
            "Delta E (kcal/mol)": f"{c.delta_dft_kcalmol:.2f}",
            "Delta E (kJ/mol)": f"{c.delta_dft_kjmol:.2f}",
            "Popolazione (%)": f"{c.boltzmann_pop:.1f}%"
        }
        if has_gibbs:
            row["Delta G (kcal/mol)"] = f"{c.delta_gibbs_kcalmol:.2f}" if c.delta_gibbs_kcalmol is not None else "-"
            row["ZPVE (Eh)"] = f"{c.zpve_hartree:.5f}" if c.zpve_hartree is not None else "-"
            row["Entropia S"] = f"{c.entropy_cal_mol_k:.1f}" if c.entropy_cal_mol_k is not None else "-"

        row.update({
            "Dipolo (Debye)": f"{c.dipole_debye:.2f}" if c.dipole_debye is not None else "-",
            "HOMO (eV)": f"{c.homo_ev:.2f}" if c.homo_ev is not None else "-",
            "LUMO (eV)": f"{c.lumo_ev:.2f}" if c.lumo_ev is not None else "-",
            "Gap (eV)": f"{c.gap_ev:.2f}" if c.gap_ev is not None else "-",
            "FF Energy": f"{c.ff_energy:.2f}",
            "Delta E FF": f"{c.delta_ff_kcalmol:.2f}",
            "File XYZ": os.path.basename(c.output_path) if c.output_path else ""
        })
        table_rows.append(row)

    df_display = pd.DataFrame(table_rows)
    st_dataframe_stretched(df_display, hide_index=True)

    c_exp1, c_exp2, c_exp3, c_exp4 = st.columns(4)

    if res.excel_path and os.path.exists(res.excel_path):
        with open(res.excel_path, "rb") as f:
            c_exp1.download_button(
                label="Scarica Excel (.xlsx)",
                data=f.read(),
                file_name=os.path.basename(res.excel_path),
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                icon=":material/download:"
            )

    if res.best_output_path and os.path.exists(res.best_output_path):
        with open(res.best_output_path, "r", encoding="utf-8") as f:
            c_exp2.download_button(
                label="Scarica Minimo (.xyz)",
                data=f.read(),
                file_name=os.path.basename(res.best_output_path),
                mime="text/plain",
                icon=":material/download:"
            )

    try:
        zip_data = create_zip_bundle(res)
        smi_clean = res.smiles.replace("/", "_").replace("\\", "_")
        c_exp3.download_button(
            label="Scarica Tutto (.zip)",
            data=zip_data,
            file_name=f"calcolo_{smi_clean}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip",
            mime="application/zip",
            icon=":material/folder_zip:"
        )
    except Exception as e:
        c_exp3.error(f"Errore ZIP: {e}")

    if c_exp4.button("Apri Cartella", icon=":material/folder_open:"):
        open_folder(res.output_dir)

