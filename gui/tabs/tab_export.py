"""Tab 4: Data Table and Computational Chemistry Interoperability Export (ORCA, Gaussian, Trajectory, Recipe, ZIP)."""

import os
import json
import subprocess
import platform
from datetime import datetime
import pandas as pd
import streamlit as st

from energyinator.models import CalculationResult
from energyinator.export_chem import (
    generate_orca_input,
    generate_gaussian_input,
    generate_multi_xyz_trajectory,
    generate_run_recipe
)
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
    """Renders the detailed results table and all download buttons including HPC exports."""
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

        if getattr(res, "used_xtb", False):
            row["Delta E xTB (kcal/mol)"] = f"{c.delta_xtb_kcalmol:.2f}" if c.delta_xtb_kcalmol is not None else "-"

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

    st.markdown("---")
    st.markdown("#### Esportazione e Interoperabilità per Chimica Computazionale (HPC & FAIR)")
    st.caption("Scarica geometrie, tabelle e template di calcolo pronti per i cluster o per software di visualizzazione molecolare:")

    col_std1, col_std2, col_std3, col_std4 = st.columns(4)

    # 1. Excel
    if res.excel_path and os.path.exists(res.excel_path):
        with open(res.excel_path, "rb") as f:
            col_std1.download_button(
                label="Tabella Excel (.xlsx)",
                data=f.read(),
                file_name=os.path.basename(res.excel_path),
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                icon=":material/table_chart:",
                use_container_width=True
            )

    # 2. Minimo XYZ
    if res.best_output_path and os.path.exists(res.best_output_path):
        with open(res.best_output_path, "r", encoding="utf-8") as f:
            col_std2.download_button(
                label="Minimo Globale (.xyz)",
                data=f.read(),
                file_name=os.path.basename(res.best_output_path),
                mime="text/plain",
                icon=":material/deployed_code:",
                use_container_width=True
            )

    # 3. Trajectory Multi-XYZ
    if res.conformers:
        traj_data = generate_multi_xyz_trajectory(res.conformers, smiles=res.smiles)
        col_std3.download_button(
            label="Insieme Multi-XYZ (.xyz)",
            data=traj_data,
            file_name=f"ensemble_trajectory_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xyz",
            mime="text/plain",
            icon=":material/movie:",
            help="Traiettoria con tutti i conformeri numerati, visualizzabile come animazione/ensemble in PyMOL, VMD o ChimeraX.",
            use_container_width=True
        )

    # 4. ZIP All-inclusive
    try:
        zip_data = create_zip_bundle(res)
        smi_clean = res.smiles.replace("/", "_").replace("\\", "_")
        col_std4.download_button(
            label="Pacchetto Completo (.zip)",
            data=zip_data,
            file_name=f"calcolo_{smi_clean}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip",
            mime="application/zip",
            icon=":material/folder_zip:",
            use_container_width=True
        )
    except Exception as e:
        col_std4.error(f"Errore ZIP: {e}")

    # Second row: ORCA, Gaussian, Run Recipe, Open Folder
    col_hpc1, col_hpc2, col_hpc3, col_hpc4 = st.columns(4)

    best_conf = res.conformers[0] if res.conformers else None
    if best_conf:
        # ORCA template
        orca_content = generate_orca_input(
            best_conf,
            method=res.method,
            basis=res.basis,
            solvent=res.solvent,
            charge=res.formal_charge,
            multiplicity=res.spin_multiplicity
        )
        col_hpc1.download_button(
            label="Input per ORCA (.inp)",
            data=orca_content,
            file_name=f"orca_opt_conf_{best_conf.conf_id}.inp",
            mime="text/plain",
            icon=":material/terminal:",
            help="Template di input pronto all'uso per il programma quantistico ORCA (ottimizzazione + frequenze).",
            use_container_width=True
        )

        # Gaussian template
        gau_content = generate_gaussian_input(
            best_conf,
            method=res.method,
            basis=res.basis,
            solvent=res.solvent,
            charge=res.formal_charge,
            multiplicity=res.spin_multiplicity
        )
        col_hpc2.download_button(
            label="Input per Gaussian (.gjf)",
            data=gau_content,
            file_name=f"gaussian_opt_conf_{best_conf.conf_id}.gjf",
            mime="text/plain",
            icon=":material/code:",
            help="Template pronto per Gaussian 09/16 con coordinate del minimo e specifiche di solvente/base.",
            use_container_width=True
        )

    # FAIR Recipe
    recipe_dict = generate_run_recipe(res)
    col_hpc3.download_button(
        label="FAIR Run Recipe (.json)",
        data=json.dumps(recipe_dict, indent=2),
        file_name=f"calculation_recipe_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
        mime="application/json",
        icon=":material/receipt_long:",
        help="Metadati scientifici completi con versioni software, hardware e parametri per la riproducibilità totale del calcolo.",
        use_container_width=True
    )

    if col_hpc4.button("Apri Cartella Output", icon=":material/folder_open:", use_container_width=True):
        open_folder(res.output_dir)
