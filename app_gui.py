"""Energyinator - Scientific Molecular Modeling & Quantum Conformer Pipeline.

Web graphical interface powered by Streamlit, RDKit, py3Dmol, Ketcher, and Psi4.
Modularized architecture delegating components to the gui/ package.
"""

import streamlit as st

from gui.styles import apply_theme
from gui.components.sidebar import render_sidebar
from gui.components.molecule_header import render_molecule_header
from gui.components.ketcher_editor import render_ketcher_editor
from gui.components.hero_summary import render_hero_summary
from gui.services.runner import run_calculation_job, stop_active_calculation

from gui.tabs.tab_3d import render_tab_3d
from gui.tabs.tab_pes import render_tab_pes
from gui.tabs.tab_thermo import render_tab_thermo
from gui.tabs.tab_export import render_tab_export
from gui.tabs.tab_logs import render_tab_logs


# ==============================================================================
# 1. Configurazione Pagina e Stile Scientifico Professionale
# ==============================================================================
st.set_page_config(
    page_title="Energyinator - Molecular Studio",
    page_icon=":material/science:",
    layout="wide",
    initial_sidebar_state="expanded"
)

apply_theme()


# ==============================================================================
# 2. Stato Applicazione & Sidebar
# ==============================================================================
if "calc_result" not in st.session_state:
    st.session_state.calc_result = None

cfg = render_sidebar()

if cfg["stop_calc"]:
    stop_active_calculation()


# ==============================================================================
# 3. Header Principale, Topologia 2D ed Editor Molecolare (Ketcher)
# ==============================================================================
res = st.session_state.calc_result

# Se è presente un risultato attivo, mostra i parametri e la struttura del risultato
# a meno che l'utente non abbia modificato esplicitamente lo SMILES nella sidebar
if res is not None and (not cfg["smiles"] or cfg["smiles"] == res.smiles):
    header_smiles = res.smiles
    header_method = res.method
    header_basis = res.basis
    header_solvent = res.solvent if res.solvent else "Gas / Vuoto"
    header_ff = res.force_field
else:
    header_smiles = cfg["smiles"]
    header_method = cfg["method"]
    header_basis = cfg["basis"]
    header_solvent = cfg["solv_choice"] if cfg["solvent_val"] else "Gas / Vuoto"
    header_ff = cfg["ff_choice"]

render_molecule_header(
    smiles=header_smiles,
    method=header_method,
    basis=header_basis,
    solvent=header_solvent,
    ff=header_ff,
    profile_name=cfg["profile_choice"]
)

render_ketcher_editor()


# ==============================================================================
# 4. Gestione Esecuzione Calcolo (Process Isolation)
# ==============================================================================
if cfg["avvia"]:
    run_calculation_job(cfg)


# ==============================================================================
# 5. Visualizzazione Risultati (Hero Summary Dashboard & 5 Schede Modulari)
# ==============================================================================
res = st.session_state.calc_result

if res is not None:
    render_hero_summary(res)

    tab3d, tab_pes, tab_stats, tab_table, tab_log = st.tabs([
        ":material/view_in_ar: Struttura 3D e Geometria",
        ":material/timeline: Scansione Torsionale (PES)",
        ":material/analytics: Analisi Energetica e Popolazioni",
        ":material/table_chart: Tabella Dati ed Esportazione",
        ":material/terminal: Log di Calcolo"
    ])

    with tab3d:
        render_tab_3d(res)

    with tab_pes:
        render_tab_pes(res)

    with tab_stats:
        render_tab_thermo(res)

    with tab_table:
        render_tab_export(res)

    with tab_log:
        render_tab_logs(res)
