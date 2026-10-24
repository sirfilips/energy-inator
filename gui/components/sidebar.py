"""Sidebar configuration panel with presets, molecular inputs, and parameters."""

import os
import streamlit as st
from gui.services.pubchem import fetch_pubchem_smiles
from gui.services.export import mol_file_to_smiles
from energyinator.models import CalculationResult

PROFILES = {
    "Profilo Legacy": {
        "n_confs": 50,
        "rmsd": 0.20,
        "post_rmsd": 0.0,
        "ff": "uff",
        "enable_window": False,
        "window": 0.0,
        "method": "b3lyp",
        "basis": "6-31G(d)",
        "solvent": "Nessuno / Gas (Vuoto)",
        "sp_only": False,
        "calc_freq": False,
        "desc": "Parametri identici alla v1.0 per garantire riproducibilità e comparabilità totale (B3LYP/6-31G(d), UFF)."
    },
    "Standard / Bilanciato": {
        "n_confs": 30,
        "rmsd": 0.20,
        "post_rmsd": 0.15,
        "ff": "mmff94",
        "enable_window": True,
        "window": 10.0,
        "method": "b3lyp",
        "basis": "6-31G(d)",
        "solvent": "Nessuno / Gas (Vuoto)",
        "sp_only": False,
        "calc_freq": False,
        "desc": "Miglior rapporto accuratezza/efficienza per composti organici (B3LYP/6-31G(d), MMFF94, finestra 10 kcal/mol)."
    },
    "Solvatazione Acquosa (PCM)": {
        "n_confs": 30,
        "rmsd": 0.20,
        "post_rmsd": 0.15,
        "ff": "mmff94",
        "enable_window": True,
        "window": 10.0,
        "method": "b3lyp",
        "basis": "6-31G(d)",
        "solvent": "Acqua (Water)",
        "sp_only": False,
        "calc_freq": False,
        "desc": "Simulazione in ambiente acquoso con modello di continuo dielettrico IEFPCM."
    },
    "Alta Precisione (Dispersion-Corrected)": {
        "n_confs": 40,
        "rmsd": 0.15,
        "post_rmsd": 0.15,
        "ff": "mmff94",
        "enable_window": True,
        "window": 12.0,
        "method": "wb97x-d",
        "basis": "def2-SVP",
        "solvent": "Nessuno / Gas (Vuoto)",
        "sp_only": False,
        "calc_freq": False,
        "desc": "Inclusione di dispersione di London per anelli e interazioni non covalenti (omegaB97X-D / def2-SVP)."
    },
    "Screening Rapido (Single-Point)": {
        "n_confs": 15,
        "rmsd": 0.25,
        "post_rmsd": 0.0,
        "ff": "mmff94",
        "enable_window": True,
        "window": 8.0,
        "method": "b3lyp",
        "basis": "3-21G",
        "solvent": "Nessuno / Gas (Vuoto)",
        "sp_only": True,
        "calc_freq": False,
        "desc": "Valutazione rapida con Single-Point e base compatta."
    },
    "Personalizzato": {
        "n_confs": 30,
        "rmsd": 0.20,
        "post_rmsd": 0.15,
        "ff": "mmff94",
        "enable_window": True,
        "window": 10.0,
        "method": "b3lyp",
        "basis": "6-31G(d)",
        "solvent": "Nessuno / Gas (Vuoto)",
        "sp_only": False,
        "calc_freq": False,
        "desc": "Configurazione manuale libera di tutti i parametri."
    }
}

EXAMPLES = {
    "Personalizzato": "",
    "Butano (Anti vs Gauche)": "CCCC",
    "Cicloesano (Sedia vs Barca)": "C1CCCCC1",
    "Bifenile (Angolo Diedro ~38°)": "c1ccccc1-c2ccccc2",
    "n-Esano (Equilibri Conformazionali)": "CCCCCC",
    "Etanolo": "CCO",
    "Caffeina": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}

SOLVENT_MAP = {
    "Nessuno / Gas (Vuoto)": None,
    "Acqua (Water)": "water",
    "Cloroformio (Chloroform)": "chloroform",
    "DMSO (Dimethylsulfoxide)": "dmso",
    "Metanolo (Methanol)": "methanol",
    "Toluene": "toluene",
    "Acetonitrile": "acetonitrile",
    "DCM (Dichloromethane)": "dcm"
}


def render_sidebar():
    """Renders the sidebar and returns configuration parameters dictionary."""
    st.sidebar.markdown("### Configurazione Calcolo")

    profile_choice = st.sidebar.selectbox("Profilo di calcolo:", list(PROFILES.keys()), index=0)
    pref = PROFILES[profile_choice]
    st.sidebar.caption(f"*{pref['desc']}*")

    st.sidebar.markdown("---")

    # Example picker
    if "pending_example_selection" in st.session_state:
        st.session_state["sel_mol_example"] = st.session_state.pop("pending_example_selection")
    elif "sel_mol_example" not in st.session_state:
        st.session_state["sel_mol_example"] = "Personalizzato"

    scelta_esempio = st.sidebar.selectbox("Esempi molecolari:", list(EXAMPLES.keys()), key="sel_mol_example")
    default_smiles = EXAMPLES[scelta_esempio] if scelta_esempio != "Personalizzato" else "CCCC"

    if "active_smiles" not in st.session_state:
        st.session_state["active_smiles"] = default_smiles

    if scelta_esempio != "Personalizzato":
        if st.session_state.get("prev_selected_example") != scelta_esempio:
            st.session_state["prev_selected_example"] = scelta_esempio
            st.session_state["active_smiles"] = default_smiles
            st.session_state["last_ketcher_applied"] = default_smiles
            st.session_state["ketcher_version"] = st.session_state.get("ketcher_version", 0) + 1

    # PubChem search
    with st.sidebar.expander("Trova Molecola per Nome (PubChem)", expanded=False):
        nome_query = st.text_input("Nome comune:", placeholder="es. aspirina, paracetamolo", key="input_pubchem_name")
        btn_search_pubchem = st.button("Cerca su PubChem", icon=":material/search:", key="btn_pubchem")
        if btn_search_pubchem and nome_query.strip():
            with st.spinner(f"Ricerca in corso per '{nome_query.strip()}'..."):
                found_smi = fetch_pubchem_smiles(nome_query.strip())
            if found_smi:
                st.session_state["active_smiles"] = found_smi
                st.session_state["last_ketcher_applied"] = found_smi
                st.session_state["pending_example_selection"] = "Personalizzato"
                st.session_state["prev_selected_example"] = "Personalizzato"
                st.session_state["ketcher_version"] = st.session_state.get("ketcher_version", 0) + 1
                st.session_state["pubchem_search_result"] = f"Trovato: {found_smi}"
                st.rerun()
            else:
                st.error(f"Nessun risultato trovato per '{nome_query}'.")

        if "pubchem_search_result" in st.session_state:
            st.success(st.session_state["pubchem_search_result"])

    # File upload
    with st.sidebar.expander("Importa da File (.mol, .sdf, .smi)", expanded=False):
        up_file = st.file_uploader("Carica struttura:", type=["mol", "sdf", "smi", "txt"], key="mol_file_uploader")
        if up_file is not None and st.session_state.get("last_uploaded_filename") != up_file.name:
            file_smi = mol_file_to_smiles(up_file.read(), up_file.name)
            if file_smi:
                st.session_state["last_uploaded_filename"] = up_file.name
                st.session_state["active_smiles"] = file_smi
                st.session_state["last_ketcher_applied"] = file_smi
                st.session_state["pending_example_selection"] = "Personalizzato"
                st.session_state["prev_selected_example"] = "Personalizzato"
                st.session_state["ketcher_version"] = st.session_state.get("ketcher_version", 0) + 1
                st.success(f"Caricato: {file_smi}")
                st.rerun()
            else:
                st.error("Impossibile convertire il file in una molecola valida.")

    # Target SMILES input box
    if "smiles_input_box" not in st.session_state or st.session_state["smiles_input_box"] != st.session_state["active_smiles"]:
        st.session_state["smiles_input_box"] = st.session_state["active_smiles"]

    def on_smiles_change():
        st.session_state["active_smiles"] = st.session_state["smiles_input_box"]
        st.session_state["last_ketcher_applied"] = st.session_state["smiles_input_box"]

    smiles_input = st.sidebar.text_input(
        "Codice SMILES target:",
        key="smiles_input_box",
        on_change=on_smiles_change
    )

    # Archived sessions loader
    output_base = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output")
    available_runs = []
    if os.path.isdir(output_base):
        for entry in sorted(os.listdir(output_base), reverse=True):
            full_dir = os.path.join(output_base, entry)
            if os.path.isdir(full_dir) and os.path.exists(os.path.join(full_dir, "result.json")):
                available_runs.append(entry)

    def on_load_archived_session():
        sel_run = st.session_state.get("sel_archived_run")
        if not sel_run:
            return
        jpath = os.path.join(output_base, sel_run, "result.json")
        try:
            loaded_res = CalculationResult.from_json(jpath)
            st.session_state.calc_result = loaded_res
            if loaded_res.smiles:
                st.session_state["active_smiles"] = loaded_res.smiles
                st.session_state["smiles_input_box"] = loaded_res.smiles
                st.session_state["last_ketcher_applied"] = loaded_res.smiles
                st.session_state["pending_example_selection"] = "Personalizzato"
                st.session_state["prev_selected_example"] = "Personalizzato"
                st.session_state["sel_mol_example"] = "Personalizzato"
                st.session_state["ketcher_version"] = st.session_state.get("ketcher_version", 0) + 1
            if "pubchem_search_result" in st.session_state:
                del st.session_state["pubchem_search_result"]
            st.session_state["archived_session_msg"] = f"Caricata sessione: {sel_run}"
        except Exception as e:
            st.session_state["archived_session_err"] = f"Errore caricamento sessione: {e}"

    if available_runs:
        with st.sidebar.expander("Sessioni Archiviate", expanded=False):
            st.selectbox("Seleziona esecuzione:", available_runs, key="sel_archived_run")
            st.button(
                "Carica nei Risultati",
                icon=":material/history:",
                key="btn_load_archived",
                on_click=on_load_archived_session
            )
            if "archived_session_msg" in st.session_state:
                st.success(st.session_state.pop("archived_session_msg"))
            if "archived_session_err" in st.session_state:
                st.error(st.session_state.pop("archived_session_err"))

    # Scientific Benchmark Suite
    with st.sidebar.expander("Validazione & Benchmark Scientifico", expanded=False):
        st.caption("Suite di validazione accademica: calcola l'accuratezza contro equilibri e barriere di riferimento di letteratura.")
        if st.button("Esegui Benchmark Rapido", icon=":material/speed:", key="btn_run_benchmark"):
            with st.spinner("Esecuzione benchmark su molecole modello (Butano, 1,2-Etandiolo, Cicloesano)..."):
                from energyinator.benchmark import run_benchmark
                bench_data = run_benchmark(method="b3lyp", basis="3-21G", sp_only=True)
                st.session_state["benchmark_results"] = bench_data
                st.success("Benchmark completato!")

        if "benchmark_results" in st.session_state:
            b_res = st.session_state["benchmark_results"]
            st.dataframe(
                b_res["df_results"][["System", "Status", "Delta_E_Calc", "Delta_E_Ref", "Error", "Elapsed_s"]],
                hide_index=True
            )
            st.caption(f"**MAE:** {b_res['mean_absolute_error_kcalmol']} kcal/mol &bull; **Tempo:** {b_res['total_elapsed_s']} s")
            csv_data = b_res["df_results"].to_csv(index=False).encode('utf-8')
            st.download_button(
                "Esporta Tabella Tesi (CSV)",
                data=csv_data,
                file_name="benchmark_energyinator.csv",
                mime="text/csv",
                icon=":material/download:",
                use_container_width=True
            )

    st.sidebar.markdown("---")

    ff_list = ["mmff94", "mmff94s", "uff"]
    ff_idx = ff_list.index(pref["ff"]) if pref["ff"] in ff_list else 0

    method_list = ["b3lyp", "wb97x-d", "pbe", "hf", "mp2"]
    method_idx = method_list.index(pref["method"]) if pref["method"] in method_list else 0

    basis_list = ["6-31G(d)", "6-311+G(d,p)", "3-21G", "sto-3g", "def2-SVP"]
    basis_idx = basis_list.index(pref["basis"]) if pref["basis"] in basis_list else 0

    solv_keys = list(SOLVENT_MAP.keys())
    solv_idx = solv_keys.index(pref.get("solvent", "Nessuno / Gas (Vuoto)")) if pref.get("solvent") in solv_keys else 0

    # Conformational Sampling Parameters
    with st.sidebar.expander("Campionamento Conformazionale (RDKit)", expanded=True):
        n_confs = st.slider("Conformeri generati (ETKDGv3):", min_value=5, max_value=200, value=pref["n_confs"], step=5, key=f"nc_{profile_choice}")
        rmsd_thresh = st.number_input("Soglia RMSD classica (Å):", min_value=0.05, max_value=2.0, value=pref["rmsd"], step=0.05, key=f"rmsd_{profile_choice}")
        post_rmsd = st.number_input("Deduplicazione post-DFT (Å):", min_value=0.0, max_value=1.0, value=pref.get("post_rmsd", 0.15), step=0.05, help="Rimuove geometrie collassate nello stesso minimo quantistico.", key=f"p_rmsd_{profile_choice}")
        ff_choice = st.selectbox("Campo di forze (pre-rilassamento):", ff_list, index=ff_idx, key=f"ff_{profile_choice}")
        enable_window = st.checkbox("Filtro Energy Window FF", value=pref["enable_window"], key=f"ew_{profile_choice}")
        energy_window = st.slider("Finestra energetica FF (kcal/mol):", 1.0, 30.0, pref["window"] if pref["window"] > 0 else 10.0, 1.0, key=f"ew_val_{profile_choice}") if enable_window else 0.0
        limit_top_n = st.checkbox("Limita conformeri inviati a DFT", value=False)
        top_n_val = st.number_input("Max conformeri per DFT:", min_value=1, max_value=50, value=10) if limit_top_n else None

    # Quantum Chemistry Parameters
    with st.sidebar.expander("Motore Quantistico (Psi4)", expanded=True):
        method = st.selectbox("Metodo Quantistico:", method_list, index=method_idx, key=f"m_{profile_choice}")
        basis = st.selectbox("Basis Set:", basis_list, index=basis_idx, key=f"b_{profile_choice}")
        solv_choice = st.selectbox("Solvatazione Implicita (PCM):", solv_keys, index=solv_idx, key=f"s_{profile_choice}")
        solvent_val = SOLVENT_MAP[solv_choice]
        threads = st.number_input("Thread CPU:", min_value=1, max_value=os.cpu_count() or 8, value=min(8, os.cpu_count() or 4))
        memory = st.selectbox("Memoria RAM:", ["4 GB", "8 GB", "16 GB", "32 GB"], index=1)
        sp_only = st.checkbox("Solo Single-Point (senza ottimizzazione geometrica)", value=pref["sp_only"], key=f"sp_{profile_choice}")

    # Electronic and Thermodynamics
    with st.sidebar.expander("Parametri Elettronici e Termodinamici", expanded=False):
        temp_k = st.number_input("Temperatura Boltzmann (K):", value=298.15, step=5.0)
        calc_freq = st.checkbox("Calcola Termochimica (Gibbs Delta G, ZPVE)", value=pref.get("calc_freq", False), help="Esegue calcolo di frequenze vibrazionali sul ground state.")
        manual_charge = st.number_input("Carica formale manuale:", value=0)
        override_charge = st.checkbox("Forza carica manuale", value=False)
        manual_mult = st.number_input("Molteplicità di spin manuale:", min_value=1, value=1)
        override_mult = st.checkbox("Forza molteplicità manuale", value=False)

    col_btn1, col_btn2 = st.sidebar.columns([2.5, 1])
    avvia = col_btn1.button("AVVIA CALCOLO", icon=":material/play_arrow:", type="primary", use_container_width=True)
    stop_calc = col_btn2.button("STOP", icon=":material/stop:", type="secondary", use_container_width=True, help="Interrompe il calcolo in corso")

    return {
        "avvia": avvia,
        "stop_calc": stop_calc,
        "profile_choice": profile_choice,
        "smiles": smiles_input.strip(),
        "n_confs": n_confs,
        "rmsd_thresh": rmsd_thresh,
        "post_rmsd": post_rmsd,
        "ff_choice": ff_choice,
        "enable_window": enable_window,
        "energy_window": energy_window,
        "limit_top_n": limit_top_n,
        "top_n_val": top_n_val,
        "method": method,
        "basis": basis,
        "solv_choice": solv_choice,
        "solvent_val": solvent_val,
        "threads": threads,
        "memory": memory,
        "sp_only": sp_only,
        "temp_k": temp_k,
        "calc_freq": calc_freq,
        "manual_charge": manual_charge,
        "override_charge": override_charge,
        "manual_mult": manual_mult,
        "override_mult": override_mult
    }

