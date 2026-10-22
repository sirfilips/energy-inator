"""Energyinator - Scientific Molecular Modeling & Quantum Conformer Pipeline.

Web graphical interface powered by Streamlit, RDKit, py3Dmol, and Psi4.
"""

import os
import io
import re
import sys
import math
import subprocess
from datetime import datetime
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import py3Dmol
import streamlit.components.v1 as components

from rdkit import Chem
from rdkit.Chem import Draw

from energyinator.config import CalculationConfig
from energyinator.models import CalculationResult, ConformerResult
from energyinator.scan import run_dihedral_scan


# ==============================================================================
# 1. Configurazione Pagina e Stile Scientifico Professionale
# ==============================================================================
st.set_page_config(
    page_title="Energyinator - Molecular Studio",
    page_icon=":material/science:",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Design scientifico moderno: tipografia pulita, palette neutra e layout arioso
st.markdown("""
<style>
    .stApp {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    }
    .main-title {
        font-size: 1.85rem;
        font-weight: 700;
        letter-spacing: -0.02em;
        color: #3b82f6;
        margin-bottom: 2px;
    }
    .sub-title {
        font-size: 0.95rem;
        color: #94a3b8;
        font-weight: 400;
        margin-bottom: 18px;
    }
    .info-card {
        background-color: rgba(255, 255, 255, 0.02);
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 8px;
        padding: 14px 18px;
        margin-bottom: 18px;
    }
    .metric-badge {
        display: inline-block;
        padding: 4px 10px;
        border-radius: 4px;
        font-size: 0.8rem;
        font-weight: 600;
        background: rgba(59, 130, 246, 0.12);
        color: #60a5fa;
        border: 1px solid rgba(59, 130, 246, 0.25);
        margin-right: 6px;
    }
    .tag-min {
        background: rgba(34, 197, 94, 0.12);
        color: #4ade80;
        border: 1px solid rgba(34, 197, 94, 0.3);
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 6px;
        border-bottom: 1px solid rgba(255, 255, 255, 0.08);
        margin-bottom: 16px;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 10px 18px;
        font-weight: 500;
        font-size: 0.92rem;
        border-radius: 6px 6px 0 0;
    }
    .stButton>button {
        border-radius: 6px;
        font-weight: 600;
        font-size: 0.9rem;
        letter-spacing: 0.01em;
        transition: all 0.2s ease;
    }
</style>
""", unsafe_allow_html=True)


# ==============================================================================
# 2. Funzioni Ausiliarie per Geometria 3D & Allineamento Kabsch
# ==============================================================================
def parse_xyz(xyz_str: str):
    """Estrae i simboli atomici e la matrice delle coordinate NumPy da un blocco XYZ."""
    lines = xyz_str.strip().splitlines()
    if len(lines) >= 2 and lines[0].strip().isdigit():
        lines = lines[2:]
    symbols = []
    coords = []
    for l in lines:
        parts = l.split()
        if len(parts) >= 4:
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords)


def calc_distance(p1, p2):
    return float(np.linalg.norm(p1 - p2))


def calc_angle(p1, p2, p3):
    v1 = p1 - p2
    v2 = p3 - p2
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom < 1e-9:
        return 0.0
    cosine = np.dot(v1, v2) / denom
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def calc_dihedral(p1, p2, p3, p4):
    b0 = -1.0 * (p2 - p1)
    b1 = p3 - p2
    b2 = p4 - p3
    norm_b1 = np.linalg.norm(b1)
    if norm_b1 < 1e-9:
        return 0.0
    b1 /= norm_b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def kabsch_align_xyz(ref_xyz_str: str, target_xyz_str: str) -> str:
    """Allinea e sovrappone le coordinate target su quelle di riferimento con l'algoritmo di Kabsch."""
    symbols_ref, P = parse_xyz(ref_xyz_str)
    symbols_tgt, Q = parse_xyz(target_xyz_str)

    if len(P) != len(Q) or len(P) == 0:
        return target_xyz_str

    p_center = np.mean(P, axis=0)
    q_center = np.mean(Q, axis=0)

    P_c = P - p_center
    Q_c = Q - q_center

    H = np.dot(Q_c.T, P_c)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = np.dot(Vt.T, U.T)

    Q_aligned = np.dot(Q_c, R.T) + p_center

    lines = [f"{len(symbols_tgt)}", "Aligned conformer"]
    for sym, (x, y, z) in zip(symbols_tgt, Q_aligned):
        lines.append(f"{sym:<2} {x:14.8f} {y:14.8f} {z:14.8f}")
    return "\n".join(lines)


def render_3d_mol(xyz_str, style='ball_and_stick', show_labels=True, spin=False, height=580):
    """Genera e renderizza il visualizzatore 3D WebGL con py3Dmol."""
    view = py3Dmol.view(width=None, height=height)
    view.addModel(xyz_str, "xyz")

    if style == 'stick':
        view.setStyle({'stick': {'radius': 0.16}})
    elif style == 'ball_and_stick':
        view.setStyle({'stick': {'radius': 0.12}, 'sphere': {'scale': 0.26}})
    elif style == 'sphere':
        view.setStyle({'sphere': {'scale': 0.8}})
    elif style == 'wireframe':
        view.setStyle({'line': {}})

    if show_labels:
        symbols, coords = parse_xyz(xyz_str)
        for idx, (sym, (x, y, z)) in enumerate(zip(symbols, coords)):
            view.addLabel(
                f"{sym}{idx}",
                {
                    'position': {'x': x, 'y': y, 'z': z},
                    'backgroundColor': 'rgba(15, 23, 42, 0.75)',
                    'fontColor': '#f8fafc',
                    'fontSize': 11
                }
            )

    view.zoomTo()
    if spin:
        view.spin(True)
    html = view._make_html()
    components.html(html, height=height + 20)


def render_3d_mol_overlay(ref_xyz_str: str, active_xyz_str: str, height=580):
    """Renderizza due conformeri sovrapposti in 3D: riferimento in ciano trasparente e attivo in CPK."""
    aligned_active = kabsch_align_xyz(ref_xyz_str, active_xyz_str)
    view = py3Dmol.view(width=None, height=height)

    # Modello 0: Minimo Globale (Riferimento ciano)
    view.addModel(ref_xyz_str, "xyz")
    view.setStyle({'model': 0}, {
        'stick': {'color': '#00bcd4', 'opacity': 0.40, 'radius': 0.12},
        'sphere': {'color': '#00bcd4', 'opacity': 0.40, 'scale': 0.22}
    })

    # Modello 1: Conformero attivo allineato (CPK standard)
    view.addModel(aligned_active, "xyz")
    view.setStyle({'model': 1}, {
        'stick': {'radius': 0.15},
        'sphere': {'scale': 0.28}
    })

    view.zoomTo()
    html = view._make_html()
    components.html(html, height=height + 20)


def generate_2d_depiction(smiles: str):
    """Disegna la struttura 2D con indici atomici."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetProp("atomLabel", f"{atom.GetSymbol()}{atom.GetIdx()}")
    img = Draw.MolToImage(mol, size=(360, 240))
    return img


# ==============================================================================
# 3. Sidebar: Parametri, Profili di Calcolo e Sessioni Archiviate
# ==============================================================================
st.sidebar.markdown("### Configurazione Calcolo")

# Profili di calcolo pre-configurati
profili = {
    "Profilo Legacy": {
        "n_confs": 50,
        "rmsd": 0.20,
        "ff": "uff",
        "enable_window": False,
        "window": 0.0,
        "method": "b3lyp",
        "basis": "6-31G(d)",
        "sp_only": False,
        "desc": "Parametri identici alla v1.0 per garantire riproducibilità e comparabilità totale (B3LYP/6-31G(d), UFF)."
    },
    "Standard / Bilanciato": {
        "n_confs": 30,
        "rmsd": 0.20,
        "ff": "mmff94",
        "enable_window": True,
        "window": 10.0,
        "method": "b3lyp",
        "basis": "6-31G(d)",
        "sp_only": False,
        "desc": "Miglior rapporto accuratezza/efficienza per composti organici (B3LYP/6-31G(d), MMFF94, finestra 10 kcal/mol)."
    },
    "Alta Precisione (Dispersion-Corrected)": {
        "n_confs": 40,
        "rmsd": 0.15,
        "ff": "mmff94",
        "enable_window": True,
        "window": 12.0,
        "method": "wb97x-d",
        "basis": "def2-SVP",
        "sp_only": False,
        "desc": "Inclusione di dispersione di London per anelli e interazioni non covalenti (ωB97X-D / def2-SVP)."
    },
    "Screening Rapido (Single-Point)": {
        "n_confs": 15,
        "rmsd": 0.25,
        "ff": "mmff94",
        "enable_window": True,
        "window": 8.0,
        "method": "b3lyp",
        "basis": "3-21G",
        "sp_only": True,
        "desc": "Valutazione rapida con Single-Point e base compatta."
    },
    "Personalizzato": {
        "n_confs": 30,
        "rmsd": 0.20,
        "ff": "mmff94",
        "enable_window": True,
        "window": 10.0,
        "method": "b3lyp",
        "basis": "6-31G(d)",
        "sp_only": False,
        "desc": "Configurazione manuale libera di tutti i parametri."
    }
}

scelta_profilo = st.sidebar.selectbox("Profilo di calcolo:", list(profili.keys()), index=0)
pref = profili[scelta_profilo]
st.sidebar.caption(f"*{pref['desc']}*")

st.sidebar.markdown("---")

# Esempi rapidi
esempi = {
    "Personalizzato": "",
    "Butano (Anti vs Gauche)": "CCCC",
    "Cicloesano (Sedia vs Barca)": "C1CCCCC1",
    "Bifenile (Angolo Diedro ~38°)": "c1ccccc1-c2ccccc2",
    "n-Esano (Equilibri Conformazionali)": "CCCCCC",
    "Etanolo": "CCO",
    "Caffeina": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}

scelta_esempio = st.sidebar.selectbox("Esempi molecolari:", list(esempi.keys()))
default_smiles = esempi[scelta_esempio] if scelta_esempio != "Personalizzato" else "CCCC"
smiles_input = st.sidebar.text_input("Codice SMILES target:", value=default_smiles)

# Caricamento sessioni archiviate da disco
output_base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
available_runs = []
if os.path.isdir(output_base):
    for entry in sorted(os.listdir(output_base), reverse=True):
        full_dir = os.path.join(output_base, entry)
        if os.path.isdir(full_dir) and os.path.exists(os.path.join(full_dir, "result.json")):
            available_runs.append(entry)

if available_runs:
    with st.sidebar.expander("Sessioni Archiviate", expanded=False):
        sel_run = st.selectbox("Seleziona esecuzione:", available_runs, key="sel_archived_run")
        if st.button("Carica nei Risultati", icon=":material/history:", key="btn_load_archived"):
            jpath = os.path.join(output_base, sel_run, "result.json")
            try:
                st.session_state.calc_result = CalculationResult.from_json(jpath)
                st.success(f"Caricata sessione: {sel_run}")
                st.rerun()
            except Exception as e:
                st.error(f"Errore di caricamento: {e}")

st.sidebar.markdown("---")

ff_list = ["mmff94", "mmff94s", "uff"]
ff_idx = ff_list.index(pref["ff"]) if pref["ff"] in ff_list else 0

method_list = ["b3lyp", "wb97x-d", "pbe", "hf", "mp2"]
method_idx = method_list.index(pref["method"]) if pref["method"] in method_list else 0

basis_list = ["6-31G(d)", "6-311+G(d,p)", "3-21G", "sto-3g", "def2-SVP"]
basis_idx = basis_list.index(pref["basis"]) if pref["basis"] in basis_list else 0

with st.sidebar.expander("Campionamento Conformazionale (RDKit)", expanded=True):
    n_confs = st.slider("Conformeri generati (ETKDGv3):", min_value=5, max_value=200, value=pref["n_confs"], step=5, key=f"nc_{scelta_profilo}")
    rmsd_thresh = st.number_input("Soglia RMSD di simmetria (Å):", min_value=0.05, max_value=2.0, value=pref["rmsd"], step=0.05, key=f"rmsd_{scelta_profilo}")
    ff_choice = st.selectbox("Campo di forze (pre-rilassamento):", ff_list, index=ff_idx, key=f"ff_{scelta_profilo}")
    enable_window = st.checkbox("Filtro Energy Window FF", value=pref["enable_window"], key=f"ew_{scelta_profilo}")
    energy_window = st.slider("Finestra energetica FF (kcal/mol):", 1.0, 30.0, pref["window"] if pref["window"] > 0 else 10.0, 1.0, key=f"ew_val_{scelta_profilo}") if enable_window else 0.0
    limit_top_n = st.checkbox("Limita conformeri inviati a DFT", value=False)
    top_n_val = st.number_input("Max conformeri per DFT:", min_value=1, max_value=50, value=10) if limit_top_n else None

with st.sidebar.expander("Motore Quantistico (Psi4)", expanded=True):
    method = st.selectbox("Metodo Quantistico:", method_list, index=method_idx, key=f"m_{scelta_profilo}")
    basis = st.selectbox("Basis Set:", basis_list, index=basis_idx, key=f"b_{scelta_profilo}")
    threads = st.number_input("Thread CPU:", min_value=1, max_value=os.cpu_count() or 8, value=min(8, os.cpu_count() or 4))
    memory = st.selectbox("Memoria RAM:", ["4 GB", "8 GB", "16 GB", "32 GB"], index=1)
    sp_only = st.checkbox("Solo Single-Point (senza ottimizzazione geometrica)", value=pref["sp_only"], key=f"sp_{scelta_profilo}")

with st.sidebar.expander("Parametri Elettronici e Termodinamici", expanded=False):
    temp_k = st.number_input("Temperatura Boltzmann (K):", value=298.15, step=5.0)
    manual_charge = st.number_input("Carica formale manuale:", value=0)
    override_charge = st.checkbox("Forza carica manuale", value=False)
    manual_mult = st.number_input("Molteplicità di spin manuale:", min_value=1, value=1)
    override_mult = st.checkbox("Forza molteplicità manuale", value=False)

avvia = st.sidebar.button("AVVIA CALCOLO", icon=":material/play_arrow:", type="primary")


# ==============================================================================
# 4. Header Principale & Topologia 2D
# ==============================================================================
st.markdown('<div class="main-title">Energyinator</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Framework ibrido per campionamento conformazionale e chimica quantistica (RDKit + Psi4)</div>', unsafe_allow_html=True)

col_top1, col_top2 = st.columns([3, 1])

with col_top1:
    st.markdown(f"""
    <div class="info-card">
        <span class="metric-badge">SMILES: {smiles_input.strip() or 'N/A'}</span>
        <span class="metric-badge">Metodo: {method.upper()} / {basis}</span>
        <span class="metric-badge">Campo di Forze: {ff_choice.upper()}</span>
        <span class="metric-badge">Profilo: {scelta_profilo}</span>
    </div>
    """, unsafe_allow_html=True)

with col_top2:
    if smiles_input.strip():
        img_2d = generate_2d_depiction(smiles_input.strip())
        if img_2d:
            st.image(img_2d, caption="Topologia 2D con indici atomici", use_container_width=True)
        else:
            st.warning("Struttura SMILES non valida.")


# ==============================================================================
# 5. Gestione Esecuzione del Calcolo (Process Isolation)
# ==============================================================================
if "calc_result" not in st.session_state:
    st.session_state.calc_result = None

if avvia:
    if not smiles_input.strip():
        st.error("Inserire una stringa SMILES valida prima di avviare il calcolo.")
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        run_output_dir = os.path.join(script_dir, "output", f"result_{timestamp}")
        os.makedirs(run_output_dir, exist_ok=True)

        cmd = [
            sys.executable, "-m", "energyinator.cli",
            smiles_input.strip(),
            "-n", str(n_confs),
            "-r", str(rmsd_thresh),
            "--ff", ff_choice,
            "-w", str(energy_window if enable_window else 0.0),
            "-m", method,
            "-b", basis,
            "-t", str(threads),
            "--memory", memory,
            "--temp", str(temp_k),
            "-o", run_output_dir
        ]
        if limit_top_n and top_n_val:
            cmd.extend(["--top-n", str(int(top_n_val))])
        if sp_only:
            cmd.append("--sp-only")
        if override_charge:
            cmd.extend(["--charge", str(int(manual_charge))])
        if override_mult:
            cmd.extend(["--multiplicity", str(int(manual_mult))])

        prog_bar = st.progress(5, text="Inizializzazione del processo quantistico isolato...")

        log_expander = st.expander("Console di Esecuzione (Output in tempo reale)", expanded=True)
        with log_expander:
            log_box = st.empty()

        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
            env=env
        )

        all_lines = []
        curr_pct = 5
        curr_text = "Inizializzazione..."

        for line in iter(proc.stdout.readline, ''):
            all_lines.append(line)
            lower_line = line.lower()

            if "generazione di" in lower_line or "generati" in lower_line:
                curr_pct = 15
                curr_text = "Generazione conformeri 3D con ETKDGv3..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "ottimizzazione conformeri" in lower_line or "mmff" in lower_line or "uff" in lower_line:
                curr_pct = 25
                curr_text = f"Pre-ottimizzazione con campo di forze ({ff_choice.upper()})..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "conformeri unici selezionati" in lower_line or "finestra energetica" in lower_line:
                curr_pct = 35
                curr_text = "Filtraggio per simmetria RMSD ed energy window..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "progresso dft:" in lower_line or "psi4 quantum" in lower_line:
                match = re.search(r'\[(\d+)/(\d+)\]', line)
                if match:
                    step = int(match.group(1))
                    total = int(match.group(2))
                    pct = 35 + int((step / total) * 55)
                    curr_pct = min(90, max(curr_pct, pct))
                    curr_text = f"Calcolo quantistico Psi4 ({method.upper()}/{basis}): conformero {step} di {total}..."
                    prog_bar.progress(curr_pct, text=curr_text)

            elif "report excel" in lower_line or "esportato con successo" in lower_line:
                curr_pct = 95
                curr_text = "Analisi termodinamica di Boltzmann ed esportazione dati..."
                prog_bar.progress(curr_pct, text=curr_text)

            log_box.code("".join(all_lines[-20:]), language="log")

        proc.stdout.close()
        ret_code = proc.wait()

        if ret_code == 0:
            prog_bar.progress(100, text="Calcolo completato con successo!")
            json_file = os.path.join(run_output_dir, "result.json")
            if os.path.exists(json_file):
                res = CalculationResult.from_json(json_file)
                st.session_state.calc_result = res
                st.success(f"Calcolo terminato con successo in {res.elapsed_time:.1f} secondi.")
                st.rerun()
            else:
                st.error("Il processo è terminato ma 'result.json' non è stato trovato.")
        else:
            prog_bar.empty()
            st.error(f"Errore durante l'esecuzione del calcolo (Codice di uscita: {ret_code}).")
            st.code("".join(all_lines), language="log")


# ==============================================================================
# 6. Visualizzazione Risultati: Tab Struttura, PES Scan, Popolazioni, Tabella, Log
# ==============================================================================
res: CalculationResult = st.session_state.calc_result

if res is not None:
    tab3d, tab_pes, tab_stats, tab_table, tab_log = st.tabs([
        ":material/view_in_ar: Struttura 3D e Geometria",
        ":material/timeline: Scansione Torsionale (PES)",
        ":material/analytics: Analisi Energetica e Popolazioni",
        ":material/table_chart: Tabella Dati ed Esportazione",
        ":material/terminal: Log di Calcolo"
    ])

    # --------------------------------------------------------------------------
    # TAB 1: Visualizzatore 3D & Ispettore Geometrico
    # --------------------------------------------------------------------------
    with tab3d:
        col_ctrl, col_view = st.columns([1, 2])

        with col_ctrl:
            st.markdown("#### Selezione Conformero")
            conf_options = {}
            for c in res.conformers:
                label = f"Conformero {c.conf_id}  |  ΔE = {c.delta_dft_kcalmol:.2f} kcal/mol  |  Pop = {c.boltzmann_pop:.1f}%"
                if c.conf_id == res.conformers[0].conf_id:
                    label += "  [Minimo Globale]"
                conf_options[label] = c

            selected_label = st.selectbox("Conformero attivo:", list(conf_options.keys()))
            active_conf: ConformerResult = conf_options[selected_label]

            # Proprietà elettroniche essenziali
            if active_conf.dipole_debye is not None or active_conf.gap_ev is not None:
                st.markdown("#### Proprietà Elettroniche")
                e1, e2 = st.columns(2)
                if active_conf.dipole_debye is not None:
                    e1.metric("Momento Dipolo (μ)", f"{active_conf.dipole_debye:.2f} D")
                if active_conf.gap_ev is not None:
                    e2.metric("Gap HOMO-LUMO", f"{active_conf.gap_ev:.2f} eV")
                if active_conf.homo_ev is not None and active_conf.lumo_ev is not None:
                    st.caption(f"**HOMO:** {active_conf.homo_ev:.2f} eV  |  **LUMO:** {active_conf.lumo_ev:.2f} eV")

            st.markdown("---")
            st.markdown("#### Parametri di Visualizzazione")
            style_choice = st.selectbox(
                "Stile grafico 3D:",
                ["ball_and_stick", "stick", "sphere", "wireframe"],
                format_func=lambda x: {
                    "ball_and_stick": "Ball & Stick (Sfere e Cilindri)",
                    "stick": "Stick (Bastoncini)",
                    "sphere": "Spacefill / CPK (Sfere vdW)",
                    "wireframe": "Wireframe (Filo di ferro)"
                }[x]
            )
            labels_on = st.checkbox("Mostra etichette atomiche (C0, C1, H2...)", value=True)
            spin_on = st.checkbox("Rotazione automatica", value=False)
            overlay_on = st.checkbox("Sovrapponi al Minimo Globale (Overlay)", value=False)

            st.markdown("---")
            # Ispettore geometrico
            st.markdown("#### Ispettore Geometrico")
            st.caption("Misura distanze, angoli di legame e diedri direttamente sulle coordinate calcolate:")

            symbols, coords = parse_xyz(active_conf.xyz_optimized)
            num_atoms = len(symbols)
            atom_indices = list(range(num_atoms))
            atom_names = [f"{sym}{i}" for i, sym in enumerate(symbols)]

            measure_type = st.radio(
                "Misura geometrica:",
                ["Distanza (2 atomi)", "Angolo di legame (3 atomi)", "Angolo Diedro (4 atomi)"]
            )

            if measure_type == "Distanza (2 atomi)":
                c1, c2 = st.columns(2)
                a1 = c1.selectbox("Atomo 1", atom_indices, format_func=lambda i: atom_names[i], index=0, key="dist_a1")
                a2 = c2.selectbox("Atomo 2", atom_indices, format_func=lambda i: atom_names[i], index=min(1, num_atoms-1), key="dist_a2")
                if a1 != a2:
                    dist = calc_distance(coords[a1], coords[a2])
                    st.metric(label=f"Distanza {atom_names[a1]} - {atom_names[a2]}", value=f"{dist:.3f} Å")
                else:
                    st.warning("Selezionare due atomi distinti.")

            elif measure_type == "Angolo di legame (3 atomi)":
                c1, c2, c3 = st.columns(3)
                a1 = c1.selectbox("Atomo 1", atom_indices, format_func=lambda i: atom_names[i], index=0, key="ang_a1")
                a2 = c2.selectbox("Vertice 2", atom_indices, format_func=lambda i: atom_names[i], index=min(1, num_atoms-1), key="ang_a2")
                a3 = c3.selectbox("Atomo 3", atom_indices, format_func=lambda i: atom_names[i], index=min(2, num_atoms-1), key="ang_a3")
                if len({a1, a2, a3}) == 3:
                    ang = calc_angle(coords[a1], coords[a2], coords[a3])
                    st.metric(label=f"Angolo {atom_names[a1]}-{atom_names[a2]}-{atom_names[a3]}", value=f"{ang:.2f}°")
                else:
                    st.warning("Selezionare tre atomi distinti.")

            else:  # Diedro
                c1, c2, c3, c4 = st.columns(4)
                a1 = c1.selectbox("Atomo 1", atom_indices, format_func=lambda i: atom_names[i], index=0, key="dih_a1")
                a2 = c2.selectbox("Atomo 2", atom_indices, format_func=lambda i: atom_names[i], index=min(1, num_atoms-1), key="dih_a2")
                a3 = c3.selectbox("Atomo 3", atom_indices, format_func=lambda i: atom_names[i], index=min(2, num_atoms-1), key="dih_a3")
                a4 = c4.selectbox("Atomo 4", atom_indices, format_func=lambda i: atom_names[i], index=min(3, num_atoms-1), key="dih_a4")
                if len({a1, a2, a3, a4}) == 4:
                    dihedral = calc_dihedral(coords[a1], coords[a2], coords[a3], coords[a4])
                    st.metric(label=f"Diedro {atom_names[a1]}-{atom_names[a2]}-{atom_names[a3]}-{atom_names[a4]}", value=f"{dihedral:.2f}°")
                    st.caption(f"Valore assoluto torsione: **{abs(dihedral):.2f}°**")
                else:
                    st.warning("Selezionare quattro atomi distinti.")

        with col_view:
            if overlay_on and active_conf.conf_id != res.conformers[0].conf_id:
                st.markdown(f"#### Struttura 3D Sovrapposta — Conformero {active_conf.conf_id} vs Minimo Globale")
                st.caption(f"Blu ciano semitrasparente: **Minimo Globale (Conf {res.conformers[0].conf_id})**  |  Colori CPK: **Conformero {active_conf.conf_id} allineato**")
                render_3d_mol_overlay(res.conformers[0].xyz_optimized, active_conf.xyz_optimized, height=620)
            else:
                st.markdown(f"#### Struttura 3D — Conformero {active_conf.conf_id}")
                st.caption("Interazione: rotazione (tasto sinistro), zoom (rotella), traslazione (tasto destro).")
                render_3d_mol(active_conf.xyz_optimized, style=style_choice, show_labels=labels_on, spin=spin_on, height=620)

    # --------------------------------------------------------------------------
    # TAB 2: Scansione Torsionale (PES Scan 1D) — Scheda Dedicata
    # --------------------------------------------------------------------------
    with tab_pes:
        st.markdown("### Scansione Torsionale 1D (Potential Energy Surface)")
        st.caption("Calcolo continuo della barriera rotazionale da 0° a 360° attorno a un angolo diedro target:")

        symbols_all, coords_all = parse_xyz(res.conformers[0].xyz_optimized)
        num_atoms_all = len(symbols_all)
        atom_indices_all = list(range(num_atoms_all))
        atom_names_all = [f"{sym}{i}" for i, sym in enumerate(symbols_all)]

        # Barra orizzontale di configurazione scansione
        cfg_col1, cfg_col2, cfg_col3, cfg_col4, cfg_col5, cfg_col6 = st.columns([1, 1, 1, 1, 1.2, 1.2])

        # Default intelligente per molecole lineari come butano (0, 1, 2, 3)
        pes_a1 = cfg_col1.selectbox("Atomo 1", atom_indices_all, format_func=lambda i: atom_names_all[i], index=0, key="tab_pes_a1")
        pes_a2 = cfg_col2.selectbox("Atomo 2", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(1, num_atoms_all-1), key="tab_pes_a2")
        pes_a3 = cfg_col3.selectbox("Atomo 3", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(2, num_atoms_all-1), key="tab_pes_a3")
        pes_a4 = cfg_col4.selectbox("Atomo 4", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(3, num_atoms_all-1), key="tab_pes_a4")

        scan_res_choice = cfg_col5.selectbox("Risoluzione:", [15, 30, 10, 5], index=0, format_func=lambda x: f"{x}° ({360//x} punti)")
        scan_engine_choice = cfg_col6.selectbox("Metodo energetico:", ["mmff94", "uff", "psi4"], format_func=lambda x: {"mmff94": "MMFF94 (Istantaneo)", "uff": "UFF (Istantaneo)", "psi4": f"Psi4 DFT ({res.method.upper()})"}[x])

        btn_run_pes = st.button("CALCOLA BARRIERA DI ROTAZIONE", icon=":material/query_stats:", type="primary", key="btn_run_pes_tab")

        pes_cache_key = f"pes_tab_{res.smiles}_{pes_a1}_{pes_a2}_{pes_a3}_{pes_a4}_{scan_engine_choice}_{scan_res_choice}"

        if btn_run_pes:
            if len({pes_a1, pes_a2, pes_a3, pes_a4}) != 4:
                st.error("Selezionare quattro atomi distinti per definire l'angolo diedro.")
            else:
                with st.spinner(f"Calcolo della superficie torsionale con {scan_engine_choice.upper()}..."):
                    try:
                        pes_data = run_dihedral_scan(
                            smiles=res.smiles,
                            xyz_str=res.conformers[0].xyz_optimized,
                            atom_indices=(pes_a1, pes_a2, pes_a3, pes_a4),
                            step_deg=int(scan_res_choice),
                            engine=scan_engine_choice,
                            psi4_method=res.method,
                            psi4_basis=res.basis,
                            charge=res.formal_charge,
                            multiplicity=res.spin_multiplicity
                        )
                        st.session_state[pes_cache_key] = pes_data
                        st.success(f"Curva generata con successo. Barriera torsionale massima: {pes_data['max_barrier_kcalmol']:.2f} kcal/mol")
                    except Exception as err:
                        st.error(f"Errore durante la scansione PES: {err}")

        # Se i dati PES sono presenti in sessione, visualizzali su layout ampio a due colonne
        if pes_cache_key in st.session_state:
            pes_result = st.session_state[pes_cache_key]
            st.markdown("---")

            pes_col_chart, pes_col_3d = st.columns([3, 2])

            with pes_col_chart:
                df_pes = pd.DataFrame({
                    "Angolo": pes_result["angles"],
                    "Delta_E_kcalmol": pes_result["delta_kcalmol"]
                })

                fig_pes = go.Figure()
                fig_pes.add_trace(go.Scatter(
                    x=df_pes["Angolo"],
                    y=df_pes["Delta_E_kcalmol"],
                    mode="lines+markers",
                    name="Energia Relativa",
                    line=dict(color="#3b82f6", width=2.5),
                    marker=dict(size=6, color="#60a5fa")
                ))

                fig_pes.update_layout(
                    title=f"Curva di Potenziale Torsionale — Diedro {atom_names_all[pes_a1]}-{atom_names_all[pes_a2]}-{atom_names_all[pes_a3]}-{atom_names_all[pes_a4]}",
                    xaxis_title="Angolo Diedro (°)",
                    yaxis_title="ΔE Relativo (kcal/mol)",
                    xaxis=dict(tickmode='linear', tick0=0, dtick=30, range=[-5, 365]),
                    hovermode="x unified",
                    height=520,
                    margin=dict(l=40, r=40, t=50, b=40)
                )

                st.plotly_chart(fig_pes, use_container_width=True)

                st.markdown(f"""
                <div class="info-card">
                    <span class="metric-badge">Minimo Energetico: {pes_result['min_angle']}°</span>
                    <span class="metric-badge">Barriera Massima (ΔE‡): {pes_result['max_barrier_kcalmol']:.2f} kcal/mol</span>
                    <span class="metric-badge">Risoluzione: {scan_res_choice}°</span>
                    <span class="metric-badge">Metodo: {scan_engine_choice.upper()}</span>
                </div>
                """, unsafe_allow_html=True)

            with pes_col_3d:
                st.markdown("#### Ispezione Geometria 3D lungo la Curva")
                scrub_angle = st.select_slider(
                    "Seleziona angolo torsionale:",
                    options=pes_result["angles"],
                    value=pes_result["min_angle"],
                    key=f"slider_{pes_cache_key}"
                )
                idx_pt = pes_result["angles"].index(scrub_angle)
                scrub_xyz = pes_result["xyz_list"][idx_pt]
                scrub_e = pes_result["delta_kcalmol"][idx_pt]

                st.caption(f"Angolo: **{scrub_angle}°**  |  ΔE relativo: **{scrub_e:.2f} kcal/mol**")
                render_3d_mol(scrub_xyz, style="ball_and_stick", show_labels=True, height=450)

    # --------------------------------------------------------------------------
    # TAB 3: Analisi Energetica & Popolazioni di Boltzmann
    # --------------------------------------------------------------------------
    with tab_stats:
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Minimo Globale", f"{res.min_energy_hartree:.6f} Eh", f"{res.min_energy_kcalmol:.2f} kcal/mol")
        m2.metric("Stato ad Alta Energia", f"{res.max_energy_hartree:.6f} Eh", f"+{res.max_energy_kcalmol - res.min_energy_kcalmol:.2f} kcal/mol")
        m3.metric("Conformeri Unici DFT", f"{res.unique_dft} / {res.calculated_dft}")
        m4.metric("Tempo di Calcolo", f"{res.elapsed_time:.1f} s")

        st.markdown("---")

        col_g1, col_g2 = st.columns([3, 2])

        df_chart = pd.DataFrame([{
            "Conformer": f"Conf {c.conf_id}",
            "Delta_E_kcalmol": c.delta_dft_kcalmol,
            "Boltzmann_Pop": c.boltzmann_pop,
            "Degeneracy": c.degeneracy
        } for c in res.conformers])

        with col_g1:
            st.markdown("#### Spettro delle Energie Relative (ΔE)")
            fig_bar = px.bar(
                df_chart,
                x="Conformer",
                y="Delta_E_kcalmol",
                color="Boltzmann_Pop",
                color_continuous_scale="Viridis",
                labels={"Delta_E_kcalmol": "ΔE (kcal/mol)", "Boltzmann_Pop": "Popolazione (%)"},
                title="Differenza Energetica rispetto al Ground State"
            )
            fig_bar.update_layout(xaxis_title="Conformero", yaxis_title="ΔE (kcal/mol)")
            st.plotly_chart(fig_bar, use_container_width=True)

        with col_g2:
            st.markdown(f"#### Distribuzione di Boltzmann (T = 298.15 K)")
            fig_pie = px.pie(
                df_chart,
                names="Conformer",
                values="Boltzmann_Pop",
                hole=0.45,
                title="Frazioni Molari all'Equilibrio"
            )
            st.plotly_chart(fig_pie, use_container_width=True)

        # Momento di dipolo se disponibile
        has_dipoles = any(c.dipole_debye is not None for c in res.conformers)
        if has_dipoles:
            st.markdown("---")
            st.markdown("#### Momenti di Dipolo Elettrico (Polarità Molecolare)")
            df_dip = pd.DataFrame([{
                "Conformer": f"Conf {c.conf_id}",
                "Dipolo_Debye": c.dipole_debye if c.dipole_debye is not None else 0.0,
                "Delta_E_kcalmol": c.delta_dft_kcalmol
            } for c in res.conformers])
            fig_dip = px.bar(
                df_dip,
                x="Conformer",
                y="Dipolo_Debye",
                color="Delta_E_kcalmol",
                color_continuous_scale="Plasma",
                labels={"Dipolo_Debye": "Momento di Dipolo μ (Debye)", "Delta_E_kcalmol": "ΔE (kcal/mol)"},
                title="Momento di Dipolo per Conformero (Debye)"
            )
            fig_dip.update_layout(xaxis_title="Conformero", yaxis_title="Momento di Dipolo (Debye)")
            st.plotly_chart(fig_dip, use_container_width=True)

    # --------------------------------------------------------------------------
    # TAB 4: Tabella Dati & Esportazione
    # --------------------------------------------------------------------------
    with tab_table:
        st.markdown("#### Tabella Risultati Dettagliata")

        table_rows = []
        for c in res.conformers:
            table_rows.append({
                "Conformer ID": c.conf_id,
                "Molteplicità (g_i)": c.degeneracy,
                "Energia (Hartree)": f"{c.dft_energy_hartree:.6f}",
                "Energia (kcal/mol)": f"{c.dft_energy_kcalmol:.2f}",
                "ΔE (kcal/mol)": f"{c.delta_dft_kcalmol:.2f}",
                "ΔE (kJ/mol)": f"{c.delta_dft_kjmol:.2f}",
                "Popolazione (%)": f"{c.boltzmann_pop:.1f}%",
                "Dipolo (Debye)": f"{c.dipole_debye:.2f}" if c.dipole_debye is not None else "-",
                "HOMO (eV)": f"{c.homo_ev:.2f}" if c.homo_ev is not None else "-",
                "LUMO (eV)": f"{c.lumo_ev:.2f}" if c.lumo_ev is not None else "-",
                "Gap (eV)": f"{c.gap_ev:.2f}" if c.gap_ev is not None else "-",
                "FF Energy": f"{c.ff_energy:.2f}",
                "ΔE FF": f"{c.delta_ff_kcalmol:.2f}",
                "File XYZ": os.path.basename(c.output_path) if c.output_path else ""
            })

        df_display = pd.DataFrame(table_rows)
        st.dataframe(df_display, use_container_width=True, hide_index=True)

        c_exp1, c_exp2, c_exp3 = st.columns(3)

        if res.excel_path and os.path.exists(res.excel_path):
            with open(res.excel_path, "rb") as f:
                c_exp1.download_button(
                    label="Scarica Tabella Excel (.xlsx)",
                    data=f.read(),
                    file_name=os.path.basename(res.excel_path),
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    icon=":material/download:"
                )

        if res.best_output_path and os.path.exists(res.best_output_path):
            with open(res.best_output_path, "r", encoding="utf-8") as f:
                c_exp2.download_button(
                    label="Scarica Geometria Minima (.xyz)",
                    data=f.read(),
                    file_name=os.path.basename(res.best_output_path),
                    mime="text/plain",
                    icon=":material/download:"
                )

        if c_exp3.button("Apri Cartella Risultati", icon=":material/folder_open:"):
            if os.path.isdir(res.output_dir):
                subprocess.Popen(f'explorer "{res.output_dir}"')

    # --------------------------------------------------------------------------
    # TAB 5: Log di Calcolo
    # --------------------------------------------------------------------------
    with tab_log:
        st.markdown("#### Log di Esecuzione")
        log_file = os.path.join(res.output_dir, "run.log")
        if os.path.exists(log_file):
            with open(log_file, "r", encoding="utf-8") as f:
                st.code(f.read(), language="log")
        else:
            st.info("Nessun file di log trovato.")

        psi4_out = os.path.join(res.output_dir, "psi4_output.dat")
        if os.path.exists(psi4_out):
            with st.expander("Output scientifico grezzo di Psi4 (ultime righe)"):
                with open(psi4_out, "r", encoding="utf-8") as f:
                    st.code(f.read()[-6000:], language="text")
