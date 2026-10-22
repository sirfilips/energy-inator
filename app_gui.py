"""Energyinator - Scientific Molecular Modeling & Quantum Conformer Pipeline.

Web graphical interface powered by Streamlit, RDKit, py3Dmol, Ketcher, and Psi4.
"""

import os
import io
import re
import sys
import math
import zipfile
import subprocess
import urllib.request
import urllib.parse
import json
from datetime import datetime
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import py3Dmol
import streamlit.components.v1 as components

from rdkit import Chem
from rdkit.Chem import Draw, AllChem

from energyinator.config import CalculationConfig
from energyinator.models import CalculationResult, ConformerResult
from energyinator.scan import run_dihedral_scan, run_2d_dihedral_scan
from energyinator.conformers import compute_pairwise_rmsd_matrix
from energyinator.quantum import generate_orbital_cube

# Importazione dell'editor molecolare 2D Ketcher
try:
    from streamlit_ketcher import st_ketcher
    HAS_KETCHER = True
except ImportError:
    HAS_KETCHER = False


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
# 2. Funzioni Ausiliarie per Rendering Adattivo, Geometria 3D, PubChem & File IO
# ==============================================================================
def st_image_stretched(image, caption=None):
    """Renderizza immagini a tutta larghezza compatibile con versioni recenti e legacy di Streamlit."""
    try:
        st.image(image, caption=caption, width="stretch")
    except (TypeError, ValueError):
        st.image(image, caption=caption, use_container_width=True)


def st_plotly_stretched(fig):
    """Renderizza grafici Plotly a tutta larghezza compatibile con versioni recenti e legacy di Streamlit."""
    try:
        st.plotly_chart(fig, width="stretch")
    except (TypeError, ValueError):
        st.plotly_chart(fig, use_container_width=True)


def st_dataframe_stretched(df, hide_index=True):
    """Renderizza tabelle DataFrame a tutta larghezza compatibile con versioni recenti e legacy di Streamlit."""
    try:
        st.dataframe(df, width="stretch", hide_index=hide_index)
    except (TypeError, ValueError):
        st.dataframe(df, use_container_width=True, hide_index=hide_index)


def fetch_pubchem_smiles(query: str):
    """Interroga PubChem, NCBI E-utilities e Cactus per convertire nomi comuni (in italiano o inglese) in SMILES."""
    query_clean = query.strip()
    if not query_clean:
        return None

    candidates = [query_clean]
    low = query_clean.lower()

    # Dizionario delle sostanze organiche comuni in italiano
    it_dict = {
        'aspirina': 'aspirin',
        'paracetamolo': 'paracetamol',
        'tachipirina': 'paracetamol',
        'glucosio': 'glucose',
        'fruttosio': 'fructose',
        'saccarosio': 'sucrose',
        'galattosio': 'galactose',
        'lattosio': 'lactose',
        'maltosio': 'maltose',
        'acqua': 'water',
        'ammoniaca': 'ammonia',
        'acido acetilsalicilico': 'acetylsalicylic acid',
        'acido ascorbico': 'ascorbic acid',
        'acido citrico': 'citric acid',
        'acido ossalico': 'oxalic acid',
        'acido tartarico': 'tartaric acid',
        'acido piruvico': 'pyruvic acid',
        'acido succinico': 'succinic acid',
        'acido malico': 'malic acid',
        'acido fumarico': 'fumaric acid',
        'acido benzoico': 'benzoic acid',
        'acido acetico': 'acetic acid',
        'acido formico': 'formic acid',
        'acido propionico': 'propionic acid',
        'acido butirrico': 'butyric acid',
        'acido lattico': 'lactic acid',
        'acido valerico': 'valeric acid',
    }
    if low in it_dict:
        candidates.append(it_dict[low])

    # Regola grammaticale 'acido ...ico' -> '...ic acid'
    m = re.match(r'^acido\s+(.+)$', low)
    if m:
        core = m.group(1)
        if core.endswith('ico'):
            candidates.append(core[:-3] + 'ic acid')
        elif core.endswith('oso'):
            candidates.append(core[:-3] + 'ous acid')

    # Suffissi chimici italiani -> inglesi
    if low.endswith('olo'):
        candidates.append(low[:-3] + 'ol')
    elif low.endswith('ano'):
        candidates.append(low[:-3] + 'ane')
    elif low.endswith('ina'):
        candidates.append(low[:-3] + 'ine')
    elif low.endswith('osio'):
        candidates.append(low[:-4] + 'ose')

    candidates = list(dict.fromkeys(candidates))
    headers = {"User-Agent": "Energyinator-App/1.0 (Quantum Molecular Modeling Studio)"}

    # Livello 1: PubChem PUG REST diretto su compound/name
    for cand in candidates:
        try:
            q = urllib.parse.quote(cand)
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/CanonicalSMILES,ConnectivitySMILES,IsomericSMILES/JSON"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=4) as response:
                data = json.loads(response.read().decode("utf-8"))
                props = data["PropertyTable"]["Properties"][0]
                smi = props.get("IsomericSMILES") or props.get("SMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
                if smi:
                    return smi
        except Exception:
            continue

    # Livello 2: NCBI E-Utilities (ricerca semantica su tutti i sinonimi e marchi)
    for cand in candidates[:2]:
        try:
            q = urllib.parse.quote(cand)
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term={q}&retmode=json"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=4) as response:
                data = json.loads(response.read().decode("utf-8"))
                id_list = data.get("esearchresult", {}).get("idlist", [])
                if id_list:
                    cid = id_list[0]
                    u2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,ConnectivitySMILES,IsomericSMILES/JSON"
                    r2 = urllib.request.Request(u2, headers=headers)
                    with urllib.request.urlopen(r2, timeout=4) as res2:
                        d2 = json.loads(res2.read().decode("utf-8"))
                        props = d2["PropertyTable"]["Properties"][0]
                        smi = props.get("IsomericSMILES") or props.get("SMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
                        if smi:
                            return smi
        except Exception:
            continue

    # Livello 3: NIH Cactus Chemical Identifier Resolver
    for cand in candidates[:2]:
        try:
            q = urllib.parse.quote(cand)
            url = f"https://cactus.nci.nih.gov/chemical/structure/{q}/smiles"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=4) as response:
                smi = response.read().decode("utf-8").strip()
                if smi and "<html" not in smi.lower() and "error" not in smi.lower():
                    return smi
        except Exception:
            continue

    return None


def mol_file_to_smiles(file_bytes: bytes, filename: str):
    """Converte il contenuto di un file molecolare (.mol, .sdf, .smi, .txt) in SMILES."""
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
    """Genera un archivio ZIP in memoria con tutte le geometrie XYZ, tabella Excel, JSON e riassunto."""
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


def render_3d_mol(
    xyz_str,
    style='ball_and_stick',
    show_labels=True,
    spin=False,
    surface_type='Nessuna',
    surface_opacity=0.65,
    orbital_cube_data=None,
    orbital_isoval=0.03,
    orbital_opacity=0.75,
    height=580
):
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

    if surface_type == 'Van der Waals (VDW)':
        view.addSurface(py3Dmol.VDW, {'opacity': surface_opacity, 'color': 'white'})
    elif surface_type == 'Accessibile al Solvente (SAS)':
        view.addSurface(py3Dmol.SAS, {'opacity': surface_opacity, 'color': '#38bdf8'})

    # Renderizzazione Isosuperfici Orbitale Molecolare (lobo positivo blu ciano, lobo negativo rosso corallo)
    if orbital_cube_data:
        view.addVolumetricData(
            orbital_cube_data,
            "cube",
            {'isoval': float(orbital_isoval), 'color': '#2563eb', 'opacity': float(orbital_opacity)}
        )
        view.addVolumetricData(
            orbital_cube_data,
            "cube",
            {'isoval': -float(orbital_isoval), 'color': '#dc2626', 'opacity': float(orbital_opacity)}
        )

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


def smiles_to_2d_molblock(smiles: str) -> str:
    """Converte una stringa SMILES in un MolBlock V2000 con coordinate 2D pulite generate da RDKit per Ketcher."""
    if not smiles or not smiles.strip():
        return ""
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return ""
        AllChem.Compute2DCoords(mol)
        return Chem.MolToMolBlock(mol)
    except Exception:
        return ""


def mol_output_to_smiles(raw_output: str) -> str:
    """Converte l'output di Ketcher (MolBlock 2D o stringa SMILES) in SMILES canonico."""
    if not raw_output or not raw_output.strip():
        return ""
    if "\n" not in raw_output:
        return raw_output.strip()
    try:
        mol = Chem.MolFromMolBlock(raw_output)
        if mol is not None:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass
    return raw_output.strip()


# ==============================================================================
# 3. Sidebar: Parametri, Profili di Calcolo e Sessioni Archiviate
# ==============================================================================
st.sidebar.markdown("### Configurazione Calcolo")

# Profili di calcolo pre-configurati
profili = {
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

if "pending_example_selection" in st.session_state:
    st.session_state["sel_mol_example"] = st.session_state.pop("pending_example_selection")
elif "sel_mol_example" not in st.session_state:
    st.session_state["sel_mol_example"] = "Personalizzato"

scelta_esempio = st.sidebar.selectbox("Esempi molecolari:", list(esempi.keys()), key="sel_mol_example")
default_smiles = esempi[scelta_esempio] if scelta_esempio != "Personalizzato" else "CCCC"

# Inizializzazione stato per sincronizzazione SMILES
if "active_smiles" not in st.session_state:
    st.session_state["active_smiles"] = default_smiles
if "last_ketcher_applied" not in st.session_state:
    st.session_state["last_ketcher_applied"] = default_smiles

# Rilevamento cambio esempio
if scelta_esempio != "Personalizzato":
    if st.session_state.get("prev_selected_example") != scelta_esempio:
        st.session_state["prev_selected_example"] = scelta_esempio
        st.session_state["active_smiles"] = default_smiles
        st.session_state["last_ketcher_applied"] = default_smiles
        st.session_state["ketcher_version"] = st.session_state.get("ketcher_version", 0) + 1

# Ricerca rapida PubChem
with st.sidebar.expander("Trova Molecola per Nome (PubChem)", expanded=False):
    nome_query = st.text_input("Nome comune:", placeholder="es. aspirina, paracetamolo, caffeina", key="input_pubchem_name")
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
            st.error(f"Nessun risultato trovato per '{nome_query}'. Prova con un sinonimo o in inglese.")

    if "pubchem_search_result" in st.session_state:
        st.success(st.session_state["pubchem_search_result"])

# Importa da file molecolare
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

# Assegnazione del session state di smiles_input_box PRIMA della creazione del widget
if "smiles_input_box" not in st.session_state or st.session_state["smiles_input_box"] != st.session_state["active_smiles"]:
    st.session_state["smiles_input_box"] = st.session_state["active_smiles"]

def on_sidebar_smiles_change():
    st.session_state["active_smiles"] = st.session_state["smiles_input_box"]
    st.session_state["last_ketcher_applied"] = st.session_state["smiles_input_box"]

smiles_input = st.sidebar.text_input(
    "Codice SMILES target:",
    key="smiles_input_box",
    on_change=on_sidebar_smiles_change
)

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

solvent_map = {
    "Nessuno / Gas (Vuoto)": None,
    "Acqua (Water)": "water",
    "Cloroformio (Chloroform)": "chloroform",
    "DMSO (Dimethylsulfoxide)": "dmso",
    "Metanolo (Methanol)": "methanol",
    "Toluene": "toluene",
    "Acetonitrile": "acetonitrile",
    "DCM (Dichloromethane)": "dcm"
}
solv_keys = list(solvent_map.keys())
solv_idx = solv_keys.index(pref.get("solvent", "Nessuno / Gas (Vuoto)")) if pref.get("solvent") in solv_keys else 0

with st.sidebar.expander("Campionamento Conformazionale (RDKit)", expanded=True):
    n_confs = st.slider("Conformeri generati (ETKDGv3):", min_value=5, max_value=200, value=pref["n_confs"], step=5, key=f"nc_{scelta_profilo}")
    rmsd_thresh = st.number_input("Soglia RMSD classica (Å):", min_value=0.05, max_value=2.0, value=pref["rmsd"], step=0.05, key=f"rmsd_{scelta_profilo}")
    post_rmsd = st.number_input("Deduplicazione post-DFT (Å):", min_value=0.0, max_value=1.0, value=pref.get("post_rmsd", 0.15), step=0.05, help="Rimuove geometrie collassate nello stesso minimo quantistico.", key=f"p_rmsd_{scelta_profilo}")
    ff_choice = st.selectbox("Campo di forze (pre-rilassamento):", ff_list, index=ff_idx, key=f"ff_{scelta_profilo}")
    enable_window = st.checkbox("Filtro Energy Window FF", value=pref["enable_window"], key=f"ew_{scelta_profilo}")
    energy_window = st.slider("Finestra energetica FF (kcal/mol):", 1.0, 30.0, pref["window"] if pref["window"] > 0 else 10.0, 1.0, key=f"ew_val_{scelta_profilo}") if enable_window else 0.0
    limit_top_n = st.checkbox("Limita conformeri inviati a DFT", value=False)
    top_n_val = st.number_input("Max conformeri per DFT:", min_value=1, max_value=50, value=10) if limit_top_n else None

with st.sidebar.expander("Motore Quantistico (Psi4)", expanded=True):
    method = st.selectbox("Metodo Quantistico:", method_list, index=method_idx, key=f"m_{scelta_profilo}")
    basis = st.selectbox("Basis Set:", basis_list, index=basis_idx, key=f"b_{scelta_profilo}")
    solv_choice = st.selectbox("Solvatazione Implicita (PCM):", solv_keys, index=solv_idx, key=f"s_{scelta_profilo}")
    solvent_val = solvent_map[solv_choice]
    threads = st.number_input("Thread CPU:", min_value=1, max_value=os.cpu_count() or 8, value=min(8, os.cpu_count() or 4))
    memory = st.selectbox("Memoria RAM:", ["4 GB", "8 GB", "16 GB", "32 GB"], index=1)
    sp_only = st.checkbox("Solo Single-Point (senza ottimizzazione geometrica)", value=pref["sp_only"], key=f"sp_{scelta_profilo}")

with st.sidebar.expander("Parametri Elettronici e Termodinamici", expanded=False):
    temp_k = st.number_input("Temperatura Boltzmann (K):", value=298.15, step=5.0)
    calc_freq = st.checkbox("Calcola Termochimica (Gibbs Delta G, ZPVE)", value=pref.get("calc_freq", False), help="Esegue calcolo di frequenze vibrazionali sul ground state.")
    manual_charge = st.number_input("Carica formale manuale:", value=0)
    override_charge = st.checkbox("Forza carica manuale", value=False)
    manual_mult = st.number_input("Molteplicità di spin manuale:", min_value=1, value=1)
    override_mult = st.checkbox("Forza molteplicità manuale", value=False)

avvia = st.sidebar.button("AVVIA CALCOLO", icon=":material/play_arrow:", type="primary")


# ==============================================================================
# 4. Header Principale, Topologia 2D ed Editor Molecolare Integrato (Ketcher)
# ==============================================================================
st.markdown('<div class="main-title">Energyinator</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Framework ibrido per campionamento conformazionale e chimica quantistica (RDKit + Psi4)</div>', unsafe_allow_html=True)

col_top1, col_top2 = st.columns([3, 1])

with col_top1:
    solv_display = solv_choice if solvent_val else "Gas / Vuoto"
    st.markdown(f"""
    <div class="info-card">
        <span class="metric-badge">SMILES: {smiles_input.strip() or 'N/A'}</span>
        <span class="metric-badge">Metodo: {method.upper()} / {basis}</span>
        <span class="metric-badge">Solvente: {solv_display}</span>
        <span class="metric-badge">Campo di Forze: {ff_choice.upper()}</span>
        <span class="metric-badge">Profilo: {scelta_profilo}</span>
    </div>
    """, unsafe_allow_html=True)

with col_top2:
    if smiles_input.strip():
        img_2d = generate_2d_depiction(smiles_input.strip())
        if img_2d:
            st_image_stretched(img_2d, caption="Topologia 2D con indici atomici")
        else:
            st.warning("Struttura SMILES non valida.")

# Editor Molecolare Visuale 2D (Ketcher) senza loop di freeze
with st.expander("Editor Molecolare 2D (Ketcher)", expanded=False):
    if HAS_KETCHER:
        st.caption("Disegna o modifica la struttura molecolare, poi premi **'Applica Struttura nel Calcolo'** per sincronizzarla:")
        
        current_target = st.session_state.get("active_smiles", "CCCC")
        ketcher_key = f"ketcher_canvas_{st.session_state.get('ketcher_version', 0)}"

        # Genera MolBlock 2D pulito con RDKit per Ketcher (risolve topologie corrotte di anelli/zuccheri)
        molblock_2d = smiles_to_2d_molblock(current_target)
        if molblock_2d:
            ketcher_val = molblock_2d
            mol_format = "MOLFILE"
        else:
            ketcher_val = current_target
            mol_format = "SMILES"

        ketcher_output = st_ketcher(value=ketcher_val, height=440, molecule_format=mol_format, key=ketcher_key)

        col_k1, col_k2 = st.columns([1.8, 3])
        btn_apply_ketcher = col_k1.button("Applica Struttura nel Calcolo", icon=":material/check:", type="primary", key="btn_apply_ketcher")

        if btn_apply_ketcher:
            new_smi = mol_output_to_smiles(ketcher_output) if ketcher_output else current_target
            if new_smi:
                st.session_state["last_ketcher_applied"] = new_smi
                st.session_state["active_smiles"] = new_smi
                st.session_state["pending_example_selection"] = "Personalizzato"
                st.session_state["prev_selected_example"] = "Personalizzato"
                st.session_state["ketcher_version"] = st.session_state.get("ketcher_version", 0) + 1
                st.rerun()
    else:
        st.info("Modulo Ketcher non trovato. Installalo con: pip install streamlit-ketcher")


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
            "--post-dft-rmsd", str(post_rmsd),
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
        if solvent_val:
            cmd.extend(["--solvent", solvent_val])
        if calc_freq:
            cmd.append("--frequencies")
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
                    pct = 35 + int((step / total) * 45)
                    curr_pct = min(85, max(curr_pct, pct))
                    curr_text = f"Calcolo quantistico Psi4 ({method.upper()}/{basis}): conformero {step} di {total}..."
                    prog_bar.progress(curr_pct, text=curr_text)

            elif "deduplicazione post-dft" in lower_line or "conformeri collassati" in lower_line:
                curr_pct = 88
                curr_text = "Deduplicazione geometrica post-DFT dei minimi collassati..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "calcolo frequenze" in lower_line or "frequenze vibrazionali" in lower_line:
                curr_pct = 92
                curr_text = "Calcolo frequenze vibrazionali e correzioni termiche (Gibbs)..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "report excel" in lower_line or "esportato con successo" in lower_line:
                curr_pct = 96
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
                label = f"Conformero {c.conf_id}  |  Delta E = {c.delta_dft_kcalmol:.2f} kcal/mol  |  Pop = {c.boltzmann_pop:.1f}%"
                if c.conf_id == res.conformers[0].conf_id:
                    label += "  [Minimo Globale]"
                conf_options[label] = c

            selected_label = st.selectbox("Conformero attivo:", list(conf_options.keys()))
            active_conf: ConformerResult = conf_options[selected_label]

            # Proprietà elettroniche e termochimiche
            st.markdown("#### Proprietà Fisiche ed Elettroniche")
            e1, e2 = st.columns(2)
            if active_conf.dipole_debye is not None:
                e1.metric("Momento Dipolo (mu)", f"{active_conf.dipole_debye:.2f} D")
            if active_conf.gap_ev is not None:
                e2.metric("Gap HOMO-LUMO", f"{active_conf.gap_ev:.2f} eV")
            if active_conf.homo_ev is not None and active_conf.lumo_ev is not None:
                st.caption(f"**HOMO:** {active_conf.homo_ev:.2f} eV  |  **LUMO:** {active_conf.lumo_ev:.2f} eV")

            if active_conf.gibbs_hartree is not None:
                st.markdown("---")
                st.markdown("#### Termochimica (Gibbs)")
                g1, g2 = st.columns(2)
                g1.metric("Energia Libera G", f"{active_conf.gibbs_hartree:.6f} Eh")
                if active_conf.delta_gibbs_kcalmol is not None:
                    g2.metric("Delta G Relativo", f"{active_conf.delta_gibbs_kcalmol:.2f} kcal/mol")
                if active_conf.zpve_hartree is not None:
                    st.caption(f"**ZPVE:** {active_conf.zpve_hartree:.6f} Eh  |  **Entropia:** {active_conf.entropy_cal_mol_k:.1f} cal/(mol*K)")

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
            surface_choice = st.selectbox(
                "Superficie molecolare:",
                ["Nessuna", "Van der Waals (VDW)", "Accessibile al Solvente (SAS)"]
            )
            labels_on = st.checkbox("Mostra etichette atomiche (C0, C1, H2...)", value=True)
            spin_on = st.checkbox("Rotazione automatica", value=False)
            overlay_on = st.checkbox("Sovrapponi al Minimo Globale (Overlay)", value=False)

            st.markdown("---")
            st.markdown("#### Orbitale Molecolare 3D (Isosuperficie)")
            st.caption("Visualizza i lobi della funzione d'onda quantistica &psi;(&mathbf{r}) calcolati con Psi4:")

            mo_choice = st.selectbox(
                "Orbitale molecolare:",
                ["Nessuno", "HOMO", "LUMO", "HOMO-1", "LUMO+1", "HOMO-2", "LUMO+2"],
                key="mo_orbital_choice"
            )

            cube_data_to_render = None
            mo_isoval = 0.03
            mo_opacity = 0.75

            if mo_choice != "Nessuno":
                cube_cache_key = f"cube_{res.smiles}_{active_conf.conf_id}_{mo_choice}_{res.method}_{res.basis}"
                
                col_btn_calc, _ = st.columns([1.5, 1])
                calc_mo_btn = col_btn_calc.button(
                    f"Calcola Isosuperficie {mo_choice}",
                    icon=":material/science:",
                    type="primary",
                    key=f"btn_calc_mo_{mo_choice}_{active_conf.conf_id}"
                )

                if calc_mo_btn and cube_cache_key not in st.session_state:
                    with st.spinner(f"Calcolo della funzione d'onda e griglia 3D per {mo_choice} (Psi4)..."):
                        cube_info = generate_orbital_cube(
                            xyz_str=active_conf.xyz_optimized,
                            orbital_type=mo_choice,
                            method=res.method,
                            basis=res.basis,
                            charge=getattr(res, "formal_charge", 0),
                            multiplicity=getattr(res, "spin_multiplicity", 1),
                            solvent=getattr(res, "solvent", None)
                        )
                        if cube_info and cube_info.get("cube_data"):
                            st.session_state[cube_cache_key] = cube_info
                            st.success(f"{mo_choice} calcolato con successo.")
                        else:
                            st.error(f"Errore durante il calcolo dell'orbitale {mo_choice} con Psi4.")

                if cube_cache_key in st.session_state:
                    cached_cube = st.session_state[cube_cache_key]
                    cube_data_to_render = cached_cube["cube_data"]

                    st.markdown(
                        f"""
                        <div style="background-color: #1e293b; border-radius: 8px; padding: 10px 14px; margin: 8px 0; border: 1px solid #334155;">
                            <div style="font-size: 0.9rem; font-weight: 600; color: #f8fafc; margin-bottom: 4px;">
                                {cached_cube['orbital_type']} (Orbitale Molecolare #{cached_cube['orbital_idx']})
                            </div>
                            <div style="font-size: 0.85rem; color: #cbd5e1;">
                                Energia: <strong style="color: #38bdf8;">{cached_cube['orbital_energy_ev']:.3f} eV</strong> &bull; Metodo: {cached_cube['method'].upper()}/{cached_cube['basis']}
                            </div>
                            <div style="font-size: 0.8rem; margin-top: 6px;">
                                <span style="display: inline-block; width: 12px; height: 12px; background: #2563eb; border-radius: 2px; vertical-align: middle; margin-right: 4px;"></span> Lobo positivo (&psi; &gt; 0)
                                <span style="display: inline-block; width: 12px; height: 12px; background: #dc2626; border-radius: 2px; vertical-align: middle; margin-left: 14px; margin-right: 4px;"></span> Lobo negativo (&psi; &lt; 0)
                            </div>
                        </div>
                        """,
                        unsafe_allow_html=True
                    )

                    c_iso1, c_iso2 = st.columns(2)
                    mo_isoval = c_iso1.slider("Soglia densità (Isoval):", min_value=0.01, max_value=0.08, value=0.03, step=0.005, key="mo_slider_isoval")
                    mo_opacity = c_iso2.slider("Opacità superficie:", min_value=0.10, max_value=1.0, value=0.75, step=0.05, key="mo_slider_opacity")

            st.markdown("---")
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
                title_suffix = f" con Isosuperficie {mo_choice}" if (cube_data_to_render and mo_choice != "Nessuno") else ""
                st.markdown(f"#### Struttura 3D — Conformero {active_conf.conf_id}{title_suffix}")
                st.caption("Interazione: rotazione (tasto sinistro), zoom (rotella), traslazione (tasto destro).")
                render_3d_mol(
                    active_conf.xyz_optimized,
                    style=style_choice,
                    show_labels=labels_on,
                    spin=spin_on,
                    surface_type=surface_choice,
                    orbital_cube_data=cube_data_to_render,
                    orbital_isoval=mo_isoval,
                    orbital_opacity=mo_opacity,
                    height=620
                )

    # --------------------------------------------------------------------------
    # TAB 2: Scansione Torsionale (PES 1D & 2D)
    # --------------------------------------------------------------------------
    with tab_pes:
        pes_mode = st.radio("Tipologia di Scansione Torsionale:", ["Scansione 1D (Diedro Singolo)", "Mappa 2D (Coppia di Diedri / Ramachandran)"], horizontal=True)

        symbols_all, coords_all = parse_xyz(res.conformers[0].xyz_optimized)
        num_atoms_all = len(symbols_all)
        atom_indices_all = list(range(num_atoms_all))
        atom_names_all = [f"{sym}{i}" for i, sym in enumerate(symbols_all)]

        if pes_mode == "Scansione 1D (Diedro Singolo)":
            st.markdown("### Scansione Torsionale 1D (Potential Energy Surface)")
            st.caption("Calcolo della barriera rotazionale da 0° a 360° con supporto per rilassamento geometrico dei gradi di libertà adiacenti:")

            cfg_col1, cfg_col2, cfg_col3, cfg_col4, cfg_col5, cfg_col6, cfg_col7 = st.columns([1, 1, 1, 1, 1, 1.2, 1.2])

            pes_a1 = cfg_col1.selectbox("Atomo 1", atom_indices_all, format_func=lambda i: atom_names_all[i], index=0, key="tab_pes_a1")
            pes_a2 = cfg_col2.selectbox("Atomo 2", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(1, num_atoms_all-1), key="tab_pes_a2")
            pes_a3 = cfg_col3.selectbox("Atomo 3", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(2, num_atoms_all-1), key="tab_pes_a3")
            pes_a4 = cfg_col4.selectbox("Atomo 4", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(3, num_atoms_all-1), key="tab_pes_a4")

            scan_res_choice = cfg_col5.selectbox("Risoluzione:", [15, 30, 10, 5], index=0, format_func=lambda x: f"{x}° ({360//x} punti)")
            scan_engine_choice = cfg_col6.selectbox("Metodo energetico:", ["mmff94", "uff", "psi4"], format_func=lambda x: {"mmff94": "MMFF94", "uff": "UFF", "psi4": f"Psi4 ({res.method.upper()})"}[x])
            relaxed_choice = cfg_col7.checkbox("Scansione Rilassata", value=True, help="Rilassa gli altri angoli e legami per eliminare picchi sterici artificiali.")

            btn_run_pes = st.button("Calcola Barriera di Rotazione 1D", icon=":material/query_stats:", type="primary", key="btn_run_pes_tab")

            pes_cache_key = f"pes_tab_{res.smiles}_{pes_a1}_{pes_a2}_{pes_a3}_{pes_a4}_{scan_engine_choice}_{scan_res_choice}_{relaxed_choice}"

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
                                relaxed=relaxed_choice,
                                psi4_method=res.method,
                                psi4_basis=res.basis,
                                charge=res.formal_charge,
                                multiplicity=res.spin_multiplicity,
                                solvent=res.solvent
                            )
                            st.session_state[pes_cache_key] = pes_data
                            st.success(f"Curva generata con successo. Barriera torsionale: {pes_data['max_barrier_kcalmol']:.2f} kcal/mol")
                        except Exception as err:
                            st.error(f"Errore durante la scansione PES: {err}")

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
                        name="Profilo Energetico",
                        line=dict(color="#38bdf8", width=2.5),
                        marker=dict(size=6, color="#0ea5e9")
                    ))

                    # Evidenziazione del minimo globale con stella smeraldo
                    min_idx = pes_result.get("min_idx", 0)
                    min_angle = pes_result["angles"][min_idx]
                    min_energy = pes_result["delta_kcalmol"][min_idx]
                    fig_pes.add_trace(go.Scatter(
                        x=[min_angle],
                        y=[min_energy],
                        mode="markers",
                        name="Minimo Assoluto",
                        marker=dict(symbol="star", size=14, color="#10b981", line=dict(color="#ffffff", width=1.5)),
                        hoverinfo="text",
                        hovertext=f"Minimo Globale: {min_angle}° (0.00 kcal/mol)"
                    ))

                    fig_pes.update_layout(
                        title=f"Curva di Potenziale Torsionale — Diedro {atom_names_all[pes_a1]}-{atom_names_all[pes_a2]}-{atom_names_all[pes_a3]}-{atom_names_all[pes_a4]}",
                        xaxis_title="Angolo Diedro (°)",
                        yaxis_title="Delta E Relativo (kcal/mol)",
                        xaxis=dict(tickmode='linear', tick0=0, dtick=30, range=[-5, 365]),
                        hovermode="x unified",
                        height=520,
                        margin=dict(l=40, r=40, t=50, b=40)
                    )
                    st_plotly_stretched(fig_pes)

                    st.markdown(f"""
                    <div class="info-card">
                        <span class="metric-badge">Minimo Energetico: {pes_result['min_angle']}°</span>
                        <span class="metric-badge">Barriera Massima: {pes_result['max_barrier_kcalmol']:.2f} kcal/mol</span>
                        <span class="metric-badge">Modalità: {'Rilassata' if pes_result.get('relaxed') else 'Rigida'}</span>
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

                    st.caption(f"Angolo: **{scrub_angle}°**  |  Delta E relativo: **{scrub_e:.2f} kcal/mol**")
                    render_3d_mol(scrub_xyz, style="ball_and_stick", show_labels=True, height=450)

        else:
            # Scansione 2D (Coppia di Diedri)
            st.markdown("### Mappa Conformazionale 2D (PES Scan Bidimensionale)")
            st.caption("Esplorazione della superficie torsionale 2D per molecole con due rotori accoppiati:")

            col_d1, col_d2 = st.columns(2)
            with col_d1:
                st.markdown("##### Diedro 1 (Asse X)")
                c1, c2, c3, c4 = st.columns(4)
                d1_a1 = c1.selectbox("A1", atom_indices_all, format_func=lambda i: atom_names_all[i], index=0, key="2d_d1_a1")
                d1_a2 = c2.selectbox("A2", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(1, num_atoms_all-1), key="2d_d1_a2")
                d1_a3 = c3.selectbox("A3", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(2, num_atoms_all-1), key="2d_d1_a3")
                d1_a4 = c4.selectbox("A4", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(3, num_atoms_all-1), key="2d_d1_a4")

            with col_d2:
                st.markdown("##### Diedro 2 (Asse Y)")
                c1, c2, c3, c4 = st.columns(4)
                d2_a1 = c1.selectbox("B1", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(1, num_atoms_all-1), key="2d_d2_a1")
                d2_a2 = c2.selectbox("B2", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(2, num_atoms_all-1), key="2d_d2_a2")
                d2_a3 = c3.selectbox("B3", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(3, num_atoms_all-1), key="2d_d2_a3")
                d2_a4 = c4.selectbox("B4", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(4, num_atoms_all-1), key="2d_d2_a4")

            c_step, c_eng, c_rel = st.columns([1, 1, 1])
            res_2d_choice = c_step.selectbox("Risoluzione angolare:", [30, 45, 20], index=0, format_func=lambda x: f"{x}° ({(360//x)+1}x{(360//x)+1} punti)")
            eng_2d_choice = c_eng.selectbox("Campo di forze:", ["mmff94", "uff"], index=0)
            rel_2d_choice = c_rel.checkbox("Rilassamento 2D", value=True)

            btn_run_2d = st.button("Calcola Mappa Conformazionale 2D", icon=":material/grid_view:", type="primary", key="btn_run_2d")

            cache_2d_key = f"pes_2d_{res.smiles}_{d1_a1}_{d1_a2}_{d1_a3}_{d1_a4}_{d2_a1}_{d2_a2}_{d2_a3}_{d2_a4}_{eng_2d_choice}_{res_2d_choice}_{rel_2d_choice}"

            if btn_run_2d:
                if len({d1_a1, d1_a2, d1_a3, d1_a4}) != 4 or len({d2_a1, d2_a2, d2_a3, d2_a4}) != 4:
                    st.error("Selezionare quattro atomi distinti per ciascun diedro.")
                else:
                    with st.spinner("Calcolo della griglia bidimensionale in corso..."):
                        try:
                            data_2d = run_2d_dihedral_scan(
                                smiles=res.smiles,
                                xyz_str=res.conformers[0].xyz_optimized,
                                dihedral_1=(d1_a1, d1_a2, d1_a3, d1_a4),
                                dihedral_2=(d2_a1, d2_a2, d2_a3, d2_a4),
                                step_deg=int(res_2d_choice),
                                engine=eng_2d_choice,
                                relaxed=rel_2d_choice
                            )
                            st.session_state[cache_2d_key] = data_2d
                            st.success(f"Mappa 2D generata. Barriera massima: {data_2d['max_barrier_kcalmol']:.2f} kcal/mol")
                        except Exception as e:
                            st.error(f"Errore calcolo scansione 2D: {e}")

            if cache_2d_key in st.session_state:
                res_2d = st.session_state[cache_2d_key]
                st.markdown("---")

                col_map1, col_map2 = st.columns([3, 2])
                with col_map1:
                    fig_contour = go.Figure(data=go.Contour(
                        z=res_2d["delta_z"],
                        x=res_2d["angles_2"],
                        y=res_2d["angles_1"],
                        colorscale="Turbo",
                        colorbar=dict(title="Delta E (kcal/mol)")
                    ))
                    fig_contour.update_layout(
                        title=f"Mappa di Potenziale 2D — Diedro 1 vs Diedro 2",
                        xaxis_title=f"Diedro 2 (°) [{atom_names_all[d2_a1]}-{atom_names_all[d2_a2]}-{atom_names_all[d2_a3]}-{atom_names_all[d2_a4]}]",
                        yaxis_title=f"Diedro 1 (°) [{atom_names_all[d1_a1]}-{atom_names_all[d1_a2]}-{atom_names_all[d1_a3]}-{atom_names_all[d1_a4]}]",
                        height=540
                    )
                    st_plotly_stretched(fig_contour)

                with col_map2:
                    fig_surf = go.Figure(data=[go.Surface(
                        z=res_2d["delta_z"],
                        x=res_2d["angles_2"],
                        y=res_2d["angles_1"],
                        colorscale="Turbo",
                        colorbar=dict(title="Delta E (kcal/mol)")
                    )])
                    fig_surf.update_layout(
                        title="Superficie di Potenziale 3D",
                        scene=dict(
                            xaxis_title="Diedro 2 (°)",
                            yaxis_title="Diedro 1 (°)",
                            zaxis_title="Delta E (kcal/mol)"
                        ),
                        height=540,
                        margin=dict(l=10, r=10, b=10, t=40)
                    )
                    st_plotly_stretched(fig_surf)

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
            "Delta_G_kcalmol": c.delta_gibbs_kcalmol if c.delta_gibbs_kcalmol is not None else c.delta_dft_kcalmol,
            "Boltzmann_Pop": c.boltzmann_pop,
            "Degeneracy": c.degeneracy
        } for c in res.conformers])

        with col_g1:
            st.markdown("#### Spettro delle Energie Relative (Delta E)")
            fig_bar = px.bar(
                df_chart,
                x="Conformer",
                y="Delta_E_kcalmol",
                color="Delta_E_kcalmol",
                color_continuous_scale="Plasma",
                labels={"Delta_E_kcalmol": "Delta E (kcal/mol)", "Boltzmann_Pop": "Popolazione (%)"},
                hover_data={"Boltzmann_Pop": ":.2f%"},
                title="Differenza Energetica rispetto al Ground State"
            )
            fig_bar.update_layout(xaxis_title="Conformero", yaxis_title="Delta E (kcal/mol)")
            st_plotly_stretched(fig_bar)

        with col_g2:
            calc_temp = getattr(res.config, 'temperature_k', 298.15) if hasattr(res, 'config') and res.config else 298.15
            st.markdown(f"#### Distribuzione di Boltzmann (T = {calc_temp:.1f} K)")
            fig_pie = px.pie(
                df_chart,
                names="Conformer",
                values="Boltzmann_Pop",
                hole=0.45,
                title="Frazioni Molari all'Equilibrio"
            )
            st_plotly_stretched(fig_pie)

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
                labels={"Dipolo_Debye": "Momento di Dipolo mu (Debye)", "Delta_E_kcalmol": "Delta E (kcal/mol)"},
                title="Momento di Dipolo per Conformero (Debye)"
            )
            fig_dip.update_layout(xaxis_title="Conformero", yaxis_title="Momento di Dipolo (Debye)")
            st_plotly_stretched(fig_dip)

        # Matrice di Diversità Conformazionale (Pairwise RMSD)
        if len(res.conformers) > 1:
            st.markdown("---")
            st.markdown("#### Matrice di Diversità Conformazionale (RMSD Pairwise)")
            st.caption("Matrice di similarità geometrica all-vs-all con permutazioni di simmetria (GetBestRMS). Valori espressi in Ångström (Å):")

            rmsd_labels, rmsd_matrix = compute_pairwise_rmsd_matrix(res.smiles, res.conformers)
            if rmsd_labels and rmsd_matrix:
                fig_rmsd = px.imshow(
                    rmsd_matrix,
                    x=rmsd_labels,
                    y=rmsd_labels,
                    color_continuous_scale="Tealgrn",
                    text_auto=".3f",
                    labels=dict(color="RMSD (Å)"),
                    title="Matrice RMSD tra Conformeri Unici (Å)"
                )
                fig_rmsd.update_layout(xaxis_title="Conformero", yaxis_title="Conformero", height=420)
                st_plotly_stretched(fig_rmsd)

        # Diagramma Orbitali di Frontiera HOMO-LUMO
        has_orbitals = any(c.homo_ev is not None and c.lumo_ev is not None for c in res.conformers)
        if has_orbitals:
            st.markdown("---")
            st.markdown("#### Struttura Elettronica di Frontiera (HOMO – LUMO)")
            st.caption("Livelli energetici orbitali di frontiera, band gap fondamentale e indici di reattività chimica:")

            conf_hl_options = [c for c in res.conformers if c.homo_ev is not None and c.lumo_ev is not None]
            c_hl1, c_hl2 = st.columns([1.2, 2])
            with c_hl1:
                sel_hl_conf = st.selectbox(
                    "Conformero da ispezionare:",
                    conf_hl_options,
                    format_func=lambda c: f"Conf {c.conf_id} (Delta E = {c.delta_dft_kcalmol:.2f} kcal/mol)",
                    key="sel_hl_conf"
                )

                homo_val = sel_hl_conf.homo_ev
                lumo_val = sel_hl_conf.lumo_ev
                gap_val = sel_hl_conf.gap_ev if sel_hl_conf.gap_ev is not None else (lumo_val - homo_val)
                ip_approx = -homo_val
                ea_approx = -lumo_val
                hardness = gap_val / 2.0

                st.markdown(f"""
                <div class="info-card">
                    <span class="metric-badge">HOMO: {homo_val:.2f} eV</span>
                    <span class="metric-badge">LUMO: {lumo_val:.2f} eV</span>
                    <span class="metric-badge">Band Gap: {gap_val:.2f} eV</span>
                    <br><br>
                    <span class="metric-badge">Potenziale Ionizzazione (IP): {ip_approx:.2f} eV</span>
                    <span class="metric-badge">Affinità Elettronica (EA): {ea_approx:.2f} eV</span>
                    <span class="metric-badge">Durezza Chimica (η): {hardness:.2f} eV</span>
                </div>
                """, unsafe_allow_html=True)

            with c_hl2:
                fig_hl = go.Figure()
                fig_hl.add_trace(go.Scatter(
                    x=[-0.6, 0.6], y=[homo_val, homo_val],
                    mode="lines",
                    name="HOMO",
                    line=dict(color="#06b6d4", width=5),
                    hoverinfo="text",
                    hovertext=f"HOMO: {homo_val:.3f} eV ({homo_val * 23.06:.1f} kcal/mol)"
                ))
                fig_hl.add_trace(go.Scatter(
                    x=[-0.6, 0.6], y=[lumo_val, lumo_val],
                    mode="lines",
                    name="LUMO",
                    line=dict(color="#f43f5e", width=5),
                    hoverinfo="text",
                    hovertext=f"LUMO: {lumo_val:.3f} eV ({lumo_val * 23.06:.1f} kcal/mol)"
                ))
                fig_hl.add_annotation(
                    x=0.0, y=(homo_val + lumo_val) / 2.0,
                    text=f"<b>ΔE Gap = {gap_val:.2f} eV</b><br>({gap_val * 23.06:.1f} kcal/mol)",
                    showarrow=False,
                    font=dict(size=12, color="#38bdf8")
                )
                fig_hl.add_annotation(
                    x=0.62, y=homo_val,
                    text=f"HOMO ({homo_val:.2f} eV)",
                    showarrow=False,
                    xanchor="left",
                    font=dict(color="#06b6d4", size=11)
                )
                fig_hl.add_annotation(
                    x=0.62, y=lumo_val,
                    text=f"LUMO ({lumo_val:.2f} eV)",
                    showarrow=False,
                    xanchor="left",
                    font=dict(color="#f43f5e", size=11)
                )
                margin_y = max(gap_val * 0.35, 1.5)
                fig_hl.update_layout(
                    title=f"Diagramma Livelli di Frontiera — Conformero {sel_hl_conf.conf_id}",
                    xaxis=dict(showticklabels=False, range=[-1.0, 1.4], showgrid=False, zeroline=False),
                    yaxis=dict(title="Energia Orbitale (eV)", range=[homo_val - margin_y, lumo_val + margin_y]),
                    height=320,
                    showlegend=False
                )
                st_plotly_stretched(fig_hl)

            st.info("Per visualizzare la forma spaziale 3D dei lobi quantistici (&psi; > 0 in blu, &psi; < 0 in rosso) per HOMO, LUMO e orbitali adiacenti, apri la scheda **Struttura 3D e Geometria** e seleziona l'orbitale nel pannello dedicato.", icon=":material/view_in_ar:")

    # --------------------------------------------------------------------------
    # TAB 4: Tabella Dati & Esportazione
    # --------------------------------------------------------------------------
    with tab_table:
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

        # Download Pacchetto Completo ZIP (tutti gli XYZ, Excel, JSON e Report)
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
