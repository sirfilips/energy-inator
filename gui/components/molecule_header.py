"""Header component displaying molecule details, configuration badges, and 2D topology."""

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors, Descriptors
from gui.styles import st_image_stretched


@st.cache_data(show_spinner=False)
def generate_2d_depiction(smiles: str):
    """Draws 2D molecular structure with atomic labels and indices."""
    if not smiles or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetProp("atomLabel", f"{atom.GetSymbol()}{atom.GetIdx()}")
    return Draw.MolToImage(mol, size=(380, 240))


def get_mol_formula_and_mw(smiles: str):
    """Calculates chemical formula and molecular weight."""
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol:
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mw = Descriptors.MolWt(mol)
            return formula, mw
    except Exception:
        pass
    return "", 0.0


def render_molecule_header(smiles: str, method: str, basis: str, solvent: str, ff: str, profile_name: str):
    """Renders the top banner containing configuration badges and 2D topology."""
    st.markdown('<div class="main-title">Energyinator Studio</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Framework per modellistica molecolare, campionamento conformazionale e chimica quantistica (RDKit + Psi4)</div>', unsafe_allow_html=True)

    formula, mw = get_mol_formula_and_mw(smiles) if smiles else ("", 0.0)
    formula_badge = f'<span class="metric-badge" style="background: rgba(16, 185, 129, 0.12); color: #059669; border-color: rgba(16, 185, 129, 0.3);">Formula: {formula} ({mw:.2f} g/mol)</span>' if formula else ""

    col_top1, col_top2 = st.columns([2.8, 1.2])

    with col_top1:
        solv_display = solvent if solvent else "Gas / Vuoto"
        st.markdown(f"""
        <div class="info-card">
            {formula_badge}
            <span class="metric-badge">SMILES: {smiles.strip() or 'N/A'}</span>
            <span class="metric-badge">Metodo: {method.upper()} / {basis}</span>
            <span class="metric-badge">Solvente: {solv_display}</span>
            <span class="metric-badge">Campo di Forze: {ff.upper()}</span>
            <span class="metric-badge">Profilo: {profile_name}</span>
        </div>
        """, unsafe_allow_html=True)

    with col_top2:
        if smiles.strip():
            img_2d = generate_2d_depiction(smiles.strip())
            if img_2d:
                st_image_stretched(img_2d, caption="Topologia 2D con indici atomici")
            else:
                st.warning("Struttura SMILES non valida.")
