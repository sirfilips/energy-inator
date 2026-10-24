"""Ketcher 2D Molecular Sketcher integration component."""

import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    from streamlit_ketcher import st_ketcher
    HAS_KETCHER = True
except ImportError:
    HAS_KETCHER = False


def smiles_to_2d_molblock(smiles: str) -> str:
    """Generates clean 2D coordinates and converts a SMILES string to MolBlock V2000 for Ketcher."""
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
    """Converts Ketcher output (MolBlock or SMILES) to canonical SMILES."""
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


def render_ketcher_editor():
    """Renders the Ketcher 2D molecular editor with synchronization controls."""
    with st.expander("Editor Molecolare 2D (Ketcher)", expanded=False):
        if HAS_KETCHER:
            st.caption("Disegna o modifica la struttura molecolare, poi premi **'Applica Struttura nel Calcolo'** per sincronizzarla:")
            
            current_target = st.session_state.get("active_smiles", "CCCC")
            ketcher_key = f"ketcher_canvas_{st.session_state.get('ketcher_version', 0)}"

            molblock_2d = smiles_to_2d_molblock(current_target)
            if molblock_2d:
                ketcher_val = molblock_2d
                mol_format = "MOLFILE"
            else:
                ketcher_val = current_target
                mol_format = "SMILES"

            ketcher_output = st_ketcher(value=ketcher_val, height=440, molecule_format=mol_format, key=ketcher_key)

            col_k1, col_k2 = st.columns([1.8, 3])
            btn_apply_ketcher = col_k1.button(
                "Applica Struttura nel Calcolo",
                icon=":material/check:",
                type="primary",
                key="btn_apply_ketcher"
            )

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

