"""Hero Summary Dashboard Component providing a clean executive overview of calculation results."""

import streamlit as st
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from energyinator.models import CalculationResult


def get_molecular_formula_and_weight(smiles: str):
    """Calculates chemical formula and molecular weight using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mw = Descriptors.MolWt(mol)
            return formula, mw
    except Exception:
        pass
    return "N/A", 0.0


def render_hero_summary(res: CalculationResult):
    """Renders the executive summary card using native Streamlit containers and metrics."""
    formula, mw = get_molecular_formula_and_weight(res.smiles)
    
    best_conf = res.conformers[0] if res.conformers else None
    ground_pop = f"{best_conf.boltzmann_pop:.1f}%" if best_conf and best_conf.boltzmann_pop is not None else "N/A"
    ground_deg = best_conf.degeneracy if best_conf else 1
    delta_span = (res.max_energy_kcalmol - res.min_energy_kcalmol)
    
    ens_dip = getattr(res, "ensemble_dipole_debye", None)
    ens_gap = getattr(res, "ensemble_gap_ev", None)
    dipole_str = f"⟨μ⟩ = {ens_dip:.2f} D" if ens_dip is not None else (f"{best_conf.dipole_debye:.2f} D" if best_conf and best_conf.dipole_debye is not None else "-")
    gap_str = f"⟨Gap⟩ = {ens_gap:.2f} eV" if ens_gap is not None else (f"{best_conf.gap_ev:.2f} eV" if best_conf and best_conf.gap_ev is not None else "-")
    solv_label = res.solvent.capitalize() if res.solvent else "Gas / Vuoto"

    with st.container(border=True):
        col_hdr1, col_hdr2 = st.columns([3, 2])
        col_hdr1.markdown(f"#### Sintesi Risultati — {formula} *(PM: {mw:.2f} g/mol)*")
        col_hdr2.caption(f"**Livello teorico:** {res.method.upper()}/{res.basis} &bull; **Ambiente:** {solv_label}")

        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric(
            label="Minimo Globale (E₀)",
            value=f"{res.min_energy_hartree:.5f} Eh",
            delta=f"{res.min_energy_kcalmol:.1f} kcal/mol",
            delta_color="normal"
        )
        c2.metric(
            label="Conformeri Unici",
            value=f"{res.unique_dft} / {res.total_generated}",
            delta=f"ΔE: +{delta_span:.2f} kcal/mol",
            delta_color="off"
        )
        c3.metric(
            label="Popolazione Ground State",
            value=ground_pop,
            delta=f"Molteplicità g₀ = {ground_deg}",
            delta_color="off"
        )
        c4.metric(
            label="Frontiera Elettronica",
            value=gap_str,
            delta=f"Dipolo: {dipole_str}",
            delta_color="off"
        )
        c5.metric(
            label="Tempo Calcolo",
            value=f"{res.elapsed_time:.1f} s",
            delta=f"FF: {res.force_field.upper()}",
            delta_color="off"
        )
