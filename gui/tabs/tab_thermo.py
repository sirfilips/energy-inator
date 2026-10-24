"""Tab 3: Energetic Spectrum, Boltzmann Distributions, Dipoles, and Frontier Orbitals."""

import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

from energyinator.models import CalculationResult
from energyinator.conformers import compute_pairwise_rmsd_matrix
from gui.styles import st_plotly_stretched


def render_tab_thermo(res: CalculationResult):
    """Renders Boltzmann population analysis, energy levels, dipole moments, and HOMO-LUMO diagrams."""
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

    # Ensemble-averaged Observables & Conformational Entropy Card
    ens_dip = getattr(res, "ensemble_dipole_debye", None)
    ens_gap = getattr(res, "ensemble_gap_ev", None)
    s_conf = getattr(res, "conformational_entropy_cal_mol_k", None)
    g_conf = getattr(res, "conformational_free_energy_kcalmol", None)

    st.markdown("---")
    st.markdown("#### Proprietà Medie d'Insieme & Termodinamica Conformazionale")
    st.caption("Proprietà macroscopiche osservabili calcolate come media ponderata di Boltzmann su tutto l'insieme:")

    c_ens1, c_ens2, c_ens3, c_ens4 = st.columns(4)
    c_ens1.metric(
        label="Dipolo Medio ⟨μ⟩",
        value=f"{ens_dip:.2f} D" if ens_dip is not None else "-",
        delta=f"Ground: {res.conformers[0].dipole_debye:.2f} D" if (res.conformers and res.conformers[0].dipole_debye is not None) else None,
        delta_color="off"
    )
    c_ens2.metric(
        label="Band Gap Medio ⟨Gap⟩",
        value=f"{ens_gap:.2f} eV" if ens_gap is not None else "-",
        delta=f"Ground: {res.conformers[0].gap_ev:.2f} eV" if (res.conformers and res.conformers[0].gap_ev is not None) else None,
        delta_color="off"
    )
    c_ens3.metric(
        label="Entropia Conformazionale S_conf",
        value=f"{s_conf:.2f} cal/(mol*K)" if s_conf is not None else "-",
        help="Contributo entropico configurazionale di miscelamento dell'insieme: S = -R * sum(p_i * ln(p_i))"
    )
    c_ens4.metric(
        label="Stabilizzazione Libera (-T·S_conf)",
        value=f"{g_conf:.3f} kcal/mol" if g_conf is not None else "-",
        help="Stabilizzazione sull'energia libera di Gibbs derivante dalla molteplicità dell'insieme conformazionale a 298.15 K."
    )
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

    # Conformational Diversity Matrix (Pairwise RMSD)
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

    # HOMO-LUMO diagram
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

