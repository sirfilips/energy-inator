"""Tab 2: 1D and 2D Dihedral Torsional Scan (Potential Energy Surface)."""

import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from rdkit import Chem

from energyinator.models import CalculationResult
from energyinator.scan import run_dihedral_scan, run_2d_dihedral_scan
from energyinator.geometry import parse_xyz, find_rotatable_dihedrals, mol_from_smiles_and_xyz
from gui.components.viewer3d import render_3d_mol
from gui.styles import st_plotly_stretched


def render_tab_pes(res: CalculationResult):
    """Renders the Potential Energy Surface (PES) 1D and 2D torsional scanning interface."""
    pes_mode = st.radio(
        "Tipologia di Scansione Torsionale:",
        ["Scansione 1D (Diedro Singolo)", "Mappa 2D (Coppia di Diedri / Ramachandran)"],
        horizontal=True
    )

    symbols_all, coords_all = parse_xyz(res.conformers[0].xyz_optimized)
    num_atoms_all = len(symbols_all)
    atom_indices_all = list(range(num_atoms_all))
    atom_names_all = [f"{sym}{i}" for i, sym in enumerate(symbols_all)]

    # Detect rotatable bonds automatically with RDKit
    suggested_dihedrals = []
    try:
        base_mol = mol_from_smiles_and_xyz(res.smiles, res.conformers[0].xyz_optimized)
        suggested_dihedrals = find_rotatable_dihedrals(base_mol)
    except Exception:
        pass

    if pes_mode == "Scansione 1D (Diedro Singolo)":
        st.markdown("### Scansione Torsionale 1D (Potential Energy Surface)")
        st.caption("Calcolo della barriera rotazionale da 0° a 360° con supporto per rilassamento geometrico dei gradi di libertà adiacenti:")

        # Quick selection of rotatable bonds if detected
        if suggested_dihedrals:
            dih_opts = ["Selezione manuale atomi"] + [d["label"] for d in suggested_dihedrals]
            sel_sugg = st.selectbox("Legami rotatabili rilevati automaticamente:", dih_opts, index=1 if len(dih_opts) > 1 else 0)
            if sel_sugg != "Selezione manuale atomi":
                matching = next(d for d in suggested_dihedrals if d["label"] == sel_sugg)
                def_a1, def_a2, def_a3, def_a4 = matching["atoms"]
            else:
                def_a1, def_a2, def_a3, def_a4 = 0, min(1, num_atoms_all-1), min(2, num_atoms_all-1), min(3, num_atoms_all-1)
        else:
            def_a1, def_a2, def_a3, def_a4 = 0, min(1, num_atoms_all-1), min(2, num_atoms_all-1), min(3, num_atoms_all-1)

        cfg_col1, cfg_col2, cfg_col3, cfg_col4, cfg_col5, cfg_col6, cfg_col7 = st.columns([1, 1, 1, 1, 1, 1.2, 1.2])

        pes_a1 = cfg_col1.selectbox("Atomo 1", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(def_a1, num_atoms_all-1), key="tab_pes_a1")
        pes_a2 = cfg_col2.selectbox("Atomo 2", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(def_a2, num_atoms_all-1), key="tab_pes_a2")
        pes_a3 = cfg_col3.selectbox("Atomo 3", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(def_a3, num_atoms_all-1), key="tab_pes_a3")
        pes_a4 = cfg_col4.selectbox("Atomo 4", atom_indices_all, format_func=lambda i: atom_names_all[i], index=min(def_a4, num_atoms_all-1), key="tab_pes_a4")

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
        # 2D Scan
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
                with st.spinner("Calcolo ottimizzato della griglia bidimensionale in corso..."):
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
                    title="Mappa di Potenziale 2D — Diedro 1 vs Diedro 2",
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

