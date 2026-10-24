"""Tab 1: 3D Molecular Structure Inspector & Geometric Measurements."""

import streamlit as st
from energyinator.models import CalculationResult, ConformerResult
from energyinator.quantum import generate_orbital_cube
from energyinator.geometry import parse_xyz, calc_distance, calc_angle, calc_dihedral
from gui.components.viewer3d import render_3d_mol, render_3d_mol_overlay


def render_tab_3d(res: CalculationResult):
    """Renders the 3D structure viewer, electronic properties, orbital isosurface, and geometry inspector."""
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

        # Electronic & Thermochemical properties
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
                    <div style="background-color: #1e293b; border-radius: 8px; padding: 12px 16px; margin: 10px 0; border: 1px solid #334155; color: #f8fafc;">
                        <div style="font-size: 0.95rem; font-weight: 600; color: #ffffff; margin-bottom: 4px;">
                            {cached_cube['orbital_type']} (Orbitale Molecolare #{cached_cube['orbital_idx']})
                        </div>
                        <div style="font-size: 0.85rem; color: #cbd5e1; margin-bottom: 8px;">
                            Energia: <strong style="color: #38bdf8;">{cached_cube['orbital_energy_ev']:.3f} eV</strong> &bull; Metodo: <span style="color: #e2e8f0;">{cached_cube['method'].upper()}/{cached_cube['basis']}</span>
                        </div>
                        <div style="font-size: 0.85rem; display: flex; align-items: center; gap: 16px; flex-wrap: wrap;">
                            <span style="display: inline-flex; align-items: center; color: #93c5fd; font-weight: 500;">
                                <span style="display: inline-block; width: 13px; height: 13px; background: #2563eb; border-radius: 3px; margin-right: 6px; box-shadow: 0 0 4px #2563eb;"></span>
                                Lobo positivo (&psi; &gt; 0)
                            </span>
                            <span style="display: inline-flex; align-items: center; color: #fca5a5; font-weight: 500;">
                                <span style="display: inline-block; width: 13px; height: 13px; background: #dc2626; border-radius: 3px; margin-right: 6px; box-shadow: 0 0 4px #dc2626;"></span>
                                Lobo negativo (&psi; &lt; 0)
                            </span>
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

