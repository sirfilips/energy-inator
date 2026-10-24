"""3D Molecular Visualization Component powered by py3Dmol and WebGL."""

import py3Dmol
import streamlit.components.v1 as components
from energyinator.geometry import parse_xyz, kabsch_align_xyz


def render_3d_mol(
    xyz_str: str,
    style: str = 'ball_and_stick',
    show_labels: bool = True,
    spin: bool = False,
    surface_type: str = 'Nessuna',
    surface_opacity: float = 0.65,
    orbital_cube_data: str = None,
    orbital_isoval: float = 0.03,
    orbital_opacity: float = 0.75,
    height: int = 580
):
    """Generates and renders an interactive WebGL 3D molecular viewer with py3Dmol."""
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

    # Molecular orbital isosurface rendering (positive lobe in cyan/blue, negative lobe in coral/red)
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
        for idx, (sym, pt) in enumerate(zip(symbols, coords)):
            view.addLabel(
                f"{sym}{idx}",
                {
                    'position': {'x': float(pt[0]), 'y': float(pt[1]), 'z': float(pt[2])},
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


def render_3d_mol_overlay(ref_xyz_str: str, active_xyz_str: str, height: int = 580):
    """Renders two overlaid conformers in 3D: reference in semitransparent cyan, active aligned in CPK."""
    aligned_active = kabsch_align_xyz(ref_xyz_str, active_xyz_str)
    view = py3Dmol.view(width=None, height=height)

    # Model 0: Reference (cyan)
    view.addModel(ref_xyz_str, "xyz")
    view.setStyle({'model': 0}, {
        'stick': {'color': '#00bcd4', 'opacity': 0.40, 'radius': 0.12},
        'sphere': {'color': '#00bcd4', 'opacity': 0.40, 'scale': 0.22}
    })

    # Model 1: Aligned active conformer (CPK standard)
    view.addModel(aligned_active, "xyz")
    view.setStyle({'model': 1}, {
        'stick': {'radius': 0.15},
        'sphere': {'scale': 0.28}
    })

    view.zoomTo()
    html = view._make_html()
    components.html(html, height=height + 20)

