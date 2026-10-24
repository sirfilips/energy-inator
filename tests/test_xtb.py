"""Unit tests for the GFN-xTB semiempirical engine integration."""

import pytest
from energyinator.xtb import (
    find_xtb_executable,
    is_xtb_available,
    get_xtb_version,
    run_xtb_calculation,
    refine_conformers_with_xtb
)
from energyinator.models import ConformerResult

SAMPLE_WATER_XYZ = """3
water
O  0.000000  0.000000  0.000000
H  0.000000 -0.757000  0.586000
H  0.000000  0.757000  0.586000
"""


def test_xtb_discovery():
    """Verifies that xtb executable is discovered on the host system."""
    assert is_xtb_available() is True
    exe = find_xtb_executable()
    assert exe is not None
    assert "xtb" in exe.lower()
    version_str = get_xtb_version()
    assert "6." in version_str or "xtb" in version_str.lower()


def test_xtb_single_point_calculation():
    """Runs a single-point GFN2-xTB calculation on water and validates energy and properties."""
    res = run_xtb_calculation(
        xyz_str=SAMPLE_WATER_XYZ,
        method="gfn2",
        opt=False,
        charge=0,
        multiplicity=1
    )
    assert res["converged"] is True
    assert res["energy_hartree"] is not None
    # Water GFN2 energy is roughly -5.07 Hartree
    assert -5.5 < res["energy_hartree"] < -4.8
    assert res["energy_kcalmol"] < 0
    assert res["gap_ev"] is not None
    assert res["gap_ev"] > 5.0  # HOMO-LUMO gap of water is large (> 10 eV)


def test_xtb_geometry_optimization():
    """Runs a geometry optimization with GFN2-xTB on water."""
    unrelaxed_water = """3
unrelaxed
O  0.000000  0.000000  0.000000
H  0.000000 -0.850000  0.650000
H  0.000000  0.850000  0.650000
"""
    res = run_xtb_calculation(
        xyz_str=unrelaxed_water,
        method="gfn2",
        opt=True,
        charge=0,
        multiplicity=1
    )
    assert res["converged"] is True
    assert res["xyz_optimized"] is not None
    assert "O" in res["xyz_optimized"]
    assert "H" in res["xyz_optimized"]


def test_refine_conformers_with_xtb():
    """Tests sorting, Delta E computation, and filtering of conformers using xTB."""
    c0 = ConformerResult(conf_id=0, ff_energy=2.0, xyz_initial=SAMPLE_WATER_XYZ)
    c1 = ConformerResult(conf_id=1, ff_energy=3.0, xyz_initial=SAMPLE_WATER_XYZ)

    refined = refine_conformers_with_xtb(
        conformers=[c0, c1],
        method="gfn2",
        opt=False,
        charge=0,
        multiplicity=1,
        energy_window_kcalmol=5.0
    )
    assert len(refined) == 2
    assert refined[0].delta_xtb_kcalmol == 0.0
    assert refined[0].xtb_energy_hartree is not None

