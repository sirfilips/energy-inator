import pytest
from energyinator.quantum import generate_orbital_cube

H2O_XYZ = """3
Water molecule
O  0.000000  0.000000 -0.079136
H  0.000000 -0.707107  0.627971
H  0.000000  0.707107  0.627971"""


def test_generate_orbital_cube_homo():
    res = generate_orbital_cube(
        xyz_str=H2O_XYZ,
        orbital_type="HOMO",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.45
    )
    assert res is not None
    assert res["orbital_type"] == "HOMO"
    assert res["orbital_idx"] == 5  # H2O has 10 electrons -> 5 doubly occupied orbitals
    assert isinstance(res["orbital_energy_ev"], float)
    assert res["orbital_energy_ev"] < 0  # Bound occupied orbital
    assert isinstance(res["cube_data"], str)
    assert len(res["cube_data"]) > 1000
    assert "Psi4" in res["cube_data"] or "Cube" in res["cube_data"]


def test_generate_orbital_cube_lumo():
    homo_res = generate_orbital_cube(
        xyz_str=H2O_XYZ,
        orbital_type="HOMO",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.5
    )
    lumo_res = generate_orbital_cube(
        xyz_str=H2O_XYZ,
        orbital_type="LUMO",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.5
    )
    assert homo_res is not None
    assert lumo_res is not None
    assert lumo_res["orbital_idx"] == 6
    assert lumo_res["orbital_energy_ev"] > homo_res["orbital_energy_ev"]


def test_generate_orbital_cube_invalid():
    # Test graceful handling with invalid orbital or broken XYZ
    bad_res = generate_orbital_cube(
        xyz_str="bad xyz",
        orbital_type="HOMO",
        method="scf",
        basis="sto-3g"
    )
    assert bad_res is None

