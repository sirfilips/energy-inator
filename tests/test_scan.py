import pytest
from energyinator.scan import run_dihedral_scan, run_2d_dihedral_scan


def test_dihedral_scan_mmff94():
    # Scan butane C0-C1-C2-C3
    res = run_dihedral_scan(
        smiles="CCCC",
        xyz_str="",
        atom_indices=(0, 1, 2, 3),
        step_deg=30,
        engine="mmff94",
        relaxed=True
    )
    
    assert len(res["angles"]) == 13  # 0 to 360 in 30 deg steps
    assert res["min_angle"] == 180  # anti is the global minimum
    assert res["max_barrier_kcalmol"] > 3.0
    assert len(res["xyz_list"]) == 13
    assert res["delta_kcalmol"][6] == 0.0  # 180 degrees corresponds to index 6
    assert res["relaxed"] is True


def test_dihedral_scan_2d():
    # 2D scan on pentane or butane
    res_2d = run_2d_dihedral_scan(
        smiles="CCCC",
        xyz_str="",
        dihedral_1=(0, 1, 2, 3),
        dihedral_2=(0, 1, 2, 3),
        step_deg=60,
        engine="mmff94",
        relaxed=False
    )
    assert len(res_2d["angles_1"]) == 7  # 0, 60, 120, 180, 240, 300, 360
    assert len(res_2d["angles_2"]) == 7
    assert len(res_2d["delta_z"]) == 7
    assert len(res_2d["delta_z"][0]) == 7
    assert res_2d["max_barrier_kcalmol"] >= 0.0
