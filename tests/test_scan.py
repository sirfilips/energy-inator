import pytest
from energyinator.scan import run_dihedral_scan


def test_dihedral_scan_mmff94():
    # Scan butane C0-C1-C2-C3
    res = run_dihedral_scan(
        smiles="CCCC",
        xyz_str="",
        atom_indices=(0, 1, 2, 3),
        step_deg=30,
        engine="mmff94"
    )
    
    assert len(res["angles"]) == 13  # 0 to 360 in 30 deg steps
    assert res["min_angle"] == 180  # anti is the global minimum
    assert res["max_barrier_kcalmol"] > 5.0  # eclipsed barrier is significant
    assert len(res["xyz_list"]) == 13
    assert res["delta_kcalmol"][6] == 0.0  # 180 degrees corresponds to index 6

