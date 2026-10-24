"""Tests for the energyinator.geometry module."""

import numpy as np
import pytest
from rdkit import Chem
from energyinator.geometry import (
    parse_xyz,
    parse_xyz_coords,
    format_xyz,
    calc_distance,
    calc_angle,
    calc_dihedral,
    kabsch_align_xyz,
    mol_from_smiles_and_xyz,
    find_rotatable_dihedrals
)

SAMPLE_XYZ = """3
Water molecule
O   0.00000000   0.00000000   0.11740000
H   0.00000000   0.75700000  -0.46960000
H   0.00000000  -0.75700000  -0.46960000"""


def test_parse_and_format_xyz():
    symbols, coords = parse_xyz(SAMPLE_XYZ)
    assert symbols == ["O", "H", "H"]
    assert coords.shape == (3, 3)
    assert np.isclose(coords[0, 2], 0.1174)

    coords_list = parse_xyz_coords(SAMPLE_XYZ)
    assert len(coords_list) == 3
    assert isinstance(coords_list[0], tuple)

    formatted = format_xyz(symbols, coords, comment="Reconstructed")
    lines = formatted.splitlines()
    assert lines[0] == "3"
    assert lines[1] == "Reconstructed"
    assert lines[2].startswith("O ")


def test_geometric_measures():
    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([1.0, 0.0, 0.0])
    p3 = np.array([1.0, 1.0, 0.0])
    p4 = np.array([1.0, 1.0, 1.0])

    dist = calc_distance(p1, p2)
    assert np.isclose(dist, 1.0)

    angle = calc_angle(p1, p2, p3)
    assert np.isclose(angle, 90.0)

    dihedral = calc_dihedral(p1, p2, p3, p4)
    assert np.isclose(abs(dihedral), 90.0)


def test_kabsch_alignment():
    # Test identical geometry has RMSD near zero
    aligned = kabsch_align_xyz(SAMPLE_XYZ, SAMPLE_XYZ)
    _, coords_orig = parse_xyz(SAMPLE_XYZ)
    _, coords_aligned = parse_xyz(aligned)
    assert np.allclose(coords_orig, coords_aligned, atol=1e-5)


def test_mol_from_smiles_and_xyz():
    mol = mol_from_smiles_and_xyz("O", SAMPLE_XYZ)
    assert mol.GetNumAtoms() == 3
    assert mol.GetNumConformers() == 1
    pos0 = mol.GetConformer(0).GetAtomPosition(0)
    assert np.isclose(pos0.z, 0.1174, atol=1e-3)


def test_find_rotatable_dihedrals():
    # Butane: CCCC has central C-C rotatable bond
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
    dihedrals = find_rotatable_dihedrals(mol)
    assert len(dihedrals) >= 1
    atoms = dihedrals[0]["atoms"]
    assert len(atoms) == 4
    assert len(set(atoms)) == 4

