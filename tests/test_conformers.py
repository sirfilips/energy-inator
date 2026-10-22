from energyinator.conformers import (
    prepare_molecule,
    generate_conformers_rdkit,
    optimize_force_field,
    filter_conformers
)


def test_prepare_molecule_neutral():
    mol, charge, mult = prepare_molecule("CCO")
    assert mol is not None
    assert charge == 0
    assert mult == 1


def test_prepare_molecule_charged():
    mol, charge, mult = prepare_molecule("CC(=O)[O-]")
    assert mol is not None
    assert charge == -1
    assert mult == 1


def test_prepare_molecule_radical():
    mol, charge, mult = prepare_molecule("[CH3]")
    assert mol is not None
    assert charge == 0
    assert mult == 2


def test_generate_and_filter_conformers():
    mol, _, _ = prepare_molecule("CCCC")  # Butane
    confs = generate_conformers_rdkit(mol, n_conformers=10, rmsd_threshold=0.2)
    assert len(confs) > 0

    energies, eff_ff = optimize_force_field(mol, confs, force_field="mmff94")
    assert len(energies) > 0
    assert eff_ff in ["mmff94", "uff"]

    filtered = filter_conformers(mol, energies, rmsd_threshold=0.2, energy_window=5.0)
    assert len(filtered) > 0
    # Minimum energy conformer must have delta_ff == 0.0
    assert abs(filtered[0][2]) < 1e-6

