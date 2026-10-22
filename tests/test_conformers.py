from energyinator.conformers import (
    prepare_molecule,
    generate_conformers_rdkit,
    optimize_force_field,
    filter_conformers,
    deduplicate_post_dft
)
from energyinator.models import ConformerResult
from rdkit import Chem


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


def test_deduplicate_post_dft():
    mol, _, _ = prepare_molecule("CCCC")
    confs = generate_conformers_rdkit(mol, n_conformers=5, rmsd_threshold=0.1)
    xyz_0 = Chem.MolToXYZBlock(mol, confId=confs[0])
    
    # Create two identical conformer results
    c0 = ConformerResult(conf_id=0, ff_energy=2.0, dft_energy_hartree=-157.0, dft_energy_kcalmol=-98000.0, xyz_optimized=xyz_0)
    c1 = ConformerResult(conf_id=1, ff_energy=2.1, dft_energy_hartree=-157.0, dft_energy_kcalmol=-98000.0, xyz_optimized=xyz_0)

    # With post-DFT deduplication, c1 should collapse into c0
    unique = deduplicate_post_dft(mol, [c0, c1], rmsd_threshold=0.15)
    assert len(unique) == 1
    assert unique[0].conf_id == 0
    assert unique[0].degeneracy == 2


def test_compute_pairwise_rmsd_matrix():
    from energyinator.conformers import compute_pairwise_rmsd_matrix
    mol, _, _ = prepare_molecule("CCCC")
    confs = generate_conformers_rdkit(mol, n_conformers=3, rmsd_threshold=0.01)
    if len(confs) >= 2:
        xyz_0 = Chem.MolToXYZBlock(mol, confId=confs[0])
        xyz_1 = Chem.MolToXYZBlock(mol, confId=confs[1])
        c0 = ConformerResult(conf_id=0, ff_energy=2.0, dft_energy_hartree=-157.0, dft_energy_kcalmol=-98000.0, xyz_optimized=xyz_0)
        c1 = ConformerResult(conf_id=1, ff_energy=2.5, dft_energy_hartree=-156.9, dft_energy_kcalmol=-97998.0, xyz_optimized=xyz_1)
        labels, matrix = compute_pairwise_rmsd_matrix("CCCC", [c0, c1])
        assert len(labels) == 2
        assert len(matrix) == 2
        assert matrix[0][0] == 0.0
        assert matrix[1][1] == 0.0
        assert matrix[0][1] == matrix[1][0]
        assert matrix[0][1] >= 0.0
