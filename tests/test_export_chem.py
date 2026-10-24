"""Tests for the energyinator.export_chem module."""

from energyinator.models import ConformerResult, CalculationResult
from energyinator.export_chem import (
    generate_orca_input,
    generate_gaussian_input,
    generate_multi_xyz_trajectory,
    generate_run_recipe
)

SAMPLE_XYZ = """3
Water
O   0.00000000   0.00000000   0.11740000
H   0.00000000   0.75700000  -0.46960000
H   0.00000000  -0.75700000  -0.46960000"""


def test_generate_orca_input():
    c = ConformerResult(conf_id=0, ff_energy=1.0, xyz_optimized=SAMPLE_XYZ)
    orca_inp = generate_orca_input(
        conf=c,
        method="b3lyp",
        basis="def2-SVP",
        solvent="water",
        charge=0,
        multiplicity=1,
        nprocs=4
    )
    assert "! B3LYP D3BJ def2-SVP Opt Freq CPCM(water)" in orca_inp
    assert "%pal nprocs 4 end" in orca_inp
    assert "* xyz 0 1" in orca_inp
    assert "O " in orca_inp
    assert orca_inp.strip().endswith("*")


def test_generate_gaussian_input():
    c = ConformerResult(conf_id=0, ff_energy=1.0, xyz_optimized=SAMPLE_XYZ)
    gau_inp = generate_gaussian_input(
        conf=c,
        method="b3lyp",
        basis="6-31G(d)",
        solvent="water",
        charge=0,
        multiplicity=1,
        nprocs=8,
        memory="16GB"
    )
    assert "%nprocshared=8" in gau_inp
    assert "%mem=16GB" in gau_inp
    assert "#p b3lyp/6-31G(d) opt freq scrf=(pcm,solvent=water)" in gau_inp
    assert "0 1" in gau_inp
    assert "O " in gau_inp


def test_generate_multi_xyz_trajectory():
    c1 = ConformerResult(conf_id=0, ff_energy=1.0, dft_energy_hartree=-76.4, delta_dft_kcalmol=0.0, boltzmann_pop=60.0, xyz_optimized=SAMPLE_XYZ)
    c2 = ConformerResult(conf_id=1, ff_energy=1.5, dft_energy_hartree=-76.399, delta_dft_kcalmol=0.63, boltzmann_pop=40.0, xyz_optimized=SAMPLE_XYZ)
    traj = generate_multi_xyz_trajectory([c1, c2], smiles="O")
    
    assert traj.count("3\n") == 2
    assert "Frame 0 | Conf 0" in traj
    assert "Frame 1 | Conf 1" in traj


def test_generate_run_recipe():
    res = CalculationResult(
        smiles="O",
        formal_charge=0,
        spin_multiplicity=1,
        method="b3lyp",
        basis="6-31G(d)",
        force_field="mmff94",
        min_energy_hartree=-76.4,
        min_energy_kcalmol=-47941.0,
        ensemble_dipole_debye=1.85,
        conformational_entropy_cal_mol_k=1.38
    )
    recipe = generate_run_recipe(res)
    assert recipe["target_molecule"]["smiles"] == "O"
    assert recipe["computational_parameters"]["method"] == "b3lyp"
    assert recipe["ensemble_results_summary"]["ensemble_dipole_debye"] == 1.85
    assert recipe["ensemble_results_summary"]["conformational_entropy_cal_mol_k"] == 1.38
    assert "rdkit_version" in recipe["platform"]

