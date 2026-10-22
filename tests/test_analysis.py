from energyinator.models import ConformerResult
from energyinator.analysis import calculate_thermodynamics


def test_calculate_thermodynamics():
    # Construct 3 mock conformers
    # Conf 0: ground state
    # Conf 1: gauche state (+0.85 kcal/mol)
    # Conf 2: degenerate gauche state (+0.85 kcal/mol)
    c0 = ConformerResult(conf_id=0, ff_energy=3.0, dft_energy_hartree=-100.0, dft_energy_kcalmol=-62750.9)
    c1 = ConformerResult(conf_id=1, ff_energy=4.0, dft_energy_hartree=-100.0 + 0.85/627.509, dft_energy_kcalmol=-62750.9 + 0.85)
    c2 = ConformerResult(conf_id=2, ff_energy=4.0, dft_energy_hartree=-100.0 + 0.85/627.509, dft_energy_kcalmol=-62750.9 + 0.85)

    res = calculate_thermodynamics([c1, c2, c0], temperature_k=298.15, energy_tol_kcalmol=0.01)

    # Sorted ascending: c0 first
    assert res[0].conf_id == 0
    assert abs(res[0].delta_dft_kcalmol) < 1e-5
    assert abs(res[1].delta_dft_kcalmol - 0.85) < 1e-4
    assert abs(res[2].delta_dft_kcalmol - 0.85) < 1e-4

    # Degeneracy of c1 and c2 must be 2
    assert res[1].degeneracy == 2
    assert res[2].degeneracy == 2
    assert res[0].degeneracy == 1

    # Sum of Boltzmann populations must equal 100%
    total_pop = sum(c.boltzmann_pop for c in res)
    assert abs(total_pop - 100.0) < 1e-3

