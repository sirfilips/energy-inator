"""Tests for ensemble observables and conformational entropy in energyinator.analysis."""

import math
from energyinator.models import ConformerResult
from energyinator.analysis import calculate_ensemble_observables


def test_calculate_ensemble_observables():
    # Construct 2 conformers: 50% and 50%
    c1 = ConformerResult(
        conf_id=0,
        ff_energy=1.0,
        boltzmann_pop=50.0,
        dipole_debye=2.0,
        homo_ev=-6.5,
        lumo_ev=-1.5,
        gap_ev=5.0
    )
    c2 = ConformerResult(
        conf_id=1,
        ff_energy=1.2,
        boltzmann_pop=50.0,
        dipole_debye=3.0,
        homo_ev=-6.1,
        lumo_ev=-1.1,
        gap_ev=5.0
    )

    obs = calculate_ensemble_observables([c1, c2], temperature_k=298.15)

    # Weighted dipole should be exactly 2.5 D
    assert math.isclose(obs["ensemble_dipole_debye"], 2.50, abs_tol=1e-2)

    # Weighted HOMO should be -6.3 eV
    assert math.isclose(obs["ensemble_homo_ev"], -6.30, abs_tol=1e-2)

    # Weighted LUMO should be -1.3 eV
    assert math.isclose(obs["ensemble_lumo_ev"], -1.30, abs_tol=1e-2)

    # Conformational entropy for equimolar binary mixture: S = R * ln(2) = 1.9872 * 0.6931 = 1.377 cal/(mol*K)
    expected_s = 1.987204 * math.log(2)
    assert math.isclose(obs["conformational_entropy_cal_mol_k"], expected_s, abs_tol=0.05)

    # Free energy term: -T * S / 1000 = -298.15 * 1.377 / 1000 = -0.41 kcal/mol
    expected_g = -(298.15 * expected_s) / 1000.0
    assert math.isclose(obs["conformational_free_energy_kcalmol"], expected_g, abs_tol=0.05)

