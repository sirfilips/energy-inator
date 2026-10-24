"""Tests for the energyinator.benchmark validation suite."""

from energyinator.benchmark import BENCHMARK_SYSTEMS, run_benchmark


def test_benchmark_systems_catalog():
    assert "butane" in BENCHMARK_SYSTEMS
    assert "ethylene_glycol" in BENCHMARK_SYSTEMS
    assert "cyclohexane" in BENCHMARK_SYSTEMS

    for key, data in BENCHMARK_SYSTEMS.items():
        assert "smiles" in data
        assert "reference_delta_e_kcalmol" in data
        assert data["reference_delta_e_kcalmol"] > 0.0
