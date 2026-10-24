"""Tests for the energyinator.benchmark validation suite."""

import pandas as pd
from energyinator.benchmark import (
    BENCHMARK_SYSTEMS,
    run_benchmark,
    generate_benchmark_csv,
    generate_benchmark_markdown
)


def test_benchmark_systems_catalog():
    assert "butane" in BENCHMARK_SYSTEMS
    assert "ethylene_glycol" in BENCHMARK_SYSTEMS
    assert "cyclohexane" in BENCHMARK_SYSTEMS

    for key, data in BENCHMARK_SYSTEMS.items():
        assert "smiles" in data
        assert "reference_delta_e_kcalmol" in data
        assert "process_description" in data
        assert "citation" in data
        assert "scientific_significance" in data
        assert data["reference_delta_e_kcalmol"] > 0.0


def test_benchmark_csv_and_markdown_generation():
    df_dummy = pd.DataFrame([{
        "Sistema": "n-Butano",
        "SMILES": "CCCC",
        "Processo_Conformazionale": "Equilibrio torsionale Anti -> Gauche",
        "Conformero_Rif_A": "Anti (180°)",
        "Conformero_Rif_B": "Gauche (60°)",
        "Metodo_Calcolo": "B3LYP/6-31G(d)",
        "Campo_Forze": "MMFF94",
        "Delta_E_Calc_kcal_mol": 0.88,
        "Delta_E_Ref_kcal_mol": 0.86,
        "Livello_Riferimento": "CCSD(T)/CBS",
        "Errore_Assoluto_kcal_mol": 0.02,
        "Tolleranza_kcal_mol": 0.40,
        "Esito_Validazione": "ACCURATO (PASS)",
        "Tempo_s": 5.4,
        "Note_Geometriche": "Anti (180.0°) vs Gauche (64.2°)",
        "Citazione_Letteratura": "Dixon et al., J. Phys. Chem. A 2000",
        "Significato_Scientifico": "Repulsione sterica metili"
    }])

    dummy_bench = {
        "df_results": df_dummy,
        "mean_absolute_error_kcalmol": 0.02,
        "total_elapsed_s": 5.4,
        "method": "B3LYP",
        "basis": "6-31G(d)",
        "force_field": "MMFF94",
        "solvent": "Fase Gassosa (Vuoto)",
        "n_systems": 1,
        "n_passed": 1,
        "timestamp": "2026-10-24 13:30:00"
    }

    csv_text = generate_benchmark_csv(dummy_bench)
    assert "REPORT DI VALIDAZIONE METODOLOGICA" in csv_text
    assert "B3LYP/6-31G(d)" in csv_text
    assert "Delta_E_Calc_kcal_mol" in csv_text
    assert "0.88" in csv_text
    assert "Dixon et al." in csv_text

    md_text = generate_benchmark_markdown(dummy_bench)
    assert "# Report di Validazione Metodologica & Benchmark Conformazionale" in md_text
    assert "MAE" in md_text
    assert "0.02 kcal/mol" in md_text
    assert "n-Butano" in md_text
    assert "CCSD(T)/CBS" in md_text
