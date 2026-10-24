"""Conformational Benchmark and Scientific Validation Suite.

Provides standard test molecular systems with theoretical / experimental literature references:
  - n-Butane: Anti vs Gauche rotational equilibrium
  - Ethylene Glycol (1,2-ethanediol): Intramolecular hydrogen bonding (gauche vs anti)
  - Cyclohexane: Chair vs Boat conformational energy difference
"""

import time
import logging
from typing import Dict, Any, List, Optional
import pandas as pd

from energyinator.config import CalculationConfig
from energyinator.pipeline import run_pipeline
from energyinator.models import CalculationResult

logger = logging.getLogger("energyinator")

# Literature-curated benchmark reference data
BENCHMARK_SYSTEMS = {
    "butane": {
        "name": "n-Butano (Equilibrio Anti vs Gauche)",
        "smiles": "CCCC",
        "description": "Barriera e differenza energetica tra conformazione estesa (Anti) e ripiegata (Gauche).",
        "expected_min_conformation": "Anti",
        "reference_delta_e_kcalmol": 0.85,  # Gauche is ~0.85 kcal/mol higher than Anti
        "tolerance_kcalmol": 0.40
    },
    "ethylene_glycol": {
        "name": "1,2-Etandiolo (Legame H Intramolecolare)",
        "smiles": "OCCO",
        "description": "Test di legame a idrogeno interno: la conformazione Gauche è stabilizzata dal legame H rispetto all'Anti.",
        "expected_min_conformation": "Gauche (H-bonded)",
        "reference_delta_e_kcalmol": 1.20,  # Gauche lower than anti in vacuum due to H-bond
        "tolerance_kcalmol": 0.60
    },
    "cyclohexane": {
        "name": "Cicloesano (Sedia vs Barca)",
        "smiles": "C1CCCCC1",
        "description": "Grande differenza conformazionale tra la forma a sedia (ground state) e forme barca/twist.",
        "expected_min_conformation": "Sedia (Chair)",
        "reference_delta_e_kcalmol": 5.70,
        "tolerance_kcalmol": 1.50
    }
}


def run_benchmark(
    systems: Optional[List[str]] = None,
    method: str = "b3lyp",
    basis: str = "3-21G",
    force_field: str = "mmff94",
    sp_only: bool = True
) -> Dict[str, Any]:
    """Executes the standard conformational benchmark suite and computes statistical accuracy metrics."""
    target_keys = systems or list(BENCHMARK_SYSTEMS.keys())
    results_rows = []
    start_total = time.time()

    for key in target_keys:
        if key not in BENCHMARK_SYSTEMS:
            continue
        item = BENCHMARK_SYSTEMS[key]
        cfg = CalculationConfig(
            smiles=item["smiles"],
            n_conformers=15,
            rmsd_threshold=0.20,
            force_field=force_field,
            psi4_method=method,
            psi4_basis=basis,
            sp_only=sp_only,
            psi4_threads=2,
            psi4_memory="2 GB"
        )
        
        t0 = time.time()
        status = "PASS"
        notes = "Calcolo completato"
        res: Optional[CalculationResult] = None
        computed_delta = 0.0

        try:
            res = run_pipeline(cfg)
            elapsed = time.time() - t0
            
            if len(res.conformers) >= 2:
                computed_delta = res.conformers[1].delta_dft_kcalmol or 0.0
                err = abs(computed_delta - item["reference_delta_e_kcalmol"])
                if err > item["tolerance_kcalmol"]:
                    status = "WARN"
                    notes = f"Delta E calcolato ({computed_delta:.2f}) scosta dal ref ({item['reference_delta_e_kcalmol']:.2f})"
            else:
                status = "WARN"
                notes = "Meno di 2 conformeri unici isolati"

            results_rows.append({
                "System": item["name"],
                "SMILES": item["smiles"],
                "Status": status,
                "Confs_DFT": res.unique_dft,
                "Delta_E_Calc": round(computed_delta, 2),
                "Delta_E_Ref": item["reference_delta_e_kcalmol"],
                "Error": round(abs(computed_delta - item["reference_delta_e_kcalmol"]), 2),
                "Elapsed_s": round(elapsed, 1),
                "Notes": notes
            })

        except Exception as e:
            results_rows.append({
                "System": item["name"],
                "SMILES": item["smiles"],
                "Status": "FAIL",
                "Confs_DFT": 0,
                "Delta_E_Calc": 0.0,
                "Delta_E_Ref": item["reference_delta_e_kcalmol"],
                "Error": None,
                "Elapsed_s": round(time.time() - t0, 1),
                "Notes": f"Errore: {e}"
            })

    total_time = round(time.time() - start_total, 2)
    df = pd.DataFrame(results_rows)
    mae = df["Error"].dropna().mean() if not df["Error"].dropna().empty else 0.0

    return {
        "df_results": df,
        "mean_absolute_error_kcalmol": round(float(mae), 3),
        "total_elapsed_s": total_time,
        "method": method,
        "basis": basis,
        "force_field": force_field
    }


def main():
    print("=" * 65)
    print(" ENERGYINATOR - Conformational Benchmark & Validation Suite")
    print("=" * 65)
    bench = run_benchmark(method="b3lyp", basis="3-21G", sp_only=True)
    print("\nRISULTATI DEL BENCHMARK:")
    print(bench["df_results"].to_string(index=False))
    print("-" * 65)
    print(f"Errore Medio Assoluto (MAE): {bench['mean_absolute_error_kcalmol']} kcal/mol")
    print(f"Tempo totale esecuzione:    {bench['total_elapsed_s']} s")
    print("=" * 65)


if __name__ == "__main__":
    main()
