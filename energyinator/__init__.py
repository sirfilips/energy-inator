"""Energyinator - Hybrid Conformer Generation and Quantum Chemistry Energy Calculator."""

from energyinator.config import CalculationConfig
from energyinator.models import ConformerResult, CalculationResult
from energyinator.pipeline import run_pipeline
from energyinator.scan import run_dihedral_scan, run_2d_dihedral_scan
from energyinator import geometry, export_chem, benchmark


def calcola_energia_psi4(
    smiles: str,
    n_conformers: int = 50,
    rmsd_threshold: float = 0.2,
    force_field: str = "mmff94",
    psi4_method: str = "b3lyp",
    psi4_basis: str = "6-31G(d)"
):
    """Legacy compatibility function returning the 11-element tuple:
    (min_e_hartree, min_e_kcal, min_e_ev, best_xyz, best_path,
     max_e_hartree, max_e_kcal, max_e_ev, worst_xyz, worst_path, elapsed_time)
    """
    config = CalculationConfig(
        smiles=smiles,
        n_conformers=n_conformers,
        rmsd_threshold=rmsd_threshold,
        force_field=force_field,
        psi4_method=psi4_method,
        psi4_basis=psi4_basis
    )
    result = run_pipeline(config)
    return result.as_legacy_tuple()


__all__ = [
    "CalculationConfig",
    "CalculationResult",
    "ConformerResult",
    "run_pipeline",
    "run_dihedral_scan",
    "run_2d_dihedral_scan",
    "calcola_energia_psi4",
    "geometry",
    "export_chem",
    "benchmark"
]
