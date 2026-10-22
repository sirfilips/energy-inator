from dataclasses import dataclass
from typing import Optional


@dataclass
class CalculationConfig:
    """Configuration settings for conformational search and quantum chemistry calculations."""
    smiles: str = ""
    n_conformers: int = 50
    rmsd_threshold: float = 0.2
    force_field: str = "mmff94"  # "mmff94", "mmff94s", or "uff"
    ff_energy_window: Optional[float] = 10.0  # kcal/mol relative to FF minimum (None to disable)
    top_n: Optional[int] = None  # Max number of conformers for Psi4 (None for all)
    
    # Psi4 quantum calculation settings
    psi4_memory: str = "8 GB"
    psi4_threads: int = 8
    psi4_method: str = "b3lyp-d3bj"  # Default with D3BJ dispersion correction
    psi4_basis: str = "6-31G(d)"
    sp_only: bool = False  # Single Point only (faster) vs full geometry optimization
    solvent: Optional[str] = None  # Implicit solvent (PCM): None / "gas", "water", "chloroform", "dmso", etc.
    calc_frequencies: bool = False  # Compute vibrational frequencies and Gibbs free energy (G)
    freq_top_n: int = 1  # How many top conformers to compute frequencies for (default: 1, global min)
    
    # Post-DFT analysis settings
    post_dft_rmsd: float = 0.15  # RMSD threshold in Angstrom to detect collapsed conformers
    energy_degeneracy_tol: float = 0.05  # kcal/mol threshold for degeneracy clustering
    temperature_k: float = 298.15  # Temperature for Boltzmann distribution calculation
    
    # Molecular properties (auto-detected if None)
    charge: Optional[int] = None
    multiplicity: Optional[int] = None
    
    # Paths & logging
    output_dir: Optional[str] = None
    verbose: bool = True
