from dataclasses import dataclass, field, asdict
import json
from typing import List, Optional


@dataclass
class ConformerResult:
    """Stores geometric and energetic data for a single conformer."""
    conf_id: int
    ff_energy: float
    delta_ff_kcalmol: float = 0.0

    # Semiempirical GFN-xTB results (Tier Intermedio)
    xtb_energy_hartree: Optional[float] = None
    xtb_energy_kcalmol: Optional[float] = None
    delta_xtb_kcalmol: Optional[float] = None
    xtb_gap_ev: Optional[float] = None
    xtb_dipole_debye: Optional[float] = None

    # Quantum chemistry results (electronic)
    dft_energy_hartree: Optional[float] = None
    dft_energy_kcalmol: Optional[float] = None
    dft_energy_ev: Optional[float] = None
    delta_dft_kcalmol: Optional[float] = None
    delta_dft_kjmol: Optional[float] = None
    
    # Thermochemistry (Vibrational / Gibbs free energy)
    gibbs_hartree: Optional[float] = None
    delta_gibbs_kcalmol: Optional[float] = None
    zpve_hartree: Optional[float] = None
    enthalpy_hartree: Optional[float] = None
    entropy_cal_mol_k: Optional[float] = None
    
    # Statistical thermodynamics
    boltzmann_pop: Optional[float] = None  # in %
    degeneracy: int = 1  # Number of equivalent/degenerate conformers
    
    # Electronic & molecular properties
    dipole_debye: Optional[float] = None
    homo_ev: Optional[float] = None
    lumo_ev: Optional[float] = None
    gap_ev: Optional[float] = None
    
    # Coordinates & files
    xyz_initial: str = ""
    xyz_optimized: Optional[str] = None
    output_path: Optional[str] = None

    @property
    def xyz(self) -> str:
        """Returns optimized XYZ string if available, otherwise initial XYZ."""
        return self.xyz_optimized or self.xyz_initial


@dataclass
class CalculationResult:
    """Comprehensive summary of the entire conformational search and quantum calculation."""
    smiles: str
    formal_charge: int
    spin_multiplicity: int
    method: str
    basis: str
    force_field: str
    solvent: Optional[str] = None
    used_xtb: bool = False
    xtb_method: Optional[str] = None
    
    total_generated: int = 0
    unique_ff: int = 0
    calculated_dft: int = 0
    unique_dft: int = 0
    
    min_energy_hartree: float = 0.0
    min_energy_kcalmol: float = 0.0
    min_energy_ev: float = 0.0
    
    max_energy_hartree: float = 0.0
    max_energy_kcalmol: float = 0.0
    max_energy_ev: float = 0.0
    
    has_frequencies: bool = False
    
    best_xyz_optimized: Optional[str] = None
    best_output_path: Optional[str] = None
    worst_xyz_optimized: Optional[str] = None
    worst_output_path: Optional[str] = None
    
    # Ensemble-averaged observables (Boltzmann weighted)
    ensemble_dipole_debye: Optional[float] = None
    ensemble_homo_ev: Optional[float] = None
    ensemble_lumo_ev: Optional[float] = None
    ensemble_gap_ev: Optional[float] = None

    # Configurational / Conformational entropy
    conformational_entropy_cal_mol_k: Optional[float] = None
    conformational_free_energy_kcalmol: Optional[float] = None

    conformers: List[ConformerResult] = field(default_factory=list)
    elapsed_time: float = 0.0
    output_dir: str = ""
    excel_path: Optional[str] = None
    json_path: Optional[str] = None

    def to_dict(self) -> dict:
        """Converts the CalculationResult into a JSON-serializable dictionary."""
        return asdict(self)

    def save_json(self, filepath: str) -> None:
        """Serializes the calculation results to a JSON file."""
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_dict(cls, data: dict) -> "CalculationResult":
        """Reconstructs a CalculationResult object from a dictionary."""
        data_copy = dict(data)
        confs_raw = data_copy.pop("conformers", [])
        
        # Filter valid ConformerResult fields
        conf_fields = ConformerResult.__dataclass_fields__.keys()
        conformers = [
            ConformerResult(**{k: v for k, v in c.items() if k in conf_fields})
            for c in confs_raw
        ]
        
        # Filter valid CalculationResult fields
        calc_fields = cls.__dataclass_fields__.keys()
        filtered_data = {k: v for k, v in data_copy.items() if k in calc_fields}
        return cls(conformers=conformers, **filtered_data)

    @classmethod
    def from_json(cls, filepath: str) -> "CalculationResult":
        """Loads and reconstructs a CalculationResult from a JSON file."""
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)
        return cls.from_dict(data)

    def as_legacy_tuple(self):
        """Maintains backward compatibility with original return signature:
        (min_e_hartree, min_e_kcal, min_e_ev, best_xyz, best_path,
         max_e_hartree, max_e_kcal, max_e_ev, worst_xyz, worst_path, elapsed_time)
        """
        return (
            self.min_energy_hartree, self.min_energy_kcalmol, self.min_energy_ev,
            self.best_xyz_optimized, self.best_output_path,
            self.max_energy_hartree, self.max_energy_kcalmol, self.max_energy_ev,
            self.worst_xyz_optimized, self.worst_output_path,
            self.elapsed_time
        )
