import os
import tempfile
from energyinator.config import CalculationConfig
from energyinator.pipeline import run_pipeline


def test_pipeline_smoke():
    with tempfile.TemporaryDirectory() as tmpdir:
        config = CalculationConfig(
            smiles="O",  # Water molecule
            n_conformers=2,
            rmsd_threshold=0.1,
            force_field="uff",
            psi4_method="b3lyp",
            psi4_basis="sto-3g",
            sp_only=True,  # Quick single point calculation
            psi4_threads=1,
            psi4_memory="1 GB",
            output_dir=tmpdir
        )
        res = run_pipeline(config)
        
        assert res.smiles == "O"
        assert res.formal_charge == 0
        assert res.spin_multiplicity == 1
        assert res.min_energy_hartree < -70.0
        assert len(res.conformers) >= 1
        assert res.excel_path is not None
        assert res.json_path is not None and os.path.exists(res.json_path)

