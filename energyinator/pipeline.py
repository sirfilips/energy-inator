import os
import time
import logging
from datetime import datetime
from tqdm import tqdm
from rdkit import Chem

from energyinator.config import CalculationConfig
from energyinator.models import ConformerResult, CalculationResult
from energyinator.conformers import (
    prepare_molecule,
    generate_conformers_rdkit,
    optimize_force_field,
    filter_conformers,
    deduplicate_post_dft
)
from energyinator.quantum import setup_psi4, close_psi4, run_psi4_calculation, cleanup_psi4_residual_files
from energyinator.analysis import calculate_thermodynamics
from energyinator.io import setup_run_logger, close_run_logger, save_conformer_xyz, export_results_excel

logger = logging.getLogger("energyinator")


def run_pipeline(config: CalculationConfig) -> CalculationResult:
    """Executes the complete hybrid conformational search and quantum chemistry workflow."""
    start_time = time.time()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # 1. Setup output directory
    if config.output_dir:
        output_dir = os.path.abspath(config.output_dir)
    else:
        script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        output_dir = os.path.join(script_dir, "output", f"result_{timestamp}")
        
    os.makedirs(output_dir, exist_ok=True)
    log_handler = setup_run_logger(output_dir)

    logger.info("=" * 65)
    logger.info("  ENERGYINATOR - Hybrid Quantum Conformer Search Pipeline")
    logger.info("=" * 65)
    logger.info(f"SMILES target: {config.smiles}")
    if config.solvent:
        logger.info(f"Solvatazione implicita (PCM): {config.solvent.upper()}")

    try:
        # 2. Prepare molecule and determine electronic state
        mol, auto_charge, auto_mult = prepare_molecule(config.smiles)
        charge = config.charge if config.charge is not None else auto_charge
        multiplicity = config.multiplicity if config.multiplicity is not None else auto_mult
        logger.info(f"Stato elettronico molecola: Carica={charge}, Molteplicità di spin={multiplicity}")

        # 3. Generate initial 3D conformers with RDKit
        initial_ids = generate_conformers_rdkit(
            mol,
            n_conformers=config.n_conformers,
            rmsd_threshold=config.rmsd_threshold
        )

        # 4. Classical Force Field optimization (MMFF94 / UFF)
        ff_results, effective_ff = optimize_force_field(mol, initial_ids, force_field=config.force_field)

        # 5. Symmetrical RMSD and Energy-window filtering
        selected_confs = filter_conformers(
            mol,
            ff_results,
            rmsd_threshold=config.rmsd_threshold,
            energy_window=config.ff_energy_window,
            top_n=config.top_n
        )

        if not selected_confs:
            raise RuntimeError("Nessun conformero ha superato la fase di filtraggio classico.")

        # 6. Setup Psi4
        scratch_dir = os.path.join(output_dir, "scratch")
        setup_psi4(
            memory=config.psi4_memory,
            threads=config.psi4_threads,
            basis=config.psi4_basis,
            log_file=os.path.join(output_dir, "psi4_output.dat"),
            scratch_dir=scratch_dir
        )

        # 7. Quantum calculations
        smiles_clean = Chem.MolToSmiles(Chem.RemoveHs(mol))
        calc_type_desc = "Single-Point" if config.sp_only else "Ottimizzazione"
        logger.info(
            f"Avvio calcoli quantistici Psi4 ({calc_type_desc}) su {len(selected_confs)} conformeri "
            f"[Metodo: {config.psi4_method}, Base: {config.psi4_basis}]..."
        )

        conformer_results = []
        effective_method_used = config.psi4_method
        total_selected = len(selected_confs)

        for idx, (conf_id, ff_energy, delta_ff) in enumerate(tqdm(selected_confs, desc="Psi4 Quantum", unit="conf"), 1):
            logger.info(f"Progresso DFT: [{idx}/{total_selected}] Conformero {conf_id} ({idx}/{total_selected})...")
            xyz_block = Chem.MolToXYZBlock(mol, confId=conf_id)
            
            calc_out = run_psi4_calculation(
                xyz_block=xyz_block,
                conf_id=conf_id,
                charge=charge,
                multiplicity=multiplicity,
                method=config.psi4_method,
                basis=config.psi4_basis,
                sp_only=config.sp_only,
                solvent=config.solvent,
                calc_frequencies=False,
                smiles_label=smiles_clean
            )

            if calc_out is not None:
                effective_method_used = calc_out.get("method_used", config.psi4_method)
                c_res = ConformerResult(
                    conf_id=conf_id,
                    ff_energy=ff_energy,
                    delta_ff_kcalmol=delta_ff,
                    dft_energy_hartree=calc_out["energy_hartree"],
                    dft_energy_kcalmol=calc_out["energy_kcalmol"],
                    dft_energy_ev=calc_out["energy_ev"],
                    dipole_debye=calc_out.get("dipole_debye"),
                    homo_ev=calc_out.get("homo_ev"),
                    lumo_ev=calc_out.get("lumo_ev"),
                    gap_ev=calc_out.get("gap_ev"),
                    xyz_initial=xyz_block,
                    xyz_optimized=calc_out["xyz_optimized"]
                )
                conformer_results.append(c_res)

        if not conformer_results:
            raise RuntimeError("Nessun conformero è stato calcolato con successo con Psi4.")

        # 8. Post-DFT RMSD Deduplication (Filter collapsed geometries)
        if config.post_dft_rmsd > 0 and len(conformer_results) > 1 and not config.sp_only:
            logger.info(f"Verifica conformeri collassati post-DFT (soglia RMSD: {config.post_dft_rmsd} Å)...")
            unique_confs = deduplicate_post_dft(
                base_mol=mol,
                conformer_results=conformer_results,
                rmsd_threshold=config.post_dft_rmsd
            )
        else:
            unique_confs = conformer_results

        # 9. Optional Vibrational Frequencies & Gibbs free energy on top conformers
        if config.calc_frequencies and not config.sp_only:
            n_freq = min(max(1, config.freq_top_n), len(unique_confs))
            logger.info(f"Calcolo frequenze vibrazionali e correzioni termiche sui primi {n_freq} conformeri...")
            for i in range(n_freq):
                c = unique_confs[i]
                freq_out = run_psi4_calculation(
                    xyz_block=c.xyz_optimized,
                    conf_id=c.conf_id,
                    charge=charge,
                    multiplicity=multiplicity,
                    method=config.psi4_method,
                    basis=config.psi4_basis,
                    sp_only=False,
                    solvent=config.solvent,
                    calc_frequencies=True,
                    smiles_label=f"{smiles_clean}_freq"
                )
                if freq_out:
                    c.gibbs_hartree = freq_out.get("gibbs_hartree")
                    c.zpve_hartree = freq_out.get("zpve_hartree")
                    c.enthalpy_hartree = freq_out.get("enthalpy_hartree")
                    c.entropy_cal_mol_k = freq_out.get("entropy_cal_mol_k")

        # Save XYZ files for unique conformers
        for c_res in unique_confs:
            save_conformer_xyz(output_dir, c_res, timestamp)

        # 10. Post-DFT Analysis (Relative energies, Degeneracy, Boltzmann populations)
        processed_confs = calculate_thermodynamics(
            unique_confs,
            temperature_k=config.temperature_k,
            energy_tol_kcalmol=config.energy_degeneracy_tol
        )

        # 11. Export Excel
        excel_path = export_results_excel(output_dir, processed_confs, timestamp)

        # 12. Summary result building
        best_conf = processed_confs[0]
        worst_conf = processed_confs[-1]
        elapsed = time.time() - start_time

        has_freqs = any(c.gibbs_hartree is not None for c in processed_confs)

        res = CalculationResult(
            smiles=config.smiles,
            formal_charge=charge,
            spin_multiplicity=multiplicity,
            method=effective_method_used,
            basis=config.psi4_basis,
            force_field=effective_ff,
            solvent=config.solvent,
            total_generated=len(initial_ids),
            unique_ff=len(selected_confs),
            calculated_dft=len(conformer_results),
            unique_dft=len(processed_confs),
            min_energy_hartree=best_conf.dft_energy_hartree,
            min_energy_kcalmol=best_conf.dft_energy_kcalmol,
            min_energy_ev=best_conf.dft_energy_ev,
            max_energy_hartree=worst_conf.dft_energy_hartree,
            max_energy_kcalmol=worst_conf.dft_energy_kcalmol,
            max_energy_ev=worst_conf.dft_energy_ev,
            has_frequencies=has_freqs,
            best_xyz_optimized=best_conf.xyz_optimized,
            best_output_path=best_conf.output_path,
            worst_xyz_optimized=worst_conf.xyz_optimized,
            worst_output_path=worst_conf.output_path,
            conformers=processed_confs,
            elapsed_time=elapsed,
            output_dir=output_dir,
            excel_path=excel_path
        )

        # 13. Save structured JSON for GUI and automated pipelines
        json_path = os.path.join(output_dir, "result.json")
        try:
            res.save_json(json_path)
            res.json_path = json_path
            logger.info(f"Riepilogo strutturato salvato in: {json_path}")
        except Exception as e:
            logger.warning(f"Impossibile salvare il riepilogo JSON: {e}")

        logger.info("=" * 65)
        logger.info(f"Calcolo completato in {elapsed:.2f} secondi!")
        logger.info(f"Minimo globale: {res.min_energy_hartree:.6f} Eh ({res.min_energy_kcalmol:.2f} kcal/mol)")
        logger.info(f"Geometria migliore: {res.best_output_path}")
        logger.info("=" * 65)

        return res

    finally:
        # Crucial for Windows: release psi4 output and file logging locks, and clean scratch
        close_psi4()
        cleanup_psi4_residual_files(output_dir)
        cleanup_psi4_residual_files(os.getcwd())
        close_run_logger(log_handler)
