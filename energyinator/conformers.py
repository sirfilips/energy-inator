import logging
from typing import List, Optional, Tuple, TYPE_CHECKING
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

if TYPE_CHECKING:
    from energyinator.models import ConformerResult

logger = logging.getLogger("energyinator")


def prepare_molecule(smiles: str) -> Tuple[Chem.Mol, int, int]:
    """Parses a SMILES string, adds explicit hydrogens, and determines formal charge and spin multiplicity."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES non valido o impossibile da convertire in molecola: '{smiles}'")
    
    mol = Chem.AddHs(mol)
    formal_charge = Chem.GetFormalCharge(mol)
    
    # Calculate spin multiplicity: 2 * S + 1 = total_radicals + 1 (assuming parallel radical spins)
    num_radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    multiplicity = num_radicals + 1
    
    return mol, formal_charge, multiplicity


def generate_conformers_rdkit(
    mol: Chem.Mol,
    n_conformers: int = 50,
    rmsd_threshold: float = 0.2
) -> List[int]:
    """Generates initial 3D conformers using RDKit's ETKDGv3 algorithm with multithreading."""
    logger.info(f"Generazione di {n_conformers} conformeri con algoritmo ETKDGv3 (filtro RMSD: {rmsd_threshold} Å)...")
    
    params = AllChem.ETKDGv3()
    params.pruneRmsThresh = float(rmsd_threshold)
    params.numThreads = 0  # Use all available CPU cores
    
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(n_conformers), params=params))
    
    # Fallback to standard EmbedMultipleConfs if ETKDGv3 fails on older/unusual geometries
    if not conf_ids:
        logger.warning("ETKDGv3 non ha generato conformeri. Tentativo con parametri standard...")
        conf_ids = list(AllChem.EmbedMultipleConfs(
            mol,
            numConfs=int(n_conformers),
            pruneRmsThresh=float(rmsd_threshold),
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
            numThreads=0
        ))
        
    if not conf_ids:
        raise ValueError("Impossibile generare conformeri 3D per la molecola fornita.")
        
    logger.info(f"Generati {len(conf_ids)} conformeri iniziali.")
    return conf_ids


def optimize_force_field(
    mol: Chem.Mol,
    conf_ids: List[int],
    force_field: str = "mmff94"
) -> Tuple[List[Tuple[int, float]], str]:
    """Optimizes all conformers using MMFF94/MMFF94s or UFF.
    
    Falls back to UFF if MMFF parameters are missing for the molecule.
    Returns:
        List of tuples: (conf_id, ff_energy_kcal_mol)
        Effective force field name used
    """
    ff_choice = force_field.lower().strip()
    effective_ff = ff_choice
    results = []

    use_mmff = ff_choice in ["mmff94", "mmff94s"]
    mmff_variant = "MMFF94s" if ff_choice == "mmff94s" else "MMFF94"

    if use_mmff:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            logger.info(f"Ottimizzazione conformeri con campo di forze {mmff_variant} (multithread)...")
            try:
                # Parallel optimization across all conformers
                AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant=mmff_variant, numThreads=0, maxIters=500)
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=mmff_variant)
                for cid in conf_ids:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=cid)
                    if ff is not None:
                        results.append((cid, ff.CalcEnergy()))
                    else:
                        logger.warning(f"Impossibile calcolare energia {mmff_variant} per conformero {cid}")
            except Exception as e:
                logger.warning(f"Errore durante ottimizzazione {mmff_variant}: {e}. Fallback su UFF.")
                results = []
                use_mmff = False
        else:
            logger.warning(f"Molecola priva di parametri completi per {mmff_variant}. Fallback automatico su UFF.")
            use_mmff = False

    if not use_mmff:
        effective_ff = "uff"
        logger.info("Ottimizzazione conformeri con Universal Force Field (UFF) (multithread)...")
        AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=500)
        for cid in conf_ids:
            try:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
                if ff is not None:
                    results.append((cid, ff.CalcEnergy()))
            except Exception as e:
                logger.warning(f"Errore calcolo UFF conformero {cid}: {e}")

    return results, effective_ff


def filter_conformers(
    mol: Chem.Mol,
    conf_energies: List[Tuple[int, float]],
    rmsd_threshold: float = 0.2,
    energy_window: Optional[float] = 10.0,
    top_n: Optional[int] = None
) -> List[Tuple[int, float, float]]:
    """Filters conformers by Force Field energy window, symmetrical RMSD, and optional top_n cutoff.
    
    Returns:
        List of tuples: (conf_id, ff_energy, delta_ff_kcalmol)
    """
    if not conf_energies:
        return []

    # Sort ascending by FF energy
    sorted_confs = sorted(conf_energies, key=lambda x: x[1])
    min_ff_energy = sorted_confs[0][1]

    # Step 1: Energy window filtering (discard conformers with high relative strain)
    filtered_by_window = []
    for conf_id, energy in sorted_confs:
        delta_ff = energy - min_ff_energy
        if energy_window is not None and energy_window > 0 and delta_ff > energy_window:
            continue
        filtered_by_window.append((conf_id, energy, delta_ff))

    logger.info(
        f"Conformeri entro la finestra energetica FF di {energy_window} kcal/mol: "
        f"{len(filtered_by_window)} / {len(conf_energies)}"
    )

    # Step 2: RMSD deduplication with molecular symmetry consideration
    unique_confs = []
    ref_ids = []
    
    for conf_id, energy, delta_ff in filtered_by_window:
        is_unique = True
        for ref_id in ref_ids:
            # GetBestRMS computes RMSD testing all automorphism permutations
            rms = AllChem.GetBestRMS(mol, mol, ref_id, conf_id)
            if rms < rmsd_threshold:
                is_unique = False
                break
        
        if is_unique:
            unique_confs.append((conf_id, energy, delta_ff))
            ref_ids.append(conf_id)
            
            if top_n is not None and len(unique_confs) >= top_n:
                logger.info(f"Raggiunto il limite top_n ({top_n}) conformeri unici da calcolare.")
                break

    logger.info(f"Conformeri unici selezionati per calcolo DFT: {len(unique_confs)}")
    return unique_confs


def parse_xyz_coords(xyz_str: str) -> List[Tuple[float, float, float]]:
    """Extracts 3D coordinates from an XYZ formatted string."""
    lines = xyz_str.strip().splitlines()
    if len(lines) >= 2 and lines[0].strip().isdigit():
        lines = lines[2:]
    coords = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
    return coords


def deduplicate_post_dft(
    base_mol: Chem.Mol,
    conformer_results: List["ConformerResult"],
    rmsd_threshold: float = 0.15
) -> List["ConformerResult"]:
    """Deduplicates conformers after quantum geometry optimization.
    
    If multiple starting structures collapsed into the same quantum potential energy well
    (RMSD < rmsd_threshold), keeps only the lowest-energy conformer and increments its degeneracy.
    """
    if rmsd_threshold <= 0 or len(conformer_results) <= 1:
        return conformer_results

    from rdkit.Geometry import Point3D

    # Sort ascending by DFT energy (Hartree)
    sorted_confs = sorted(conformer_results, key=lambda c: c.dft_energy_hartree if c.dft_energy_hartree is not None else 1e9)

    # Build reference molecule containing conformers
    mol_work = Chem.Mol(base_mol)
    mol_work.RemoveAllConformers()

    valid_confs = []
    for c in sorted_confs:
        if not c.xyz_optimized:
            continue
        coords = parse_xyz_coords(c.xyz_optimized)
        if len(coords) != mol_work.GetNumAtoms():
            continue
        conf = Chem.Conformer(mol_work.GetNumAtoms())
        for idx, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(idx, Point3D(x, y, z))
        cid = mol_work.AddConformer(conf, assignId=True)
        valid_confs.append((c, cid))

    if not valid_confs:
        return conformer_results

    unique_results: List["ConformerResult"] = []
    accepted_cids: List[int] = []

    for c, cid in valid_confs:
        collapsed_match = None
        for ref_conf, ref_cid in zip(unique_results, accepted_cids):
            rms = AllChem.GetBestRMS(mol_work, mol_work, ref_cid, cid)
            if rms < rmsd_threshold:
                collapsed_match = ref_conf
                break

        if collapsed_match is not None:
            collapsed_match.degeneracy += 1
            logger.info(
                f"Conformero {c.conf_id} collassato su {collapsed_match.conf_id} "
                f"(RMSD: {rms:.3f} < {rmsd_threshold} Å). Degenerazione: {collapsed_match.degeneracy}."
            )
        else:
            unique_results.append(c)
            accepted_cids.append(cid)

    logger.info(f"Deduplicazione post-DFT: {len(unique_results)} conformeri unici su {len(valid_confs)} calcolati.")
    return unique_results


def compute_pairwise_rmsd_matrix(
    smiles: str,
    conformer_results: List["ConformerResult"]
) -> Tuple[List[str], List[List[float]]]:
    """Calcola la matrice di diversità RMSD simmetrico-allineata (tutti contro tutti) tra conformeri unici."""
    mol_base = Chem.MolFromSmiles(smiles)
    if mol_base is None:
        return [], []
    mol_work = Chem.AddHs(mol_base)
    mol_work.RemoveAllConformers()

    labels: List[str] = []
    cids: List[int] = []

    for c in conformer_results:
        if not c.xyz_optimized:
            continue
        coords = parse_xyz_coords(c.xyz_optimized)
        if len(coords) != mol_work.GetNumAtoms():
            continue
        conf = Chem.Conformer(mol_work.GetNumAtoms())
        for idx, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(idx, Point3D(x, y, z))
        cid = mol_work.AddConformer(conf, assignId=True)
        cids.append(cid)
        labels.append(f"Conf {c.conf_id}")

    n = len(cids)
    if n == 0:
        return [], []

    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            rms = AllChem.GetBestRMS(mol_work, mol_work, cids[i], cids[j])
            matrix[i][j] = round(rms, 3)
            matrix[j][i] = round(rms, 3)

    return labels, matrix
