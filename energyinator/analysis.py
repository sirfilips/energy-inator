import math
import logging
from typing import List, Dict, Optional
from energyinator.models import ConformerResult

logger = logging.getLogger("energyinator")

# Universal gas constant in kcal / (mol * K)
R_GAS_KCAL = 1.98720425864083e-3
R_GAS_CAL = 1.98720425864083
HARTREE_TO_KCALMOL = 627.5094740631


def calculate_thermodynamics(
    conformers: List[ConformerResult],
    temperature_k: float = 298.15,
    energy_tol_kcalmol: float = 0.05
) -> List[ConformerResult]:
    """Calculates relative energies (Delta E in kcal/mol and kJ/mol),
    statistical degeneracy handling, and Boltzmann populations at specified temperature.
    Supports Gibbs free energies (Delta G) when vibrational frequency data is available.
    
    Boltzmann populations correctly incorporate conformer degeneracies (g_i):
        P_i = (g_i * exp(-Delta E_i / RT)) / Q
        Q = sum_j (g_j * exp(-Delta E_j / RT))
    """
    valid_confs = [c for c in conformers if c.dft_energy_hartree is not None]
    if not valid_confs:
        return conformers

    # Sort ascending by quantum energy (global minimum first)
    valid_confs.sort(key=lambda c: c.dft_energy_hartree)
    min_energy_kcal = valid_confs[0].dft_energy_kcalmol

    # 1. Compute Delta E (electronic)
    for conf in valid_confs:
        delta_kcal = conf.dft_energy_kcalmol - min_energy_kcal
        # Clamp tiny numerical precision artifacts to 0.0
        if abs(delta_kcal) < 1e-6:
            delta_kcal = 0.0
        conf.delta_dft_kcalmol = delta_kcal
        conf.delta_dft_kjmol = delta_kcal * 4.184

    # 1b. Compute Delta Gibbs (free energy) if available
    valid_gibbs = [c.gibbs_hartree for c in valid_confs if c.gibbs_hartree is not None]
    if valid_gibbs:
        min_gibbs = min(valid_gibbs)
        for conf in valid_confs:
            if conf.gibbs_hartree is not None:
                d_g = (conf.gibbs_hartree - min_gibbs) * HARTREE_TO_KCALMOL
                conf.delta_gibbs_kcalmol = 0.0 if abs(d_g) < 1e-6 else d_g

    # 2. Boltzmann population calculation with statistical degeneracy (g_i)
    use_gibbs = len(valid_gibbs) == len(valid_confs)
    rt = R_GAS_KCAL * temperature_k
    boltzmann_factors = []

    for conf in valid_confs:
        energy_term = conf.delta_gibbs_kcalmol if (use_gibbs and conf.delta_gibbs_kcalmol is not None) else conf.delta_dft_kcalmol
        exponent = -energy_term / rt
        boltz_exp = math.exp(exponent) if exponent > -700 else 0.0
        # Incorporate intrinsic degeneracy (collapsed microstates / symmetry)
        g_i = max(1, getattr(conf, "degeneracy", 1))
        boltzmann_factors.append(g_i * boltz_exp)

    q_partition = sum(boltzmann_factors)
    if q_partition > 0:
        for conf, factor in zip(valid_confs, boltzmann_factors):
            conf.boltzmann_pop = (factor / q_partition) * 100.0
    else:
        for conf in valid_confs:
            conf.boltzmann_pop = 0.0

    # 3. Group into degenerate states if within energy_tol_kcalmol (reporting/clustering)
    clusters = []
    for conf in valid_confs:
        matched_cluster = None
        for cl in clusters:
            if abs(conf.delta_dft_kcalmol - cl[0].delta_dft_kcalmol) <= energy_tol_kcalmol:
                matched_cluster = cl
                break
        if matched_cluster is not None:
            matched_cluster.append(conf)
        else:
            clusters.append([conf])

    for cl in clusters:
        deg = len(cl)
        for conf in cl:
            conf.degeneracy = max(conf.degeneracy, deg)

    return valid_confs


def calculate_ensemble_observables(
    conformers: List[ConformerResult],
    temperature_k: float = 298.15
) -> Dict[str, Optional[float]]:
    """Calculates Boltzmann ensemble-averaged properties and conformational entropy.
    
    Observables computed:
      - ensemble_dipole_debye: <mu> = sum(p_i * mu_i)
      - ensemble_homo_ev: <HOMO> = sum(p_i * HOMO_i)
      - ensemble_lumo_ev: <LUMO> = sum(p_i * LUMO_i)
      - ensemble_gap_ev: <Gap> = sum(p_i * Gap_i)
      - conformational_entropy_cal_mol_k: S_conf = -R * sum(p_i * ln(p_i))
      - conformational_free_energy_kcalmol: -T * S_conf
    """
    valid_confs = [c for c in conformers if c.boltzmann_pop is not None and c.boltzmann_pop > 0]
    if not valid_confs:
        return {
            "ensemble_dipole_debye": None,
            "ensemble_homo_ev": None,
            "ensemble_lumo_ev": None,
            "ensemble_gap_ev": None,
            "conformational_entropy_cal_mol_k": None,
            "conformational_free_energy_kcalmol": None
        }

    total_pop = sum(c.boltzmann_pop for c in valid_confs)
    weights = [c.boltzmann_pop / total_pop for c in valid_confs]

    # 1. Conformational entropy of mixing: S_conf = -R * sum(p_i * ln(p_i))
    s_conf = 0.0
    for p in weights:
        if p > 1e-12:
            s_conf -= R_GAS_CAL * p * math.log(p)

    # Free energy stabilization from conformational multiplicity: -T * S_conf (kcal/mol)
    g_conf = -(temperature_k * s_conf) / 1000.0

    # 2. Ensemble-averaged dipole moment
    dipoles = [c.dipole_debye for c in valid_confs if c.dipole_debye is not None]
    ensemble_dipole = None
    if len(dipoles) == len(valid_confs):
        ensemble_dipole = round(float(sum(w * d for w, d in zip(weights, dipoles))), 2)

    # 3. Ensemble-averaged HOMO, LUMO, Gap
    homos = [c.homo_ev for c in valid_confs if c.homo_ev is not None]
    lumos = [c.lumo_ev for c in valid_confs if c.lumo_ev is not None]
    gaps = [c.gap_ev for c in valid_confs if c.gap_ev is not None]

    ensemble_homo = None
    ensemble_lumo = None
    ensemble_gap = None
    if len(homos) == len(valid_confs):
        ensemble_homo = round(float(sum(w * h for w, h in zip(weights, homos))), 2)
    if len(lumos) == len(valid_confs):
        ensemble_lumo = round(float(sum(w * l for w, l in zip(weights, lumos))), 2)
    if len(gaps) == len(valid_confs):
        ensemble_gap = round(float(sum(w * g for w, g in zip(weights, gaps))), 2)
    elif ensemble_homo is not None and ensemble_lumo is not None:
        ensemble_gap = round(ensemble_lumo - ensemble_homo, 2)

    return {
        "ensemble_dipole_debye": ensemble_dipole,
        "ensemble_homo_ev": ensemble_homo,
        "ensemble_lumo_ev": ensemble_lumo,
        "ensemble_gap_ev": ensemble_gap,
        "conformational_entropy_cal_mol_k": round(float(s_conf), 2),
        "conformational_free_energy_kcalmol": round(float(g_conf), 3)
    }
