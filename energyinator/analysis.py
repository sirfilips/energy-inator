import math
import logging
from typing import List
from energyinator.models import ConformerResult

logger = logging.getLogger("energyinator")

# Universal gas constant in kcal / (mol * K)
R_GAS_KCAL = 1.98720425864083e-3
HARTREE_TO_KCALMOL = 627.5094740631


def calculate_thermodynamics(
    conformers: List[ConformerResult],
    temperature_k: float = 298.15,
    energy_tol_kcalmol: float = 0.05
) -> List[ConformerResult]:
    """Calculates relative energies (Delta E in kcal/mol and kJ/mol),
    statistical degeneracy for isoenergetic states, and Boltzmann populations at specified temperature.
    Supports Gibbs free energies (Delta G) when vibrational frequency data is available.
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
        # Clamp tiny negative numerical precision artifacts to 0.0
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

    # 2. Group into degenerate states if within energy_tol_kcalmol
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

    # Assign degeneracy count (or preserve existing multiplicity if already higher)
    for cl in clusters:
        deg = len(cl)
        for conf in cl:
            conf.degeneracy = max(conf.degeneracy, deg)

    # 3. Calculate Boltzmann populations
    # P_i = exp(-Delta E_i / (R * T)) / Partition_Function
    # If all conformers have Gibbs free energy, prioritize Gibbs
    use_gibbs = len(valid_gibbs) == len(valid_confs)
    rt = R_GAS_KCAL * temperature_k
    boltzmann_factors = []

    for conf in valid_confs:
        energy_term = conf.delta_gibbs_kcalmol if (use_gibbs and conf.delta_gibbs_kcalmol is not None) else conf.delta_dft_kcalmol
        exponent = -energy_term / rt
        factor = math.exp(exponent) if exponent > -700 else 0.0
        boltzmann_factors.append(factor)

    q_partition = sum(boltzmann_factors)
    if q_partition > 0:
        for conf, factor in zip(valid_confs, boltzmann_factors):
            conf.boltzmann_pop = (factor / q_partition) * 100.0
    else:
        for conf in valid_confs:
            conf.boltzmann_pop = 0.0

    return valid_confs
