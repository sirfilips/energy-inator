import math
import logging
from typing import List
from energyinator.models import ConformerResult

logger = logging.getLogger("energyinator")

# Universal gas constant in kcal / (mol * K)
R_GAS_KCAL = 1.98720425864083e-3


def calculate_thermodynamics(
    conformers: List[ConformerResult],
    temperature_k: float = 298.15,
    energy_tol_kcalmol: float = 0.05
) -> List[ConformerResult]:
    """Calculates relative energies (Delta E in kcal/mol and kJ/mol),
    statistical degeneracy for isoenergetic states, and Boltzmann populations at specified temperature.
    """
    valid_confs = [c for c in conformers if c.dft_energy_hartree is not None]
    if not valid_confs:
        return conformers

    # Sort ascending by quantum energy (global minimum first)
    valid_confs.sort(key=lambda c: c.dft_energy_hartree)
    min_energy_hartree = valid_confs[0].dft_energy_hartree
    min_energy_kcal = valid_confs[0].dft_energy_kcalmol

    # 1. Compute Delta E
    for conf in valid_confs:
        delta_kcal = conf.dft_energy_kcalmol - min_energy_kcal
        # Clamp tiny negative numerical precision artifacts to 0.0
        if abs(delta_kcal) < 1e-6:
            delta_kcal = 0.0
        conf.delta_dft_kcalmol = delta_kcal
        conf.delta_dft_kjmol = delta_kcal * 4.184

    # 2. Group into degenerate states if within energy_tol_kcalmol
    # Group tracking
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

    # Assign degeneracy count
    for cl in clusters:
        deg = len(cl)
        for conf in cl:
            conf.degeneracy = deg

    # 3. Calculate Boltzmann populations
    # P_i = exp(-Delta E_i / (R * T)) / Partition_Function
    rt = R_GAS_KCAL * temperature_k
    boltzmann_factors = []
    for conf in valid_confs:
        # Exponent protection against underflow
        exponent = -conf.delta_dft_kcalmol / rt
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

