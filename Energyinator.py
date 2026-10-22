"""Energyinator: Conformer Generation and Quantum Chemistry Energy Calculator.

Backward-compatible entry point delegating to the modular energyinator package.
"""

import sys
from energyinator.config import CalculationConfig
from energyinator.pipeline import run_pipeline
from energyinator.cli import main
from energyinator import calcola_energia_psi4

# Maintain legacy module-level variables for backward compatibility
PSI4_MEMORY = '8 GB'
PSI4_THREADS = 8
PSI4_METHOD = 'b3lyp'
PSI4_BASIS = '6-31G(d)'

if __name__ == "__main__":
    main()
