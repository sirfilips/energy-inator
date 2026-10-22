"""1D Dihedral Potential Energy Surface (PES) Scanner for conformational analysis."""

import numpy as np
from typing import List, Dict, Any, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, AllChem
from rdkit.Geometry import Point3D

from energyinator.conformers import prepare_molecule
from energyinator.quantum import run_psi4_calculation


def parse_xyz_string(xyz_str: str) -> Tuple[List[str], np.ndarray]:
    """Parses an XYZ coordinate string into symbols and numpy coordinates."""
    lines = xyz_str.strip().splitlines()
    if len(lines) >= 2 and lines[0].strip().isdigit():
        lines = lines[2:]
    symbols = []
    coords = []
    for l in lines:
        parts = l.split()
        if len(parts) >= 4:
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords)


def create_mol_from_smiles_and_xyz(smiles: str, xyz_str: str) -> Chem.Mol:
    """Creates a hydrated RDKit molecule with 3D conformer coordinates set from an XYZ block."""
    mol, _, _ = prepare_molecule(smiles)
    symbols, coords = parse_xyz_string(xyz_str)
    
    if len(symbols) != mol.GetNumAtoms():
        # Fallback: embed directly
        AllChem.EmbedMolecule(mol)
        return mol

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, pt in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(float(pt[0]), float(pt[1]), float(pt[2])))
        
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


def run_dihedral_scan(
    smiles: str,
    xyz_str: str,
    atom_indices: Tuple[int, int, int, int],
    step_deg: int = 15,
    engine: str = "mmff94",
    psi4_method: str = "b3lyp",
    psi4_basis: str = "sto-3g",
    charge: int = 0,
    multiplicity: int = 1
) -> Dict[str, Any]:
    """Performs a 1D torsional potential energy surface (PES) scan over 0° to 360°.
    
    Args:
        smiles: Target molecule SMILES string.
        xyz_str: Reference 3D optimized coordinates.
        atom_indices: Tuple of 4 atom indices (a1, a2, a3, a4) defining the dihedral angle.
        step_deg: Step size in degrees (e.g. 15° or 30°).
        engine: "mmff94", "uff", or "psi4".
        psi4_method: DFT method if engine is "psi4".
        psi4_basis: Basis set if engine is "psi4".
        charge: Formal charge.
        multiplicity: Spin multiplicity.

    Returns:
        Dictionary containing scan angles, absolute energies, relative energies (kcal/mol),
        XYZ coordinate strings for each point, and summary metrics.
    """
    a1, a2, a3, a4 = atom_indices
    mol = create_mol_from_smiles_and_xyz(smiles, xyz_str)
    conf = mol.GetConformer(0)

    # Prepare force field properties if classical
    mp = None
    if engine.lower() == "mmff94":
        mp = AllChem.MMFFGetMoleculeProperties(mol)

    angles = list(range(0, 361, step_deg))
    raw_energies = []
    xyz_list = []

    for deg in angles:
        rdMolTransforms.SetDihedralDeg(conf, a1, a2, a3, a4, float(deg))
        current_xyz = Chem.MolToXYZBlock(mol)
        xyz_list.append(current_xyz)

        if engine.lower() == "mmff94" and mp is not None:
            ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=0)
            if ff is not None:
                energy_kcal = ff.CalcEnergy()
            else:
                energy_kcal = 0.0
            raw_energies.append(energy_kcal)

        elif engine.lower() == "uff":
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
            if ff is not None:
                energy_kcal = ff.CalcEnergy()
            else:
                energy_kcal = 0.0
            raw_energies.append(energy_kcal)

        elif engine.lower() == "psi4":
            calc_out = run_psi4_calculation(
                xyz_block=current_xyz,
                conf_id=deg,
                charge=charge,
                multiplicity=multiplicity,
                method=psi4_method,
                basis=psi4_basis,
                sp_only=True,
                smiles_label=f"PES_Scan_{deg}deg"
            )
            if calc_out is not None:
                raw_energies.append(calc_out["energy_kcalmol"])
            else:
                raw_energies.append(np.nan)
        else:
            raw_energies.append(0.0)

    raw_energies = np.array(raw_energies)
    valid_mask = ~np.isnan(raw_energies)
    if not np.any(valid_mask):
        raise RuntimeError("Impossibile calcolare le energie lungo la scansione del diedro.")

    min_val = np.min(raw_energies[valid_mask])
    delta_kcalmol = [float(e - min_val) if not np.isnan(e) else None for e in raw_energies]
    max_barrier = float(np.max(raw_energies[valid_mask]) - min_val)

    min_idx = int(np.argmin(raw_energies[valid_mask]))
    min_angle = angles[min_idx]

    return {
        "angles": angles,
        "energies_raw": [float(e) if not np.isnan(e) else None for e in raw_energies],
        "delta_kcalmol": delta_kcalmol,
        "xyz_list": xyz_list,
        "min_angle": min_angle,
        "max_barrier_kcalmol": max_barrier,
        "atom_names": (a1, a2, a3, a4),
        "engine": engine
    }

