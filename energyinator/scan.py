"""1D and 2D Dihedral Potential Energy Surface (PES) Scanner for conformational analysis."""

import numpy as np
from typing import List, Dict, Any, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, AllChem

from energyinator.geometry import (
    parse_xyz,
    mol_from_smiles_and_xyz,
    find_rotatable_dihedrals
)
from energyinator.quantum import run_psi4_calculation


def run_dihedral_scan(
    smiles: str,
    xyz_str: str,
    atom_indices: Tuple[int, int, int, int],
    step_deg: int = 15,
    engine: str = "mmff94",
    relaxed: bool = True,
    psi4_method: str = "b3lyp",
    psi4_basis: str = "sto-3g",
    charge: int = 0,
    multiplicity: int = 1,
    solvent: Optional[str] = None
) -> Dict[str, Any]:
    """Performs a 1D torsional potential energy surface (PES) scan over 0° to 360°.
    
    Supports both relaxed (constrained optimization) and rigid torsional scans.
    Highly optimized: pre-computes force-field parameters outside the loop.
    """
    a1, a2, a3, a4 = atom_indices
    base_mol = mol_from_smiles_and_xyz(smiles, xyz_str)

    angles = list(range(0, 361, step_deg))
    raw_energies = []
    xyz_list = []

    eng = engine.lower()
    
    # Pre-calculate MMFF molecule properties once if applicable
    mmff_props = None
    if eng == "mmff94":
        mmff_props = AllChem.MMFFGetMoleculeProperties(base_mol)

    for deg in angles:
        mol_step = Chem.Mol(base_mol)
        conf = mol_step.GetConformer(0)
        rdMolTransforms.SetDihedralDeg(conf, a1, a2, a3, a4, float(deg))

        if eng == "mmff94":
            if mmff_props is not None:
                ff = AllChem.MMFFGetMoleculeForceField(mol_step, mmff_props, confId=0)
                if ff is not None:
                    if relaxed:
                        ff.MMFFAddTorsionConstraint(a1, a2, a3, a4, False, float(deg), float(deg), 1000.0)
                        ff.Minimize(maxIts=250)
                    energy_kcal = ff.CalcEnergy()
                else:
                    energy_kcal = np.nan
            else:
                energy_kcal = np.nan
            raw_energies.append(energy_kcal)
            xyz_list.append(Chem.MolToXYZBlock(mol_step))

        elif eng == "uff":
            ff = AllChem.UFFGetMoleculeForceField(mol_step, confId=0)
            if ff is not None:
                if relaxed:
                    ff.UFFAddTorsionConstraint(a1, a2, a3, a4, False, float(deg), float(deg), 1000.0)
                    ff.Minimize(maxIts=250)
                energy_kcal = ff.CalcEnergy()
            else:
                energy_kcal = np.nan
            raw_energies.append(energy_kcal)
            xyz_list.append(Chem.MolToXYZBlock(mol_step))

        elif eng == "psi4":
            current_xyz = Chem.MolToXYZBlock(mol_step)
            xyz_list.append(current_xyz)
            calc_out = run_psi4_calculation(
                xyz_block=current_xyz,
                conf_id=deg,
                charge=charge,
                multiplicity=multiplicity,
                method=psi4_method,
                basis=psi4_basis,
                sp_only=True,
                solvent=solvent,
                calc_frequencies=False,
                smiles_label=f"PES_Scan_{deg}deg"
            )
            if calc_out is not None:
                raw_energies.append(calc_out["energy_kcalmol"])
            else:
                raw_energies.append(np.nan)
        else:
            raw_energies.append(0.0)
            xyz_list.append(Chem.MolToXYZBlock(mol_step))

    raw_energies = np.array(raw_energies)
    valid_mask = ~np.isnan(raw_energies)
    if not np.any(valid_mask):
        raise RuntimeError("Impossibile calcolare le energie lungo la scansione del diedro.")

    valid_indices = np.where(valid_mask)[0]
    min_sub_idx = int(np.argmin(raw_energies[valid_mask]))
    min_idx = int(valid_indices[min_sub_idx])

    min_val = float(raw_energies[min_idx])
    delta_kcalmol = [float(e - min_val) if not np.isnan(e) else None for e in raw_energies]
    max_barrier = float(np.max(raw_energies[valid_mask]) - min_val)
    min_angle = angles[min_idx]

    return {
        "angles": angles,
        "energies_raw": [float(e) if not np.isnan(e) else None for e in raw_energies],
        "delta_kcalmol": delta_kcalmol,
        "xyz_list": xyz_list,
        "min_angle": min_angle,
        "max_barrier_kcalmol": max_barrier,
        "atom_names": (a1, a2, a3, a4),
        "engine": engine,
        "relaxed": relaxed
    }


def run_2d_dihedral_scan(
    smiles: str,
    xyz_str: str,
    dihedral_1: Tuple[int, int, int, int],
    dihedral_2: Tuple[int, int, int, int],
    step_deg: int = 30,
    engine: str = "mmff94",
    relaxed: bool = True
) -> Dict[str, Any]:
    """Performs a 2D torsional potential energy surface scan across two coupled dihedral angles.
    
    Highly optimized for performance: pre-computes force field parameters once outside the double loop.
    """
    a1, a2, a3, a4 = dihedral_1
    b1, b2, b3, b4 = dihedral_2

    base_mol = mol_from_smiles_and_xyz(smiles, xyz_str)
    angles_1 = list(range(0, 361, step_deg))
    angles_2 = list(range(0, 361, step_deg))

    eng = engine.lower()
    mmff_props = None
    if eng == "mmff94":
        mmff_props = AllChem.MMFFGetMoleculeProperties(base_mol)

    z_matrix = []
    min_e = float("inf")
    min_coords = (0, 0)

    for deg1 in angles_1:
        row = []
        for deg2 in angles_2:
            mol_step = Chem.Mol(base_mol)
            conf = mol_step.GetConformer(0)
            rdMolTransforms.SetDihedralDeg(conf, a1, a2, a3, a4, float(deg1))
            rdMolTransforms.SetDihedralDeg(conf, b1, b2, b3, b4, float(deg2))

            if eng == "mmff94":
                if mmff_props is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(mol_step, mmff_props, confId=0)
                    if ff is not None:
                        if relaxed:
                            ff.MMFFAddTorsionConstraint(a1, a2, a3, a4, False, float(deg1), float(deg1), 1000.0)
                            ff.MMFFAddTorsionConstraint(b1, b2, b3, b4, False, float(deg2), float(deg2), 1000.0)
                            ff.Minimize(maxIts=200)
                        e_val = ff.CalcEnergy()
                    else:
                        e_val = np.nan
                else:
                    e_val = np.nan
            elif eng == "uff":
                ff = AllChem.UFFGetMoleculeForceField(mol_step, confId=0)
                if ff is not None:
                    if relaxed:
                        ff.UFFAddTorsionConstraint(a1, a2, a3, a4, False, float(deg1), float(deg1), 1000.0)
                        ff.UFFAddTorsionConstraint(b1, b2, b3, b4, False, float(deg2), float(deg2), 1000.0)
                        ff.Minimize(maxIts=200)
                    e_val = ff.CalcEnergy()
                else:
                    e_val = np.nan
            else:
                e_val = 0.0

            if not np.isnan(e_val) and e_val < min_e:
                min_e = e_val
                min_coords = (deg1, deg2)
            row.append(e_val)
        z_matrix.append(row)

    # Compute relative delta energies
    delta_z = []
    max_val = 0.0
    for row in z_matrix:
        d_row = []
        for val in row:
            if np.isnan(val) or np.isinf(min_e):
                d_row.append(None)
            else:
                diff = float(val - min_e)
                d_row.append(diff)
                if diff > max_val:
                    max_val = diff
        delta_z.append(d_row)

    return {
        "angles_1": angles_1,
        "angles_2": angles_2,
        "delta_z": delta_z,
        "min_deg1": min_coords[0],
        "min_deg2": min_coords[1],
        "max_barrier_kcalmol": max_val,
        "dihedral_1": dihedral_1,
        "dihedral_2": dihedral_2,
        "engine": engine,
        "relaxed": relaxed
    }
