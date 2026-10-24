"""Geometric utilities, coordinate transformations, and XYZ handling for molecular modeling."""

import re
import numpy as np
from typing import List, Tuple, Optional, Dict
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, AllChem
from rdkit.Geometry import Point3D


def parse_xyz(xyz_str: str) -> Tuple[List[str], np.ndarray]:
    """Parses an XYZ coordinate string into atomic symbols and an Nx3 NumPy coordinate array.
    
    Correctly handles standard 2-line XYZ headers (atom count and comment line)
    as well as raw coordinate lines.
    """
    lines = xyz_str.strip().splitlines()
    if len(lines) >= 2 and lines[0].strip().isdigit():
        lines = lines[2:]
        
    symbols: List[str] = []
    coords: List[List[float]] = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            
    return symbols, np.array(coords, dtype=np.float64)


def parse_xyz_coords(xyz_str: str) -> List[Tuple[float, float, float]]:
    """Convenience helper returning coordinates as a list of 3-tuples (x, y, z)."""
    _, coords_np = parse_xyz(xyz_str)
    return [tuple(pt) for pt in coords_np]


def format_xyz(symbols: List[str], coords: np.ndarray, comment: str = "") -> str:
    """Formats atomic symbols and coordinates into a standard XYZ formatted string."""
    natoms = len(symbols)
    lines = [str(natoms), comment]
    for sym, (x, y, z) in zip(symbols, coords):
        lines.append(f"{sym:<2} {x:14.8f} {y:14.8f} {z:14.8f}")
    return "\n".join(lines)


def calc_distance(p1: np.ndarray, p2: np.ndarray) -> float:
    """Calculates the Euclidean distance in Angstroms between two points."""
    return float(np.linalg.norm(np.asarray(p1) - np.asarray(p2)))


def calc_angle(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """Calculates the valence angle in degrees formed by p1-p2-p3 (vertex at p2)."""
    v1 = np.asarray(p1) - np.asarray(p2)
    v2 = np.asarray(p3) - np.asarray(p2)
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom < 1e-9:
        return 0.0
    cosine = np.dot(v1, v2) / denom
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def calc_dihedral(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """Calculates the dihedral (torsional) angle in degrees defined by four points p1-p2-p3-p4."""
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)
    p3 = np.asarray(p3)
    p4 = np.asarray(p4)
    
    b0 = -1.0 * (p2 - p1)
    b1 = p3 - p2
    b2 = p4 - p3
    norm_b1 = np.linalg.norm(b1)
    if norm_b1 < 1e-9:
        return 0.0
    b1 /= norm_b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def kabsch_align_xyz(ref_xyz_str: str, target_xyz_str: str) -> str:
    """Aligns and overlays target coordinates onto reference coordinates using the Kabsch algorithm."""
    symbols_ref, P = parse_xyz(ref_xyz_str)
    symbols_tgt, Q = parse_xyz(target_xyz_str)

    if len(P) != len(Q) or len(P) == 0:
        return target_xyz_str

    p_center = np.mean(P, axis=0)
    q_center = np.mean(Q, axis=0)

    P_c = P - p_center
    Q_c = Q - q_center

    H = np.dot(Q_c.T, P_c)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = np.dot(Vt.T, U.T)

    Q_aligned = np.dot(Q_c, R.T) + p_center
    return format_xyz(symbols_tgt, Q_aligned, comment="Aligned conformer")


def mol_from_smiles_and_xyz(smiles: str, xyz_str: str) -> Chem.Mol:
    """Builds an RDKit Chem.Mol with explicit hydrogens and conformer 0 populated from an XYZ string."""
    mol_base = Chem.MolFromSmiles(smiles)
    if mol_base is None:
        raise ValueError(f"SMILES non valido: '{smiles}'")
    mol = Chem.AddHs(mol_base)
    
    symbols, coords = parse_xyz(xyz_str)
    if len(symbols) != mol.GetNumAtoms():
        AllChem.EmbedMolecule(mol)
        return mol

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, pt in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(float(pt[0]), float(pt[1]), float(pt[2])))

    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


def find_rotatable_dihedrals(mol: Chem.Mol) -> List[Dict[str, any]]:
    """Detects rotatable bonds and suggests 4-atom dihedral definitions (a1-a2-a3-a4).
    
    Returns a list of dicts with:
      - 'atoms': (a1, a2, a3, a4)
      - 'bond': (a2, a3)
      - 'label': human-readable description (e.g. 'C1 - C2')
    """
    # Smarts query for single bonds between non-terminal heavy atoms or key polar rotors
    pattern = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")
    if not pattern:
        return []
    
    matches = mol.GetSubstructMatches(pattern)
    dihedrals = []
    
    for (a2, a3) in matches:
        atom2 = mol.GetAtomWithIdx(a2)
        atom3 = mol.GetAtomWithIdx(a3)
        
        # Pick a neighbor for a2 distinct from a3
        a1 = None
        for n in atom2.GetNeighbors():
            idx = n.GetIdx()
            if idx != a3:
                # Prefer heavy atoms over H if possible
                if a1 is None or (n.GetAtomicNum() > 1 and mol.GetAtomWithIdx(a1).GetAtomicNum() == 1):
                    a1 = idx
                    
        # Pick a neighbor for a3 distinct from a2
        a4 = None
        for n in atom3.GetNeighbors():
            idx = n.GetIdx()
            if idx != a2:
                if a4 is None or (n.GetAtomicNum() > 1 and mol.GetAtomWithIdx(a4).GetAtomicNum() == 1):
                    a4 = idx
                    
        if a1 is not None and a4 is not None and len({a1, a2, a3, a4}) == 4:
            sym1 = mol.GetAtomWithIdx(a1).GetSymbol()
            sym2 = atom2.GetSymbol()
            sym3 = atom3.GetSymbol()
            sym4 = mol.GetAtomWithIdx(a4).GetSymbol()
            label = f"{sym2}{a2} - {sym3}{a3} ({sym1}{a1}-{sym2}{a2}-{sym3}{a3}-{sym4}{a4})"
            dihedrals.append({
                "atoms": (a1, a2, a3, a4),
                "bond": (a2, a3),
                "label": label
            })
            
    return dihedrals

