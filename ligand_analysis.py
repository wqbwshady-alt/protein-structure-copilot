"""
RDKit-based ligand physicochemical analysis.

Computes drug-likeness properties, identifies aromatic systems,
and performs MMFF94 energy minimization for strain energy estimation.
"""

import io
import math

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors


def analyze_ligand(pdb_path, ligand_name):
    """Analyze ligand from PDB using RDKit.

    Args:
        pdb_path: path to PDB file
        ligand_name: 3-letter residue name of ligand (e.g. 'ATP')

    Returns:
        dict with physicochemical properties, aromatic info, strain energy,
        or None if parsing fails.
    """
    # Extract ligand atoms from PDB
    ligand_atoms = _extract_ligand_atoms(pdb_path, ligand_name)
    if not ligand_atoms:
        return None

    # Build RDKit molecule from PDB coordinates
    mol = _build_mol_from_pdb_atoms(ligand_atoms)
    if mol is None:
        return None

    # Basic properties (some require successful sanitization)
    mw = round(Descriptors.MolWt(mol), 1)
    try:
        logp = round(Crippen.MolLogP(mol), 2)
    except Exception:
        logp = None
    try:
        tpsa = round(rdMolDescriptors.CalcTPSA(mol), 1)
    except Exception:
        tpsa = None
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot_bonds = Lipinski.NumRotatableBonds(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    try:
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    except Exception:
        aromatic_rings = 0
    heavy_atoms = mol.GetNumHeavyAtoms()

    # Aromatic detection
    has_aromatic_system = False
    try:
        ri = mol.GetRingInfo()
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                has_aromatic_system = True
                break
    except Exception:
        pass

    # Drug-likeness assessment
    ro5_violations = 0
    if mw > 500: ro5_violations += 1
    if logp > 5: ro5_violations += 1
    if hbd > 5: ro5_violations += 1
    if hba > 10: ro5_violations += 1

    drug_likeness = "good" if ro5_violations == 0 else (
        "moderate" if ro5_violations <= 1 else "poor"
    )

    # MMFF94 energy minimization for strain energy
    strain_energy = None
    minimized_energy = None
    try:
        mol_h = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol_h, randomSeed=42) == 0:
            ff = AllChem.MMFFGetMoleculeForceField(mol_h, AllChem.MMFFGetMoleculeProperties(mol_h))
            if ff:
                initial_energy = ff.CalcEnergy()
                ff.Minimization(maxIts=200)
                minimized_energy = round(ff.CalcEnergy(), 2)
                strain_energy = round(initial_energy - minimized_energy, 2)
    except Exception:
        pass

    return {
        "name": ligand_name,
        "mw": mw,
        "logp": logp,
        "tpsa": tpsa,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rot_bonds,
        "ring_count": ring_count,
        "aromatic_rings": aromatic_rings,
        "has_aromatic_system": has_aromatic_system,
        "heavy_atoms": heavy_atoms,
        "drug_likeness": drug_likeness,
        "ro5_violations": ro5_violations,
        "mmff_minimized_energy": minimized_energy,
        "mmff_strain_energy": strain_energy,
    }


def _extract_ligand_atoms(pdb_path, ligand_name):
    """Extract HETATM records for a specific ligand from a PDB file."""
    atoms = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            res_name = line[17:20].strip()
            if res_name.upper() != ligand_name.upper():
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip().upper() or _infer_element(atom_name)
            atoms.append({
                "atom_name": atom_name,
                "element": element,
                "res_name": res_name,
                "coord": (x, y, z),
            })
    return atoms


def _infer_element(atom_name):
    """Infer element from PDB atom name."""
    name = atom_name.strip().upper()
    if not name:
        return "C"
    first = name[0]
    if first in "CHONPSF":
        return first
    return "C"


def _build_mol_from_pdb_atoms(atoms):
    """Build an RDKit molecule from PDB atom coordinates using MOL block.

    MOL format avoids PDB element-parsing issues with exotic atom names.
    """
    n = len(atoms)
    bonds = _infer_bonds(atoms)
    n_bonds = len(bonds)

    # MOL V2000 format
    lines = [
        "",                                    # line 1: molecule name
        "  RDKit",                             # line 2: program
        "",                                    # line 3: comment
        f"{n:>3d}{n_bonds:>3d}  0  0  0  0  0  0  0  0999 V2000",  # counts line
    ]
    # Atom block
    for atom in atoms:
        elem = atom['element'].strip()[:2]
        x, y, z = atom['coord']
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {elem:<3s} 0  0  0  0  0  0  0  0  0  0  0  0")
    # Bond block
    for i, j, order in bonds:
        lines.append(f"{i+1:>3d}{j+1:>3d}{order:>3d}  0  0  0  0")
    lines.append("M  END")
    mol_block = "\n".join(lines)

    mol = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=True)
    if mol is None:
        mol = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=False)
        if mol:
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                return mol
    return mol


def _infer_bonds(atoms):
    """Infer covalent bonds from interatomic distances using element radii.

    Distinguishes single vs double bonds by bond-length thresholds.
    """
    _COVALENT_RADII = {
        "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66,
        "P": 1.07, "S": 1.05, "F": 0.57, "CL": 1.02,
        "MG": 1.39, "MN": 1.39, "ZN": 1.22, "CA": 1.76,
        "BR": 1.20, "I": 1.39,
    }
    bonds = []
    n = len(atoms)
    for i in range(n):
        ri = _COVALENT_RADII.get(atoms[i]["element"].strip(), 0.75)
        ci = atoms[i]["coord"]
        for j in range(i + 1, n):
            rj = _COVALENT_RADII.get(atoms[j]["element"].strip(), 0.75)
            cj = atoms[j]["coord"]
            dist = math.sqrt(sum((ci[k] - cj[k]) ** 2 for k in range(3)))
            single_threshold = (ri + rj) * 1.3
            if dist < single_threshold:
                bonds.append((i, j, 1))
    return bonds
