"""
Per-residue interaction energy decomposition.

Computes approximate Lennard-Jones (vdW) and Coulomb (electrostatic) energies
between protein pocket residues and the ligand using simplified AMBER-like
force field parameters. Pure Python — no external dependencies.

Physics:
  V_LJ    = 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]
  V_Coul  = 332.0637 * q1 * q2 / (epsilon_r * r)    [kcal/mol]
  epsilon_r = 4*r  (distance-dependent dielectric for protein interior)
"""

import math

# ---- Physical constants ----
_COULOMB_PREFACTOR = 332.0637  # kcal·A / (mol·e^2)
_COULOMB_SCALE = 0.12  # scale factor for coarse partial charges (calibrated for relative ranking)
_CUTOFF = 6.0  # Angstrom — interaction cutoff
_DIELECTRIC_SCALE = 4.0  # ε(r) = scale * r (Sigworth-Honig model)

# ---- Simplified AMBER-like atom type parameters ----
# (sigma in Angstrom, epsilon in kcal/mol, charge in e)

_ATOM_PARAMS = {
    # Carbon types
    "C":   (1.908, 0.086,  0.00),   # aliphatic C (CH, CH2, CH3)
    "CA":  (1.908, 0.086,  0.00),   # aromatic C (PHE/TYR/TRP/HIS ring)
    "CO":  (1.908, 0.086,  0.50),   # carbonyl/amide C (backbone C=O)
    "CC":  (1.908, 0.086,  0.70),   # carboxylate C (ASP/GLU side chain)
    # Nitrogen types
    "N":   (1.824, 0.170, -0.50),   # amide N (backbone NH, ASN/GLN side chain)
    "N3":  (1.824, 0.170,  1.00),   # ammonium N (LYS NZ, NH3+)
    "NR":  (1.824, 0.170,  0.00),   # aromatic N (HIS ring, neutral)
    "NG":  (1.824, 0.170,  0.80),   # guanidinium N (ARG NH1/NH2, positive)
    # Oxygen types
    "O":   (1.661, 0.210, -0.50),   # carbonyl O (backbone C=O, ASN/GLN side chain)
    "O2":  (1.661, 0.210, -0.80),   # carboxylate O (ASP/GLU side chain)
    "OH":  (1.721, 0.210, -0.60),   # hydroxyl O (SER/THR/TYR)
    "OW":  (1.768, 0.152, -0.83),   # water O (TIP3P-like, for HOH)
    # Sulfur types
    "S":   (2.000, 0.250, -0.15),   # thioether S (MET)
    "SH":  (2.000, 0.250, -0.25),   # thiol S (CYS)
    # Metal ions
    "MG":  (1.180, 0.015,  2.00),   # Mg2+
    "MN":  (1.340, 0.015,  2.00),   # Mn2+
    "ZN":  (1.100, 0.015,  2.00),   # Zn2+
    "CA":  (1.780, 0.015,  2.00),   # Ca2+
}

# ---- Atom type classification from element + atom name + residue context ----

def _classify_atom(element, atom_name, res_name):
    """Determine simplified force-field atom type from PDB atom name and residue."""
    name = atom_name.strip().upper()

    if element == "C":
        # Aromatic carbons
        if name in ("CG", "CD1", "CD2", "CE1", "CE2", "CZ",
                     "CD", "CE", "CG1", "CG2"):
            if res_name in ("PHE", "TYR", "TRP", "HIS"):
                return "CA"
        # PHE/TYR ring carbons specifically
        if name in ("CD1", "CD2", "CE1", "CE2", "CG", "CZ"):
            return "CA"
        # Backbone carbonyl C
        if name == "C":
            return "CO"
        # ASP/GLU side chain carboxylate C
        if name == "CD" and res_name in ("ASP", "GLU"):
            return "CC"
        return "C"

    elif element == "N":
        # LYS side chain amine
        if name == "NZ":
            return "N3"
        # ARG guanidinium
        if name in ("NH1", "NH2", "NE"):
            return "NG"
        # HIS ring
        if name in ("ND1", "NE2") and res_name == "HIS":
            return "NR"
        # Backbone N or amide N
        return "N"

    elif element == "O":
        # ASP/GLU carboxylate O
        if name in ("OD1", "OD2", "OE1", "OE2"):
            return "O2"
        # Hydroxyl O
        if name in ("OG", "OG1", "OH"):
            return "OH"
        # Water
        if name == "O" and res_name in ("HOH", "WAT"):
            return "OW"
        return "O"

    elif element == "S":
        if name == "SG" and res_name == "CYS":
            return "SH"
        return "S"

    elif element == "MG":
        return "MG"
    elif element == "MN":
        return "MN"
    elif element == "ZN":
        return "ZN"
    elif element == "CA":
        return "CA"

    # Fallback — use generic parameters
    return element


def compute_interaction_energy(protein_atoms, ligand_atoms, cutoff=_CUTOFF):
    """Compute per-residue LJ + Coulomb interaction energy with ligand.

    Args:
        protein_atoms: list of ATOM dicts (atom_type=="ATOM") from parse_pdb_atoms
        ligand_atoms:  list of HETATM dicts from parse_pdb_atoms
        cutoff:        distance cutoff in Angstrom (default 8.0)

    Returns:
        dict with:
            per_residue: {residue_key: {vdw, coulomb, total, atom_pairs}}
            total_vdw: float (kcal/mol)
            total_coulomb: float (kcal/mol)
            total_energy: float (kcal/mol)
    """
    if not protein_atoms or not ligand_atoms:
        return _empty_result()

    per_residue = {}
    total_vdw = 0.0
    total_coul = 0.0

    for patom in protein_atoms:
        p_coord = patom.get("coord")
        if not p_coord:
            continue
        p_elem = patom.get("element", "")
        p_type = _classify_atom(p_elem, patom.get("atom_name", ""),
                                patom.get("res_name", ""))
        p_params = _ATOM_PARAMS.get(p_type)
        if not p_params:
            continue

        residue_key = f"{patom['chain_id']}:{patom['res_name']}{patom['res_id']}"
        if residue_key not in per_residue:
            per_residue[residue_key] = {
                "vdw": 0.0, "coulomb": 0.0, "total": 0.0,
                "atom_pairs": 0, "res_name": patom["res_name"],
                "chain_id": patom["chain_id"], "res_id": patom["res_id"],
            }

        p_sigma, p_eps, p_q = p_params

        for latom in ligand_atoms:
            l_coord = latom.get("coord")
            if not l_coord:
                continue
            l_elem = latom.get("element", "")
            l_type = _classify_atom(l_elem, latom.get("atom_name", ""),
                                    latom.get("res_name", ""))
            l_params = _ATOM_PARAMS.get(l_type)
            if not l_params:
                continue

            # Distance
            dx = p_coord[0] - l_coord[0]
            dy = p_coord[1] - l_coord[1]
            dz = p_coord[2] - l_coord[2]
            r2 = dx*dx + dy*dy + dz*dz
            if r2 < 0.01 or r2 > cutoff * cutoff:
                continue

            r = math.sqrt(r2)
            l_sigma, l_eps, l_q = l_params

            # Lorentz-Berthelot combining rules
            sigma = 0.5 * (p_sigma + l_sigma)
            epsilon = math.sqrt(p_eps * l_eps)

            # LJ 12-6
            sr = sigma / r
            sr6 = sr ** 6
            vdw = 4.0 * epsilon * (sr6 * sr6 - sr6)

            # Coulomb with distance-dependent dielectric + coarse-charge scaling
            eps_r = _DIELECTRIC_SCALE * r
            coulomb = _COULOMB_SCALE * _COULOMB_PREFACTOR * p_q * l_q / eps_r

            per_residue[residue_key]["vdw"] += vdw
            per_residue[residue_key]["coulomb"] += coulomb
            per_residue[residue_key]["atom_pairs"] += 1
            total_vdw += vdw
            total_coul += coulomb

    # Compute total and round
    for key in per_residue:
        per_residue[key]["vdw"] = round(per_residue[key]["vdw"], 2)
        per_residue[key]["coulomb"] = round(per_residue[key]["coulomb"], 2)
        per_residue[key]["total"] = round(
            per_residue[key]["vdw"] + per_residue[key]["coulomb"], 2
        )

    return {
        "per_residue": per_residue,
        "total_vdw": round(total_vdw, 2),
        "total_coulomb": round(total_coul, 2),
        "total_energy": round(total_vdw + total_coul, 2),
    }


def _empty_result():
    return {
        "per_residue": {},
        "total_vdw": 0.0,
        "total_coulomb": 0.0,
        "total_energy": 0.0,
    }
