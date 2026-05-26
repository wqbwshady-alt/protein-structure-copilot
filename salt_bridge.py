"""
Salt bridge detection — geometric identification of protein-protein
and protein-ligand salt bridges.

A salt bridge is a close-range electrostatic interaction between
oppositely charged side chains. Criteria:
  - Acidic: ASP (OD1/OD2), GLU (OE1/OE2) — carboxylate O
  - Basic:  LYS (NZ), ARG (NH1/NH2), HIS (ND1/NE2) — ammonium/guanidinium N
  - Distance between any acidic O and basic N < 4.0 A
"""

import math

_SALT_BRIDGE_DIST_MAX = 4.0  # A — close-range electrostatic

_ACIDIC_RESIDUES = {"ASP", "GLU"}
_BASIC_RESIDUES = {"LYS", "ARG", "HIS"}

_ACIDIC_ATOMS = {"OD1", "OD2", "OE1", "OE2"}
_BASIC_ATOMS = {
    "LYS": {"NZ"},
    "ARG": {"NH1", "NH2", "NE"},
    "HIS": {"ND1", "NE2"},
}


def detect_salt_bridges(protein_atoms):
    """Detect salt bridges among protein residues.

    Args:
        protein_atoms: list of ATOM dicts from parse_pdb_atoms

    Returns:
        dict with bridges list and summary counts
    """
    if not protein_atoms:
        return _empty_result()

    # Collect acidic and basic atoms
    acidic_atoms = []
    basic_atoms = []

    for atom in protein_atoms:
        if atom["atom_type"] != "ATOM":
            continue
        name = atom["atom_name"].strip()
        res = atom["res_name"]

        if res in _ACIDIC_RESIDUES and name in _ACIDIC_ATOMS:
            acidic_atoms.append({
                "residue_key": f"{atom['chain_id']}:{atom['res_name']}{atom['res_id']}",
                "res_name": res,
                "chain_id": atom["chain_id"],
                "res_id": atom["res_id"],
                "atom_name": name,
                "coord": atom["coord"],
                "charge": "negative",
            })
        elif res in _BASIC_RESIDUES and name in _BASIC_ATOMS.get(res, set()):
            basic_atoms.append({
                "residue_key": f"{atom['chain_id']}:{atom['res_name']}{atom['res_id']}",
                "res_name": res,
                "chain_id": atom["chain_id"],
                "res_id": atom["res_id"],
                "atom_name": name,
                "coord": atom["coord"],
                "charge": "positive",
            })

    # Find close-range pairs
    bridges = []
    seen_pairs = set()

    for acid in acidic_atoms:
        for base in basic_atoms:
            # Skip same residue
            if acid["residue_key"] == base["residue_key"]:
                continue

            pair = tuple(sorted([acid["residue_key"], base["residue_key"]]))
            if pair in seen_pairs:
                continue

            dist = _distance(acid["coord"], base["coord"])
            if dist <= _SALT_BRIDGE_DIST_MAX:
                seen_pairs.add(pair)
                bridges.append({
                    "acidic_residue": acid["residue_key"],
                    "basic_residue": base["residue_key"],
                    "acidic_atom": acid["atom_name"],
                    "basic_atom": base["atom_name"],
                    "acidic_res": acid["res_name"],
                    "basic_res": base["res_name"],
                    "distance": round(dist, 2),
                    "strength": "strong" if dist < 3.2 else "moderate" if dist < 3.6 else "weak",
                })

    # Count unique residues involved
    acidic_involved = set(b["acidic_residue"] for b in bridges)
    basic_involved = set(b["basic_residue"] for b in bridges)

    return {
        "bridges": bridges,
        "summary": {
            "total": len(bridges),
            "strong": sum(1 for b in bridges if b["strength"] == "strong"),
            "moderate": sum(1 for b in bridges if b["strength"] == "moderate"),
            "weak": sum(1 for b in bridges if b["strength"] == "weak"),
            "acidic_residues_involved": len(acidic_involved),
            "basic_residues_involved": len(basic_involved),
        },
    }


def _distance(c1, c2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


def _empty_result():
    return {
        "bridges": [],
        "summary": {
            "total": 0, "strong": 0, "moderate": 0, "weak": 0,
            "acidic_residues_involved": 0, "basic_residues_involved": 0,
        },
    }
