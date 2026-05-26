"""
Water-mediated protein-ligand contact detection.

Identifies water molecules (HOH) that bridge between protein residues
and ligand atoms — a common feature in binding interfaces that
contributes to binding affinity and specificity.

Criteria: water O atom within 3.5A of both a protein atom and a ligand atom.
"""

import math

_WATER_BRIDGE_DIST = 3.5  # A — water must be within this distance of both partners


def detect_water_bridges(protein_atoms, ligand_atoms):
    """Detect water molecules bridging protein and ligand.

    Args:
        protein_atoms: list of atom dicts (includes HOH HETATM records)
        ligand_atoms:  list of ligand HETATM dicts

    Returns:
        dict with bridges list and summary counts
    """
    if not protein_atoms or not ligand_atoms:
        return _empty_result()

    # Separate water molecules from protein atoms
    water_atoms = []
    actual_protein = []
    for a in protein_atoms:
        if a.get("res_name", "").strip() in ("HOH", "WAT"):
            water_atoms.append(a)
        elif a.get("atom_type") == "ATOM":
            actual_protein.append(a)

    if not water_atoms:
        return _empty_result()

    bridges = []
    seen_waters = set()

    for water in water_atoms:
        w_coord = water.get("coord")
        if not w_coord:
            continue

        water_id = f"{water['chain_id']}:{water['res_name']}{water['res_id']}"

        # Find protein contacts for this water
        protein_contacts = []
        for patom in actual_protein:
            p_coord = patom.get("coord")
            if not p_coord:
                continue
            dist = _distance(w_coord, p_coord)
            if dist <= _WATER_BRIDGE_DIST:
                protein_contacts.append({
                    "residue_key": f"{patom['chain_id']}:{patom['res_name']}{patom['res_id']}",
                    "atom_name": patom["atom_name"],
                    "distance": round(dist, 2),
                })

        # Find ligand contacts for this water
        ligand_contacts = []
        for latom in ligand_atoms:
            l_coord = latom.get("coord")
            if not l_coord:
                continue
            dist = _distance(w_coord, l_coord)
            if dist <= _WATER_BRIDGE_DIST:
                ligand_contacts.append({
                    "atom_name": latom.get("atom_name", ""),
                    "res_name": latom.get("res_name", ""),
                    "distance": round(dist, 2),
                })

        # Bridge = water contacts both protein AND ligand
        if protein_contacts and ligand_contacts:
            for pc in protein_contacts:
                bridges.append({
                    "water_id": water_id,
                    "protein_residue": pc["residue_key"],
                    "protein_atom": pc["atom_name"],
                    "protein_water_dist": pc["distance"],
                    "ligand_atom": ligand_contacts[0]["atom_name"],
                    "ligand_water_dist": ligand_contacts[0]["distance"],
                    "category": "water_bridge",
                })
            seen_waters.add(water_id)

    # Count unique residues involved
    protein_residues = set(b["protein_residue"] for b in bridges)

    return {
        "bridges": bridges,
        "summary": {
            "total": len(bridges),
            "water_molecules_involved": len(seen_waters),
            "protein_residues_bridged": len(protein_residues),
        },
    }


def _distance(c1, c2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


def _empty_result():
    return {
        "bridges": [],
        "summary": {
            "total": 0,
            "water_molecules_involved": 0,
            "protein_residues_bridged": 0,
        },
    }
