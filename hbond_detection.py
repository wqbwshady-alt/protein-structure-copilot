"""
Hydrogen bond geometry detection using Baker-Hubbard / IUPAC criteria.

PDB files typically lack hydrogen atoms, so H positions are estimated
from heavy-atom geometry (idealized bond lengths and angles).

Criteria:
  - H···A distance < 2.5 A
  - D···A distance < 3.5 A
  - D-H···A angle > 120 deg
"""

import math

# ---- Thresholds ----
_H_A_DIST_MAX = 2.5     # H···A distance (A)
_D_A_DIST_MAX = 3.5     # D···A distance (A)
_ANGLE_MIN = 120.0      # D-H···A angle (degrees)
_ANGLE_MODERATE = 90.0  # Minimum for "possible H-bond"


def detect_hbonds(protein_atoms, ligand_atoms, contact_residues):
    """Detect hydrogen bonds among protein residues and protein-ligand.

    Args:
        protein_atoms: list of ATOM dicts from parse_pdb_atoms
        ligand_atoms:  list of HETATM dicts from parse_pdb_atoms
        contact_residues: {(chain, res_name, res_id): atom_dict}

    Returns:
        dict with:
            hbonds: list of {type, donor_key, acceptor_key, donor_atom,
                            acceptor_atom, d_a_dist, h_a_dist, angle, category}
            protein_protein: list of protein-protein H-bonds
            protein_ligand: list of protein-ligand H-bonds
            summary: {total, validated, possible, protein_protein, protein_ligand}
    """
    if not protein_atoms:
        return _empty_result()

    # Build per-residue atom lookup
    residue_atoms = {}
    for atom in protein_atoms:
        rkey = (atom["chain_id"], atom["res_name"], atom["res_id"])
        if rkey not in residue_atoms:
            residue_atoms[rkey] = []
        residue_atoms[rkey].append(atom)

    # Identify donors and acceptors
    donors = _find_donors(protein_atoms, residue_atoms)
    acceptors = _find_acceptors(protein_atoms)

    # Also add ligand acceptors
    ligand_acceptors = _find_ligand_acceptors(ligand_atoms)
    all_acceptors = acceptors + ligand_acceptors

    all_hbonds = []

    # Protein-protein (or protein-ligand) H-bonds
    for donor in donors:
        d_coord = donor["coord"]
        h_coord = donor["h_coord"]
        donor_rkey = donor["residue_key"]

        for acceptor in all_acceptors:
            a_coord = acceptor["coord"]
            acceptor_rkey = acceptor["residue_key"]

            # Skip self
            if donor_rkey == acceptor_rkey:
                continue

            # D···A distance
            d_a_dist = _distance(d_coord, a_coord)
            if d_a_dist > _D_A_DIST_MAX:
                continue

            # H···A distance
            h_a_dist = _distance(h_coord, a_coord)
            if h_a_dist > _H_A_DIST_MAX:
                continue

            # D-H···A angle
            angle = _angle_dha(d_coord, h_coord, a_coord)

            if angle >= _ANGLE_MIN:
                category = "validated"
            elif angle >= _ANGLE_MODERATE:
                category = "possible"
            else:
                continue

            a_source = "ligand" if "ligand" in acceptor else "protein"
            hb_type = f"protein_{a_source}"

            all_hbonds.append({
                "type": hb_type,
                "donor_key": donor_rkey,
                "acceptor_key": acceptor_rkey,
                "donor_atom": donor["atom_name"],
                "acceptor_atom": acceptor["atom_name"],
                "donor_res_name": donor["res_name"],
                "acceptor_res_name": acceptor["res_name"],
                "d_a_dist": round(d_a_dist, 2),
                "h_a_dist": round(h_a_dist, 2),
                "angle": round(angle, 1),
                "category": category,
            })

    # Also detect intra-protein H-bonds (backbone-backbone, sidechain-backbone)
    pp_hbonds = [h for h in all_hbonds if h["type"] == "protein_protein"]
    pl_hbonds = [h for h in all_hbonds if h["type"] == "protein_ligand"]
    validated = [h for h in all_hbonds if h["category"] == "validated"]
    possible = [h for h in all_hbonds if h["category"] == "possible"]

    return {
        "hbonds": all_hbonds,
        "protein_protein": pp_hbonds,
        "protein_ligand": pl_hbonds,
        "summary": {
            "total": len(all_hbonds),
            "validated": len(validated),
            "possible": len(possible),
            "protein_protein": len(pp_hbonds),
            "protein_ligand": len(pl_hbonds),
        },
    }


# ---- Donor/acceptor identification ----

def _find_donors(protein_atoms, residue_atoms):
    """Find H-bond donor atoms with estimated H positions."""
    donors = []
    for atom in protein_atoms:
        rkey = (atom["chain_id"], atom["res_name"], atom["res_id"])
        h_coord = _estimate_h_position(atom, residue_atoms.get(rkey, []))
        if h_coord is not None:
            donors.append({
                "residue_key": f"{atom['chain_id']}:{atom['res_name']}{atom['res_id']}",
                "res_name": atom["res_name"],
                "chain_id": atom["chain_id"],
                "res_id": atom["res_id"],
                "atom_name": atom["atom_name"],
                "coord": atom["coord"],
                "h_coord": h_coord,
            })
    return donors


def _find_acceptors(protein_atoms):
    """Find H-bond acceptor atoms (O, N with lone pairs)."""
    acceptors = []
    for atom in protein_atoms:
        name = atom["atom_name"].strip()
        res = atom["res_name"]
        # Backbone carbonyl O
        if name == "O":
            acceptors.append(_acceptor_entry(atom, "protein"))
        # Side chain O acceptors
        elif name in ("OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH"):
            acceptors.append(_acceptor_entry(atom, "protein"))
        # HIS ring N acceptors (when deprotonated)
        elif name in ("ND1", "NE2") and res == "HIS":
            acceptors.append(_acceptor_entry(atom, "protein"))
        # Water oxygen
        elif name == "O" and res in ("HOH", "WAT"):
            acceptors.append(_acceptor_entry(atom, "protein"))
    return acceptors


def _find_ligand_acceptors(ligand_atoms):
    """Find ligand acceptor atoms (O, N)."""
    acceptors = []
    for atom in (ligand_atoms or []):
        elem = atom.get("element", "").strip().upper()
        if elem in ("O", "N"):
            acceptors.append(_acceptor_entry(atom, "ligand"))
    return acceptors


def _acceptor_entry(atom, source):
    return {
        "residue_key": f"{atom['chain_id']}:{atom['res_name']}{atom['res_id']}"
                       if source == "protein" else f"ligand:{atom.get('res_name','')}",
        "res_name": atom.get("res_name", ""),
        "chain_id": atom.get("chain_id", ""),
        "res_id": atom.get("res_id", ""),
        "atom_name": atom.get("atom_name", ""),
        "coord": atom.get("coord"),
        "source": source,
    }


# ---- H position estimation from heavy atoms ----

def _estimate_h_position(atom, residue_atoms_list):
    """Estimate hydrogen position for a donor heavy atom from local geometry.

    Returns (x, y, z) or None if geometry cannot be reliably estimated.
    """
    name = atom["atom_name"].strip()
    res = atom["res_name"]
    coord = atom["coord"]
    neighbors = _get_bonded_neighbors(atom, residue_atoms_list)

    # Backbone amide N-H
    if name == "N" and res not in ("PRO", "HYP"):
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=1.01, out=True)
        if h_pos:
            return h_pos

    # Side chain donors
    if name in ("ND2",) and res in ("ASN",):  # ASN side chain NH2
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=1.01, out=True)
        if h_pos:
            return h_pos

    if name in ("NE2",) and res in ("GLN",):  # GLN side chain NH2
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=1.01, out=True)
        if h_pos:
            return h_pos

    if name == "NZ" and res == "LYS":  # LYS NH3+
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=1.01, out=True)
        if h_pos:
            return h_pos

    if name in ("NH1", "NH2") and res == "ARG":  # ARG guanidinium
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=1.01, out=True)
        if h_pos:
            return h_pos

    if name in ("ND1", "NE2") and res == "HIS":  # HIS ring NH
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=1.01, out=True)
        if h_pos:
            return h_pos

    if name == "NE1" and res == "TRP":  # TRP indole NH
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=1.01, out=True)
        if h_pos:
            return h_pos

    # Hydroxyl donors (SER, THR, TYR)
    if name in ("OG",) and res == "SER":
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=0.96, out=True)
        if h_pos:
            return h_pos
    if name in ("OG1",) and res == "THR":
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=0.96, out=True)
        if h_pos:
            return h_pos
    if name in ("OH",) and res == "TYR":
        h_pos = _extend_from_geometry(coord, neighbors, bond_len=0.96, out=True)
        if h_pos:
            return h_pos

    return None


def _get_bonded_neighbors(atom, residue_atoms):
    """Get coordinates of atoms bonded to the given atom (same residue, within 1.8 A)."""
    neighbors = []
    for other in residue_atoms:
        if other is atom:
            continue
        dist = _distance(atom["coord"], other["coord"])
        if dist < 1.8:
            neighbors.append(other["coord"])
    return neighbors


def _extend_from_geometry(coord, neighbors, bond_len=1.01, out=True):
    """Estimate H position by extending outward from bonded neighbors.

    For a donor atom D with bonded neighbor B:
      H = D + bond_len * unit_vector(D - B)   [if one neighbor]
      H = D + bond_len * unit_vector(D - centroid(neighbors))   [if multiple]
    """
    if not neighbors:
        return None

    # Centroid of bonded neighbors
    cx = sum(n[0] for n in neighbors) / len(neighbors)
    cy = sum(n[1] for n in neighbors) / len(neighbors)
    cz = sum(n[2] for n in neighbors) / len(neighbors)

    # Direction from neighbors toward donor (outward for H)
    dx = coord[0] - cx
    dy = coord[1] - cy
    dz = coord[2] - cz
    norm = math.sqrt(dx * dx + dy * dy + dz * dz)
    if norm < 1e-8:
        return None

    # H position
    scale = bond_len / norm
    return (coord[0] + dx * scale, coord[1] + dy * scale, coord[2] + dz * scale)


# ---- Geometry helpers ----

def _distance(c1, c2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


def _angle_dha(d, h, a):
    """Compute D-H···A angle in degrees."""
    # Vectors: D→H and H→A
    dh_x = h[0] - d[0]
    dh_y = h[1] - d[1]
    dh_z = h[2] - d[2]
    ha_x = a[0] - h[0]
    ha_y = a[1] - h[1]
    ha_z = a[2] - h[2]

    dh_norm = math.sqrt(dh_x * dh_x + dh_y * dh_y + dh_z * dh_z)
    ha_norm = math.sqrt(ha_x * ha_x + ha_y * ha_y + ha_z * ha_z)
    if dh_norm < 1e-8 or ha_norm < 1e-8:
        return 0.0

    dot = (dh_x * ha_x + dh_y * ha_y + dh_z * ha_z) / (dh_norm * ha_norm)
    # We want the D-H···A angle, which is the angle between vectors H→D and H→A
    # = angle between -(D→H) and H→A
    # So we negate the dot product
    dot = -dot
    dot = max(-1.0, min(1.0, dot))
    return math.degrees(math.acos(dot))


def _empty_result():
    return {
        "hbonds": [],
        "protein_protein": [],
        "protein_ligand": [],
        "summary": {
            "total": 0, "validated": 0, "possible": 0,
            "protein_protein": 0, "protein_ligand": 0,
        },
    }
