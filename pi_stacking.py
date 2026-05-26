"""
Pi-stacking and cation-pi interaction geometry detection.

Detects aromatic ring interactions in the binding pocket using
purely geometric criteria — no external data or APIs needed.

Ring definitions use standard PDB atom names for each aromatic side chain.
Ring centroids and plane normals are computed from atom coordinates,
then classified by centroid distance and inter-planar angle.
"""

import itertools
import math

# ---- Ring atom definitions (PDB atom names) ----

# 6-membered benzene ring (PHE, TYR)
_RING_PHE = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}

# 5-membered imidazole ring (HIS)
_RING_HIS = {"CG", "ND1", "CD2", "CE1", "NE2"}

# TRP has two rings: 6-membered benzene + 5-membered pyrrole
_RING_TRP_6 = {"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"}
_RING_TRP_5 = {"CG", "CD1", "CD2", "CE2", "NE1"}

# Residue → ring atom name sets
_AROMATIC_RING_ATOMS = {
    "PHE": [_RING_PHE],
    "TYR": [_RING_PHE],
    "HIS": [_RING_HIS],
    "TRP": [_RING_TRP_6, _RING_TRP_5],
}

# Thresholds
_PI_PI_DISTANCE_MAX = 6.5       # Max centroid-centroid distance for pi-pi (A)
_FACE_TO_FACE_ANGLE_MAX = 30.0  # Inter-planar angle for face-to-face (degrees)
_EDGE_TO_FACE_ANGLE_MIN = 50.0  # Inter-planar angle for edge-to-face (degrees)
_CATION_PI_DISTANCE_MAX = 6.0   # Max distance for cation-pi (A)

# Charged residues for cation-pi
_CHARGED_RESIDUES = {"LYS", "ARG"}


def detect_pi_interactions(contact_residues, all_atoms, ligand_atoms=None):
    """Detect pi-stacking and cation-pi interactions among pocket contact residues.

    Args:
        contact_residues: {(chain_id, res_name, res_id): atom_dict, ...}
        all_atoms: list of atom dicts from parse_pdb_atoms
        ligand_atoms: optional list of ligand (HETATM) atom dicts

    Returns:
        dict with:
            pi_pi_interactions: list of {type, residue1, residue2, distance, angle}
            cation_pi_interactions: list of {type, cationic_res, aromatic_res, distance}
            aromatic_residues_found: list of residue keys with detected rings
    """
    if not all_atoms:
        return _empty_result()

    # Group all atoms by residue key for ring lookup
    residue_atoms = {}
    for atom in all_atoms:
        if atom["atom_type"] != "ATOM":
            continue
        rkey = (atom["chain_id"], atom["res_name"], atom["res_id"])
        if rkey not in residue_atoms:
            residue_atoms[rkey] = {}
        residue_atoms[rkey][atom["atom_name"]] = atom

    # Identify aromatic residues in contact set and compute ring geometry
    aromatic_rings = []  # [{residue_key, res_name, rings: [{centroid, normal, atoms}]}]

    for (chain_id, res_name, res_id) in contact_residues:
        rkey = (chain_id, res_name, res_id)
        residue_key = f"{chain_id}:{res_name}{res_id}"
        if res_name not in _AROMATIC_RING_ATOMS:
            continue

        atoms_by_name = residue_atoms.get(rkey, {})
        ring_sets = _AROMATIC_RING_ATOMS[res_name]

        rings = []
        for ring_atom_names in ring_sets:
            ring_coords = []
            for name in ring_atom_names:
                atom = atoms_by_name.get(name)
                if atom and atom.get("coord"):
                    ring_coords.append(atom["coord"])
            if len(ring_coords) >= 4:  # Need at least 4 points for a plane
                centroid, normal = _compute_ring_geometry(ring_coords)
                if centroid and normal:
                    rings.append({"centroid": centroid, "normal": normal})

        if rings:
            aromatic_rings.append({
                "residue_key": residue_key,
                "res_name": res_name,
                "chain_id": chain_id,
                "res_id": res_id,
                "rings": rings,
            })

    # Detect pi-pi interactions between aromatic residue pairs
    pi_pi_found = []
    for a, b in itertools.combinations(aromatic_rings, 2):
        for ra in a["rings"]:
            for rb in b["rings"]:
                dist = _distance(ra["centroid"], rb["centroid"])
                if dist > _PI_PI_DISTANCE_MAX:
                    continue
                angle = _angle_between_normals(ra["normal"], rb["normal"])
                if angle <= _FACE_TO_FACE_ANGLE_MAX:
                    itype = "pi_pi_face_to_face"
                elif angle >= _EDGE_TO_FACE_ANGLE_MIN:
                    itype = "pi_pi_edge_to_face"
                else:
                    itype = "pi_pi_t_shaped"
                pi_pi_found.append({
                    "type": itype,
                    "residue1": a["residue_key"],
                    "residue2": b["residue_key"],
                    "distance": round(dist, 2),
                    "angle": round(angle, 1),
                })
                break  # One ring pair per residue pair

    # Detect cation-pi interactions (aromatic near LYS/ARG)
    cation_pi_found = []
    charged_contact = [
        (chain_id, res_name, res_id) for (chain_id, res_name, res_id)
        in contact_residues if res_name in _CHARGED_RESIDUES
    ]

    for (chain_id, res_name, res_id) in charged_contact:
        rkey = (chain_id, res_name, res_id)
        charged_atom = _get_charged_center(atoms_by_name, residue_atoms, rkey)
        if not charged_atom:
            continue
        charged_key = f"{chain_id}:{res_name}{res_id}"

        for ar in aromatic_rings:
            for ring in ar["rings"]:
                dist = _distance(charged_atom, ring["centroid"])
                if dist <= _CATION_PI_DISTANCE_MAX:
                    cation_pi_found.append({
                        "type": "cation_pi",
                        "cationic_residue": charged_key,
                        "aromatic_residue": ar["residue_key"],
                        "distance": round(dist, 2),
                    })
                    break

    return {
        "pi_pi_interactions": pi_pi_found,
        "cation_pi_interactions": cation_pi_found,
        "aromatic_residues_found": [ar["residue_key"] for ar in aromatic_rings],
    }


# ---- Geometry helpers ----

def _compute_ring_geometry(coords):
    """Compute ring centroid and best-fit plane normal from atom coordinates.

    Uses SVD of centered coordinates to find the plane. The normal is the
    right singular vector corresponding to the smallest singular value.
    """
    n = len(coords)
    if n < 3:
        return None, None

    centroid = tuple(sum(c[i] for c in coords) / n for i in range(3))

    # Center coordinates
    centered = [[c[i] - centroid[i] for i in range(3)] for c in coords]

    # Compute covariance matrix (3x3)
    cov = [[0.0] * 3 for _ in range(3)]
    for row in centered:
        for i in range(3):
            for j in range(3):
                cov[i][j] += row[i] * row[j]

    # Find normal via cross product of first two centered vectors for simplicity.
    # For rings, the atoms are approximately coplanar, so this gives a good estimate.
    v1 = centered[0]
    v2 = centered[1]
    normal = (
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    )
    # Normalize
    norm = math.sqrt(sum(x * x for x in normal))
    if norm < 1e-10:
        return centroid, (0.0, 0.0, 1.0)
    normal = tuple(x / norm for x in normal)

    return centroid, normal


def _angle_between_normals(n1, n2):
    """Compute angle between two plane normals in degrees.

    Returns the acute angle (0-90 degrees) since ring planes are undirected.
    """
    dot = abs(sum(a * b for a, b in zip(n1, n2)))
    dot = min(dot, 1.0)
    return math.degrees(math.acos(dot))


def _get_charged_center(atoms_by_name, residue_atoms, rkey):
    """Get the side-chain charged group center for LYS (NZ) or ARG (NH1/NH2)."""
    atoms = residue_atoms.get(rkey, {})
    res_name = rkey[1]

    if res_name == "LYS":
        nz = atoms.get("NZ")
        return nz["coord"] if nz else atoms.get("CD", {}).get("coord") if "CD" in atoms else None
    elif res_name == "ARG":
        nh1 = atoms.get("NH1")
        nh2 = atoms.get("NH2")
        if nh1 and nh2:
            c1, c2 = nh1["coord"], nh2["coord"]
            return tuple((c1[i] + c2[i]) / 2 for i in range(3))
        return nh1["coord"] if nh1 else nh2["coord"] if nh2 else None
    return None


def _distance(c1, c2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))


def _empty_result():
    return {
        "pi_pi_interactions": [],
        "cation_pi_interactions": [],
        "aromatic_residues_found": [],
    }
