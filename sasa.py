"""
Shrake-Rupley solvent accessible surface area (SASA) computation.

Computes per-atom and per-residue SASA using the rolling-probe algorithm.
Identifies buried vs exposed residues in the binding pocket.

Algorithm: Shrake & Rupley (1973) — place probe sphere points on each
atom's surface, count unoccluded points.
"""

import math


# Van der Waals radii (A) from Bondi / CHARMM
_VDW_RADII = {
    "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85,
    "I": 1.98, "H": 1.20,
    "MG": 1.40, "MN": 1.45, "ZN": 1.30, "CA": 1.80,
}

_PROBE_RADIUS = 1.4  # A — water probe radius

# Relative SASA thresholds for burial classification
_BURIED_THRESHOLD = 0.25     # < 25% of max = buried
_EXPOSED_THRESHOLD = 0.50    # > 50% of max = exposed


def _SPHERE_POINTS(n):
    """Generate n approximately uniformly distributed points on a unit sphere.

    Uses the Fibonacci sphere algorithm for even distribution.
    Returns list of (x, y, z) tuples.
    """
    points = []
    phi = math.pi * (3.0 - math.sqrt(5.0))  # Golden angle
    for i in range(n):
        y = 1.0 - (i / float(n - 1)) * 2.0 if n > 1 else 0.0
        radius = math.sqrt(1.0 - y * y)
        theta = phi * i
        points.append((math.cos(theta) * radius, y, math.sin(theta) * radius))
    return points


def _shrake_rupley(atoms, radii, n_points=200):
    """Compute per-atom SASA using Shrake-Rupley algorithm.

    Args:
        atoms: list of dicts with "coord" (x, y, z) and "element" keys
        radii: dict mapping element → vdW radius
        n_points: number of probe points per atom

    Returns:
        list of per-atom SASA values (A^2)
    """
    sphere = _SPHERE_POINTS(n_points)
    n = len(atoms)
    sasa = [0.0] * n

    for i in range(n):
        elem_i = atoms[i].get("element", "C").strip().upper()
        r_i = radii.get(elem_i, 1.70) + _PROBE_RADIUS
        coord_i = atoms[i]["coord"]

        # Neighbors within possible overlap distance
        neighbors = []
        for j in range(n):
            if i == j:
                continue
            elem_j = atoms[j].get("element", "C").strip().upper()
            r_j = radii.get(elem_j, 1.70) + _PROBE_RADIUS
            dist = _distance(coord_i, atoms[j]["coord"])
            if dist < r_i + r_j:
                neighbors.append((atoms[j]["coord"], r_j))

        area_per_point = 4.0 * math.pi * r_i * r_i / n_points
        accessible = 0

        for sp in sphere:
            # Probe point position on atom i's surface
            px = coord_i[0] + r_i * sp[0]
            py = coord_i[1] + r_i * sp[1]
            pz = coord_i[2] + r_i * sp[2]

            occluded = False
            for n_coord, n_r in neighbors:
                dx = px - n_coord[0]
                dy = py - n_coord[1]
                dz = pz - n_coord[2]
                if dx * dx + dy * dy + dz * dz < n_r * n_r:
                    occluded = True
                    break

            if not occluded:
                accessible += 1

        sasa[i] = accessible * area_per_point

    return sasa


def compute_sasa(all_atoms, contact_residues):
    """Compute per-residue SASA for pocket contact residues.

    Args:
        all_atoms: list of atom dicts from parse_pdb_atoms
        contact_residues: {(chain, res_name, res_id): atom_dict}

    Returns:
        dict with per_residue SASA, total SASA, and burial classification
    """
    if not all_atoms:
        return {"per_residue": {}, "total_sasa": 0.0}

    # Compute per-atom SASA for ALL atoms
    atom_sasa = _shrake_rupley(all_atoms, _VDW_RADII)

    # Group by residue for contact residues
    residue_atoms = {}
    residue_indices = {}
    for i, atom in enumerate(all_atoms):
        rkey = (atom.get("chain_id", ""), atom.get("res_name", ""), atom.get("res_id", ""))
        if rkey not in residue_atoms:
            residue_atoms[rkey] = []
            residue_indices[rkey] = []
        residue_atoms[rkey].append(atom)
        residue_indices[rkey].append(i)

    # Compute per-residue SASA for contact residues only
    per_residue = {}
    total_sasa = 0.0

    for (chain_id, res_name, res_id) in contact_residues:
        rkey = (chain_id, res_name, res_id)
        residue_key = f"{chain_id}:{res_name}{res_id}"

        indices = residue_indices.get(rkey, [])
        if not indices:
            continue

        res_sasa = sum(atom_sasa[i] for i in indices)
        res_sasa = round(res_sasa, 1)

        # Estimate max SASA for this residue type (rough reference values)
        max_sasa = _estimate_max_sasa(res_name, residue_atoms.get(rkey, []))
        if max_sasa > 0:
            relative = res_sasa / max_sasa
        else:
            relative = 0.0

        if relative < _BURIED_THRESHOLD:
            classification = "buried"
        elif relative > _EXPOSED_THRESHOLD:
            classification = "exposed"
        else:
            classification = "partially_buried"

        per_residue[residue_key] = {
            "sasa": res_sasa,
            "relative_sasa": round(relative, 3),
            "classification": classification,
        }
        total_sasa += res_sasa

    return {
        "per_residue": per_residue,
        "total_sasa": round(total_sasa, 1),
    }


def _estimate_max_sasa(res_name, atoms):
    """Estimate maximum possible SASA for a residue type (Gly-X-Gly tripeptide reference)."""
    ref = {
        "ALA": 115, "ARG": 240, "ASN": 160, "ASP": 155, "CYS": 140,
        "GLN": 190, "GLU": 190, "GLY": 80,  "HIS": 195, "ILE": 180,
        "LEU": 180, "LYS": 215, "MET": 205, "PHE": 215, "PRO": 145,
        "SER": 120, "THR": 145, "TRP": 260, "TYR": 230, "VAL": 160,
    }
    return ref.get(res_name, sum(_VDW_RADII.get(a.get("element", "C"), 1.7) * 8 for a in atoms))


def _distance(c1, c2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))
