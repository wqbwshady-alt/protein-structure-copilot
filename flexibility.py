"""
B-factor pocket flexibility analysis.

PDB ATOM records contain temperature factors (B-factors) that reflect
atomic displacement — a proxy for local structural flexibility.
High B-factor = flexible/disordered; low B-factor = rigid/well-ordered.

This module compares pocket residue B-factors against the whole protein
to classify pocket regions by flexibility.
"""

import statistics


def compute_pocket_flexibility(contact_residues, all_atoms):
    """Compute B-factor flexibility for pocket contact residues.

    Args:
        contact_residues: {(chain_id, res_name, res_id): atom_dict, ...}
        all_atoms: list of atom dicts from parse_pdb_atoms (with b_factor field)

    Returns:
        dict with:
            per_residue: {residue_key: {mean_b, classification, z_score}}
            pocket_summary: {mean_b, global_mean_b, global_std_b,
                            rigid_count, flexible_count, highly_flexible_count,
                            flexibility_ratio, classification}
    """
    if not all_atoms:
        return _empty_result()

    # Group atoms by residue for per-residue B-factor averaging
    all_bfactors = [a["b_factor"] for a in all_atoms if a.get("b_factor", 0) > 0]
    if not all_bfactors:
        return _empty_result()

    global_mean = statistics.mean(all_bfactors)
    global_std = statistics.stdev(all_bfactors) if len(all_bfactors) > 1 else 1.0

    # Compute per-residue mean B-factor for all residues in the protein
    residue_bfactors = {}
    for atom in all_atoms:
        b = atom.get("b_factor", 0)
        if b <= 0:
            continue
        key = (atom["chain_id"], atom["res_name"], atom["res_id"])
        if key not in residue_bfactors:
            residue_bfactors[key] = []
        residue_bfactors[key].append(b)

    # Classify pocket contact residues
    per_residue = {}
    pocket_bfactors = []

    for (chain_id, res_name, res_id) in contact_residues:
        rkey = (chain_id, res_name, res_id)
        residue_key = f"{chain_id}:{res_name}{res_id}"
        bfactors = residue_bfactors.get(rkey, [])

        if not bfactors:
            per_residue[residue_key] = {
                "mean_b": None,
                "z_score": None,
                "classification": "unknown",
                "label": "no B-factor data",
            }
            continue

        mean_b = statistics.mean(bfactors)
        pocket_bfactors.append(mean_b)
        z_score = (mean_b - global_mean) / global_std if global_std > 0 else 0

        if mean_b < global_mean * 0.5:
            classification = "rigid"
            label = "Rigid"
        elif mean_b > global_mean * 2.0:
            classification = "highly_flexible"
            label = "Highly flexible"
        elif mean_b > global_mean * 1.5:
            classification = "flexible"
            label = "Flexible"
        else:
            classification = "normal"
            label = "Normal"

        per_residue[residue_key] = {
            "mean_b": round(mean_b, 2),
            "z_score": round(z_score, 2),
            "classification": classification,
            "label": label,
        }

    # Pocket-level summary
    if pocket_bfactors:
        pocket_mean = statistics.mean(pocket_bfactors)
        rigid_count = sum(1 for v in per_residue.values() if v["classification"] == "rigid")
        flexible_count = sum(1 for v in per_residue.values() if v["classification"] == "flexible")
        highly_flexible_count = sum(1 for v in per_residue.values() if v["classification"] == "highly_flexible")

        if rigid_count > len(pocket_bfactors) * 0.4:
            pocket_class = "rigid_pocket"
            pocket_label = "Rigid binding site — well-ordered, suitable for structure-based design"
        elif highly_flexible_count > 0:
            pocket_class = "flexible_pocket"
            pocket_label = "Contains highly flexible regions — may undergo induced fit"
        elif flexible_count > len(pocket_bfactors) * 0.3:
            pocket_class = "partially_flexible"
            pocket_label = "Partially flexible pocket — some conformational adaptability"
        else:
            pocket_class = "normal_pocket"
            pocket_label = "Normal flexibility — typical crystalline binding site"
    else:
        pocket_mean = global_mean
        rigid_count = flexible_count = highly_flexible_count = 0
        pocket_class = "unknown"
        pocket_label = "Insufficient B-factor data for pocket residues"

    return {
        "per_residue": per_residue,
        "pocket_summary": {
            "mean_b": round(pocket_mean, 2),
            "global_mean_b": round(global_mean, 2),
            "global_std_b": round(global_std, 2),
            "rigid_count": rigid_count,
            "flexible_count": flexible_count,
            "highly_flexible_count": highly_flexible_count,
            "total_in_pocket": len(pocket_bfactors),
            "flexibility_ratio": round(pocket_mean / global_mean, 2) if global_mean > 0 else 1.0,
            "classification": pocket_class,
            "label": pocket_label,
        },
    }


def _empty_result():
    return {
        "per_residue": {},
        "pocket_summary": {
            "mean_b": None,
            "global_mean_b": None,
            "global_std_b": None,
            "rigid_count": 0,
            "flexible_count": 0,
            "highly_flexible_count": 0,
            "total_in_pocket": 0,
            "flexibility_ratio": None,
            "classification": "no_data",
            "label": "No B-factor data available",
        },
    }
