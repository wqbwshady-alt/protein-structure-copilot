import re

from analysis_core import analyze_ligand_pocket, classify_residue, parse_pdb_atoms


ONE_TO_THREE = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

AA_PROPERTIES = {
    "ALA": {"charge": "neutral", "polarity": "nonpolar", "hydrophobicity": "hydrophobic", "aromaticity": "non-aromatic", "sidechain_size": "small"},
    "ARG": {"charge": "positive", "polarity": "charged polar", "hydrophobicity": "hydrophilic", "aromaticity": "non-aromatic", "sidechain_size": "large"},
    "ASN": {"charge": "neutral", "polarity": "polar", "hydrophobicity": "hydrophilic", "aromaticity": "non-aromatic", "sidechain_size": "medium"},
    "ASP": {"charge": "negative", "polarity": "charged polar", "hydrophobicity": "hydrophilic", "aromaticity": "non-aromatic", "sidechain_size": "medium"},
    "CYS": {"charge": "neutral", "polarity": "polarizable", "hydrophobicity": "weakly hydrophobic", "aromaticity": "non-aromatic", "sidechain_size": "small"},
    "GLN": {"charge": "neutral", "polarity": "polar", "hydrophobicity": "hydrophilic", "aromaticity": "non-aromatic", "sidechain_size": "medium"},
    "GLU": {"charge": "negative", "polarity": "charged polar", "hydrophobicity": "hydrophilic", "aromaticity": "non-aromatic", "sidechain_size": "large"},
    "GLY": {"charge": "neutral", "polarity": "weakly polar", "hydrophobicity": "neutral", "aromaticity": "non-aromatic", "sidechain_size": "tiny"},
    "HIS": {"charge": "weak positive", "polarity": "polar", "hydrophobicity": "mixed", "aromaticity": "aromatic", "sidechain_size": "medium"},
    "ILE": {"charge": "neutral", "polarity": "nonpolar", "hydrophobicity": "hydrophobic", "aromaticity": "non-aromatic", "sidechain_size": "large"},
    "LEU": {"charge": "neutral", "polarity": "nonpolar", "hydrophobicity": "hydrophobic", "aromaticity": "non-aromatic", "sidechain_size": "large"},
    "LYS": {"charge": "positive", "polarity": "charged polar", "hydrophobicity": "hydrophilic", "aromaticity": "non-aromatic", "sidechain_size": "large"},
    "MET": {"charge": "neutral", "polarity": "nonpolar", "hydrophobicity": "hydrophobic", "aromaticity": "non-aromatic", "sidechain_size": "large"},
    "PHE": {"charge": "neutral", "polarity": "nonpolar", "hydrophobicity": "hydrophobic", "aromaticity": "aromatic", "sidechain_size": "large"},
    "PRO": {"charge": "neutral", "polarity": "nonpolar", "hydrophobicity": "hydrophobic", "aromaticity": "non-aromatic", "sidechain_size": "medium"},
    "SER": {"charge": "neutral", "polarity": "polar", "hydrophobicity": "hydrophilic", "aromaticity": "non-aromatic", "sidechain_size": "small"},
    "THR": {"charge": "neutral", "polarity": "polar", "hydrophobicity": "hydrophilic", "aromaticity": "non-aromatic", "sidechain_size": "medium"},
    "TRP": {"charge": "neutral", "polarity": "nonpolar", "hydrophobicity": "hydrophobic", "aromaticity": "aromatic", "sidechain_size": "large"},
    "TYR": {"charge": "neutral", "polarity": "polar aromatic", "hydrophobicity": "mixed", "aromaticity": "aromatic", "sidechain_size": "large"},
    "VAL": {"charge": "neutral", "polarity": "nonpolar", "hydrophobicity": "hydrophobic", "aromaticity": "non-aromatic", "sidechain_size": "medium"},
}


class MutationScanError(ValueError):
    pass


def parse_mutation(mutation_text):
    match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", mutation_text.strip().upper())

    if not match:
        raise MutationScanError("Mutation must look like R273H, K120A, or Y220C.")

    wild_type_letter, residue_id, mutant_letter = match.groups()

    if wild_type_letter not in ONE_TO_THREE or mutant_letter not in ONE_TO_THREE:
        raise MutationScanError("Mutation contains an unsupported amino acid code.")

    return {
        "mutation": f"{wild_type_letter}{residue_id}{mutant_letter}",
        "wild_type_residue": ONE_TO_THREE[wild_type_letter],
        "residue_id": residue_id,
        "mutant_residue": ONE_TO_THREE[mutant_letter],
    }


def find_residue(pdb_path, mutation, chain_id=None):
    matches = {}

    for atom in parse_pdb_atoms(pdb_path):
        if atom["atom_type"] != "ATOM":
            continue
        if atom["res_id"] != mutation["residue_id"]:
            continue
        if atom["res_name"] != mutation["wild_type_residue"]:
            continue
        if chain_id and atom["chain_id"] != chain_id:
            continue

        key = (atom["chain_id"], atom["res_name"], atom["res_id"])
        matches[key] = {
            "chain_id": atom["chain_id"],
            "res_name": atom["res_name"],
            "res_id": atom["res_id"],
            "atom_name": "",
            "element": "",
        }

    if not matches:
        raise MutationScanError(
            f"Could not find {mutation['wild_type_residue']}{mutation['residue_id']} in the uploaded PDB."
        )

    sorted_matches = sorted(matches.values(), key=lambda item: item["chain_id"])
    selected = sorted_matches[0]
    note = ""

    if len(sorted_matches) > 1 and not chain_id:
        chains = ", ".join(item["chain_id"] or "(blank)" for item in sorted_matches)
        note = f"Multiple matching chains found ({chains}); chain {selected['chain_id']} was selected."

    return selected, note


def compare_properties(wild_type_residue, mutant_residue):
    wild_type_properties = AA_PROPERTIES[wild_type_residue]
    mutant_properties = AA_PROPERTIES[mutant_residue]
    changes = {}

    for key in ["charge", "polarity", "hydrophobicity", "aromaticity", "sidechain_size"]:
        old = wild_type_properties[key]
        new = mutant_properties[key]
        changes[key] = "unchanged" if old == new else f"{old} -> {new}"

    return wild_type_properties, mutant_properties, changes


def filter_residue_interactions(interactions, residue):
    return [
        interaction for interaction in interactions
        if interaction["chain_id"] == residue["chain_id"] and
        interaction["res_name"] == residue["res_name"] and
        interaction["res_id"] == residue["res_id"]
    ]


def classify_interaction_impact(original_interactions, property_changes, mutant_properties):
    possible_loss = []
    possible_gain = []

    for interaction in original_interactions:
        interaction_type = interaction["interaction_type"]

        if interaction_type == "charged / electrostatic":
            if property_changes["charge"] != "unchanged":
                possible_loss.append({
                    "interaction_type": interaction_type,
                    "reason": f"Charge changes from {property_changes['charge']}, so the original electrostatic contact may weaken.",
                    "distance": interaction["distance"],
                    "chain_id": interaction["chain_id"],
                    "res_name": interaction["res_name"],
                    "res_id": interaction["res_id"],
                    "atom_name": interaction["atom_name"],
                    "element": interaction["element"],
                })
            continue

        if interaction_type == "polar / possible H-bond":
            if "nonpolar" in mutant_properties["polarity"]:
                possible_loss.append({
                    "interaction_type": interaction_type,
                    "reason": "The mutant residue is less polar, so the original polar contact may be lost.",
                    "distance": interaction["distance"],
                    "chain_id": interaction["chain_id"],
                    "res_name": interaction["res_name"],
                    "res_id": interaction["res_id"],
                    "atom_name": interaction["atom_name"],
                    "element": interaction["element"],
                })
            continue

        if interaction_type == "hydrophobic contact":
            if mutant_properties["hydrophobicity"] not in {"hydrophobic", "weakly hydrophobic", "mixed"}:
                possible_loss.append({
                    "interaction_type": interaction_type,
                    "reason": "The mutant residue is less hydrophobic, so hydrophobic packing may weaken.",
                    "distance": interaction["distance"],
                    "chain_id": interaction["chain_id"],
                    "res_name": interaction["res_name"],
                    "res_id": interaction["res_id"],
                    "atom_name": interaction["atom_name"],
                    "element": interaction["element"],
                })
            continue

        if property_changes["sidechain_size"] != "unchanged":
            possible_loss.append({
                "interaction_type": interaction_type,
                "reason": "Sidechain size changes, so local van der Waals packing may be altered.",
                "distance": interaction["distance"],
                "chain_id": interaction["chain_id"],
                "res_name": interaction["res_name"],
                "res_id": interaction["res_id"],
                "atom_name": interaction["atom_name"],
                "element": interaction["element"],
            })

    if property_changes["aromaticity"].endswith("-> aromatic"):
        possible_gain.append({
            "interaction_type": "aromatic / packing",
            "reason": "The mutant introduces an aromatic sidechain, which could create new packing contacts if geometry permits.",
        })

    if property_changes["hydrophobicity"].endswith("-> hydrophobic"):
        possible_gain.append({
            "interaction_type": "hydrophobic contact",
            "reason": "The mutant is more hydrophobic, which could improve local nonpolar packing if positioned near the ligand.",
        })

    if property_changes["charge"] != "unchanged" and mutant_properties["charge"] != "neutral":
        possible_gain.append({
            "interaction_type": "charged / electrostatic",
            "reason": "The mutant changes charge state, which could create new electrostatic contacts if compatible ligand atoms are nearby.",
        })

    return {
        "possible_loss": possible_loss,
        "possible_gain": possible_gain,
    }


def build_mutation_interpretation(mutation_result):
    mutation = mutation_result["mutation"]
    original = mutation_result["original_residue"]
    mutant = mutation_result["mutant_residue"]
    property_changes = mutation_result["property_changes"]
    impact = mutation_result["interaction_impact"]
    original_interactions = mutation_result["original_interactions"]

    interaction_count = len(original_interactions)
    loss_count = len(impact["possible_loss"])
    gain_count = len(impact["possible_gain"])

    if interaction_count == 0:
        interaction_sentence = (
            "No direct ligand interaction from this residue was detected in the current interaction set, "
            "so the predicted binding effect is mainly property-based."
        )
    elif loss_count > 0:
        interaction_sentence = (
            f"The wild-type residue contributes {interaction_count} close ligand interaction(s), "
            f"and {loss_count} may weaken after mutation."
        )
    else:
        interaction_sentence = (
            f"The wild-type residue contributes {interaction_count} close ligand interaction(s), "
            "with no strong rule-based loss predicted."
        )

    return (
        f"Mutation Impact Summary\n\n"
        f"{mutation} replaces {original['res_name']}{original['res_id']} with {mutant['res_name']} "
        f"on chain {mutation_result['chain_id']}. Charge change: {property_changes['charge']}; "
        f"hydrophobicity change: {property_changes['hydrophobicity']}; polarity change: {property_changes['polarity']}; "
        f"aromaticity change: {property_changes['aromaticity']}; sidechain size change: {property_changes['sidechain_size']}.\n\n"
        f"{interaction_sentence} Possible gained interaction classes: {gain_count}. "
        "Because no sidechain remodeling was performed, this is a heuristic risk estimate rather than a geometry-refined mutant model."
    )


def analyze_mutation_scan(pdb_path, ligand_name, mutation_text, chain_id=None):
    parsed_mutation = parse_mutation(mutation_text)
    selected_residue, selection_note = find_residue(pdb_path, parsed_mutation, chain_id=chain_id)

    contact_residues, counts, primary_interpretation, interactions = analyze_ligand_pocket(
        pdb_path,
        ligand_name,
        max_interactions=None
    )

    if contact_residues is None:
        raise MutationScanError(f"No ligand named {ligand_name} found in this PDB file.")

    wild_type_properties, mutant_properties, property_changes = compare_properties(
        parsed_mutation["wild_type_residue"],
        parsed_mutation["mutant_residue"]
    )
    original_interactions = filter_residue_interactions(interactions, selected_residue)
    interaction_impact = classify_interaction_impact(
        original_interactions,
        property_changes,
        mutant_properties
    )

    mutation_result = {
        "mutation": parsed_mutation["mutation"],
        "chain_id": selected_residue["chain_id"],
        "selection_note": selection_note,
        "original_residue": {
            "chain_id": selected_residue["chain_id"],
            "res_name": parsed_mutation["wild_type_residue"],
            "res_id": parsed_mutation["residue_id"],
            "atom_name": "",
            "element": "",
        },
        "mutant_residue": {
            "chain_id": selected_residue["chain_id"],
            "res_name": parsed_mutation["mutant_residue"],
            "res_id": parsed_mutation["residue_id"],
            "atom_name": "",
            "element": "",
        },
        "wild_type_properties": wild_type_properties,
        "mutant_properties": mutant_properties,
        "property_changes": property_changes,
        "original_interactions": original_interactions,
        "interaction_impact": interaction_impact,
    }
    mutation_result["ai_interpretation"] = build_mutation_interpretation(mutation_result)

    return mutation_result
