import os

from ai_client import generate_ai_interpretation
from analysis_core import get_hotspot_residues, residue_key_to_text, sort_residue_keys


def count_interaction_types(interactions):
    counts = {}

    for interaction in interactions:
        interaction_type = interaction.get("interaction_type", "contact")
        counts[interaction_type] = counts.get(interaction_type, 0) + 1

    return counts


def generate_pymol_script(pdb_load_path, ligand_name, contact_residues, result_folder, output_prefix=None):
    base_name = output_prefix or os.path.splitext(os.path.basename(pdb_load_path))[0]
    pml_filename = f"highlight_{base_name}_{ligand_name}_pocket.pml"
    pml_path = os.path.join(result_folder, pml_filename)

    selections = []

    for chain_id, res_name, res_id in contact_residues.keys():
        if chain_id:
            selections.append(f"(chain {chain_id} and resi {res_id})")
        else:
            selections.append(f"(resi {res_id})")

    pocket_selection = " or ".join(selections) if selections else "none"

    script = f"""
load {pdb_load_path}

hide everything
show cartoon
color cyan, polymer

select ligand, resn {ligand_name}
show sticks, ligand
color red, ligand

select pocket, {pocket_selection}
show sticks, pocket
color yellow, pocket

zoom ligand, 8
"""

    with open(pml_path, "w", encoding="utf-8") as f:
        f.write(script)

    return pml_filename


def build_report(ligand_name, contact_residues, counts, primary_interpretation, pymol_filename, interactions):
    residue_lines = []

    for chain_id, res_name, res_id in sort_residue_keys(contact_residues.keys()):
        residue_lines.append(f"- Chain {chain_id}: {res_name}{res_id}")

    interaction_lines = []

    for interaction in interactions:
        interaction_lines.append(
            f"- Chain {interaction['chain_id']}: ligand atom {interaction.get('ligand_atom', '?')} -> "
            f"{interaction['res_name']}{interaction['res_id']} atom {interaction['atom_name']} | "
            f"{interaction.get('interaction_type', 'contact')} | {interaction['distance']} A"
        )

    hotspot_lines = []

    for hotspot in get_hotspot_residues(interactions):
        hotspot_lines.append(
            f"- Chain {hotspot['chain_id']}: {hotspot['res_name']}{hotspot['res_id']} "
            f"nearest contact {hotspot['distance']} A"
        )

    interaction_summary_lines = [
        f"- {interaction_type}: {count}"
        for interaction_type, count in sorted(count_interaction_types(interactions).items())
    ] or ["- No close interaction lines detected"]

    ai_text = generate_ai_interpretation(
        contact_residues,
        counts,
        interactions,
        ligand_name
    )

    report_text = f"""Target ligand: {ligand_name}

[1] Ligand Pocket Summary
Contact residues within 5A: {len(contact_residues)}
Hydrophobic residues: {counts["hydrophobic"]}
Polar residues: {counts["polar"]}
Positive residues: {counts["positive"]}
Negative residues: {counts["negative"]}

[2] Interaction Type Summary
{chr(10).join(interaction_summary_lines)}

[3] Contact Residues
{chr(10).join(residue_lines)}

[4] Interaction Lines
Showing nearest ligand contact per residue, max 15 lines:
{chr(10).join(interaction_lines)}

[5] Hotspot Residues
Top 5 closest ligand-contact residues:
{chr(10).join(hotspot_lines)}

[6] DeepSeek AI Structural Interpretation
{ai_text}

Primary pocket interpretation:
{primary_interpretation}

[7] Visualization
PyMOL script:
results/{pymol_filename}
"""

    return report_text, ai_text


def build_comparison_report(ligand_name, wt_contacts, mut_contacts, wt_counts, mut_counts):
    wt_set = set(wt_contacts.keys())
    mut_set = set(mut_contacts.keys())

    lost = sort_residue_keys(wt_set - mut_set)
    gained = sort_residue_keys(mut_set - wt_set)
    shared = sort_residue_keys(wt_set & mut_set)

    lost_lines = [f"- {residue_key_to_text(x)}" for x in lost] or ["- None"]
    gained_lines = [f"- {residue_key_to_text(x)}" for x in gained] or ["- None"]

    hydrophobic_change = mut_counts["hydrophobic"] - wt_counts["hydrophobic"]
    polar_change = mut_counts["polar"] - wt_counts["polar"]
    positive_change = mut_counts["positive"] - wt_counts["positive"]
    negative_change = mut_counts["negative"] - wt_counts["negative"]

    if len(lost) == 0 and len(gained) == 0:
        interpretation = (
            "No major pocket residue change detected. The mutant pocket appears similar to WT based on 5A ligand contact residues."
        )
    elif len(lost) > len(gained):
        interpretation = (
            "The mutant lost more ligand-contact residues than it gained, suggesting possible weakening or rearrangement of the binding pocket."
        )
    elif hydrophobic_change < 0:
        interpretation = (
            "The mutant shows reduced hydrophobic pocket composition, suggesting possible loss of hydrophobic packing around the ligand."
        )
    elif len(gained) > len(lost):
        interpretation = (
            "The mutant gained additional ligand-contact residues, suggesting possible pocket reorganization or altered ligand environment."
        )
    else:
        interpretation = (
            "The mutant pocket shows moderate residue-level changes, suggesting local binding environment remodeling."
        )

    report_text = f"""WT vs Mutant Pocket Comparison

Target ligand: {ligand_name}

[1] Pocket Size
WT contact residues within 5A: {len(wt_contacts)}
Mutant contact residues within 5A: {len(mut_contacts)}
Shared contact residues: {len(shared)}
Lost contact residues: {len(lost)}
Gained contact residues: {len(gained)}

[2] Residue Property Change
Hydrophobic: WT {wt_counts["hydrophobic"]} -> Mutant {mut_counts["hydrophobic"]} ({hydrophobic_change:+d})
Polar: WT {wt_counts["polar"]} -> Mutant {mut_counts["polar"]} ({polar_change:+d})
Positive: WT {wt_counts["positive"]} -> Mutant {mut_counts["positive"]} ({positive_change:+d})
Negative: WT {wt_counts["negative"]} -> Mutant {mut_counts["negative"]} ({negative_change:+d})

[3] Lost Residues
{chr(10).join(lost_lines)}

[4] Gained Residues
{chr(10).join(gained_lines)}

[5] Interpretation
{interpretation}
"""

    return report_text, lost, gained


def build_mutation_scan_report(mutation_result):
    property_lines = [
        f"- {name}: {value}"
        for name, value in mutation_result["property_changes"].items()
    ]
    loss_lines = [
        f"- {item['interaction_type']}: {item['reason']}"
        for item in mutation_result["interaction_impact"]["possible_loss"]
    ] or ["- None predicted by current rules"]
    gain_lines = [
        f"- {item['interaction_type']}: {item['reason']}"
        for item in mutation_result["interaction_impact"]["possible_gain"]
    ] or ["- None predicted by current rules"]
    interaction_lines = [
        (
            f"- Chain {item['chain_id']}: {item['res_name']}{item['res_id']} "
            f"atom {item['atom_name']} | {item['interaction_type']} | {item['distance']} A"
        )
        for item in mutation_result["original_interactions"]
    ] or ["- No direct ligand interaction from this residue was detected in the current interaction set"]

    note = ""

    if mutation_result.get("selection_note"):
        note = f"\nSelection note: {mutation_result['selection_note']}\n"

    return f"""Mutation Impact Summary

Mutation: {mutation_result["mutation"]}
Chain: {mutation_result["chain_id"]}
Original residue: {mutation_result["original_residue"]["res_name"]}{mutation_result["original_residue"]["res_id"]}
Mutant residue: {mutation_result["mutant_residue"]["res_name"]}{mutation_result["mutant_residue"]["res_id"]}
{note}
[1] Physicochemical Property Changes
{chr(10).join(property_lines)}

[2] Original Residue Interactions
{chr(10).join(interaction_lines)}

[3] Possible Interaction Loss
{chr(10).join(loss_lines)}

[4] Possible Interaction Gain
{chr(10).join(gain_lines)}

[5] AI-style Mutation Interpretation
{mutation_result["ai_interpretation"]}
"""
