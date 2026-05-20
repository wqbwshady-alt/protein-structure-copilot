import os

from openai import OpenAI


def create_deepseek_client():
    return OpenAI(
        api_key=os.getenv("DEEPSEEK_API_KEY"),
        base_url="https://api.deepseek.com"
    )


def local_ai_interpretation(contact_residues, counts, interactions, ligand_name):
    total_contacts = len(contact_residues)
    h_count = counts["hydrophobic"]
    p_count = counts["polar"]
    pos_count = counts["positive"]
    neg_count = counts["negative"]

    # Section 1: Pocket type classification
    if h_count > max(p_count, pos_count, neg_count):
        pocket_type = "Primarily hydrophobic pocket"
        pocket_detail = (
            f"Hydrophobic residues dominate ({h_count} of {total_contacts} contacts), "
            "suggesting nonpolar packing is the main stabilization mechanism. "
            "This is characteristic of small-molecule binding sites in enzymes and receptors."
        )
    elif pos_count > neg_count and pos_count >= 3:
        pocket_type = "Positively charged pocket"
        pocket_detail = (
            f"Enrichment in basic residues (Lys/Arg/His: {pos_count}) suggests affinity for "
            "negatively charged ligands, phosphate groups, or nucleic acid substrates."
        )
    elif neg_count > pos_count and neg_count >= 3:
        pocket_type = "Negatively charged pocket"
        pocket_detail = (
            f"Enrichment in acidic residues (Asp/Glu: {neg_count}) creates a locally negative "
            "electrostatic environment, potentially favoring cationic ligands or metal ion coordination."
        )
    elif p_count >= h_count:
        pocket_type = "Polar/amphipathic pocket"
        pocket_detail = (
            f"Mixed polar ({p_count}) and hydrophobic ({h_count}) character suggests "
            "a versatile binding site capable of accommodating diverse ligand chemotypes."
        )
    else:
        pocket_type = "Mixed-composition pocket"
        pocket_detail = "Multiple residue types contribute to a composite binding interface."

    # Section 2: Hydrophobic packing
    if h_count >= 8:
        hydrophobic_comment = (
            f"A dense hydrophobic core ({h_count} residues) forms extensive van der Waals contacts "
            "with the ligand. This degree of hydrophobic burial typically contributes significantly "
            "to binding affinity and suggests the ligand occupies a well-defined cavity rather than "
            "a surface-exposed cleft."
        )
    elif h_count >= 4:
        hydrophobic_comment = (
            f"Moderate hydrophobic packing ({h_count} residues) provides nonpolar contacts "
            "that likely contribute to ligand orientation and shape complementarity."
        )
    else:
        hydrophobic_comment = (
            f"Limited hydrophobic contacts ({h_count} residues) suggest the binding site may be "
            "relatively solvent-exposed or dominated by polar interactions."
        )

    # Section 3: Polar/H-bond network
    polar_interactions = [x for x in interactions if "H-bond" in x.get("interaction_type", "") or "polar" in x.get("interaction_type", "")]
    if len(polar_interactions) >= 4:
        hbond_comment = (
            f"{len(polar_interactions)} potential hydrogen-bonding contacts were identified, "
            "suggesting a directional H-bond network that may confer binding specificity. "
            "H-bonds often discriminate between chemically similar ligands by enforcing precise "
            "geometric complementarity."
        )
    elif len(polar_interactions) >= 1:
        hbond_comment = (
            f"{len(polar_interactions)} potential H-bond contacts detected. "
            "Specific polar interactions, even when few in number, can be critical for ligand "
            "orientation and binding specificity."
        )
    else:
        hbond_comment = (
            "No strong polar/H-bond contacts were identified within the 4.0 Å cutoff, suggesting "
            "binding may rely primarily on shape complementarity and hydrophobic effects."
        )

    # Section 4: Charge interactions
    charged = [x for x in interactions if "charged" in x.get("interaction_type", "") or "electrostatic" in x.get("interaction_type", "")]
    if charged:
        charge_comment = (
            f"{len(charged)} charged/electrostatic contacts indicate long-range electrostatic "
            "steering may contribute to ligand recognition. These interactions are sensitive to "
            "pH and ionic strength, which could modulate binding under physiological conditions."
        )
    else:
        charge_comment = "No dominant charge-charge interactions detected at close range."

    # Section 5: Loop/flexibility assessment
    glycine_count = sum(1 for (_, res_name, _) in contact_residues if res_name == "GLY")
    proline_count = sum(1 for (_, res_name, _) in contact_residues if res_name == "PRO")
    if glycine_count >= 2 or proline_count >= 1:
        flex_comment = (
            f"Presence of Gly ({glycine_count}) and Pro ({proline_count}) in the pocket suggests "
            "local backbone flexibility may influence binding. Glycine-rich loops can undergo "
            "conformational changes upon ligand binding, while prolines may restrict backbone mobility."
        )
    else:
        flex_comment = (
            "The pocket shows typical backbone composition without strong flexibility signatures. "
            "This may indicate a relatively rigid binding site."
        )

    # Section 6: Binding insights
    if total_contacts >= 25:
        insight = (
            f"The extensive contact network ({total_contacts} residues within 5 Å) indicates a "
            "deep, well-defined binding cavity with high shape complementarity. This architecture "
            "is typical of high-affinity drug targets and enzyme active sites. The ligand is likely "
            "tightly bound with limited solvent accessibility."
        )
    elif total_contacts >= 15:
        insight = (
            f"A moderate contact network ({total_contacts} residues) suggests a semi-exposed "
            "binding site. The pocket provides sufficient encasement for ligand recognition while "
            "maintaining some solvent accessibility, characteristic of druggable binding sites."
        )
    else:
        insight = (
            f"A relatively small contact set ({total_contacts} residues) suggests a shallow or "
            "surface-exposed binding site. This may indicate a low-affinity interaction, a cofactor "
            "binding site, or a solvent-exposed allosteric pocket."
        )

    sections = [
        f"## Pocket Classification\n{pocket_type}\n{pocket_detail}",
        f"## Hydrophobic Packing\n{hydrophobic_comment}",
        f"## Polar & H-Bond Network\n{hbond_comment}",
        f"## Charge Interactions\n{charge_comment}",
        f"## Backbone Flexibility\n{flex_comment}",
        f"## Binding & Functional Insights\n{insight}",
    ]

    return "\n\n".join(sections)


def deepseek_ai_interpretation(ligand_name, contact_residues, counts, interactions):
    residue_list = [
        f"{chain_id}:{res_name}{res_id}"
        for chain_id, res_name, res_id in sorted(contact_residues.keys())
    ]

    interaction_list = [
        (
            f"{x['chain_id']}:{x['res_name']}{x['res_id']} "
            f"{x['atom_name']} {x['interaction_type']} {x['distance']} A"
        )
        for x in interactions
    ]

    prompt = f"""You are a structural biology AI scientist analyzing a protein-ligand binding pocket.

Ligand: {ligand_name}

Pocket residue composition (within 5 Å):
Hydrophobic: {counts["hydrophobic"]} | Polar: {counts["polar"]} | Positive: {counts["positive"]} | Negative: {counts["negative"]}

Contact residues: {residue_list}

Nearest ligand-residue contacts (≤4.0 Å, max 15):
{interaction_list}

Write a scientific interpretation with these 6 sections using Markdown headers (##):

## Pocket Classification
Classify the pocket type (hydrophobic, charged, polar, amphipathic, mixed) and explain what this implies about ligand preference and binding mode. Compare to known pocket archetypes (e.g. kinase hinge, GPCR orthosteric, enzyme active site).

## Hydrophobic Packing Analysis
Assess the hydrophobic packing quality. Is there a hydrophobic core? Are aromatic residues positioned for π-stacking? Comment on shape complementarity and potential desolvation effects. Note any suboptimal packing that could be improved.

## Polar & Hydrogen Bond Network
Analyze the H-bond donor/acceptor landscape. Which residues likely form directional H-bonds with the ligand? Are there unsatisfied H-bond partners that could be exploited for ligand design? Comment on the specificity conferred by the polar network.

## Electrostatic Environment
Describe the local electrostatic potential. Is the pocket complementary to the ligand's expected charge state? Would pH changes modulate binding? Note any long-range electrostatic effects.

## Backbone Flexibility & Dynamics
Assess likely flexibility based on Gly/Pro content and loop vs helix composition of contact residues. Would the pocket undergo induced fit upon binding? Are there glycine-rich loops that could accommodate diverse ligands?

## Functional & Druggability Assessment
Synthesize the above into a functional assessment. Is this pocket likely druggable? What ligand properties would optimize binding? What biological role might this site play (catalytic, allosteric, protein-protein interface)? End with one concrete suggestion for a follow-up experiment or computational analysis.

Rules:
- Write in the style of a Nature Structural Biology brief communication
- Be quantitative where possible (cite specific residue counts and distances)
- When data is insufficient for a conclusion, state so explicitly
- Use structural biology terminology correctly (e.g. π-cation, salt bridge, β-hairpin)
- Each section should be 2-4 sentences, substantive and insightful
- Keep the total response under 500 words
"""

    response = create_deepseek_client().chat.completions.create(
        model="deepseek-chat",
        messages=[
            {
                "role": "system",
                "content": (
                    "You are an expert structural biology AI copilot. "
                    "You analyze protein-ligand interactions with the depth of a seasoned structural biologist. "
                    "You write concisely, quantitatively, and insightfully. "
                    "You cite specific residue-level details. "
                    "When uncertain, you state limitations clearly rather than over-interpreting."
                )
            },
            {
                "role": "user",
                "content": prompt
            }
        ]
    )

    return response.choices[0].message.content


def generate_ai_interpretation(contact_residues, counts, interactions, ligand_name):
    try:
        return deepseek_ai_interpretation(
            ligand_name,
            contact_residues,
            counts,
            interactions
        )
    except Exception:
        fallback = local_ai_interpretation(
            contact_residues,
            counts,
            interactions,
            ligand_name
        )

        return (
            "[DeepSeek AI unavailable, using local rule-based interpretation]\n\n"
            + fallback
            + "\n\nAI service error details were hidden. Check server logs for debugging."
        )
