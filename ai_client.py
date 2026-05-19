import os

from openai import OpenAI


def create_deepseek_client():
    return OpenAI(
        api_key=os.getenv("DEEPSEEK_API_KEY"),
        base_url="https://api.deepseek.com"
    )


def local_ai_interpretation(contact_residues, counts, interactions, ligand_name):
    comments = []
    total_contacts = len(contact_residues)

    if counts["hydrophobic"] >= 8:
        comments.append(
            "The ligand is surrounded by a dense hydrophobic environment, suggesting nonpolar packing may strongly stabilize binding."
        )
    elif counts["hydrophobic"] >= 4:
        comments.append(
            "Several hydrophobic residues are present around the ligand, indicating hydrophobic packing may contribute to pocket stability."
        )

    if counts["polar"] >= 6:
        comments.append(
            "Multiple polar residues are located near the ligand, suggesting possible hydrogen-bond-mediated stabilization."
        )

    if counts["positive"] >= 4:
        comments.append(
            "A positively charged residue cluster is present near the ligand, which may support interaction with negatively charged groups or nucleic-acid-like substrates."
        )

    if counts["negative"] >= 4:
        comments.append(
            "The pocket contains multiple acidic residues, suggesting a locally negative electrostatic environment."
        )

    if total_contacts >= 20:
        comments.append(
            "The ligand forms an extensive contact network, indicating a geometrically constrained binding pocket."
        )

    if not comments:
        comments.append(
            "No dominant structural interaction pattern was detected from the current pocket features."
        )

    return "\n\n".join([f"- {x}" for x in comments])


def deepseek_ai_interpretation(ligand_name, contact_residues, counts, interactions):
    residue_list = [
        f"{chain_id}:{res_name}{res_id}"
        for chain_id, res_name, res_id in contact_residues.keys()
    ]

    interaction_list = [
        (
            f"{x['chain_id']}:{x['res_name']}{x['res_id']} "
            f"{x['atom_name']} {x['interaction_type']} {x['distance']} A"
        )
        for x in interactions
    ]

    prompt = f"""
You are a structural biology copilot.

Analyze this protein-ligand binding pocket based on extracted structural features.

Ligand: {ligand_name}

Pocket residue composition:
Hydrophobic residues: {counts["hydrophobic"]}
Polar residues: {counts["polar"]}
Positive residues: {counts["positive"]}
Negative residues: {counts["negative"]}

Contact residues:
{residue_list}

Closest ligand-residue contacts:
{interaction_list}

Please provide:
1. Pocket type classification
2. Structural interpretation
3. Possible biological implication
4. One short limitation note

Keep the answer concise, scientific, and useful for a structural biology report.
"""

    response = create_deepseek_client().chat.completions.create(
        model="deepseek-chat",
        messages=[
            {
                "role": "system",
                "content": "You are an expert structural biology assistant."
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
