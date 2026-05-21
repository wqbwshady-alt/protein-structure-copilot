import os
import traceback

from openai import OpenAI


def create_deepseek_client():
    return OpenAI(
        api_key=os.getenv("DEEPSEEK_API_KEY"),
        base_url="https://api.deepseek.com"
    )


def _ai_error_detail(exc):
    msg = str(exc)
    low = msg.lower()
    if "timeout" in low:
        return "DeepSeek API timed out. The server may be under heavy load."
    if "rate" in low or "429" in msg:
        return "DeepSeek API rate limit reached. Try again in a few seconds."
    if "auth" in low or "401" in msg or "403" in msg or "key" in low:
        return "DeepSeek API authentication failed. Check DEEPSEEK_API_KEY."
    if "connect" in low or "refused" in low or "network" in low:
        return "Could not reach DeepSeek API. Check your network connection."
    return f"DeepSeek API error: {msg[:200]}"


# ---- ligand-pocket helpers ----

def _top_residue_mentions(contact_residues, interactions, limit=5):
    names = []
    for ix in interactions[:limit]:
        names.append(f"{ix['res_name']}{ix['res_id']}")
    if not names:
        for key in list(contact_residues.keys())[:limit]:
            _, rn, ri = key
            names.append(f"{rn}{ri}")
    return ", ".join(names) if names else "(none)"


def _aromatic_mentions(contact_residues):
    aromatic = {"PHE", "TYR", "TRP", "HIS"}
    found = []
    for _, res_name, res_id in contact_residues:
        if res_name in aromatic:
            found.append(f"{res_name}{res_id}")
    return found


# ---- Local (rule-based) fallback for single ligand-pocket analysis ----

def local_ai_interpretation(contact_residues, counts, interactions, ligand_name):
    total_contacts = len(contact_residues)
    h_count = counts["hydrophobic"]
    p_count = counts["polar"]
    pos_count = counts["positive"]
    neg_count = counts["negative"]
    top_names = _top_residue_mentions(contact_residues, interactions)
    aromatics = _aromatic_mentions(contact_residues)

    # Section 1: Pocket type classification
    if h_count > max(p_count, pos_count, neg_count):
        pocket_type = "Primarily hydrophobic pocket"
        pocket_detail = (
            f"Hydrophobic residues dominate ({h_count} of {total_contacts} contacts, e.g. {top_names}), "
            "suggesting nonpolar packing is the main stabilization mechanism. "
            "This is characteristic of small-molecule binding sites in enzymes and receptors."
        )
    elif pos_count > neg_count and pos_count >= 3:
        pocket_type = "Positively charged pocket"
        pocket_detail = (
            f"Enrichment in basic residues (Lys/Arg/His: {pos_count}) suggests affinity for "
            "negatively charged ligands, phosphate groups, or nucleic acid substrates. "
            f"Key contacts include {top_names}."
        )
    elif neg_count > pos_count and neg_count >= 3:
        pocket_type = "Negatively charged pocket"
        pocket_detail = (
            f"Enrichment in acidic residues (Asp/Glu: {neg_count}) creates a locally negative "
            "electrostatic environment, potentially favoring cationic ligands or metal ion coordination. "
            f"Key contacts include {top_names}."
        )
    elif p_count >= h_count:
        pocket_type = "Polar/amphipathic pocket"
        pocket_detail = (
            f"Mixed polar ({p_count}) and hydrophobic ({h_count}) character suggests "
            "a versatile binding site capable of accommodating diverse ligand chemotypes. "
            f"Representative contacts: {top_names}."
        )
    else:
        pocket_type = "Mixed-composition pocket"
        pocket_detail = (
            f"Multiple residue types contribute to a composite binding interface. "
            f"Notable contacts: {top_names}."
        )

    # Section 2: Hydrophobic packing
    aromatic_count = len(aromatics)
    if h_count >= 8:
        hydrophobic_comment = (
            f"A dense hydrophobic core ({h_count} residues, e.g. {top_names}) forms extensive "
            "van der Waals contacts with the ligand. "
            + (
                f"Aromatic residues ({', '.join(aromatics[:5])}) may contribute pi-stacking "
                "interactions with the ligand ring system. "
                if aromatics else ""
            ) +
            "This degree of hydrophobic burial typically contributes significantly "
            "to binding affinity and suggests the ligand occupies a well-defined cavity."
        )
    elif h_count >= 4:
        hydrophobic_comment = (
            f"Moderate hydrophobic packing ({h_count} residues including {top_names}) provides "
            "nonpolar contacts that likely contribute to ligand orientation and shape complementarity."
            + (
                f" Aromatic sidechains ({', '.join(aromatics[:3])}) may form additional pi-interactions."
                if aromatics else ""
            )
        )
    else:
        hydrophobic_comment = (
            f"Limited hydrophobic contacts ({h_count} residues, mainly {top_names}) suggest "
            "the binding site may be relatively solvent-exposed or dominated by polar interactions."
        )

    # Section 3: Polar/H-bond network
    polar_interactions = [
        x for x in interactions
        if "H-bond" in x.get("interaction_type", "") or "polar" in x.get("interaction_type", "")
    ]
    polar_residues = sorted(set(
        f"{x['res_name']}{x['res_id']}" for x in polar_interactions
    ))
    if len(polar_interactions) >= 4:
        hbond_comment = (
            f"{len(polar_interactions)} potential H-bond contacts were identified "
            f"({', '.join(polar_residues[:6])}), suggesting a directional H-bond network "
            "that may confer binding specificity. H-bonds often discriminate between chemically "
            "similar ligands by enforcing precise geometric complementarity."
        )
    elif len(polar_interactions) >= 1:
        hbond_comment = (
            f"{len(polar_interactions)} potential H-bond contact(s) detected "
            f"({', '.join(polar_residues)}). "
            "Specific polar interactions, even when few in number, can be critical for ligand "
            "orientation and binding specificity."
        )
    else:
        hbond_comment = (
            "No strong polar/H-bond contacts were identified within the 4.0 A cutoff, suggesting "
            "binding may rely primarily on shape complementarity and hydrophobic effects."
        )

    # Section 4: Charge interactions
    charged = [
        x for x in interactions
        if "charged" in x.get("interaction_type", "") or "electrostatic" in x.get("interaction_type", "")
    ]
    charged_names = sorted(set(
        f"{x['res_name']}{x['res_id']}" for x in charged
    ))
    if charged:
        charge_comment = (
            f"{len(charged)} charged/electrostatic contact(s) "
            f"({', '.join(charged_names)}) indicate long-range electrostatic "
            "steering may contribute to ligand recognition. These interactions are sensitive to "
            "pH and ionic strength, which could modulate binding under physiological conditions."
        )
    else:
        charge_comment = "No dominant charge-charge interactions detected at close range."

    # Section 5: Loop/flexibility assessment
    glycine_count = sum(1 for (_, res_name, _) in contact_residues if res_name == "GLY")
    proline_count = sum(1 for (_, res_name, _) in contact_residues if res_name == "PRO")
    gly_positions = sorted(set(
        f"{r_id}" for _, rn, r_id in contact_residues if rn == "GLY"
    ))
    pro_positions = sorted(set(
        f"{r_id}" for _, rn, r_id in contact_residues if rn == "PRO"
    ))

    if glycine_count >= 2 or proline_count >= 1:
        parts = []
        if glycine_count:
            parts.append(f"Gly at positions {', '.join(gly_positions[:4])} ({glycine_count} total)")
        if proline_count:
            parts.append(f"Pro at positions {', '.join(pro_positions[:4])} ({proline_count} total)")
        flex_comment = (
            "; ".join(parts) + ". "
            "Glycine-rich loops can undergo conformational changes upon ligand binding, "
            "while prolines may restrict backbone mobility via kinked geometry."
        )
    else:
        flex_comment = (
            "The pocket shows typical backbone composition without strong flexibility signatures "
            "(no enrichment in Gly or Pro). This may indicate a relatively rigid binding site."
        )

    # Section 6: Binding insights
    vdw_count = sum(1 for x in interactions if "van der Waals" in x.get("interaction_type", ""))
    if total_contacts >= 25:
        insight = (
            f"The extensive contact network ({total_contacts} residues within 5 A, "
            f"closest: {top_names}) indicates a deep, well-defined binding cavity with high "
            "shape complementarity. This architecture is typical of high-affinity drug targets "
            "and enzyme active sites. The ligand is likely tightly bound."
        )
    elif total_contacts >= 15:
        insight = (
            f"A moderate contact network ({total_contacts} residues, key: {top_names}) "
            "suggests a semi-exposed binding site. The pocket provides sufficient encasement "
            "for ligand recognition while maintaining some solvent accessibility — "
            "characteristic of druggable binding sites."
        )
    else:
        insight = (
            f"A relatively small contact set ({total_contacts} residues, {top_names}) "
            "suggests a shallow or surface-exposed binding site. This may indicate a "
            "low-affinity interaction, a cofactor binding site, or a solvent-exposed "
            "allosteric pocket."
        )
    if vdw_count > 0:
        insight += (
            f" {vdw_count} van der Waals-only contacts fill out the remaining interaction surface."
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


# ---- DeepSeek: ligand-pocket analysis ----

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
            ligand_name, contact_residues, counts, interactions
        )
    except Exception as exc:
        traceback.print_exc()
        fallback = local_ai_interpretation(
            contact_residues, counts, interactions, ligand_name
        )
        return (
            f"[DeepSeek AI unavailable — {_ai_error_detail(exc)}. "
            "Using local rule-based interpretation.]\n\n"
            + fallback
        )


# ---- Mutation Scan AI ----

def deepseek_mutation_interpretation(mutation_result):
    m = mutation_result
    prompt = f"""You are a structural biology AI scientist analyzing a point mutation in a protein-ligand binding pocket.

Mutation: {m['mutation']}
Chain: {m['chain_id']}
Original residue: {m['original_residue']['res_name']}{m['original_residue']['res_id']}
Mutant residue: {m['mutant_residue']['res_name']}{m['mutant_residue']['res_id']}

Physicochemical property changes:
{m['property_changes']}

Original residue interactions with ligand:
{m['original_interactions']}

Possible interaction loss:
{m['interaction_impact']['possible_loss']}

Possible interaction gain:
{m['interaction_impact']['possible_gain']}

Write a scientific interpretation with these 4 sections using Markdown headers (##):

## Mutation Overview
Summarize the mutation in one sentence. What is the chemical nature of the change (e.g. charged→hydrophobic, small→large, aliphatic→aromatic)?

## Binding Impact Assessment
Based on the original interactions and property changes, assess whether binding is likely weakened, disrupted, or tolerated. Cite specific interaction types at risk. If the residue does not directly contact the ligand, discuss potential allosteric or stability effects.

## Structural Context
Speculate on how this position fits into the larger pocket architecture. Is it a core contact, peripheral, or solvent-exposed? What secondary structure element likely hosts it?

## Recommendation
One concrete follow-up: MD simulation, alanine scanning, SPR binding assay, or a specific complementary mutation to test. End with a caution that this is a heuristic assessment, not a validated ΔΔG prediction.

Rules:
- 2-4 sentences per section
- Cite specific numbers and interaction types
- Be honest about uncertainty
- Use structural biology terminology
- Total under 400 words
"""

    response = create_deepseek_client().chat.completions.create(
        model="deepseek-chat",
        messages=[
            {
                "role": "system",
                "content": (
                    "You are an expert structural biology AI copilot specializing in "
                    "mutational analysis of protein-ligand interfaces."
                )
            },
            {"role": "user", "content": prompt}
        ]
    )

    return response.choices[0].message.content


def generate_mutation_interpretation(mutation_result):
    try:
        return deepseek_mutation_interpretation(mutation_result)
    except Exception as exc:
        traceback.print_exc()
        from services.mutation_scan import build_mutation_interpretation
        fallback = build_mutation_interpretation(mutation_result)
        return (
            f"[DeepSeek AI unavailable — {_ai_error_detail(exc)}. "
            "Using local rule-based interpretation.]\n\n"
            + fallback
        )


# ---- WT vs Mutant Comparison AI ----

def deepseek_comparison_interpretation(ligand_name, wt_contacts, mut_contacts, wt_counts, mut_counts, lost, gained, shared):
    lost_str = ", ".join(f"{c}:{r}{i}" for c, r, i in lost[:10]) or "none"
    gained_str = ", ".join(f"{c}:{r}{i}" for c, r, i in gained[:10]) or "none"
    shared_str = ", ".join(f"{c}:{r}{i}" for c, r, i in shared[:10]) or "none"

    prompt = f"""You are a structural biology AI scientist comparing WT vs mutant protein-ligand binding pockets.

Ligand: {ligand_name}

WT pocket: {len(wt_contacts)} contact residues (within 5 Å)
Mutant pocket: {len(mut_contacts)} contact residues (within 5 Å)
Shared contacts: {len(shared)}
Lost contacts: {len(lost)}
Gained contacts: {len(gained)}

WT residue composition — Hydrophobic: {wt_counts['hydrophobic']}, Polar: {wt_counts['polar']}, Positive: {wt_counts['positive']}, Negative: {wt_counts['negative']}
Mutant residue composition — Hydrophobic: {mut_counts['hydrophobic']}, Polar: {mut_counts['polar']}, Positive: {mut_counts['positive']}, Negative: {mut_counts['negative']}

Lost contact residues: {lost_str}
Gained contact residues: {gained_str}
Shared contact residues (first 10): {shared_str}

Write a scientific interpretation with these 4 sections using Markdown headers (##):

## Pocket Remodeling Summary
Quantify the net change in one sentence. Is the pocket expanding, contracting, or being remodeled? What is the net effect on residue composition?

## Key Contact Changes
Highlight the most significant lost and gained contacts. Which specific residue changes are likely to alter binding the most? Mention residue identity and likely interaction type.

## Binding Affinity Implications
Based on the pattern of contact changes, predict whether affinity is likely increased, decreased, or maintained. Consider both the number and quality (type) of contacts. What does the hydrophobic/charged balance shift suggest?

## Experimental Recommendation
One concrete suggestion: which biophysical assay (SPR, ITC, DSF) or computational method (MM-GBSA, alchemical FEP) would best validate this comparison? End with a disclaimer about heuristic limitations.

Rules:
- 2-4 sentences per section
- Cite specific numbers
- Be honest when data is insufficient
- Total under 400 words
"""

    response = create_deepseek_client().chat.completions.create(
        model="deepseek-chat",
        messages=[
            {
                "role": "system",
                "content": (
                    "You are an expert structural biology AI copilot specializing in "
                    "comparative analysis of protein-ligand interfaces."
                )
            },
            {"role": "user", "content": prompt}
        ]
    )

    return response.choices[0].message.content


def generate_comparison_interpretation(ligand_name, wt_contacts, mut_contacts, wt_counts, mut_counts, lost, gained, shared):
    try:
        return deepseek_comparison_interpretation(
            ligand_name, wt_contacts, mut_contacts, wt_counts, mut_counts, lost, gained, shared
        )
    except Exception as exc:
        traceback.print_exc()
        lost_count = len(lost)
        gained_count = len(gained)
        shared_count = len(shared)

        if lost_count == 0 and gained_count == 0:
            fallback = "No major pocket residue change detected. The mutant pocket appears identical to WT within the 5 A contact cutoff."
        elif lost_count > gained_count:
            fallback = f"The mutant lost {lost_count} contact residue(s) and gained {gained_count}, suggesting possible pocket weakening or contraction. {shared_count} contacts are preserved."
        elif gained_count > lost_count:
            fallback = f"The mutant gained {gained_count} contact residue(s) while losing {lost_count}, suggesting pocket expansion or reorganization. {shared_count} contacts are preserved."
        else:
            fallback = f"The mutant shows balanced contact changes (lost {lost_count}, gained {gained_count}) with {shared_count} shared contacts, indicating local remodeling rather than net gain or loss."

        return (
            f"[DeepSeek AI unavailable — {_ai_error_detail(exc)}. "
            "Using local rule-based interpretation.]\n\n"
            + fallback
        )


# ---- Protein-Only AI ----

def deepseek_protein_only_interpretation(summary):
    chains = summary.get("chains", [])
    chains_str = "\n".join(
        f"Chain {c['chain_id']}: {c['residue_count']} residues, {c['atom_count']} atoms, "
        f"hydrophobic={c['hydrophobic']}, polar={c['polar']}, positive={c['positive']}, negative={c['negative']}"
        for c in chains
    )
    ligand_candidates = summary.get("ligand_candidates", [])

    prompt = f"""You are a structural biology AI scientist analyzing a protein structure.

Structure overview:
{summary['chain_count']} chain(s), {summary['residue_count']} protein residues, {summary['atom_count']} protein atoms

Chain composition:
{chains_str}

Residue chemistry — Hydrophobic: {summary['residue_type_counts']['hydrophobic']}, Polar: {summary['residue_type_counts']['polar']}, Positive: {summary['residue_type_counts']['positive']}, Negative: {summary['residue_type_counts']['negative']}, Other: {summary['residue_type_counts']['other']}

Detected non-water HETATM ligands: {ligand_candidates if ligand_candidates else 'None detected'}

Write a scientific structural overview with these 4 sections using Markdown headers (##):

## Architecture Summary
Describe the overall architecture in 2 sentences. Multi-domain? Multi-chain complex? What class of protein might this be based on residue composition?

## Fold & Stability Assessment
Based on hydrophobic vs polar residue distribution, comment on likely fold class and stability. Is the hydrophobic core proportion typical for globular proteins? Any signs of intrinsic disorder?

## Functional Implications
If ligands are detected, what might their role be? If no ligands, what does the surface residue composition suggest about possible binding interfaces or catalytic sites?

## Research Context
One sentence on what experimental or computational follow-up would be valuable: cryo-EM, HDX-MS, docking screen, or homology modeling. End with a note that this is a sequence/composition-based overview without 3D structural context.

Rules:
- 2-3 sentences per section
- Cite specific numbers
- Be honest about limitations
- Total under 300 words
"""

    response = create_deepseek_client().chat.completions.create(
        model="deepseek-chat",
        messages=[
            {
                "role": "system",
                "content": (
                    "You are an expert structural biology AI copilot. "
                    "You provide concise, quantitative structural overviews."
                )
            },
            {"role": "user", "content": prompt}
        ]
    )

    return response.choices[0].message.content


def generate_protein_only_interpretation(summary):
    try:
        return deepseek_protein_only_interpretation(summary)
    except Exception as exc:
        traceback.print_exc()
        counts = summary.get("residue_type_counts", {})
        chains = summary.get("chains", [])
        ligands = summary.get("ligand_candidates", [])

        fallback = (
            f"Protein-only structural overview.\n\n"
            f"This structure contains {summary.get('chain_count', 0)} chain(s) with "
            f"{summary.get('residue_count', 0)} residues and {summary.get('atom_count', 0)} atoms. "
            f"Residue composition: {counts.get('hydrophobic', 0)} hydrophobic, "
            f"{counts.get('polar', 0)} polar, {counts.get('positive', 0)} positively charged, "
            f"{counts.get('negative', 0)} negatively charged, {counts.get('other', 0)} other.\n\n"
        )

        if chains:
            fallback += "Chain breakdown:\n"
            for c in chains[:5]:
                fallback += (
                    f"- Chain {c.get('chain_id', '?')}: {c.get('residue_count', 0)} residues, "
                    f"hydrophobic={c.get('hydrophobic', 0)}, polar={c.get('polar', 0)}, "
                    f"charged={c.get('positive', 0) + c.get('negative', 0)}\n"
                )

        if ligands:
            fallback += f"\nDetected non-water ligands: {', '.join(ligands[:8])}. "
            fallback += "These may represent cofactors, substrates, or crystallization additives."

        fallback += (
            "\n\nLigand pocket contact analysis was not performed. "
            "To analyze a specific binding pocket, enter a ligand name and re-run the analysis."
        )

        return (
            f"[DeepSeek AI unavailable — {_ai_error_detail(exc)}. "
            "Using local rule-based interpretation.]\n\n"
            + fallback
        )


# ---- Structured 7-Section AI Interpretation (new, coexists with legacy functions) ----

SYSTEM_PROMPT = (
    "You are a structural biology research assistant. "
    "Your role is evidence-based interpretation of protein-ligand binding data, not speculation. "
    "You write in the style of a scientific research report: concise, cautious, professional. "
    "Every claim must be traceable to specific structural evidence provided in the data. "
    "When evidence is absent, you explicitly state so. "
    "You do not invent protein functions, disease mechanisms, drug names, or literature conclusions."
)

_STRUCTURED_SECTIONS = [
    ("A. Executive Summary", "2-3 sentence overview: pocket type, key residues, dominant interaction modes, confidence level"),
    ("B. Pocket Physicochemical Profile", "hydrophobic/polar/charged residue composition, what the chemical environment implies for ligand binding"),
    ("C. Key Interaction Network", "which residues are likely key anchors, which interaction types dominate, what each type suggests about stabilization mechanism"),
    ("D. Residue-Level Evidence", "list each key contact with: chain ID, residue name+number, atom name, interaction type, distance (A). Format each as bullet point"),
    ("E. Mutation Impact Assessment", "if mutation data present: confidence level LOW/MEDIUM/HIGH with reasoning; if no mutation data: explicitly state N/A"),
    ("F. Biological Interpretation", "what the pocket architecture suggests about binding mode; use cautious language: suggests/may indicate/is consistent with"),
    ("G. Limitations", "list inherent limitations: static crystal structure, no MD/energy calculation, 5 A cutoff, geometric interaction classification only, no solvent/ions/pH effects"),
]

_STRUCTURED_RULES = (
    "CRITICAL RULES — every response must follow these:\n"
    "1. Every factual claim MUST cite: residue name, chain ID, residue number, interaction type, distance (A)\n"
    "2. If structural data cannot support a claim, write exactly: Insufficient structural evidence to determine this.\n"
    "3. Do NOT mention: specific protein function, disease mechanism, known drug names, literature references, biological pathways\n"
    "4. Use cautious language: suggests, may indicate, is consistent with, possibly, potentially, appears to\n"
    "5. For Section E: state confidence level LOW / MEDIUM / HIGH with explicit reasoning based on available data\n"
    "6. Each section: 2-5 sentences, concise and substantive\n"
    "7. Use correct structural biology terminology: salt bridge, pi-stacking, H-bond network, hydrophobic core, van der Waals packing\n"
    "8. Total response under 700 words"
)


def parse_ai_sections(text):
    """Split AI response into an ordered dict keyed by section letter+title."""
    from collections import OrderedDict

    if not text or not isinstance(text, str):
        return OrderedDict()

    sections = OrderedDict()
    current_key = None
    current_body = []

    for line in text.split("\n"):
        stripped = line.strip()
        if stripped.startswith("## ") and len(stripped) > 3:
            if current_key:
                sections[current_key] = "<br>".join(current_body) if current_body else ""
            current_key = stripped[3:].strip()
            current_body = []
        elif current_key:
            current_body.append(stripped)

    if current_key:
        sections[current_key] = "<br>".join(current_body) if current_body else ""

    if not sections:
        sections["A. Executive Summary"] = text.replace("\n", "<br>")

    return sections


def build_structured_prompt(mode, data):
    """Build a DeepSeek prompt for the 7-section structured interpretation."""

    section_descriptions = "\n".join(
        f"## {title}\n{desc}" for title, desc in _STRUCTURED_SECTIONS
    )

    if mode == "ligand":
        ligand_name = data.get("ligand_name", "unknown")
        counts = data.get("counts", {})
        interactions = data.get("interactions", [])

        residue_list = data.get("residue_list") or []
        if not residue_list:
            contact_residues = data.get("contact_residues")
            if contact_residues and hasattr(contact_residues, "keys"):
                residue_list = sorted(
                    f"{chain}:{rn}{ri}"
                    for chain, rn, ri in contact_residues.keys()
                )

        interaction_lines = []
        for x in (interactions or [])[:15]:
            interaction_lines.append(
                f"Chain {x.get('chain_id','?')}: {x.get('res_name','?')}{x.get('res_id','?')} "
                f"atom {x.get('atom_name','?')} | {x.get('interaction_type','contact')} | "
                f"{x.get('distance','?')} A | ligand atom {x.get('ligand_atom','?')}"
            )

        data_block = f"""Ligand: {ligand_name}

Pocket residue composition (within 5 A):
Hydrophobic: {counts.get('hydrophobic', 0)} | Polar: {counts.get('polar', 0)} | Positive: {counts.get('positive', 0)} | Negative: {counts.get('negative', 0)}

Contact residues: {residue_list}

Nearest ligand-residue contacts (<=4.0 A, max 15):
{chr(10).join(interaction_lines) if interaction_lines else 'No close contacts detected within 4.0 A.'}"""

    elif mode == "mutation":
        m = data.get("mutation_result", {})
        original = m.get("original_residue", {})
        mutant = m.get("mutant_residue", {})
        property_changes = m.get("property_changes", {})
        impact = m.get("interaction_impact", {})
        original_interactions = m.get("original_interactions", [])

        interaction_lines = []
        for x in original_interactions[:15]:
            interaction_lines.append(
                f"Chain {x.get('chain_id','?')}: {x.get('res_name','?')}{x.get('res_id','?')} "
                f"atom {x.get('atom_name','?')} | {x.get('interaction_type','contact')} | "
                f"{x.get('distance','?')} A"
            )

        loss_lines = []
        for x in impact.get("possible_loss", []):
            loss_lines.append(f"- {x.get('interaction_type','?')}: {x.get('reason','?')}")
        gain_lines = []
        for x in impact.get("possible_gain", []):
            gain_lines.append(f"- {x.get('interaction_type','?')}: {x.get('reason','?')}")

        data_block = f"""Mutation: {m.get('mutation','?')}
Chain: {m.get('chain_id','?')}
Original residue: {original.get('res_name','?')}{original.get('res_id','?')}
Mutant residue: {mutant.get('res_name','?')}{mutant.get('res_id','?')}

Physicochemical property changes:
{chr(10).join(f'- {k}: {v}' for k,v in (property_changes or {}).items())}

Original residue interactions with ligand:
{chr(10).join(interaction_lines) if interaction_lines else 'No direct ligand interactions detected for this residue.'}

Possible interaction loss:
{chr(10).join(loss_lines) if loss_lines else 'None predicted by current rules.'}

Possible interaction gain:
{chr(10).join(gain_lines) if gain_lines else 'None predicted by current rules.'}"""

    elif mode == "comparison":
        ligand_name = data.get("ligand_name", "unknown")
        wt_counts = data.get("wt_counts", {})
        mut_counts = data.get("mut_counts", {})
        wt_contact_count = data.get("wt_contact_count", 0)
        mut_contact_count = data.get("mut_contact_count", 0)
        lost = data.get("lost", [])
        gained = data.get("gained", [])
        shared = data.get("shared", [])

        lost_str = ", ".join(f"{c}:{r}{i}" for c, r, i in (lost or [])[:10]) or "none"
        gained_str = ", ".join(f"{c}:{r}{i}" for c, r, i in (gained or [])[:10]) or "none"
        shared_str = ", ".join(f"{c}:{r}{i}" for c, r, i in (shared or [])[:10]) or "none"

        data_block = f"""Ligand: {ligand_name}

WT pocket: {wt_contact_count} contact residues (within 5 A)
Mutant pocket: {mut_contact_count} contact residues (within 5 A)
Shared contacts: {len(shared or [])}
Lost contacts: {len(lost or [])}
Gained contacts: {len(gained or [])}

WT residue composition — Hydrophobic: {wt_counts.get('hydrophobic',0)}, Polar: {wt_counts.get('polar',0)}, Positive: {wt_counts.get('positive',0)}, Negative: {wt_counts.get('negative',0)}
Mutant residue composition — Hydrophobic: {mut_counts.get('hydrophobic',0)}, Polar: {mut_counts.get('polar',0)}, Positive: {mut_counts.get('positive',0)}, Negative: {mut_counts.get('negative',0)}

Lost contact residues: {lost_str}
Gained contact residues: {gained_str}
Shared contact residues (first 10): {shared_str}"""

    elif mode == "protein_only":
        summary = data.get("summary", {})
        chains = summary.get("chains", [])
        counts = summary.get("residue_type_counts", {})
        ligands = summary.get("ligand_candidates", [])

        chains_str = "\n".join(
            f"Chain {c.get('chain_id','?')}: {c.get('residue_count',0)} residues, "
            f"hydrophobic={c.get('hydrophobic',0)}, polar={c.get('polar',0)}, "
            f"charged={c.get('positive',0)+c.get('negative',0)}"
            for c in (chains or [])
        )

        data_block = f"""Structure overview:
{summary.get('chain_count',0)} chain(s), {summary.get('residue_count',0)} protein residues, {summary.get('atom_count',0)} protein atoms

Chain composition:
{chains_str if chains_str else 'No chain data available.'}

Residue chemistry — Hydrophobic: {counts.get('hydrophobic',0)}, Polar: {counts.get('polar',0)}, Positive: {counts.get('positive',0)}, Negative: {counts.get('negative',0)}, Other: {counts.get('other',0)}

Detected non-water HETATM ligands: {ligands if ligands else 'None detected'}

Note: This is a protein-only analysis. No ligand pocket or mutation data is available.
Sections C, D, E should explicitly state this limitation where applicable."""

    else:
        data_block = f"Unknown mode: {mode}\nData: {data}"

    return f"""Write a scientific interpretation with these EXACT 7 sections using Markdown headers (## ):

{section_descriptions}

{_STRUCTURED_RULES}

---ANALYSIS DATA---
{data_block}"""


def _structured_local_fallback(mode, data):
    """Rule-based 7-section interpretation when DeepSeek is unavailable."""

    def _section(title, body):
        return f"## {title}\n{body}"

    if mode == "ligand":
        ligand_name = data.get("ligand_name", "unknown")
        counts = data.get("counts", {})
        interactions = data.get("interactions", [])
        contact_residues = data.get("contact_residues", {})

        total = len(contact_residues) if isinstance(contact_residues, dict) else 0
        h = counts.get("hydrophobic", 0)
        p = counts.get("polar", 0)
        pos = counts.get("positive", 0)
        neg = counts.get("negative", 0)

        top_names = ", ".join(
            f"{x.get('res_name','')}{x.get('res_id','')}"
            for x in (interactions or [])[:5]
        ) or "(no close contacts)"

        charged = [x for x in (interactions or []) if "charged" in x.get("interaction_type", "")]
        polar_ints = [x for x in (interactions or []) if "H-bond" in x.get("interaction_type", "") or "polar" in x.get("interaction_type", "")]
        vdw = [x for x in (interactions or []) if "van der Waals" in x.get("interaction_type", "")]

        if h > max(p, pos, neg):
            pocket_type = f"Primarily hydrophobic pocket ({h}/{total} contacts hydrophobic)"
        elif pos > neg and pos >= 3:
            pocket_type = f"Positively charged pocket (Lys/Arg/His: {pos})"
        elif neg > pos and neg >= 3:
            pocket_type = f"Negatively charged pocket (Asp/Glu: {neg})"
        else:
            pocket_type = "Mixed-composition pocket"

        sections = [
            _section("A. Executive Summary",
                f"{pocket_type}. {total} contact residues within 5 A of {ligand_name}. "
                f"{len(polar_ints)} polar/H-bond contacts, {len(charged)} charged contacts, "
                f"{len(vdw)} van der Waals contacts. "
                f"Top contacts: {top_names}."
            ),
            _section("B. Pocket Physicochemical Profile",
                f"Hydrophobic: {h} residues | Polar: {p} | Positively charged: {pos} | Negatively charged: {neg}. "
                + (f"A hydrophobic-dominant pocket suggests nonpolar packing is the main stabilization mechanism. "
                   f"This is characteristic of small-molecule binding sites."
                   if h > max(p, pos, neg) else
                   f"A charged/polar pocket environment suggests electrostatic complementarity "
                   f"may play a significant role in ligand recognition.")
            ),
            _section("C. Key Interaction Network",
                f"{len(polar_ints)} potential H-bond/polar contacts identified. "
                + (f"Key polar contacts: {', '.join(x['res_name'] + x['res_id'] for x in polar_ints[:5])}. "
                   if polar_ints else "No strong polar contacts detected. ")
                + f"{len(charged)} charged/electrostatic contacts. "
                + (f"Charged contacts involve: {', '.join(x['res_name'] + x['res_id'] for x in charged[:5])}. "
                   if charged else "No charge-charge interactions at close range. ")
                + f"{len(vdw)} van der Waals contacts fill out the remaining interaction surface."
            ),
            _section("D. Residue-Level Evidence",
                "\n".join(
                    f"- Chain {x.get('chain_id','?')}: {x.get('res_name','?')}{x.get('res_id','?')} "
                    f"({x.get('atom_name','?')}) — {x.get('interaction_type','contact')} — "
                    f"{x.get('distance','?')} A (ligand atom: {x.get('ligand_atom','?')})"
                    for x in (interactions or [])[:10]
                ) or "Insufficient structural evidence — no interaction data available within cutoff distance."
            ),
            _section("E. Mutation Impact Assessment",
                "N/A — single-structure ligand pocket analysis. No mutation data was provided. "
                "To assess mutation impact, use the Mutation Scan mode with a specific point mutation."
            ),
            _section("F. Biological Interpretation",
                f"The {ligand_name} binding pocket contains {total} contact residues within 5 A. "
                + (f"The extensive contact network suggests a well-defined binding cavity with high shape complementarity. "
                   f"This architecture is consistent with a high-affinity binding site."
                   if total >= 20 else
                   f"A moderate contact network ({total} residues) suggests a semi-exposed binding site "
                   f"that may accommodate ligand recognition while maintaining some solvent accessibility."
                   if total >= 10 else
                   f"A small contact set ({total} residues) suggests a shallow or surface-exposed binding site, "
                   f"possibly a low-affinity interaction or cofactor binding region.")
            ),
            _section("G. Limitations",
                "This interpretation is based on a static crystal structure. "
                "No molecular dynamics simulation, binding free energy calculation, or quantum mechanical analysis was performed. "
                "Interaction classification uses geometric criteria only (5 A contact cutoff, 4.0 A close-contact threshold). "
                "Solvent effects, protonation states, ionic strength, and protein dynamics are not accounted for. "
                "This is a heuristic structural analysis, not a validated binding prediction. "
                "For quantitative assessment, consider MD simulation, MM-GBSA, or experimental binding assays."
            ),
        ]
        return "\n\n".join(sections)

    if mode == "mutation":
        m = data.get("mutation_result", {})
        original = m.get("original_residue", {})
        mutant = m.get("mutant_residue", {})
        property_changes = m.get("property_changes", {})
        impact = m.get("interaction_impact", {})
        loss = impact.get("possible_loss", [])
        gain = impact.get("possible_gain", [])

        loss_count = len(loss)
        gain_count = len(gain)

        if loss_count > 1:
            confidence = "HIGH"
            conf_reason = "multiple interaction types at risk of disruption"
        elif loss_count == 1:
            confidence = "MEDIUM"
            conf_reason = "one interaction type at risk, but geometry may partially compensate"
        elif gain_count > 0:
            confidence = "MEDIUM"
            conf_reason = "possible gain of new interactions, but depends on sidechain geometry"
        else:
            confidence = "LOW"
            conf_reason = "no strong evidence of interaction disruption based on geometric analysis alone"

        sections = [
            _section("A. Executive Summary",
                f"Mutation {m.get('mutation','?')} on chain {m.get('chain_id','?')} "
                f"replaces {original.get('res_name','?')}{original.get('res_id','?')} with "
                f"{mutant.get('res_name','?')}{mutant.get('res_id','?')}. "
                f"Confidence: {confidence} — {conf_reason}."
            ),
            _section("B. Pocket Physicochemical Profile",
                "Property changes: " +
                "; ".join(f"{k}: {v}" for k, v in (property_changes or {}).items()) +
                ". " +
                (f"The {original.get('res_name','?')}->{mutant.get('res_name','?')} substitution significantly alters "
                 f"the local physicochemical environment."
                 if loss_count > 0 or gain_count > 0 else
                 "The substitution does not strongly alter the local physicochemical profile based on the property analysis.")
            ),
            _section("C. Key Interaction Network",
                f"{loss_count} possible interaction loss(es), {gain_count} possible interaction gain(s). "
                + (f"Loss: {'; '.join(x.get('reason','?') for x in loss)}. " if loss else "")
                + (f"Gain: {'; '.join(x.get('reason','?') for x in gain)}. " if gain else "")
            ),
            _section("D. Residue-Level Evidence",
                "\n".join(
                    f"- Chain {x.get('chain_id','?')}: {x.get('res_name','?')}{x.get('res_id','?')} "
                    f"({x.get('atom_name','?')}) — {x.get('interaction_type','?')} — "
                    f"{x.get('distance','?')} A — {x.get('reason','?')}"
                    for x in loss
                ) if loss else (
                    "\n".join(
                        f"- Chain {x.get('chain_id','?')}: {x.get('res_name','?')}{x.get('res_id','?')} "
                        f"({x.get('atom_name','?')}) — {x.get('interaction_type','?')} — "
                        f"{x.get('distance','?')} A"
                        for x in m.get("original_interactions", [])[:10]
                    ) if m.get("original_interactions") else "Insufficient structural evidence — no direct ligand interactions detected for this residue."
                )
            ),
            _section("E. Mutation Impact Assessment",
                f"Confidence: {confidence}. {conf_reason}. "
                "This assessment is based on geometric interaction classification and physicochemical property comparison. "
                "It does not account for sidechain remodeling, backbone rearrangement, or changes in solvation."
            ),
            _section("F. Biological Interpretation",
                f"The {m.get('mutation','?')} mutation "
                + ("may significantly weaken ligand binding due to loss of key interactions with the pocket."
                   if confidence == "HIGH" else
                   "may moderately affect ligand binding, though some interactions may be preserved or compensated."
                   if confidence == "MEDIUM" else
                   "is unlikely to strongly disrupt ligand binding based on the current structural analysis.")
            ),
            _section("G. Limitations",
                "This is a heuristic mutation impact assessment. "
                "No mutant sidechain remodeling, molecular dynamics, or free energy calculation was performed. "
                "The analysis compares physicochemical properties and geometric interaction patterns only. "
                "Solvent effects, protein dynamics, and allosteric effects are not modeled. "
                "For quantitative prediction, consider alchemical free energy perturbation (FEP) or MD-based MM-GBSA."
            ),
        ]
        return "\n\n".join(sections)

    if mode == "comparison":
        ligand_name = data.get("ligand_name", "unknown")
        wt_cc = data.get("wt_contact_count", 0)
        mut_cc = data.get("mut_contact_count", 0)
        lost = data.get("lost", [])
        gained = data.get("gained", [])
        shared = data.get("shared", [])
        wt_counts = data.get("wt_counts", {})
        mut_counts = data.get("mut_counts", {})

        lost_count = len(lost or [])
        gained_count = len(gained or [])
        shared_count = len(shared or [])

        h_change = mut_counts.get("hydrophobic", 0) - wt_counts.get("hydrophobic", 0)

        if lost_count > gained_count:
            direction = "contraction / weakening"
        elif gained_count > lost_count:
            direction = "expansion / reorganization"
        else:
            direction = "local remodeling"

        sections = [
            _section("A. Executive Summary",
                f"WT vs mutant comparison for {ligand_name} binding pocket. "
                f"WT: {wt_cc} contacts, Mutant: {mut_cc} contacts. "
                f"Lost: {lost_count}, Gained: {gained_count}, Shared: {shared_count}. "
                f"Net effect: {direction}."
            ),
            _section("B. Pocket Physicochemical Profile",
                f"WT — Hydrophobic: {wt_counts.get('hydrophobic',0)}, Polar: {wt_counts.get('polar',0)}, "
                f"Positive: {wt_counts.get('positive',0)}, Negative: {wt_counts.get('negative',0)}. "
                f"Mutant — Hydrophobic: {mut_counts.get('hydrophobic',0)}, Polar: {mut_counts.get('polar',0)}, "
                f"Positive: {mut_counts.get('positive',0)}, Negative: {mut_counts.get('negative',0)}. "
                f"Hydrophobic change: {h_change:+d}"
                + (", indicating increased hydrophobic character."
                   if h_change > 0 else ", indicating decreased hydrophobic character."
                   if h_change < 0 else ", no net hydrophobic change.")
            ),
            _section("C. Key Interaction Network",
                f"{shared_count} contacts preserved between WT and mutant. "
                + (f"Lost contacts: {', '.join(f'{c}:{r}{i}' for c,r,i in (lost or [])[:10])}. "
                   if lost_count else "No contacts lost. ")
                + (f"Gained contacts: {', '.join(f'{c}:{r}{i}' for c,r,i in (gained or [])[:10])}. "
                   if gained_count else "No new contacts gained. ")
            ),
            _section("D. Residue-Level Evidence",
                "\n".join(
                    [f"Lost residues:"] +
                    [f"- Chain {c}: {r}{i}" for c, r, i in (lost or [])[:10]] +
                    [f"Gained residues:"] +
                    [f"- Chain {c}: {r}{i}" for c, r, i in (gained or [])[:10]] +
                    [f"Shared residues (first 10):"] +
                    [f"- Chain {c}: {r}{i}" for c, r, i in (shared or [])[:10]]
                ) if (lost_count or gained_count or shared_count) else "Insufficient structural evidence — no residue-level differences detected."
            ),
            _section("E. Mutation Impact Assessment",
                ("HIGH" if abs(lost_count - gained_count) >= 3 else
                 "MEDIUM" if abs(lost_count - gained_count) >= 1 else "LOW") +
                f" — the mutant shows {direction} of the {ligand_name} binding pocket "
                f"with {lost_count} contacts lost and {gained_count} gained. "
                "This assessment is based on contact comparison within 5 A cutoff only."
            ),
            _section("F. Biological Interpretation",
                f"The mutant pocket {direction} suggests "
                + ("possible weakening of ligand binding affinity. The lost contacts may reduce shape complementarity and interaction energy."
                   if lost_count > gained_count else
                   "possible expansion of the binding site. The gained contacts may create new interaction opportunities for the ligand."
                   if gained_count > lost_count else
                   "the overall binding geometry is maintained, with local residue-level remodeling rather than net gain or loss of contacts.")
            ),
            _section("G. Limitations",
                "Comparison based on static crystal structures. No molecular dynamics or energy minimization was performed. "
                "Contact comparison uses 5 A cutoff; subtle conformational changes beyond this radius are not captured. "
                "Solvent, pH, and dynamics effects are not modeled. "
                "For quantitative comparison, consider MD simulation or alchemical free energy calculations."
            ),
        ]
        return "\n\n".join(sections)

    if mode == "protein_only":
        summary = data.get("summary", {})
        counts = summary.get("residue_type_counts", {})
        chains = summary.get("chains", [])
        ligands = summary.get("ligand_candidates", [])

        sections = [
            _section("A. Executive Summary",
                f"Protein-only structural overview. {summary.get('chain_count',0)} chain(s), "
                f"{summary.get('residue_count',0)} residues, {summary.get('atom_count',0)} atoms. "
                + (f"Ligands detected: {', '.join(ligands[:8])}."
                   if ligands else "No non-water ligands detected.")
            ),
            _section("B. Pocket Physicochemical Profile",
                f"Residue composition — Hydrophobic: {counts.get('hydrophobic',0)}, "
                f"Polar: {counts.get('polar',0)}, Positive: {counts.get('positive',0)}, "
                f"Negative: {counts.get('negative',0)}, Other: {counts.get('other',0)}. "
                + ("The hydrophobic/polar ratio is consistent with a typical globular protein fold."
                   if counts.get('hydrophobic', 0) > counts.get('polar', 0) else
                   "The elevated polar content may indicate significant solvent exposure or a membrane/interface region.")
            ),
            _section("C. Key Interaction Network",
                "N/A — protein-only analysis. No ligand was specified for pocket interaction analysis. "
                + (f"Detected HETATM ligands ({', '.join(ligands[:8])}) may represent cofactors, substrates, "
                   f"or crystallization additives. To analyze a binding pocket, specify a ligand name and re-run."
                   if ligands else "To analyze a specific binding pocket, specify a ligand name and re-run the analysis.")
            ),
            _section("D. Residue-Level Evidence",
                "\n".join(
                    f"- Chain {c.get('chain_id','?')}: {c.get('residue_count',0)} residues, "
                    f"{c.get('atom_count',0)} atoms, hydrophobic={c.get('hydrophobic',0)}, "
                    f"polar={c.get('polar',0)}, charged={c.get('positive',0)+c.get('negative',0)}"
                    for c in (chains or [])[:10]
                ) if chains else "Insufficient structural evidence — no chain-level data available."
            ),
            _section("E. Mutation Impact Assessment",
                "N/A — protein-only analysis. No mutation data was provided. "
                "To assess mutation impact, use the Mutation Scan mode."
            ),
            _section("F. Biological Interpretation",
                "Insufficient structural evidence to determine specific biological function. "
                "Residue composition analysis provides compositional clues but cannot determine function "
                "without ligand interaction data, evolutionary conservation, or experimental context."
            ),
            _section("G. Limitations",
                "Protein-only analysis based on primary sequence composition and static crystal structure. "
                "No ligand interaction analysis was performed. No secondary structure assignment (DSSP). "
                "No homology modeling, functional site prediction, or evolutionary conservation analysis. "
                "This is a compositional overview only. For functional insights, analyze a specific binding pocket "
                "or perform comparative structural analysis."
            ),
        ]
        return "\n\n".join(sections)

    return f"## A. Executive Summary\nUnknown mode: {mode}\n\n## G. Limitations\nUnsupported analysis mode."


def generate_structured_interpretation(mode, data):
    """Main entry: given mode and data, return an ordered dict of 7 AI sections.

    Modes: 'ligand' | 'mutation' | 'comparison' | 'protein_only'

    Returns OrderedDict with keys like 'A. Executive Summary', etc.
    Falls back to local rule-based interpretation when DeepSeek is unavailable.
    """
    try:
        prompt = build_structured_prompt(mode, data)
        client = create_deepseek_client()
        response = client.chat.completions.create(
            model="deepseek-chat",
            messages=[
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": prompt},
            ],
        )
        text = response.choices[0].message.content
    except Exception:
        import traceback
        traceback.print_exc()
        text = _structured_local_fallback(mode, data)

    return parse_ai_sections(text)
