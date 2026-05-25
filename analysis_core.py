import math


HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
POLAR = {"SER", "THR", "ASN", "GLN", "CYS", "TYR", "GLY"}
POSITIVE = {"LYS", "ARG", "HIS"}
NEGATIVE = {"ASP", "GLU"}
POLAR_ELEMENTS = {"N", "O", "S"}
CHARGED_LIGAND_ELEMENTS = {"N", "O", "P", "S"}
IGNORED_HETATM_RESIDUES = {"HOH", "WAT"}
PDB_RECORD_TYPES = ("ATOM", "HETATM")


def distance(a, b):
    return math.sqrt(
        (a[0] - b[0]) ** 2 +
        (a[1] - b[1]) ** 2 +
        (a[2] - b[2]) ** 2
    )


def parse_pdb_atoms(pdb_path):
    atoms = []

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(PDB_RECORD_TYPES):
                atom_type = line[0:6].strip()
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_id = line[22:26].strip()
                element = line[76:78].strip().upper() or infer_element(atom_name)

                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue

                atoms.append({
                    "atom_type": atom_type,
                    "atom_name": atom_name,
                    "element": element,
                    "res_name": res_name,
                    "chain_id": chain_id,
                    "res_id": res_id,
                    "coord": (x, y, z)
                })

    return atoms


def is_pdb_file(pdb_path):
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(PDB_RECORD_TYPES):
                return True

    return False


def list_ligands(pdb_path):
    atoms = parse_pdb_atoms(pdb_path)
    ligands = {}

    for atom in atoms:
        if atom["atom_type"] != "HETATM":
            continue

        if atom["res_name"] in IGNORED_HETATM_RESIDUES:
            continue

        key = (atom["res_name"], atom["chain_id"], atom["res_id"])
        ligands[key] = {
            "res_name": atom["res_name"],
            "chain_id": atom["chain_id"],
            "res_id": atom["res_id"],
            "atom_name": atom["atom_name"],
            "element": atom["element"]
        }

    return sorted(
        ligands.values(),
        key=lambda item: (item["res_name"], item["chain_id"], item["res_id"])
    )


def format_ligand_suggestions(ligands, limit=8):
    if not ligands:
        return "No non-water HETATM ligands were detected in this PDB file."

    suggestions = []

    for ligand in ligands[:limit]:
        location = ligand["res_id"]

        if ligand["chain_id"]:
            location = f'chain {ligand["chain_id"]} residue {ligand["res_id"]}'

        suggestions.append(f'{ligand["res_name"]} ({location})')

    suffix = ""

    if len(ligands) > limit:
        suffix = f" and {len(ligands) - limit} more"

    return "Available ligand candidates: " + ", ".join(suggestions) + suffix + "."


def classify_residue(res_name):
    if res_name in HYDROPHOBIC:
        return "hydrophobic"
    if res_name in POLAR:
        return "polar"
    if res_name in POSITIVE:
        return "positive"
    if res_name in NEGATIVE:
        return "negative"
    return "other"


def infer_element(atom_name):
    cleaned_name = "".join(char for char in atom_name.strip().upper() if char.isalpha())

    if not cleaned_name:
        return ""

    return cleaned_name[0]


def classify_interaction(ligand_atom, protein_atom, contact_distance):
    ligand_element = ligand_atom.get("element") or infer_element(ligand_atom["atom_name"])
    protein_element = protein_atom.get("element") or infer_element(protein_atom["atom_name"])
    residue_class = classify_residue(protein_atom["res_name"])

    if (
        contact_distance <= 4.0 and
        residue_class in {"positive", "negative"} and
        ligand_element in CHARGED_LIGAND_ELEMENTS
    ):
        return "charged / electrostatic", "#ef4444"

    if (
        contact_distance <= 3.5 and
        ligand_element in POLAR_ELEMENTS and
        protein_element in POLAR_ELEMENTS
    ):
        return "polar / possible H-bond", "#38bdf8"

    if (
        contact_distance <= 4.2 and
        residue_class == "hydrophobic" and
        ligand_element == "C" and
        protein_element == "C"
    ):
        return "hydrophobic contact", "#facc15"

    return "van der Waals contact", "#a3e635"


def analyze_ligand_pocket(pdb_path, ligand_name, cutoff=5.0, max_interactions=15):
    atoms = parse_pdb_atoms(pdb_path)

    ligand_atoms = [
        atom for atom in atoms
        if atom["atom_type"] == "HETATM" and atom["res_name"] == ligand_name
    ]

    protein_atoms = [
        atom for atom in atoms
        if atom["atom_type"] == "ATOM"
    ]

    if not ligand_atoms:
        return None, None, None, None

    contact_residues = {}
    nearest_interactions = {}

    for latom in ligand_atoms:
        for patom in protein_atoms:
            d = distance(latom["coord"], patom["coord"])

            if d <= cutoff:
                key = (patom["chain_id"], patom["res_name"], patom["res_id"])
                contact_residues[key] = patom

                if d <= 4.0:
                    current = nearest_interactions.get(key)

                    if current is None or d < current["distance"]:
                        interaction_type, color = classify_interaction(latom, patom, d)

                        nearest_interactions[key] = {
                            "start": list(latom["coord"]),
                            "end": list(patom["coord"]),
                            "distance": round(d, 2),
                            "chain_id": patom["chain_id"],
                            "res_name": patom["res_name"],
                            "res_id": patom["res_id"],
                            "atom_name": patom["atom_name"],
                            "element": patom["element"],
                            "ligand_atom": latom["atom_name"],
                            "ligand_element": latom["element"],
                            "interaction_type": interaction_type,
                            "color": color
                        }

    interactions = sorted(
        nearest_interactions.values(),
        key=lambda x: x["distance"]
    )

    if max_interactions is not None:
        interactions = interactions[:max_interactions]

    counts = {
        "hydrophobic": 0,
        "polar": 0,
        "positive": 0,
        "negative": 0,
        "other": 0
    }

    for chain_id, res_name, res_id in contact_residues.keys():
        counts[classify_residue(res_name)] += 1

    if counts["hydrophobic"] > max(counts["polar"], counts["positive"], counts["negative"]):
        primary_interpretation = (
            f"{ligand_name} binding pocket is mainly composed of hydrophobic residues, "
            "suggesting hydrophobic interaction may dominate ligand stabilization."
        )
    elif counts["positive"] > counts["negative"]:
        primary_interpretation = (
            f"{ligand_name} binding pocket is enriched in positively charged residues, "
            "suggesting possible electrostatic interaction with negatively charged ligands or nucleic acids."
        )
    elif counts["negative"] > counts["positive"]:
        primary_interpretation = (
            f"{ligand_name} binding pocket contains more negatively charged residues, "
            "suggesting a negatively charged local electrostatic environment."
        )
    else:
        primary_interpretation = (
            f"{ligand_name} binding pocket contains mixed residue types, suggesting multiple interaction modes."
        )

    return contact_residues, counts, primary_interpretation, interactions


def get_hotspot_residues(interactions):
    hotspots = []

    for item in interactions[:5]:
        hotspots.append({
            "chain_id": item["chain_id"],
            "res_id": item["res_id"],
            "res_name": item["res_name"],
            "atom_name": item["atom_name"],
            "element": item["element"],
            "distance": item["distance"],
            "interaction_type": item.get("interaction_type", "contact"),
            "color": item.get("color", "#facc15")
        })

    return hotspots


def residue_key_to_text(key):
    chain_id, res_name, res_id = key
    return f"Chain {chain_id}: {res_name}{res_id}"


def residue_keys_to_json(residue_keys):
    return [
        {
            "chain_id": chain_id,
            "res_name": res_name,
            "res_id": res_id,
            "atom_name": "",
            "element": ""
        }
        for chain_id, res_name, res_id in residue_keys
    ]


def sort_residue_keys(residue_keys):
    return sorted(
        residue_keys,
        key=lambda x: (x[0], int(x[2]) if x[2].isdigit() else 9999)
    )


# ============================================================
#  Phase 1: ResidueRanker — deterministic multi-factor scoring
# ============================================================

INTERACTION_WEIGHTS = {
    "charged / electrostatic": 1.0,
    "polar / possible H-bond": 0.8,
    "hydrophobic contact": 0.6,
    "van der Waals contact": 0.4,
}

AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}
CHARGED_RESIDUES = POSITIVE | NEGATIVE
D_MIN = 1.5
D_MAX = 4.0
MAX_WEIGHTED_CONTACTS = 5.0
MAX_UNIQUE_TYPES = 4


class ResidueRanker:
    """Deterministic multi-factor residue importance scoring.

    Formula:
      residue_score = distance_score + interaction_score + hotspot_score
                    + diversity_score + chemistry_score

    All components in [0, 1]; total max = 1.0.
    No LLM dependency — purely geometric + rule-based.
    """

    def __init__(self, contact_residues, interactions):
        self._contacts = contact_residues or {}
        self._interactions = interactions or []

        self._hotspot_keys = set()
        hotspots = get_hotspot_residues(self._interactions)
        for h in hotspots:
            self._hotspot_keys.add((h["chain_id"], h["res_name"], h["res_id"]))

        self._residue_interactions = {}
        for ix in self._interactions:
            key = (ix["chain_id"], ix["res_name"], ix["res_id"])
            self._residue_interactions.setdefault(key, []).append(ix)

    def score_all(self):
        """Return list of dicts sorted by score descending."""
        results = []
        for (chain_id, res_name, res_id), _atom in self._contacts.items():
            residue_ixs = self._residue_interactions.get(
                (chain_id, res_name, res_id), []
            )
            if not residue_ixs:
                continue

            min_d = min(ix["distance"] for ix in residue_ixs)
            components = {
                "distance_score": self._distance_score(min_d),
                "interaction_score": self._interaction_score(residue_ixs),
                "hotspot_score": self._hotspot_score(chain_id, res_name, res_id),
                "diversity_score": self._diversity_score(residue_ixs),
                "chemistry_score": self._chemistry_score(res_name),
            }
            total = sum(components.values())
            total = round(min(total, 1.0), 4)

            results.append({
                "chain_id": chain_id,
                "res_name": res_name,
                "res_id": res_id,
                "residue_key": f"{chain_id}:{res_name}{res_id}",
                "score": total,
                "score_components": components,
                "interaction_evidence": self._build_evidence(residue_ixs),
                "contact_count": len(residue_ixs),
                "min_distance": min_d,
                "residue_confidence": self._residue_confidence(min_d, residue_ixs),
                "pocket_location": self._pocket_location(min_d, len(residue_ixs)),
                "why_matters": self._why_matters(
                    chain_id, res_name, res_id, total, min_d, residue_ixs
                ),
            })

        results.sort(key=lambda r: (-r["score"], r["min_distance"]))
        for i, r in enumerate(results):
            r["rank"] = i + 1
        return results

    # ---- Scoring components ----

    def _distance_score(self, d):
        if d <= D_MIN:
            return 0.30
        if d >= D_MAX:
            return 0.0
        return round(0.30 * (1.0 - (d - D_MIN) / (D_MAX - D_MIN)), 4)

    def _interaction_score(self, residue_ixs):
        weighted = sum(
            INTERACTION_WEIGHTS.get(ix.get("interaction_type", ""), 0.4)
            for ix in residue_ixs
        )
        return round(0.30 * min(weighted, MAX_WEIGHTED_CONTACTS) / MAX_WEIGHTED_CONTACTS, 4)

    def _hotspot_score(self, chain_id, res_name, res_id):
        return 0.15 if (chain_id, res_name, res_id) in self._hotspot_keys else 0.0

    def _diversity_score(self, residue_ixs):
        unique = len(set(ix.get("interaction_type", "") for ix in residue_ixs))
        return round(0.15 * min(unique, MAX_UNIQUE_TYPES) / MAX_UNIQUE_TYPES, 4)

    def _chemistry_score(self, res_name):
        score = 0.0
        if res_name in CHARGED_RESIDUES:
            score += 0.05
        if res_name in AROMATIC_RESIDUES:
            score += 0.05
        if res_name == "GLY":
            score += 0.05
        if res_name == "PRO":
            score += 0.05
        return round(min(score, 0.10), 4)

    # ---- Evidence ----

    def _build_evidence(self, residue_ixs):
        return [
            {
                "type": ix.get("interaction_type", "contact"),
                "distance": ix["distance"],
                "ligand_atom": ix.get("ligand_atom", ""),
                "protein_atom": ix.get("atom_name", ""),
                "confidence": self._interaction_confidence(
                    ix.get("interaction_type", ""), ix["distance"],
                    ix.get("ligand_element", ""), ix.get("element", "")
                ),
            }
            for ix in residue_ixs
        ]

    # ---- Confidence ----

    def _residue_confidence(self, min_d, residue_ixs):
        has_strong = any(
            ix.get("interaction_type", "") in ("charged / electrostatic", "polar / possible H-bond")
            for ix in residue_ixs
        )
        if min_d <= 3.2 and len(residue_ixs) >= 2 and has_strong:
            return "high"
        if min_d > 3.8 or (len(residue_ixs) == 1 and not has_strong):
            return "low"
        return "medium"

    def _interaction_confidence(self, itype, distance, lig_el, prot_el):
        if itype == "charged / electrostatic" and distance <= 3.2:
            return "high"
        if itype == "polar / possible H-bond" and distance <= 3.0:
            return "high"
        if itype == "hydrophobic contact" and distance <= 3.5 and lig_el == "C" and prot_el == "C":
            return "high"
        if itype == "van der Waals contact" and distance > 3.5:
            return "low"
        return "medium"

    # ---- Pocket location ----

    def _pocket_location(self, min_d, contact_count):
        if min_d <= 3.0 and contact_count >= 2:
            return "core"
        if min_d <= 3.8:
            return "inner"
        return "peripheral"

    # ---- Explanation ----

    def _why_matters(self, chain_id, res_name, res_id, score, min_d, residue_ixs):
        types = sorted(set(ix.get("interaction_type", "contact") for ix in residue_ixs))
        type_str = ", ".join(types)
        count = len(residue_ixs)
        is_hotspot = (chain_id, res_name, res_id) in self._hotspot_keys

        if score >= 0.70:
            level = "key binding anchor"
        elif score >= 0.45:
            level = "significant contributor"
        elif score >= 0.20:
            level = "moderate contact"
        else:
            level = "peripheral contact"

        parts = [
            f"{res_name}{res_id} ranks as a {level} (score {score:.2f})",
            f"because it forms {count} close contact(s) at minimum {min_d:.1f}A",
            f"involving {type_str}.",
        ]
        if is_hotspot:
            parts.append("It is among the top-5 hotspot residues.")
        if score >= 0.50:
            parts.append("This residue likely plays a significant role in ligand stabilization.")
        else:
            parts.append("Its contribution to binding may be secondary to the core anchor residues.")

        return " ".join(parts)


# ---- Statistical Enrichment Analysis ----


def compute_full_protein_composition(atoms):
    """Compute residue-type composition of the entire protein (ATOM records only).

    Args:
        atoms: list of atom dicts from parse_pdb_atoms()

    Returns:
        dict: {hydrophobic: N, polar: N, positive: N, negative: N, other: N, total: N,
               residues: {(chain_id, res_name, res_id), ...}}
    """
    protein_atoms = [a for a in atoms if a["atom_type"] == "ATOM"]
    unique_residues = set()
    for a in protein_atoms:
        unique_residues.add((a["chain_id"], a["res_name"], a["res_id"]))

    composition = {"hydrophobic": 0, "polar": 0, "positive": 0,
                   "negative": 0, "other": 0, "total": len(unique_residues),
                   "residues": unique_residues}

    for _, res_name, _ in unique_residues:
        cat = classify_residue(res_name)
        composition[cat] += 1

    return composition


def compute_pocket_enrichment(pocket_counts, full_composition):
    """Fisher exact test for enrichment/depletion of residue types in pocket vs whole protein.

    Compares pocket contact-residue composition against the full protein background.

    Args:
        pocket_counts: dict from analyze_ligand_pocket() — {hydrophobic, polar, positive, negative}
        full_composition: dict from compute_full_protein_composition()

    Returns:
        dict with per-type enrichment metrics + significant_types + baseline info
    """
    try:
        from scipy.stats import fisher_exact
    except ImportError:
        # scipy unavailable: return structure without p-values
        return _enrichment_fallback(pocket_counts, full_composition)

    total_pocket = sum(pocket_counts.get(k, 0) for k in
                       ["hydrophobic", "polar", "positive", "negative"])
    total_protein = full_composition["total"]

    if total_pocket == 0 or total_protein == 0:
        return _enrichment_fallback(pocket_counts, full_composition)

    result = {}
    significant_types = []

    for cat in ["hydrophobic", "polar", "positive", "negative"]:
        pocket_n = pocket_counts.get(cat, 0)
        protein_n = full_composition.get(cat, 0)
        pocket_frac = pocket_n / total_pocket if total_pocket > 0 else 0.0
        protein_frac = protein_n / total_protein if total_protein > 0 else 0.0

        fold = pocket_frac / protein_frac if protein_frac > 0 else float("inf") if pocket_frac > 0 else 0.0
        fold = round(fold, 2)

        # Fisher exact contingency table:
        #           Pocket   Rest_of_protein
        # Cat       a        b
        # Not cat   c        d
        a = pocket_n
        b = protein_n - pocket_n
        c = total_pocket - pocket_n
        d = (total_protein - total_pocket) - (protein_n - pocket_n)

        if b < 0 or c < 0 or d < 0:
            p_value = 1.0
        else:
            try:
                _, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")
                p_value = round(float(p_value), 4)
            except (ValueError, TypeError):
                p_value = 1.0

        result[cat] = {
            "pocket_count": pocket_n,
            "protein_count": protein_n,
            "pocket_fraction": round(pocket_frac, 4),
            "protein_fraction": round(protein_frac, 4),
            "fold_enrichment": fold,
            "p_value": p_value,
            "significant": p_value < 0.05,
        }

        if p_value < 0.05 and fold > 1.0:
            significant_types.append(f"{cat} ({fold}x enriched)")
        elif p_value < 0.05 and fold < 1.0:
            significant_types.append(f"{cat} ({fold}x depleted)")

    result["significant_types"] = significant_types
    result["baseline"] = "whole_protein"
    result["test"] = "fisher_exact"

    return result


def _enrichment_fallback(pocket_counts, full_composition):
    """Fallback enrichment without p-values (when scipy is unavailable)."""
    total_pocket = sum(pocket_counts.get(k, 0) for k in
                       ["hydrophobic", "polar", "positive", "negative"])
    total_protein = full_composition.get("total", 1)

    result = {}
    for cat in ["hydrophobic", "polar", "positive", "negative"]:
        pocket_n = pocket_counts.get(cat, 0)
        protein_n = full_composition.get(cat, 0)
        pocket_frac = pocket_n / total_pocket if total_pocket > 0 else 0.0
        protein_frac = protein_n / total_protein if total_protein > 0 else 0.0
        fold = pocket_frac / protein_frac if protein_frac > 0 else 0.0

        result[cat] = {
            "pocket_count": pocket_n,
            "protein_count": protein_n,
            "pocket_fraction": round(pocket_frac, 4),
            "protein_fraction": round(protein_frac, 4),
            "fold_enrichment": round(fold, 2),
            "p_value": None,
            "significant": None,
        }

    result["significant_types"] = []
    result["baseline"] = "whole_protein"
    result["test"] = "fisher_exact (unavailable — scipy not installed, p-values not computed)"
    return result


def merge_enhancement(ranked_residues, enrichment, conservation_annotations):
    """Attach enrichment + conservation annotation to each ranked residue.

    Does NOT modify ranking order or scores. Only adds new fields.

    Args:
        ranked_residues: list of dicts from ResidueRanker.score_all()
        enrichment: dict from compute_pocket_enrichment()
        conservation_annotations: dict from compute_conservation_annotation()

    Returns:
        list of dicts — same residues, with added 'enrichment' and 'conservation' fields
    """
    enhanced = []
    for r in ranked_residues:
        key = r.get("residue_key", "")
        entry = dict(r)

        # Per-residue enrichment: add the category-level enrichment for this residue type
        res_name = r.get("res_name", "")
        cat = classify_residue(res_name)
        entry["enrichment"] = {
            "category": cat,
            "category_enrichment": enrichment.get(cat) if enrichment else None,
            "overall_enrichment": {
                "significant_types": enrichment.get("significant_types", []) if enrichment else [],
                "baseline": enrichment.get("baseline", "") if enrichment else "",
                "test": enrichment.get("test", "") if enrichment else "",
            },
        }

        # Conservation + functional annotation
        entry["conservation"] = {}
        entry["functional_annotations"] = {}
        entry["evidence_tags"] = {}
        entry["residue_limitations"] = []

        if conservation_annotations and key in conservation_annotations:
            ca = conservation_annotations[key]
            entry["conservation"] = ca.get("conservation", {})
            entry["functional_annotations"] = ca.get("functional_annotations", {})
            entry["evidence_tags"] = ca.get("evidence_tags", {})
            entry["residue_limitations"] = ca.get("limitations", [])
        else:
            # No-data fallback for this residue
            blosum = 0.5
            entry["conservation"] = {
                "score": 0.5,
                "available": False,
                "source": "blosum62_proxy",
                "source_detail": "BLOSUM62 self-substitution proxy. No conservation data available.",
            }
            entry["functional_annotations"] = {
                "available": False,
                "source": "none",
                "mapping_confidence": "low",
                "features": [],
            }
            entry["evidence_tags"] = {
                "structural": True,
                "enrichment": True,
                "functional": False,
                "conservation": False,
                "proxy_only": True,
            }
            entry["residue_limitations"] = [
                "No conservation evidence available",
                "No functional annotation available",
                "No experimental validation of binding contribution",
            ]

        enhanced.append(entry)

    return enhanced


# ============================================================
#  Phase 2: ConfidenceAssessor + LimitationsBuilder
# ============================================================

class ConfidenceAssessor:
    """Evidence-based confidence classification — no LLM dependency.

    Produces: overall_analysis_confidence, per-residue distribution,
    per-interaction distribution, ai_interpretation_confidence,
    data_quality_flags.
    """

    def __init__(self, contact_residues, interactions, ligand_detected,
                 ai_available=True, ai_truncated=False):
        self._contacts = contact_residues or {}
        self._interactions = interactions or []
        self._ligand_detected = ligand_detected
        self._ai_available = ai_available
        self._ai_truncated = ai_truncated

    def assess(self, ranked_residues=None):
        overall = self._overall_confidence()
        res_dist = self._residue_distribution(ranked_residues or [])
        int_dist = self._interaction_distribution()

        return {
            "overall_analysis_confidence": overall,
            "confidence_reason": self._confidence_reason(overall),
            "residue_confidence_distribution": res_dist,
            "interaction_confidence_distribution": int_dist,
            "ai_interpretation_confidence": self._ai_confidence(overall),
            "data_quality_flags": self._quality_flags(),
        }

    def _overall_confidence(self):
        n_contacts = len(self._contacts)
        types = set(ix.get("interaction_type", "") for ix in self._interactions)
        has_close_strong = any(
            ix.get("interaction_type", "") in ("charged / electrostatic", "polar / possible H-bond")
            and ix["distance"] <= 3.5
            for ix in self._interactions
        )
        if (
            self._ligand_detected
            and n_contacts >= 15
            and len(types) >= 2
            and has_close_strong
        ):
            return "high"
        if (
            not self._ligand_detected
            or n_contacts < 10
            or len(types) <= 1
            or (not has_close_strong and n_contacts < 15)
        ):
            return "low"
        return "medium"

    def _confidence_reason(self, overall):
        n = len(self._contacts)
        types = sorted(set(ix.get("interaction_type", "") for ix in self._interactions))
        reasons = []
        if self._ligand_detected:
            reasons.append(f"ligand detected with {n} contact residues")
        else:
            reasons.append("no ligand detected")
        reasons.append(f"interaction types: {', '.join(types) if types else 'none'}")
        if n >= 15:
            reasons.append("sufficient contact count (≥15)")
        else:
            reasons.append("low contact count (<15)")
        reasons.append("static crystal structure without MD/energy validation")
        if overall == "high":
            reasons.append("multiple strong interaction types with close contacts")
        elif overall == "low":
            reasons.append("insufficient structural evidence for high-confidence interpretation")
        return "; ".join(reasons)

    def _residue_distribution(self, ranked):
        dist = {"high": 0, "medium": 0, "low": 0}
        for r in ranked:
            conf = r.get("residue_confidence", "low")
            if conf in dist:
                dist[conf] += 1
        return dist

    def _interaction_distribution(self):
        dist = {"high": 0, "medium": 0, "low": 0}
        for ix in self._interactions:
            itype = ix.get("interaction_type", "")
            d = ix["distance"]
            lig_el = ix.get("ligand_element", "")
            prot_el = ix.get("element", "")
            if itype == "charged / electrostatic" and d <= 3.2:
                dist["high"] += 1
            elif itype == "polar / possible H-bond" and d <= 3.0:
                dist["high"] += 1
            elif itype == "hydrophobic contact" and d <= 3.5 and lig_el == "C" and prot_el == "C":
                dist["high"] += 1
            elif itype == "van der Waals contact" and d > 3.5:
                dist["low"] += 1
            else:
                dist["medium"] += 1
        return dist

    def _ai_confidence(self, overall):
        if not self._ai_available:
            return "low"
        if self._ai_truncated or overall == "low":
            return "low"
        if overall == "high":
            return "medium"
        return "medium"

    def _quality_flags(self):
        n_contacts = len(self._contacts)
        types = set(ix.get("interaction_type", "") for ix in self._interactions)
        return {
            "has_valid_pdb": True,
            "ligand_detected": self._ligand_detected,
            "contact_count_sufficient": n_contacts >= 10,
            "has_multiple_evidence_types": len(types) >= 2,
            "has_close_strong_contacts": any(
                ix["distance"] <= 3.5
                and ix.get("interaction_type", "") in ("charged / electrostatic", "polar / possible H-bond")
                for ix in self._interactions
            ),
            "has_energetic_validation": False,
            "has_conservation_data": False,
            "has_md_simulation": False,
            "has_multiple_structures": False,
        }


class LimitationsBuilder:
    """Structured limitations — always present in output."""

    def build(self, mode="ligand", ligand_ambiguous=False):
        base = {
            "static_structure_only": True,
            "no_energetic_validation": True,
            "no_dynamics": True,
            "no_hydrogens": True,
            "no_solvent_modeling": True,
            "no_conservation_analysis": True,
            "no_mutation_validation": True,
            "no_docking": True,
            "distance_cutoff_5A": True,
            "interaction_classification_geometric_only": True,
            "ligand_ambiguous": ligand_ambiguous,
            "disclaimer": (
                "This analysis is based on geometric distance criteria from a "
                "single static PDB structure. It does not constitute a validated "
                "energetic prediction. For quantitative binding assessment, "
                "MD simulation, MM-GBSA, ITC, or SPR experiments are recommended."
            ),
        }
        if mode == "mutation":
            base["disclaimer"] += (
                " Mutation impact predictions are heuristic estimates, "
                "not validated ΔΔG calculations. Sidechain remodeling, "
                "backbone rearrangement, and solvent reorganization are not modeled."
            )
        if mode == "comparison":
            base["disclaimer"] += (
                " Structural comparison is based on static crystallographic "
                "conformations. Conformational changes between WT and mutant "
                "may not be fully captured by contact comparison alone."
            )
        return base


# ============================================================
#  Phase 3: SafetyGuardrails — scientific overclaiming prevention
# ============================================================

FORBIDDEN_PATTERNS = [
    (r"\bKd\s*[=~<>≈]\s*\d+", "numeric Kd value"),
    (r"\bKi\s*[=~<>≈]\s*\d+", "numeric Ki value"),
    (r"\bIC50\s*[=~<>≈]\s*\d+", "numeric IC50 value"),
    (r"\bΔG\s*[=~<>≈]\s*-?\d+", "numeric ΔG value"),
    (r"\bbinding affinity\s*[=~<>≈]\s*\d+", "numeric binding affinity"),
    (r"\b(inhibits?|treats?|cures?|therapy for)\b.*\b(disease|cancer|disorder|syndrome)\b",
     "disease mechanism claim"),
    (r"\b(ritonavir|darunavir|lopinavir|atazanavir|saquinavir|tipranavir|indinavir|nelfinavir|amprenavir)\b",
     "specific drug name"),
    (r"\bknown\s+(to|as)\s+(a\s+)?(drug|inhibitor|therapeutic|treatment)\b",
     "unverified drug classification"),
]

EVIDENCE_TAG_RULES = {
    "charged / electrostatic": "[S]",
    "polar / possible H-bond": "[S]",
    "hydrophobic contact": "[S]",
    "van der Waals contact": "[S]",
    "suggest": "[I]",
    "indicate": "[I]",
    "consistent with": "[I]",
    "may": "[H]",
    "possibly": "[H]",
    "likely": "[H]",
    "heuristic": "[H]",
    "requires experimental": "[E]",
    "validated": "[E]",
}

MANDATORY_DISCLAIMER = (
    "[E] This is a structural hypothesis based on geometric distance criteria "
    "from a single static PDB structure. It is not a validated energetic or "
    "functional prediction."
)


class SafetyGuardrails:
    """Post-process AI output for scientific safety."""

    @staticmethod
    def validate(text):
        """Check text for forbidden claims. Returns (is_safe, violations)."""
        if not text:
            return True, []
        violations = []
        import re
        for pattern, description in FORBIDDEN_PATTERNS:
            matches = re.findall(pattern, text, re.IGNORECASE)
            if matches:
                violations.append({
                    "pattern": pattern,
                    "description": description,
                    "matches": matches[:3],
                })
        return len(violations) == 0, violations

    @staticmethod
    def has_disclaimer(text):
        """Check if text contains an appropriate limitations disclaimer."""
        if not text:
            return False
        keywords = ["static structure", "structural hypothesis", "not a validated",
                     "structural inference", "heuristic", "requires experimental"]
        return any(kw in text.lower() for kw in keywords)

    @staticmethod
    def tag_evidence(text, mode="ligand"):
        """Apply [S][I][H][E] evidence tags to fallback output paragraphs."""
        if not text:
            return text
        lines = text.split("\n")
        tagged = []
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith("## "):
                tagged.append(line)
                continue
            # Determine evidence level
            if any(w in stripped.lower() for w in
                   ["distance", "contact", "atom ", "residue", "chain ", "angstrom", "a)",
                    "bullet", "- ", "charged", "polar", "hydrophobic", "van der waals"]):
                tagged.append("[S] " + line)
            elif any(w in stripped.lower() for w in
                     ["suggest", "indicate", "consistent with", "appears",
                      "characteristic of", "typical"]):
                tagged.append("[I] " + line)
            elif any(w in stripped.lower() for w in
                     ["may", "possibly", "likely", "could", "might", "heuristic",
                      "mutation", "impact"]):
                tagged.append("[H] " + line)
            else:
                tagged.append("[I] " + line)
        result = "\n".join(tagged)
        if SafetyGuardrails.has_disclaimer(result):
            return result
        return result + "\n\n" + MANDATORY_DISCLAIMER

    @staticmethod
    def safe_output_template(mode, data):
        """Build a guaranteed-safe fallback that never overclaims."""
        ranked = data.get("important_residues", [])[:5]
        conf = data.get("confidence", {})
        lim = data.get("limitations", {})

        sections = []
        sections.append("## A. Executive Summary")
        if ranked:
            top = ranked[0]
            sections.append(
                f"[S] {top['residue_key']} ranks #1 (score {top['score']:.2f}) "
                f"with {top['contact_count']} contact(s) at minimum {top['min_distance']}A."
            )
            sections.append(
                f"[I] The pocket contains {len(data.get('important_residues',[]))} "
                f"ranked residues with {conf.get('overall_analysis_confidence','unknown')} "
                f"overall confidence."
            )
        else:
            sections.append("[S] No ranked residues available for this analysis.")

        sections.append("## B. Pocket Profile")
        sections.append("[S] Pocket residue composition and interaction types "
                        "are based on geometric classification from PDB coordinates.")

        sections.append("## C. Key Interactions")
        if ranked:
            for r in ranked[:3]:
                ev = r.get("interaction_evidence", [])
                types = ", ".join(sorted(set(e["type"] for e in ev)))
                sections.append(
                    f"[S] {r['residue_key']}: {len(ev)} contact(s) at "
                    f"{r['min_distance']}A — {types}"
                )
        else:
            sections.append("[S] Insufficient structural evidence — no close "
                            "contacts detected within 4.0A cutoff.")

        sections.append("## D. Residue-Level Evidence")
        if ranked:
            for r in ranked[:5]:
                sections.append(f"- [S] {r['residue_key']}: score {r['score']:.2f}, "
                                f"min distance {r['min_distance']}A, "
                                f"{r['contact_count']} contact(s)")
        else:
            sections.append("[S] No residue-level interaction data within cutoff distance.")

        sections.append("## E. Mutation Impact Assessment")
        if mode == "mutation":
            sections.append("[H] Mutation impact is a heuristic estimate based on "
                            "physicochemical property comparison and geometric interaction "
                            "patterns. Confidence: LOW — no energetic validation performed.")
        else:
            sections.append("[S] N/A — single-structure analysis. No mutation data provided.")

        sections.append("## F. Biological Interpretation")
        sections.append(
            "[I] The structural evidence suggests a binding interface consistent "
            "with geometric pocket-ligand complementarity. Specific functional roles "
            "cannot be assigned from structural data alone."
        )

        sections.append("## G. Limitations")
        sections.append(lim.get("disclaimer", MANDATORY_DISCLAIMER))

        return "\n\n".join(sections)

    @staticmethod
    def recommended_next_steps(mode, conf):
        """Generate evidence-based recommendations — no hallucination."""
        steps = []
        overall = conf.get("overall_analysis_confidence", "medium") if conf else "medium"

        steps.append(
            "Perform 100ns molecular dynamics simulation to assess "
            "pocket flexibility and contact stability over time."
        )
        if overall != "high":
            steps.append(
                "Validate key residue contacts with alanine scanning "
                "mutagenesis or site-directed mutagenesis experiments."
            )
        steps.append(
            "Run MM-GBSA or MM-PBSA to estimate per-residue binding "
            "free energy contributions (currently distance-based only)."
        )
        if conf and conf.get("data_quality_flags", {}).get("has_multiple_evidence_types"):
            steps.append(
                "Compare with homologous structures to assess "
                "conservation of key contact residues across species."
            )
        steps.append(
            "Experimental validation: SPR or ITC to measure binding "
            "affinity and compare with structural predictions."
        )
        return steps
