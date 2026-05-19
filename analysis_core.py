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
