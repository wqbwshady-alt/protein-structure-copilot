"""
Conservation and functional annotation analysis for protein residues.

Provides:
- UniProt REST API integration (ACT_SITE, BINDING, DOMAIN, VARIANT, MUTAGEN)
- BLOSUM62 substitution proxy (NOT true conservation — clearly labeled)
- DBREF-based PDB-to-UniProt residue mapping
- Graceful degradation: all external data is optional, failures never crash

Design constraints:
- conservation.score is ONLY for real conservation data or BLOSUM62 proxy.
  UniProt functional annotations go in functional_annotations.features.
- All external dependencies (UniProt API, DBREF) are optional.
  Missing data → No-data fallback + limitations.
- Mapping confidence explicitly flagged: high / medium / low.
"""

import json
import math
import os
import re

import requests

# ---- Paths ----

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache", "uniprot")

# ---- Residue property sets ----

THREE_LETTER_CODES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
}

# Residue to index mapping for BLOSUM62
_BLOSUM_ORDER = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
]
_BLOSUM_INDEX = {aa: i for i, aa in enumerate(_BLOSUM_ORDER)}

# BLOSUM62 substitution matrix (20 standard amino acids)
# Row/column order: A R N D C Q E G H I L K M F P S T W Y V
BLOSUM62_MATRIX = [
    [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0],  # A
    [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3],  # R
    [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3],  # N
    [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3],  # D
    [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],  # C
    [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2],  # Q
    [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2],  # E
    [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3],  # G
    [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3],  # H
    [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3],  # I
    [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1],  # L
    [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2],  # K
    [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1],  # M
    [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1],  # F
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2],  # P
    [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2],  # S
    [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0],  # T
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3],  # W
    [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1],  # Y
    [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4],  # V
]

# Max self-substitution score = W (11). Used for normalization to 0–1.
_MAX_BLOSUM_SELF = max(BLOSUM62_MATRIX[i][i] for i in range(len(_BLOSUM_ORDER)))


# ---- BLOSUM62 proxy (NOT conservation) ----


def blosum62_proxy_score(res_name):
    """Return BLOSUM62 self-substitution score normalized to [0, 1].

    This is a substitution-tolerance proxy — NOT evolutionary conservation.
    High score = this amino acid type is rarely substituted in homologous proteins.
    Low score = this amino acid type is frequently substituted.
    """
    idx = _BLOSUM_INDEX.get(res_name)
    if idx is None:
        return 0.5
    return round(BLOSUM62_MATRIX[idx][idx] / _MAX_BLOSUM_SELF, 4)


# ---- UniProt REST API ----


def _ensure_cache_dir():
    os.makedirs(CACHE_DIR, exist_ok=True)


def _cache_path(accession):
    safe = re.sub(r"[^A-Za-z0-9_\-]", "_", accession)
    return os.path.join(CACHE_DIR, f"{safe}.json")


def _load_cache(accession):
    path = _cache_path(accession)
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            return None
    return None


def _save_cache(accession, data):
    _ensure_cache_dir()
    path = _cache_path(accession)
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f)
    except IOError:
        pass


def fetch_uniprot_features(accession):
    """Fetch functional annotations from UniProt REST API.

    Returns dict with features list keyed by UniProt position, or None on failure.
    Cached to disk by accession.
    """
    cached = _load_cache(accession)
    if cached is not None:
        return cached

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code != 200:
            return None
        data = resp.json()
    except (requests.RequestException, json.JSONDecodeError, ValueError):
        return None

    features = _extract_features(data)
    _save_cache(accession, features)
    return features


# UniProt feature types we care about
_FEATURE_TYPES_OF_INTEREST = {
    "ACT_SITE",
    "BINDING",
    "DOMAIN",
    "MUTAGEN",
    "DISULFID",
    "VARIANT",
    "METAL",
    "NP_BIND",
    "ZN_FING",
    "MOTIF",
    "SITE",
    "MOD_RES",
    "LIPID",
    "CARBOHYD",
    "CROSSLNK",
}


def _extract_features(uniprot_data):
    """Extract functional feature annotations from UniProt JSON response.

    Returns dict: {uniprot_position_int: [feature_dict, ...]}
    """
    features_by_pos = {}

    raw_features = uniprot_data.get("features", [])
    for feat in raw_features:
        ftype = feat.get("type", "")
        if ftype not in _FEATURE_TYPES_OF_INTEREST:
            continue

        desc = feat.get("description", "") or ""
        location = feat.get("location", {})
        start = location.get("start", {})
        end = location.get("end", {})

        # Only handle single-residue or short-range features
        pos = start.get("value") if isinstance(start, dict) else None
        end_pos = end.get("value") if isinstance(end, dict) else None

        if pos is None:
            continue

        if end_pos is None:
            end_pos = pos

        for p in range(pos, end_pos + 1):
            if p not in features_by_pos:
                features_by_pos[p] = []
            features_by_pos[p].append({
                "type": ftype,
                "description": desc[:200] if desc else "",
                "uniprot_position": p,
            })

    return features_by_pos


# ---- DBREF parsing ----


def parse_dbref(pdb_path):
    """Parse DBREF records from a PDB file.

    Uses regex-based field extraction for robustness against column-alignment variations.

    Returns list of dicts:
        {chain_id, pdb_begin, pdb_end, insert_begin, insert_end,
         database, accession, db_begin, db_end}

    Returns empty list if no DBREF records found or file can't be read.
    """
    mappings = []
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith("DBREF "):
                    continue

                # Try fixed-column first, fall back to split-based
                m = _parse_dbref_columns(line)
                if m is None:
                    m = _parse_dbref_split(line)
                if m is None:
                    continue
                mappings.append(m)

    except (IOError, OSError):
        pass

    return mappings


def _parse_dbref_columns(line):
    """Parse DBREF using standard PDB column positions (v3.30)."""
    chain_id = line[12:13].strip()
    if not chain_id:
        return None

    try:
        pdb_begin = int(line[14:18].strip())
    except ValueError:
        return None
    insert_begin = line[18:19].strip() or None

    try:
        pdb_end = int(line[20:24].strip())
    except ValueError:
        pdb_end = pdb_begin
    insert_end = line[24:25].strip() or None

    database = line[26:32].strip()
    accession = line[33:43].strip()
    if not accession:
        return None
    # If column misalignment included dbIdCode characters, take first token only
    accession = accession.split()[0]

    # db_begin and db_end: use trailing-number scan only (more robust than fixed columns)
    trailing = line[43:].strip()
    db_begin, db_end = _extract_db_positions(trailing, pdb_begin, pdb_end)

    return {
        "chain_id": chain_id,
        "pdb_begin": pdb_begin,
        "pdb_end": pdb_end,
        "insert_begin": insert_begin,
        "insert_end": insert_end,
        "database": database,
        "accession": accession,
        "db_begin": db_begin,
        "db_end": db_end,
    }


def _parse_dbref_split(line):
    """Parse DBREF by splitting on whitespace (fallback)."""
    parts = line.split()
    if len(parts) < 9 or parts[0] != "DBREF":
        return None

    chain_id = parts[2]
    if not chain_id or len(chain_id) > 1:
        return None

    pdb_begin_str = parts[3]
    pdb_end_str = parts[4]
    pdb_begin, insert_begin = _parse_res_id(pdb_begin_str)
    pdb_end, insert_end = _parse_res_id(pdb_end_str)
    if pdb_begin is None or pdb_end is None:
        return None
    insert_begin = insert_begin or None
    insert_end = insert_end or None

    database = parts[5]
    accession = parts[6]
    if not accession:
        return None

    # dbseqBegin and dbseqEnd are typically at the end
    db_begin = None
    db_end = None
    for part in parts[7:]:
        if re.match(r"^-?\d+[A-Za-z]?$", part):
            try:
                num = int(re.sub(r"[A-Za-z]", "", part))
                if db_begin is None:
                    db_begin = num
                elif db_end is None:
                    db_end = num
                    break
            except ValueError:
                continue

    if db_begin is None:
        db_begin = pdb_begin
    if db_end is None:
        db_end = db_begin

    return {
        "chain_id": chain_id,
        "pdb_begin": pdb_begin,
        "pdb_end": pdb_end,
        "insert_begin": insert_begin,
        "insert_end": insert_end,
        "database": database,
        "accession": accession,
        "db_begin": db_begin,
        "db_end": db_end,
    }


def _extract_db_positions(trailing, pdb_default_begin, pdb_default_end):
    """Extract dbseqBegin and dbseqEnd from trailing text after dbAccession.

    Scans for the last two numeric tokens in the trailing portion of the DBREF line.
    """
    tokens = re.findall(r"(-?\d+)", trailing)
    tokens = [int(t) for t in tokens]
    if len(tokens) >= 2:
        db_begin = tokens[-2]
        db_end = tokens[-1]
        # Sanity: db_begin and db_end should be positive
        if db_begin > 0 and db_end > 0 and db_begin <= db_end:
            return db_begin, db_end
        # Try first two
        if tokens[0] > 0 and tokens[1] > 0 and tokens[0] <= tokens[1]:
            return tokens[0], tokens[1]
    elif len(tokens) == 1 and tokens[0] > 0:
        return tokens[0], tokens[0]

    return pdb_default_begin, pdb_default_end


def _parse_res_id(res_id_str):
    """Split PDB residue ID into integer + insertion code.

    '22' → (22, '')
    '22A' → (22, 'A')
    '100B' → (100, 'B')
    """
    m = re.match(r"^(-?\d+)([A-Za-z]?)$", str(res_id_str))
    if not m:
        return None, None
    return int(m.group(1)), m.group(2).upper() if m.group(2) else ""


def _linear_residue_index(pdb_res_id, pdb_begin, pdb_end, insert_begin, insert_end):
    """Compute a 0-based linear index for a PDB residue within a DBREF range.

    Handles insertion codes by counting them as additional positions.
    Returns (linear_index, confidence) or (None, None) if can't map.
    """
    res_num, ins_code = _parse_res_id(pdb_res_id)
    if res_num is None:
        return None, None

    if res_num < pdb_begin or res_num > pdb_end:
        return None, None

    if not ins_code:
        # Simple case: regular residue number, no insertion
        confidence = "high"
        return (res_num - pdb_begin, confidence)

    # With insertion code
    if ins_code == "A":
        confidence = "medium"
        return (res_num - pdb_begin + 1, confidence)
    elif ins_code == "B":
        confidence = "low"
        return (res_num - pdb_begin + 2, confidence)
    else:
        confidence = "low"
        offset = ord(ins_code) - ord("A") + 1
        return (res_num - pdb_begin + offset, confidence)


def map_pdb_to_uniprot(res_id_str, chain_id, dbref_mappings):
    """Map a PDB residue ID to (accession, uniprot_position, mapping_confidence).

    Returns (None, None, "low") if no mapping found.
    """
    res_num, ins_code = _parse_res_id(res_id_str)
    if res_num is None:
        return None, None, "low"

    for m in dbref_mappings:
        if m["chain_id"] != chain_id:
            continue
        if res_num < m["pdb_begin"] or res_num > m["pdb_end"]:
            continue

        linear_idx, conf = _linear_residue_index(
            res_id_str,
            m["pdb_begin"], m["pdb_end"],
            m["insert_begin"], m["insert_end"],
        )
        if linear_idx is None:
            continue

        uniprot_pos = m["db_begin"] + linear_idx
        if uniprot_pos > m["db_end"] and not ins_code:
            continue

        if ins_code:
            conf = "low"
        elif conf is None:
            conf = "high"

        return m["accession"], uniprot_pos, conf

    return None, None, "low"


# ---- Enhancement computation ----


def compute_conservation_annotation(contact_residues, pdb_path, consurf_scores=None):
    """Main entry point: compute conservation + functional annotation for contact residues.

    Args:
        contact_residues: {(chain_id, res_name, res_id): atom_dict, ...}
        pdb_path: path to PDB file
        consurf_scores: optional dict from consurf.map_consurf_to_residues()
            {residue_key: {score, color, confidence, mapped_position, source}}

    Returns:
        dict: {residue_key: {
            conservation: {score, available, source, source_detail},
            functional_annotations: {available, source, mapping_confidence, features},
            evidence_tags: {structural, enrichment, functional, conservation, proxy_only},
            limitations: [str, ...],
        }}

    Never raises — returns no-data fallback for all residues on failure.
    """
    residue_keys = list(contact_residues.keys())
    fallback = {}
    for chain_id, res_name, res_id in residue_keys:
        key = f"{chain_id}:{res_name}{res_id}"
        fallback[key] = _no_data_fallback(res_name, consurf_scores.get(key) if consurf_scores else None)

    if not residue_keys:
        return fallback

    # Step 1: Parse DBREF
    dbref_mappings = parse_dbref(pdb_path)

    # Step 2: Determine which accessions to fetch
    accession_features = {}
    accessions_to_fetch = set()
    residue_mapping = {}  # {residue_key: (accession, uniprot_pos, confidence)}

    for chain_id, res_name, res_id in residue_keys:
        key = f"{chain_id}:{res_name}{res_id}"
        if dbref_mappings:
            accession, uniprot_pos, conf = map_pdb_to_uniprot(res_id, chain_id, dbref_mappings)
            if accession:
                residue_mapping[key] = (accession, uniprot_pos, conf)
                if accession not in accession_features:
                    accessions_to_fetch.add(accession)

    # Step 3: Fetch UniProt features (cached + network)
    for accession in accessions_to_fetch:
        feat = fetch_uniprot_features(accession)
        if feat is not None:
            accession_features[accession] = feat

    # Step 4: Build per-residue annotations
    result = {}
    for chain_id, res_name, res_id in residue_keys:
        key = f"{chain_id}:{res_name}{res_id}"
        cs = consurf_scores.get(key) if consurf_scores else None
        result[key] = _build_residue_annotation(
            res_name, key, residue_mapping, accession_features, consurf_entry=cs
        )

    return result


def _no_data_fallback(res_name, consurf_entry=None):
    """Return a no-data enhancement dict for a single residue.

    If consurf_entry is provided, uses ConSurf-DB score as primary conservation
    and marks evidence accordingly.
    """
    blosum = blosum62_proxy_score(res_name)

    if consurf_entry:
        # ConSurf-DB data available — use as primary conservation
        consurf_color = consurf_entry.get("color")
        consurf_raw = consurf_entry.get("score")
        if consurf_color is not None:
            # Normalize ConSurf color (1-9 scale → 0-1)
            consurf_normalized = round((int(consurf_color) - 1) / 8.0, 4)
        else:
            consurf_normalized = blosum

        mapping_conf = consurf_entry.get("confidence", "low")
        limitations = [
            "No experimental validation of binding contribution",
            "No energetic calculation (distance-based geometric classification only)",
            f"ConSurf-DB mapping confidence: {mapping_conf}"
            if mapping_conf != "high"
            else "ConSurf-DB conservation data available",
            "No UniProt functional annotation available (no DBREF mapping or API data)",
            "Mapping confidence: low (no PDB-to-UniProt mapping)",
        ]
        if mapping_conf == "high":
            limitations = [
                "No experimental validation of binding contribution",
                "No energetic calculation (distance-based geometric classification only)",
                "No UniProt functional annotation available (no DBREF mapping or API data)",
                "Mapping confidence: low (no PDB-to-UniProt mapping)",
            ]

        return {
            "conservation": {
                "score": consurf_normalized,
                "available": True,
                "source": "consurf_db",
                "source_detail": (
                    f"ConSurf-DB evolutionary conservation score (color={consurf_color}/9, "
                    f"raw={consurf_raw}, normalized={consurf_normalized}). "
                    "True evolutionary conservation from MSA-based calculation. "
                    f"Mapping confidence: {mapping_conf}."
                ),
            },
            "functional_annotations": {
                "available": False,
                "source": "none",
                "mapping_confidence": "low",
                "features": [],
            },
            "evidence_tags": {
                "structural": True,
                "enrichment": True,
                "functional": False,
                "conservation": True,
                "proxy_only": False,
            },
            "limitations": limitations,
        }

    # No ConSurf data — BLOSUM62 proxy only
    return {
        "conservation": {
            "score": blosum,
            "available": False,
            "source": "blosum62_proxy",
            "source_detail": (
                "BLOSUM62 self-substitution score — substitution tolerance proxy, "
                "NOT evolutionary conservation. No homology-based conservation "
                "data available for this residue."
            ),
        },
        "functional_annotations": {
            "available": False,
            "source": "none",
            "mapping_confidence": "low",
            "features": [],
        },
        "evidence_tags": {
            "structural": True,
            "enrichment": True,
            "functional": False,
            "conservation": False,
            "proxy_only": True,
        },
        "limitations": [
            "No experimental validation of binding contribution",
            "No energetic calculation (distance-based geometric classification only)",
            "Conservation evidence unavailable — using BLOSUM62 substitution proxy",
            "No UniProt functional annotation available (no DBREF mapping or API data)",
            "Mapping confidence: low (no PDB-to-UniProt mapping)",
        ],
    }


def _build_residue_annotation(res_name, residue_key, residue_mapping, accession_features, consurf_entry=None):
    """Build annotation for one residue, combining ConSurf/UniProt + BLOSUM62 fallback.

    Priority:
    1. ConSurf-DB → true evolutionary conservation
    2. BLOSUM62 → substitution tolerance proxy (fallback)
    """
    blosum = blosum62_proxy_score(res_name)
    limitations = [
        "No experimental validation of binding contribution",
        "No energetic calculation (distance-based geometric classification only)",
    ]

    # --- Conservation: ConSurf-DB first, BLOSUM62 fallback ---
    if consurf_entry:
        consurf_color = consurf_entry.get("color")
        consurf_raw = consurf_entry.get("score")
        if consurf_color is not None:
            consurf_normalized = round((int(consurf_color) - 1) / 8.0, 4)
        else:
            consurf_normalized = blosum

        consurf_mapping_conf = consurf_entry.get("confidence", "low")
        conservation_data = {
            "score": consurf_normalized,
            "available": True,
            "source": "consurf_db",
            "source_detail": (
                f"ConSurf-DB evolutionary conservation score (color={consurf_color}/9, "
                f"raw={consurf_raw}, normalized={consurf_normalized}). "
                "True evolutionary conservation from MSA-based calculation. "
                f"Mapping confidence: {consurf_mapping_conf}."
            ),
        }
        cons_tags = {
            "conservation": True,
            "proxy_only": False,
        }
        if consurf_mapping_conf != "high":
            limitations.append(
                f"ConSurf-DB position mapping confidence is {consurf_mapping_conf} "
                f"for residue {residue_key}"
            )
    else:
        conservation_data = {
            "score": blosum,
            "available": False,
            "source": "blosum62_proxy",
            "source_detail": (
                "BLOSUM62 self-substitution score — substitution tolerance proxy. "
                "No evolutionary conservation data available."
            ),
        }
        cons_tags = {
            "conservation": False,
            "proxy_only": True,
        }
        limitations.append(
            "No true evolutionary conservation data available — using BLOSUM62 substitution proxy"
        )

    # --- Functional annotations from UniProt ---
    features = []
    feat_available = False
    feat_source = "none"
    mapping_conf = "low"

    if residue_key not in residue_mapping:
        limitations.append(
            "No UniProt functional annotation available (no DBREF mapping for this residue)"
        )
    else:
        accession, uniprot_pos, mapping_conf = residue_mapping[residue_key]

        if mapping_conf == "low":
            limitations.append(
                f"PDB-to-UniProt mapping confidence is low for residue {residue_key}"
            )
        elif mapping_conf == "medium":
            limitations.append(
                f"PDB-to-UniProt mapping confidence is medium for residue {residue_key} "
                "(insertion code involved)"
            )

        if accession in accession_features:
            pos_features = accession_features[accession].get(uniprot_pos, [])
            if pos_features:
                feat_available = True
                feat_source = "uniprot"
                features = list(pos_features)
            else:
                limitations.append(
                    f"No functional annotation at UniProt position {uniprot_pos} "
                    f"(accession {accession})"
                )
        else:
            limitations.append(
                f"UniProt API data unavailable for accession {accession}"
            )

    if not feat_available and residue_key in residue_mapping:
        limitations.append("No UniProt functional annotation for this residue position")

    evidence_tags = {
        "structural": True,
        "enrichment": True,
        "functional": feat_available,
        "conservation": cons_tags["conservation"],
        "proxy_only": cons_tags["proxy_only"],
    }

    return {
        "conservation": conservation_data,
        "functional_annotations": {
            "available": feat_available,
            "source": feat_source,
            "mapping_confidence": mapping_conf,
            "features": features,
        },
        "evidence_tags": evidence_tags,
        "limitations": limitations,
    }
