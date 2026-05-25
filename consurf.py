"""
ConSurf-DB integration for real evolutionary conservation scores.

API flow (reverse-engineered from consurfDB client library):
  1. GET chain_selection.php?pdb_ID=XXXX  → extract {chain: 5-char final_id}
  2. GET DB/{final_id}/consurf_summary.txt → parse TSV grades file

Provides:
- Disk cache with 30-day TTL
- Conservative position mapping (3LATOM → PDB residue ID)
- Graceful degradation: failures return None, never crash
"""

import json
import os
import re
import time

import requests

# ---- Configuration ----

CONSURF_BASE_URL = os.getenv(
    "CONSURF_BASE_URL",
    "https://consurfdb.tau.ac.il"
)

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache", "consurf")

CACHE_TTL = 30 * 24 * 3600

REQUEST_TIMEOUT = 10

# SSL verify: ConSurf server certificate may fail validation
REQUEST_VERIFY = os.getenv("CONSURF_VERIFY_SSL", "1").strip().lower() in {"1", "true", "yes"}

# ---- Cache helpers ----

def _ensure_cache_dir():
    os.makedirs(CACHE_DIR, exist_ok=True)


def _cache_path(pdb_id):
    safe = re.sub(r"[^A-Za-z0-9_]", "_", pdb_id)
    return os.path.join(CACHE_DIR, f"{safe}.json")


def _load_cache(pdb_id):
    path = _cache_path(pdb_id)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except (json.JSONDecodeError, IOError):
        return None
    cached_at = data.get("_cached_at", 0)
    if time.time() - cached_at > CACHE_TTL:
        return None
    return data.get("_payload", None)


def _save_cache(pdb_id, payload):
    _ensure_cache_dir()
    path = _cache_path(pdb_id)
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump({"_cached_at": int(time.time()), "_payload": payload}, f)
    except IOError:
        pass


# ---- PDB ID extraction ----

def extract_pdb_id(pdb_path, filename=""):
    """Extract PDB ID from filename or PDB HEADER record."""
    if filename:
        m = re.match(r"RCSB_([A-Za-z0-9]{4})_", filename, re.IGNORECASE)
        if m:
            return m.group(1).upper()
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("HEADER") and len(line) >= 66:
                    id_code = line[62:66].strip().upper()
                    if re.fullmatch(r"[A-Z0-9]{4}", id_code):
                        return id_code
                    break
    except (IOError, OSError):
        pass
    return None


# ---- ConSurf-DB API query (2-step) ----

def query_consurf_db(pdb_id):
    """Query ConSurf-DB for precomputed conservation scores.

    Two-step API:
      1. chain_selection.php → {chain_id: 5-char final_id}
      2. DB/{final_id}/consurf_summary.txt → TSV grades

    Args:
        pdb_id: 4-character PDB ID (e.g. "1HSG")

    Returns:
        dict: {chain_id: [{pos: int, score: float, color: int, residue_name: str}, ...]}
        None: if unavailable, failed, or no data

    Never raises.
    """
    if not pdb_id or not re.fullmatch(r"[A-Za-z0-9]{4}", pdb_id):
        return None

    pdb_id = pdb_id.upper()

    cached = _load_cache(pdb_id)
    if cached is not None:
        return cached

    try:
        # Step 1: get chain → final_id mapping
        chain_map = _fetch_chain_mapping(pdb_id)
        if not chain_map:
            return None

        # Step 2: fetch grades for each chain
        all_data = {}
        for chain, final_id in sorted(chain_map.items()):
            grades = _fetch_grades(final_id)
            if grades:
                all_data[chain] = grades

        if not all_data:
            return None

    except Exception:
        return None

    _save_cache(pdb_id, all_data)
    return all_data


def _fetch_chain_mapping(pdb_id):
    """Step 1: GET chain_selection.php → extract {chain: final_id} mapping.

    Returns dict or None.
    """
    url = f"{CONSURF_BASE_URL}/scripts/chain_selection.php"
    try:
        resp = requests.get(
            url,
            params={"pdb_ID": pdb_id},
            timeout=REQUEST_TIMEOUT,
            verify=REQUEST_VERIFY,
        )
    except requests.RequestException:
        return None

    if resp.status_code != 200:
        return None

    text = resp.text
    if "No chains found for" in text:
        return None
    if "Temporary Downtime" in text:
        return None

    # Parse: <option value="A ABCDE">chain A</option>
    mapping = dict(re.findall(r'option value="(\w) (\w{5})"', text))
    return mapping if mapping else None


def _fetch_grades(final_id):
    """Step 2: GET DB/{final_id}/consurf_summary.txt → parse TSV.

    TSV columns (from ConsurfDB ResidueDataType):
        POS, SEQ, 3LATOM, SCORE, COLOR,
        CONFIDENCE INTERVAL, CONFIDENCE COLORS,
        MSA DATA, RESIDUE VARIETY

    Returns list of dicts: [{pos, score, color, residue_name, chain_id}, ...]
    """
    url = f"{CONSURF_BASE_URL}/DB/{final_id}/consurf_summary.txt"
    try:
        resp = requests.get(url, timeout=REQUEST_TIMEOUT, verify=REQUEST_VERIFY)
    except requests.RequestException:
        return None

    if resp.status_code != 200:
        return None

    return _parse_grades_text(resp.text)


def _parse_grades_text(text):
    """Parse ConSurf grades TSV text.

    Extracts per-residue: PDB position (from 3LATOM), conservation score, color.
    """
    if not text or not text.strip():
        return None

    # Column names in order (from ConsurfDB ResidueDataType definition)
    columns = [
        "POS", "SEQ", "3LATOM", "SCORE", "COLOR",
        "CONFIDENCE_INTERVAL", "CONFIDENCE_COLORS",
        "MSA_DATA_IDX", "RESIDUE_VARIETY",
    ]

    results = []
    for line in text.split("\n"):
        line = line.strip()
        if not line or not line[0].isdigit():
            continue

        parts = [v.strip() for v in line.split("\t") if v]
        if len(parts) < 5:
            continue

        row = dict(zip(columns, parts))

        # Parse 3LATOM → residue_name, position, chain_id
        # Format: "MET1:A" or "ASP25A:B" or "GLY-1:A"
        latom = row.get("3LATOM", "")
        if latom == "-":
            # Missing density — skip
            continue

        parsed = _parse_3latom(latom)
        if not parsed:
            continue

        try:
            score = float(row.get("SCORE", ""))
            color = int(row.get("COLOR", ""))
        except (ValueError, TypeError):
            continue

        results.append({
            "pos": parsed["pos"],
            "score": score,
            "color": color,
            "chain_id": parsed["chain"],
            "residue_name": parsed["res_name"],
            "insertion_code": parsed["ins_code"],
        })

    return results if results else None


def _parse_3latom(latom):
    """Parse a 3LATOM string like 'MET1:A' or 'ASP25A:B' or 'GLY-1:A'.

    Returns {res_name, pos, ins_code, chain} or None.
    """
    # Match: 3-letter code, optional negative sign, number, optional insertion code, colon, chain
    m = re.match(r"^([A-Z]{3})(-?\d+)([A-Za-z]?):([A-Za-z0-9])$", latom)
    if not m:
        return None
    return {
        "res_name": m.group(1).upper(),
        "pos": int(m.group(2)),
        "ins_code": m.group(3).upper() if m.group(3) else "",
        "chain": m.group(4).upper(),
    }


# ---- Mapping: ConSurf positions → PDB contact residues ----

def map_consurf_to_residues(consurf_data, contact_residues):
    """Map ConSurf-DB per-position scores to PDB contact residues.

    Conservative mapping:
      - Exact position + chain match → confidence "high"
      - Position match with insertion code → confidence "medium"
      - No match → skip

    Args:
        consurf_data: from query_consurf_db() — {chain_id: [{pos, score, ...}]}
        contact_residues: {(chain_id, res_name, res_id): atom_dict, ...}

    Returns:
        dict: {residue_key: {score, color, confidence, mapped_position, source}}
    """
    if not consurf_data or not contact_residues:
        return {}

    # Build lookup: {chain_id: {position_int: entry}}
    lookup = {}
    for chain_id, entries in consurf_data.items():
        chain_lookup = {}
        for entry in entries:
            chain_lookup[entry["pos"]] = entry
        lookup[chain_id] = chain_lookup

    result = {}
    for (chain_id, res_name, res_id) in contact_residues:
        if chain_id not in lookup:
            continue

        chain_lookup = lookup[chain_id]
        m = re.match(r"^(-?\d+)([A-Za-z]?)$", str(res_id).strip())
        if not m:
            continue

        res_num = int(m.group(1))
        ins_code = m.group(2).upper() if m.group(2) else ""

        if res_num in chain_lookup:
            entry = chain_lookup[res_num]
            key = f"{chain_id}:{res_name}{res_id}"
            confidence = "medium" if ins_code else "high"
            result[key] = {
                "score": entry["score"],
                "color": entry["color"],
                "confidence": confidence,
                "mapped_position": res_num,
                "source": "consurf_db",
            }

    return result
