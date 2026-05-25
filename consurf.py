"""
ConSurf-DB integration for evolutionary conservation scores.

Provides:
- ConSurf-DB REST API query (GET only, no submission)
- Disk cache with 30-day TTL
- Conservative PDB position → ConSurf position mapping
- Graceful degradation: failures return None, never crash

Design constraints:
- API endpoint is configurable, response parsing is schema-tolerant
- Position mapping is conservative: exact match → high, insertion → medium, no match → skip
- ConSurf unavailable → return None → caller falls back to BLOSUM62 proxy
"""

import json
import os
import re
import time

import requests

# ---- Configuration ----

CONSURF_DB_URL = os.getenv(
    "CONSURF_DB_URL",
    "https://consurfdb.tau.ac.il/api/"
)

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache", "consurf")

# 30 days in seconds
CACHE_TTL = 30 * 24 * 3600

# HTTP timeout in seconds
REQUEST_TIMEOUT = 10


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
            json.dump({
                "_cached_at": int(time.time()),
                "_payload": payload,
            }, f)
    except IOError:
        pass


# ---- PDB ID extraction ----

def extract_pdb_id(pdb_path, filename=""):
    """Extract PDB ID from filename or PDB HEADER record.

    Returns 4-character PDB ID string, or None.
    """
    # Strategy 1: RCSB_ prefix in filename
    if filename:
        m = re.match(r"RCSB_([A-Za-z0-9]{4})_", filename, re.IGNORECASE)
        if m:
            return m.group(1).upper()

    # Strategy 2: HEADER record in PDB file
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("HEADER"):
                    # Standard PDB format: HEADER ... cols 63-66 = idCode
                    if len(line) >= 66:
                        id_code = line[62:66].strip().upper()
                        if re.fullmatch(r"[A-Z0-9]{4}", id_code):
                            return id_code
                    break
    except (IOError, OSError):
        pass

    return None


# ---- ConSurf-DB API query ----

def query_consurf_db(pdb_id):
    """Query ConSurf-DB for precomputed conservation scores.

    Args:
        pdb_id: 4-character PDB ID (e.g. "1HSG")

    Returns:
        dict: {chain_id: [{pos: int, score: int, color: int, ...}, ...]}
        None: if unavailable, failed, or no data

    Never raises.
    """
    if not pdb_id or not re.fullmatch(r"[A-Za-z0-9]{4}", pdb_id):
        return None

    pdb_id = pdb_id.upper()

    # Check cache first
    cached = _load_cache(pdb_id)
    if cached is not None:
        return cached

    url = f"{CONSURF_DB_URL.rstrip('/')}/{pdb_id}"
    try:
        resp = requests.get(url, timeout=REQUEST_TIMEOUT)
        if resp.status_code != 200:
            return None
        raw = resp.json()
    except (requests.RequestException, json.JSONDecodeError, ValueError):
        return None

    parsed = _parse_consurf_response(raw)
    if parsed:
        _save_cache(pdb_id, parsed)
    return parsed


def _parse_consurf_response(raw):
    """Parse ConSurf-DB API response with schema tolerance.

    Expected structure (flexible):
      {"results": {"A": [{"pos": 1, "score": 7, ...}, ...], "B": [...]}}

    Also handles:
      - Wrapped in top-level array or extra keys
      - Missing "results" key → tries top-level chain keys
      - Score values as int or float
      - Position as int or str

    Returns dict {chain_id: [per_position_dicts]} or None.
    """
    if not isinstance(raw, dict):
        return None

    results = raw.get("results")
    if not isinstance(results, dict):
        # Try: maybe top-level keys ARE chain IDs
        results = {
            k: v for k, v in raw.items()
            if isinstance(v, list) and len(k) <= 3 and k.strip()
        }
        if not results:
            return None

    parsed = {}
    for chain_id, positions in results.items():
        if not isinstance(positions, list):
            continue
        cleaned = []
        for entry in positions:
            if not isinstance(entry, dict):
                continue
            pos = entry.get("pos") or entry.get("position") or entry.get("residue_number")
            if pos is None:
                continue
            try:
                pos = int(pos)
            except (ValueError, TypeError):
                continue

            score = entry.get("score")
            if score is not None:
                try:
                    score = float(score)
                except (ValueError, TypeError):
                    score = None

            color = entry.get("color")
            if color is not None:
                try:
                    color = int(color)
                except (ValueError, TypeError):
                    color = None

            cleaned.append({
                "pos": pos,
                "score": score,
                "color": color,
                # Carry through any extra fields for future use
                "_raw": {k: v for k, v in entry.items()
                        if k not in ("pos", "position", "residue_number", "score", "color")},
            })

        if cleaned:
            parsed[chain_id] = cleaned

    return parsed if parsed else None


# ---- Mapping: ConSurf positions → PDB contact residues ----

def map_consurf_to_residues(consurf_data, contact_residues):
    """Map ConSurf-DB per-position scores to PDB contact residues.

    Conservative mapping strategy:
      - Exact position match (no insertion code) → confidence "high"
      - Position match with insertion code  → confidence "medium"
      - No match in ConSurf data → skip (returns None for that residue)

    Args:
        consurf_data: dict from query_consurf_db() — {chain_id: [{pos, score, ...}]}
        contact_residues: {(chain_id, res_name, res_id): atom_dict, ...}

    Returns:
        dict: {residue_key: {score, color, confidence, mapped_position}}
              Only includes residues that were successfully mapped.

    Never raises.
    """
    if not consurf_data or not contact_residues:
        return {}

    result = {}

    # Build lookup: {chain_id: {position_int: entry_dict}}
    lookup = {}
    for chain_id, positions in consurf_data.items():
        chain_lookup = {}
        for entry in positions:
            chain_lookup[entry["pos"]] = entry
        lookup[chain_id] = chain_lookup

    for (chain_id, res_name, res_id) in contact_residues:
        if chain_id not in lookup:
            continue

        chain_lookup = lookup[chain_id]

        # Parse PDB residue ID
        m = re.match(r"^(-?\d+)([A-Za-z]?)$", str(res_id).strip())
        if not m:
            continue
        res_num = int(m.group(1))
        ins_code = m.group(2).upper() if m.group(2) else ""

        # Try exact match first
        if res_num in chain_lookup and not ins_code:
            entry = chain_lookup[res_num]
            key = f"{chain_id}:{res_name}{res_id}"
            result[key] = {
                "score": entry["score"],
                "color": entry["color"],
                "confidence": "high",
                "mapped_position": res_num,
                "source": "consurf_db",
            }
            continue

        # With insertion code: map to same position number but mark confidence lower
        if res_num in chain_lookup and ins_code:
            entry = chain_lookup[res_num]
            key = f"{chain_id}:{res_name}{res_id}"
            result[key] = {
                "score": entry["score"],
                "color": entry["color"],
                "confidence": "medium",
                "mapped_position": res_num,
                "source": "consurf_db",
            }

    return result
