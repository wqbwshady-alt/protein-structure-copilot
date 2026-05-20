import json
import os
import re
import threading
from datetime import datetime, timezone
from uuid import uuid4


ROOT_DIR = os.path.dirname(os.path.dirname(__file__))
DEFAULT_STATS_FILE = os.path.join(ROOT_DIR, "stats.json")
MAX_RECENT = 20
_lock = threading.Lock()
LEGACY_LOCAL_DEV_TOTAL = 19
LEGACY_LOCAL_DEV_PDB_NAMES = {
    "MUTSCAN_42cbe176c1ff4ad7b88bbe1378b0438d_1HSG.pdb",
    "RCSB_1HSG_test.pdb",
    "RCSB_1HSG_test.pdb vs RCSB_1HSG_V82A_test.pdb",
    "c0d96e1d82b4421094e717703556c2b6_1HSG.pdb",
    "9c66e9093b6c420c8fc99d10c8ab5b66_protein_only.pdb",
}


def _stats_file_path():
    explicit_file = os.getenv("PSC_STATS_FILE") or os.getenv("STATS_FILE")
    if explicit_file:
        return explicit_file

    persistent_dir = os.getenv("PSC_STATS_DIR") or os.getenv("RENDER_DISK_PATH")
    if persistent_dir:
        return os.path.join(persistent_dir, "stats.json")

    render_disk = "/var/data"
    if os.path.isdir(render_disk) and os.access(render_disk, os.W_OK):
        return os.path.join(render_disk, "protein_structure_copilot_stats.json")

    return DEFAULT_STATS_FILE


def _defaults():
    return {
        "total_analyses": 0,
        "last_updated": None,
        "recent_analyses": [],
    }


def _display_structure_name(raw_name):
    if not raw_name:
        return ""

    if " vs " in raw_name:
        return " vs ".join(_display_structure_name(part) for part in raw_name.split(" vs "))

    name = os.path.basename(str(raw_name))
    rcsb_match = re.search(r"RCSB_([A-Za-z0-9]{4})", name)
    if rcsb_match:
        return rcsb_match.group(1).upper()

    name = re.sub(r"^(WT_|MUT_|MUTSCAN_)", "", name)
    name = re.sub(r"^[0-9a-f]{32}_", "", name, flags=re.IGNORECASE)
    name = re.sub(r"^[0-9a-f]{8}_", "", name, flags=re.IGNORECASE)
    return name


def _normalize_entry(entry):
    analysis_type = entry.get("analysis_type") or entry.get("type") or "single"
    pdb_name = entry.get("pdb_name") or entry.get("pdb_id") or ""
    pdb_id = entry.get("pdb_id") or _display_structure_name(pdb_name)

    normalized = {
        "timestamp": entry.get("timestamp"),
        "analysis_type": analysis_type,
        "type": analysis_type,
        "pdb_id": pdb_id,
        "pdb_name": pdb_name,
        "ligand_name": entry.get("ligand_name", ""),
        "mutation": entry.get("mutation", ""),
        "source": entry.get("source", "local"),
        "mode": entry.get("mode", "ligand"),
    }

    for optional_key in ("wt_pdb_id", "mutant_pdb_id"):
        if entry.get(optional_key):
            normalized[optional_key] = entry[optional_key]

    return normalized


def _normalize_stats(data):
    if not isinstance(data, dict):
        data = {}

    recent = data.get("recent_analyses")
    if recent is None:
        recent = data.get("recent", [])
    if not isinstance(recent, list):
        recent = []

    normalized = {
        "total_analyses": int(data.get("total_analyses") or 0),
        "last_updated": data.get("last_updated"),
        "recent_analyses": [_normalize_entry(item) for item in recent[:MAX_RECENT] if isinstance(item, dict)],
    }
    return _remove_legacy_local_dev_artifacts(normalized)


def _is_legacy_local_dev_entry(entry):
    timestamp = str(entry.get("timestamp") or "")
    pdb_name = entry.get("pdb_name") or ""
    return timestamp.startswith("2026-05-20T13:28:") and pdb_name in LEGACY_LOCAL_DEV_PDB_NAMES


def _remove_legacy_local_dev_artifacts(stats):
    legacy_count = sum(1 for entry in stats["recent_analyses"] if _is_legacy_local_dev_entry(entry))
    if legacy_count == 0:
        return stats

    cleaned_recent = [
        entry for entry in stats["recent_analyses"]
        if not _is_legacy_local_dev_entry(entry)
    ][:MAX_RECENT]

    cleaned_total = max(0, stats["total_analyses"] - LEGACY_LOCAL_DEV_TOTAL)
    return {
        "total_analyses": cleaned_total,
        "last_updated": cleaned_recent[0]["timestamp"] if cleaned_recent else None,
        "recent_analyses": cleaned_recent,
    }


def _read_file(path):
    if not os.path.exists(path):
        return _defaults()
    try:
        with open(path, "r", encoding="utf-8") as f:
            return _normalize_stats(json.load(f))
    except (json.JSONDecodeError, OSError, ValueError):
        return _defaults()


def _read():
    return _read_file(_stats_file_path())


def _write(stats):
    path = _stats_file_path()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    normalized = _normalize_stats(stats)
    temp_path = f"{path}.{uuid4().hex}.tmp"
    with open(temp_path, "w", encoding="utf-8") as f:
        json.dump(normalized, f, ensure_ascii=False, indent=2)
        f.write("\n")
    os.replace(temp_path, path)


def ensure_stats_store():
    with _lock:
        stats = _read()
        _write(stats)
        return stats


def get_stats():
    with _lock:
        s = _read()
        return {
            "total_analyses": s["total_analyses"],
            "last_updated": s["last_updated"],
        }


def get_recent():
    with _lock:
        return _read()["recent_analyses"][:MAX_RECENT]


def record_analysis(metadata):
    with _lock:
        stats = _read()
        stats["total_analyses"] += 1
        now = datetime.now(timezone.utc).isoformat()
        stats["last_updated"] = now

        entry = _normalize_entry({
            "timestamp": now,
            "analysis_type": metadata.get("analysis_type") or metadata.get("type") or "single",
            "pdb_id": metadata.get("pdb_id"),
            "pdb_name": metadata.get("pdb_name", ""),
            "ligand_name": metadata.get("ligand_name", ""),
            "mutation": metadata.get("mutation", ""),
            "source": metadata.get("source", "local"),
            "mode": metadata.get("mode", "ligand"),
            "wt_pdb_id": metadata.get("wt_pdb_id"),
            "mutant_pdb_id": metadata.get("mutant_pdb_id"),
        })
        stats["recent_analyses"].insert(0, entry)
        stats["recent_analyses"] = stats["recent_analyses"][:MAX_RECENT]

        _write(stats)
        return stats
