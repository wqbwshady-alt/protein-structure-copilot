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
_db_schema_ready = False
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


def _database_url():
    return os.getenv("PSC_DATABASE_URL") or os.getenv("DATABASE_URL")


def _use_database():
    backend = os.getenv("PSC_STATS_BACKEND", "").strip().lower()
    if backend == "file":
        return False
    if backend in {"postgres", "postgresql"}:
        return bool(_database_url())
    return bool(_database_url()) and not os.getenv("PSC_STATS_FILE")


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
    has_legacy_counter_only = (
        stats["total_analyses"] == LEGACY_LOCAL_DEV_TOTAL and
        str(stats.get("last_updated") or "").startswith("2026-05-20T13:28:") and
        not stats["recent_analyses"]
    )

    if legacy_count == 0 and not has_legacy_counter_only:
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


def _connect_db():
    import psycopg

    return psycopg.connect(_database_url())


def _ensure_database_store():
    global _db_schema_ready
    if _db_schema_ready:
        return

    with _connect_db() as conn:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS usage_stats (
                    id INTEGER PRIMARY KEY CHECK (id = 1),
                    total_analyses INTEGER NOT NULL DEFAULT 0,
                    last_updated TIMESTAMPTZ
                )
            """)
            cur.execute("""
                CREATE TABLE IF NOT EXISTS recent_analyses (
                    id BIGSERIAL PRIMARY KEY,
                    timestamp TIMESTAMPTZ NOT NULL,
                    analysis_type TEXT NOT NULL,
                    pdb_id TEXT NOT NULL DEFAULT '',
                    pdb_name TEXT NOT NULL DEFAULT '',
                    ligand_name TEXT NOT NULL DEFAULT '',
                    mutation TEXT NOT NULL DEFAULT '',
                    source TEXT NOT NULL DEFAULT 'local',
                    mode TEXT NOT NULL DEFAULT 'ligand',
                    wt_pdb_id TEXT NOT NULL DEFAULT '',
                    mutant_pdb_id TEXT NOT NULL DEFAULT ''
                )
            """)
            cur.execute("""
                INSERT INTO usage_stats (id, total_analyses, last_updated)
                VALUES (1, 0, NULL)
                ON CONFLICT (id) DO NOTHING
            """)
    _db_schema_ready = True


def _isoformat(value):
    if value is None:
        return None
    if hasattr(value, "isoformat"):
        return value.isoformat()
    return str(value)


def _entry_from_db_row(row):
    return _normalize_entry({
        "timestamp": _isoformat(row["timestamp"]),
        "analysis_type": row["analysis_type"],
        "pdb_id": row["pdb_id"],
        "pdb_name": row["pdb_name"],
        "ligand_name": row["ligand_name"],
        "mutation": row["mutation"],
        "source": row["source"],
        "mode": row["mode"],
        "wt_pdb_id": row["wt_pdb_id"],
        "mutant_pdb_id": row["mutant_pdb_id"],
    })


def _get_database_stats():
    from psycopg.rows import dict_row

    _ensure_database_store()
    with _connect_db() as conn:
        with conn.cursor(row_factory=dict_row) as cur:
            cur.execute("""
                SELECT total_analyses, last_updated
                FROM usage_stats
                WHERE id = 1
            """)
            row = cur.fetchone()

    return {
        "total_analyses": int(row["total_analyses"] or 0),
        "last_updated": _isoformat(row["last_updated"]),
    }


def _get_database_recent():
    from psycopg.rows import dict_row

    _ensure_database_store()
    with _connect_db() as conn:
        with conn.cursor(row_factory=dict_row) as cur:
            cur.execute("""
                SELECT timestamp, analysis_type, pdb_id, pdb_name, ligand_name,
                       mutation, source, mode, wt_pdb_id, mutant_pdb_id
                FROM recent_analyses
                ORDER BY timestamp DESC, id DESC
                LIMIT %s
            """, (MAX_RECENT,))
            rows = cur.fetchall()

    return [_entry_from_db_row(row) for row in rows]


def _record_database_analysis(metadata):
    _ensure_database_store()
    now = datetime.now(timezone.utc)
    entry = _normalize_entry({
        "timestamp": now.isoformat(),
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

    with _connect_db() as conn:
        with conn.cursor() as cur:
            cur.execute("""
                UPDATE usage_stats
                SET total_analyses = total_analyses + 1,
                    last_updated = %s
                WHERE id = 1
            """, (now,))
            cur.execute("""
                INSERT INTO recent_analyses (
                    timestamp, analysis_type, pdb_id, pdb_name, ligand_name,
                    mutation, source, mode, wt_pdb_id, mutant_pdb_id
                )
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            """, (
                now,
                entry["analysis_type"],
                entry["pdb_id"],
                entry["pdb_name"],
                entry["ligand_name"],
                entry["mutation"],
                entry["source"],
                entry["mode"],
                entry.get("wt_pdb_id", ""),
                entry.get("mutant_pdb_id", ""),
            ))
            cur.execute("""
                DELETE FROM recent_analyses
                WHERE id NOT IN (
                    SELECT id
                    FROM recent_analyses
                    ORDER BY timestamp DESC, id DESC
                    LIMIT %s
                )
            """, (MAX_RECENT,))

    return {
        "total_analyses": _get_database_stats()["total_analyses"],
        "last_updated": now.isoformat(),
        "recent_analyses": _get_database_recent(),
    }


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
        if _use_database():
            _ensure_database_store()
            return {
                **_get_database_stats(),
                "recent_analyses": _get_database_recent(),
            }

        stats = _read()
        _write(stats)
        return stats


def get_stats():
    with _lock:
        if _use_database():
            return _get_database_stats()

        s = _read()
        return {
            "total_analyses": s["total_analyses"],
            "last_updated": s["last_updated"],
        }


def get_recent():
    with _lock:
        if _use_database():
            return _get_database_recent()

        return _read()["recent_analyses"][:MAX_RECENT]


def record_analysis(metadata):
    with _lock:
        if _use_database():
            return _record_database_analysis(metadata)

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
