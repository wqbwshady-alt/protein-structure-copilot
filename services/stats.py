import json
import os
import threading
from datetime import datetime, timezone

STATS_FILE = os.path.join(os.path.dirname(os.path.dirname(__file__)), "stats.json")
MAX_RECENT = 20
_lock = threading.Lock()


def _defaults():
    return {
        "total_analyses": 0,
        "last_updated": None,
        "recent": [],
    }


def _read():
    if not os.path.exists(STATS_FILE):
        return _defaults()
    try:
        with open(STATS_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
    except (json.JSONDecodeError, IOError):
        return _defaults()

    return {
        "total_analyses": data.get("total_analyses", 0),
        "last_updated": data.get("last_updated"),
        "recent": data.get("recent", [])[:MAX_RECENT],
    }


def _write(stats):
    with open(STATS_FILE, "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False)


def get_stats():
    with _lock:
        s = _read()
        return {
            "total_analyses": s["total_analyses"],
            "last_updated": s["last_updated"],
        }


def get_recent():
    with _lock:
        return _read()["recent"]


def record_analysis(metadata):
    with _lock:
        stats = _read()
        stats["total_analyses"] += 1
        now = datetime.now(timezone.utc).isoformat()
        stats["last_updated"] = now

        entry = {
            "timestamp": now,
            "type": metadata.get("type", "single"),
            "pdb_name": metadata.get("pdb_name", ""),
            "ligand_name": metadata.get("ligand_name", ""),
            "mutation": metadata.get("mutation", ""),
            "source": metadata.get("source", "local"),
            "mode": metadata.get("mode", "ligand"),
        }
        stats["recent"].insert(0, entry)
        stats["recent"] = stats["recent"][:MAX_RECENT]

        _write(stats)
        return stats
