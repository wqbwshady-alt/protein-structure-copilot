import json
import os
import threading

STATS_FILE = os.path.join(os.path.dirname(os.path.dirname(__file__)), "stats.json")
_lock = threading.Lock()


def _defaults():
    return {
        "total_analyses": 0,
        "structures_analyzed": 0,
        "mutation_scans": 0,
        "comparisons": 0,
    }


def _read():
    if not os.path.exists(STATS_FILE):
        return _defaults()
    try:
        with open(STATS_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
        defaults = _defaults()
        for key in defaults:
            defaults[key] = data.get(key, 0)
        return defaults
    except (json.JSONDecodeError, IOError):
        return _defaults()


def _write(stats):
    with open(STATS_FILE, "w", encoding="utf-8") as f:
        json.dump(stats, f)


def get_stats():
    with _lock:
        return _read()


def increment_analysis():
    with _lock:
        stats = _read()
        stats["total_analyses"] += 1
        stats["structures_analyzed"] += 1
        _write(stats)
        return stats


def increment_mutation_scan():
    with _lock:
        stats = _read()
        stats["mutation_scans"] += 1
        _write(stats)
        return stats


def increment_comparison():
    with _lock:
        stats = _read()
        stats["comparisons"] += 1
        _write(stats)
        return stats
