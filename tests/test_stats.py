import json
import os
import tempfile
import unittest
from unittest.mock import patch

from services.stats import ensure_stats_store, get_recent, get_stats, record_analysis


class StatsServiceTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.stats_path = os.path.join(self.tmpdir.name, "stats.json")
        self.seed_path = os.path.join(self.tmpdir.name, "seed.json")
        self.env = patch.dict(os.environ, {
            "PSC_STATS_FILE": self.stats_path,
            "PSC_STATS_SEED_FILE": self.seed_path,
        })
        self.env.start()
        self.addCleanup(self.env.stop)

    def test_initializes_and_persists_lifetime_counter(self):
        ensure_stats_store()
        self.assertTrue(os.path.exists(self.stats_path))
        self.assertEqual(get_stats()["total_analyses"], 0)

        record_analysis({
            "analysis_type": "single",
            "pdb_name": "RCSB_7VV4_test.pdb",
            "ligand_name": "CLR",
            "source": "rcsb",
            "mode": "ligand",
        })

        with open(self.stats_path, "r", encoding="utf-8") as f:
            raw = json.load(f)

        self.assertEqual(raw["total_analyses"], 1)
        self.assertEqual(get_stats()["total_analyses"], 1)
        self.assertEqual(get_recent()[0]["pdb_id"], "7VV4")
        self.assertEqual(get_recent()[0]["analysis_type"], "single")

    def test_recent_analyses_are_limited_to_twenty(self):
        ensure_stats_store()

        for idx in range(25):
            record_analysis({
                "analysis_type": "mutation",
                "pdb_name": f"local_{idx}.pdb",
                "mutation": "R273H",
                "source": "local",
                "mode": "mutation",
            })

        self.assertEqual(get_stats()["total_analyses"], 25)
        self.assertEqual(len(get_recent()), 20)
        self.assertEqual(get_recent()[0]["pdb_name"], "local_24.pdb")

    def test_missing_runtime_store_starts_from_seed(self):
        with open(self.seed_path, "w", encoding="utf-8") as f:
            json.dump({
                "total_analyses": 19,
                "last_updated": "2026-05-20T13:28:37+00:00",
                "recent_analyses": [{
                    "timestamp": "2026-05-20T13:28:37+00:00",
                    "analysis_type": "single",
                    "pdb_id": "7VV4",
                    "pdb_name": "RCSB_7VV4_seed.pdb",
                    "ligand_name": "CLR",
                    "source": "rcsb",
                    "mode": "ligand",
                }],
            }, f)

        ensure_stats_store()

        self.assertEqual(get_stats()["total_analyses"], 19)
        self.assertEqual(get_recent()[0]["pdb_id"], "7VV4")

        record_analysis({
            "analysis_type": "mutation",
            "pdb_name": "RCSB_1TSR_test.pdb",
            "ligand_name": "ZN",
            "mutation": "R273H",
            "source": "rcsb",
            "mode": "mutation",
        })

        self.assertEqual(get_stats()["total_analyses"], 20)
        self.assertEqual(get_recent()[0]["analysis_type"], "mutation")


if __name__ == "__main__":
    unittest.main()
