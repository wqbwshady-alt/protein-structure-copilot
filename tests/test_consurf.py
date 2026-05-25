import json
import os
import time
import unittest
from unittest import mock

from consurf import (
    CACHE_DIR,
    extract_pdb_id,
    map_consurf_to_residues,
    query_consurf_db,
    _load_cache,
    _save_cache,
    _parse_consurf_response,
)


ROOT_DIR = os.path.dirname(os.path.dirname(__file__))
TEST_PDB = os.path.join(ROOT_DIR, "uploads", "RCSB_1HSG_abc12345.pdb")


class ConsurfCacheTest(unittest.TestCase):
    def setUp(self):
        self.test_pdb_id = "TEST"
        os.makedirs(CACHE_DIR, exist_ok=True)

    def tearDown(self):
        cache_path = os.path.join(CACHE_DIR, "TEST.json")
        if os.path.exists(cache_path):
            os.remove(cache_path)

    def test_cache_roundtrip(self):
        payload = {"A": [{"pos": 1, "score": 7}]}
        _save_cache("TEST", payload)
        loaded = _load_cache("TEST")
        self.assertIsNotNone(loaded)
        self.assertEqual(loaded["A"][0]["score"], 7)

    def test_cache_ttl_expired(self):
        payload = {"A": []}
        _save_cache("TEST", payload)
        cache_path = os.path.join(CACHE_DIR, "TEST.json")
        with open(cache_path, "r") as f:
            data = json.load(f)
        data["_cached_at"] = int(time.time()) - 31 * 24 * 3600
        with open(cache_path, "w") as f:
            json.dump(data, f)
        loaded = _load_cache("TEST")
        self.assertIsNone(loaded)

    def test_load_cache_missing(self):
        loaded = _load_cache("NONEXISTENT_PDB_ID_9999")
        self.assertIsNone(loaded)


class ExtractPdbIdTest(unittest.TestCase):
    def test_rcsb_filename(self):
        pdb_id = extract_pdb_id("", filename="RCSB_1HSG_abc12345.pdb")
        self.assertEqual(pdb_id, "1HSG")

    def test_rcsb_filename_case(self):
        pdb_id = extract_pdb_id("", filename="rcsb_4hhb_ffffffff.pdb")
        self.assertEqual(pdb_id, "4HHB")

    def test_no_rcsb_prefix(self):
        pdb_id = extract_pdb_id("", filename="my_upload.pdb")
        self.assertIsNone(pdb_id)

    def test_empty(self):
        pdb_id = extract_pdb_id("", filename="")
        self.assertIsNone(pdb_id)


class ParseConsurfResponseTest(unittest.TestCase):
    def test_standard_format(self):
        raw = {
            "results": {
                "A": [
                    {"pos": 1, "score": 8, "color": 9},
                    {"pos": 2, "score": 4, "color": 5},
                ],
                "B": [
                    {"pos": 1, "score": 7, "color": 8},
                ],
            }
        }
        parsed = _parse_consurf_response(raw)
        self.assertIsNotNone(parsed)
        self.assertEqual(len(parsed["A"]), 2)
        self.assertEqual(parsed["A"][0]["pos"], 1)
        self.assertEqual(parsed["A"][0]["score"], 8.0)
        self.assertEqual(parsed["A"][0]["color"], 9)

    def test_float_score(self):
        raw = {"results": {"A": [{"pos": 10, "score": 6.5, "color": 7}]}}
        parsed = _parse_consurf_response(raw)
        self.assertEqual(parsed["A"][0]["score"], 6.5)

    def test_string_position(self):
        raw = {"results": {"A": [{"pos": "42", "score": 3}]}}
        parsed = _parse_consurf_response(raw)
        self.assertEqual(parsed["A"][0]["pos"], 42)

    def test_alternative_keys(self):
        raw = {
            "results": {
                "A": [
                    {"position": 5, "score": 9},
                    {"residue_number": 6, "score": 8},
                ]
            }
        }
        parsed = _parse_consurf_response(raw)
        self.assertEqual(len(parsed["A"]), 2)
        self.assertEqual(parsed["A"][0]["pos"], 5)
        self.assertEqual(parsed["A"][1]["pos"], 6)

    def test_no_results_key(self):
        raw = {
            "A": [{"pos": 1, "score": 7}],
            "B": [{"pos": 2, "score": 5}],
        }
        parsed = _parse_consurf_response(raw)
        self.assertIsNotNone(parsed)
        self.assertEqual(len(parsed["A"]), 1)
        self.assertEqual(len(parsed["B"]), 1)

    def test_missing_position(self):
        raw = {"results": {"A": [{"score": 5}]}}
        parsed = _parse_consurf_response(raw)
        self.assertIsNone(parsed)

    def test_empty_results(self):
        raw = {"results": {}}
        parsed = _parse_consurf_response(raw)
        self.assertIsNone(parsed)

    def test_non_dict(self):
        parsed = _parse_consurf_response([1, 2, 3])
        self.assertIsNone(parsed)

    def test_none(self):
        parsed = _parse_consurf_response(None)
        self.assertIsNone(parsed)


class MapConsurfToResiduesTest(unittest.TestCase):
    def setUp(self):
        self.consurf_data = {
            "A": [
                {"pos": 1, "score": 8, "color": 9},
                {"pos": 2, "score": 4, "color": 5},
                {"pos": 25, "score": 7, "color": 8},
            ],
            "B": [
                {"pos": 10, "score": 3, "color": 4},
            ],
        }
        self.contact_residues = {
            ("A", "ASP", "1"): {},
            ("A", "GLY", "2"): {},
            ("A", "ILE", "25"): {},
            ("A", "LYS", "99"): {},          # not in ConSurf
            ("B", "PHE", "10"): {},
            ("C", "VAL", "1"): {},           # chain not in ConSurf
        }

    def test_exact_match_high_confidence(self):
        result = map_consurf_to_residues(self.consurf_data, self.contact_residues)
        self.assertIn("A:ASP1", result)
        self.assertEqual(result["A:ASP1"]["score"], 8)
        self.assertEqual(result["A:ASP1"]["confidence"], "high")
        self.assertEqual(result["A:ASP1"]["source"], "consurf_db")

    def test_no_match_for_missing(self):
        result = map_consurf_to_residues(self.consurf_data, self.contact_residues)
        self.assertNotIn("A:LYS99", result)
        self.assertNotIn("C:VAL1", result)

    def test_insertion_code_medium_confidence(self):
        contact_residues = {("A", "ASP", "25A"): {}}
        result = map_consurf_to_residues(self.consurf_data, contact_residues)
        self.assertIn("A:ASP25A", result)
        self.assertEqual(result["A:ASP25A"]["confidence"], "medium")

    def test_empty_inputs(self):
        self.assertEqual(map_consurf_to_residues(None, {}), {})
        self.assertEqual(map_consurf_to_residues({}, None), {})
        self.assertEqual(map_consurf_to_residues({}, {}), {})

    def test_cross_chain_mapping(self):
        result = map_consurf_to_residues(self.consurf_data, self.contact_residues)
        self.assertIn("B:PHE10", result)
        self.assertEqual(result["B:PHE10"]["score"], 3)


class QueryConsurfDbTest(unittest.TestCase):
    def setUp(self):
        os.makedirs(CACHE_DIR, exist_ok=True)

    def tearDown(self):
        for pdb_id in ["1XYZ", "2ABC", "3DEF"]:
            cache_path = os.path.join(CACHE_DIR, f"{pdb_id}.json")
            if os.path.exists(cache_path):
                os.remove(cache_path)

    def test_invalid_pdb_id(self):
        self.assertIsNone(query_consurf_db(""))
        self.assertIsNone(query_consurf_db("INVALID"))
        self.assertIsNone(query_consurf_db("12345"))

    def test_cache_hit(self):
        pdb_id = "CACH"
        _save_cache(pdb_id, {"A": [{"pos": 1, "score": 9}]})
        result = query_consurf_db(pdb_id)
        self.assertIsNotNone(result)
        self.assertEqual(result["A"][0]["score"], 9)
        os.remove(os.path.join(CACHE_DIR, "CACH.json"))

    @mock.patch("consurf.requests.get")
    def test_api_success(self, mock_get):
        mock_resp = mock.MagicMock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {
            "results": {"A": [{"pos": 1, "score": 7}]}
        }
        mock_get.return_value = mock_resp

        result = query_consurf_db("3DEF")
        self.assertIsNotNone(result)
        self.assertEqual(result["A"][0]["score"], 7)

    @mock.patch("consurf.requests.get")
    def test_api_404(self, mock_get):
        mock_resp = mock.MagicMock()
        mock_resp.status_code = 404
        mock_get.return_value = mock_resp
        result = query_consurf_db("2ABC")
        self.assertIsNone(result)

    @mock.patch("consurf.requests.get")
    def test_api_timeout(self, mock_get):
        import requests as req_mod
        mock_get.side_effect = req_mod.Timeout
        result = query_consurf_db("1XYZ")
        self.assertIsNone(result)

    @mock.patch("consurf.requests.get")
    def test_api_connection_error(self, mock_get):
        import requests as req_mod
        mock_get.side_effect = req_mod.ConnectionError
        result = query_consurf_db("1XYZ")
        self.assertIsNone(result)
