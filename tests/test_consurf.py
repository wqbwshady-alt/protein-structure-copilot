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
    _fetch_chain_mapping,
    _fetch_grades,
    _parse_grades_text,
    _parse_3latom,
)


ROOT_DIR = os.path.dirname(os.path.dirname(__file__))


class ConsurfCacheTest(unittest.TestCase):
    def setUp(self):
        os.makedirs(CACHE_DIR, exist_ok=True)

    def tearDown(self):
        cache_path = os.path.join(CACHE_DIR, "TEST.json")
        if os.path.exists(cache_path):
            os.remove(cache_path)

    def test_cache_roundtrip(self):
        payload = {"A": [{"pos": 1, "score": 7, "color": 8, "chain_id": "A", "residue_name": "MET", "insertion_code": ""}]}
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
        self.assertEqual(extract_pdb_id("", filename="RCSB_1HSG_abc12345.pdb"), "1HSG")

    def test_rcsb_filename_case(self):
        self.assertEqual(extract_pdb_id("", filename="rcsb_4hhb_ffffffff.pdb"), "4HHB")

    def test_no_rcsb_prefix(self):
        self.assertIsNone(extract_pdb_id("", filename="my_upload.pdb"))

    def test_empty(self):
        self.assertIsNone(extract_pdb_id("", filename=""))


class Parse3latomTest(unittest.TestCase):
    def test_standard(self):
        result = _parse_3latom("MET1:A")
        self.assertEqual(result["res_name"], "MET")
        self.assertEqual(result["pos"], 1)
        self.assertEqual(result["chain"], "A")
        self.assertEqual(result["ins_code"], "")

    def test_with_insertion_code(self):
        result = _parse_3latom("ASP25A:B")
        self.assertEqual(result["res_name"], "ASP")
        self.assertEqual(result["pos"], 25)
        self.assertEqual(result["chain"], "B")
        self.assertEqual(result["ins_code"], "A")

    def test_negative_position(self):
        result = _parse_3latom("GLY-1:A")
        self.assertEqual(result["pos"], -1)

    def test_chain_as_digit(self):
        result = _parse_3latom("ALA5:0")
        self.assertEqual(result["chain"], "0")

    def test_invalid(self):
        self.assertIsNone(_parse_3latom(""))
        self.assertIsNone(_parse_3latom("INVALID"))


class ParseGradesTextTest(unittest.TestCase):
    def test_valid_grades(self):
        text = (
            "1\tM\tMET1:A\t-1.610\t9\t-1.838,-1.537\t9,9\t94/300\tM,L,V,I\n"
            "2\tG\tGLY2:A\t0.523\t4\t0.201,0.845\t4,4\t85/300\tG,A\n"
            "3\tA\tALA3:A\t-0.120\t6\t-0.401,0.161\t6,6\t90/300\tA,S,T\n"
        )
        results = _parse_grades_text(text)
        self.assertEqual(len(results), 3)
        self.assertEqual(results[0]["pos"], 1)
        self.assertEqual(results[0]["score"], -1.610)
        self.assertEqual(results[0]["color"], 9)
        self.assertEqual(results[0]["residue_name"], "MET")
        self.assertEqual(results[0]["chain_id"], "A")

    def test_missing_density_skipped(self):
        text = "1\tM\t-\t0.5\t5\t0.2,0.8\t5,5\t50/100\tA,L\n"
        results = _parse_grades_text(text)
        self.assertIsNone(results)

    def test_empty_text(self):
        self.assertIsNone(_parse_grades_text(""))
        self.assertIsNone(_parse_grades_text("   \n"))

    def test_header_lines_ignored(self):
        text = "# Header line\n1\tM\tMET1:A\t-1.610\t9\t\n"
        results = _parse_grades_text(text)
        self.assertEqual(len(results), 1)


class MapConsurfToResiduesTest(unittest.TestCase):
    def setUp(self):
        self.consurf_data = {
            "A": [
                {"pos": 1, "score": -1.610, "color": 9, "chain_id": "A", "residue_name": "ASP", "insertion_code": ""},
                {"pos": 2, "score": 0.523, "color": 4, "chain_id": "A", "residue_name": "GLY", "insertion_code": ""},
                {"pos": 25, "score": -0.120, "color": 7, "chain_id": "A", "residue_name": "ILE", "insertion_code": ""},
            ],
            "B": [
                {"pos": 10, "score": 1.2, "color": 3, "chain_id": "B", "residue_name": "PHE", "insertion_code": ""},
            ],
        }
        self.contact_residues = {
            ("A", "ASP", "1"): {},
            ("A", "GLY", "2"): {},
            ("A", "ILE", "25"): {},
            ("A", "LYS", "99"): {},
            ("B", "PHE", "10"): {},
            ("C", "VAL", "1"): {},
        }

    def test_exact_match_high_confidence(self):
        result = map_consurf_to_residues(self.consurf_data, self.contact_residues)
        self.assertIn("A:ASP1", result)
        self.assertEqual(result["A:ASP1"]["score"], -1.610)
        self.assertEqual(result["A:ASP1"]["confidence"], "high")

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
        payload = {"A": [{"pos": 1, "score": -1.610, "color": 9, "chain_id": "A"}]}
        _save_cache("CACH", payload)
        result = query_consurf_db("CACH")
        self.assertIsNotNone(result)
        self.assertEqual(result["A"][0]["score"], -1.610)
        os.remove(os.path.join(CACHE_DIR, "CACH.json"))

    @mock.patch("consurf.requests.get")
    def test_full_query_success(self, mock_get):
        # Mock Step 1: chain_selection
        # Mock Step 2: grades file
        def side_effect(url, **kwargs):
            m = mock.MagicMock()
            if "chain_selection" in url:
                m.status_code = 200
                m.text = '<option value="A ABCDE">chain A</option>'
            elif "consurf_summary" in url:
                m.status_code = 200
                m.text = "1\tM\tMET1:A\t-1.610\t9\t-1.838,-1.537\t9,9\t94/300\tM,L,V,I\n"
            else:
                m.status_code = 404
                m.text = ""
            return m

        mock_get.side_effect = side_effect
        result = query_consurf_db("3DEF")
        self.assertIsNotNone(result)
        self.assertIn("A", result)
        self.assertEqual(result["A"][0]["score"], -1.610)

    @mock.patch("consurf.requests.get")
    def test_chain_selection_no_chains(self, mock_get):
        m = mock.MagicMock()
        m.status_code = 200
        m.text = "No chains found for 1XYZ"
        mock_get.return_value = m
        result = query_consurf_db("1XYZ")
        self.assertIsNone(result)

    @mock.patch("consurf.requests.get")
    def test_downtime_page(self, mock_get):
        m = mock.MagicMock()
        m.status_code = 200
        m.text = "<html><h1>Temporary Downtime</h1></html>"
        mock_get.return_value = m
        result = query_consurf_db("1XYZ")
        self.assertIsNone(result)

    @mock.patch("consurf.requests.get")
    def test_connection_error(self, mock_get):
        import requests as req_mod
        mock_get.side_effect = req_mod.ConnectionError
        result = query_consurf_db("1XYZ")
        self.assertIsNone(result)

    @mock.patch("consurf.requests.get")
    def test_http_404(self, mock_get):
        m = mock.MagicMock()
        m.status_code = 404
        m.text = ""
        mock_get.return_value = m
        result = query_consurf_db("1XYZ")
        self.assertIsNone(result)
