from io import BytesIO
import os
import tempfile
import unittest
from unittest.mock import patch

import app as app_module


ROOT_DIR = os.path.dirname(os.path.dirname(__file__))


class AppTest(unittest.TestCase):
    def setUp(self):
        app_module.app.config["TESTING"] = True
        self.client = app_module.app.test_client()

    def test_index_loads(self):
        response = self.client.get("/")

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Protein Structure Copilot", response.data)
        self.assertIn(b"Select an analysis mode above", response.data)

    def test_health_route(self):
        response = self.client.get("/health")

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.get_json(), {"status": "ok"})

    def test_analyze_rejects_non_pdb_upload(self):
        response = self.client.post(
            "/analyze",
            data={
                "ligand_name": "MK1",
                "pdb_file": (BytesIO(b"not a pdb"), "notes.txt")
            },
            content_type="multipart/form-data"
        )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Only .pdb files are supported.", response.data)

    def test_analyze_rejects_pdb_without_atom_records(self):
        response = self.client.post(
            "/analyze",
            data={
                "ligand_name": "MK1",
                "pdb_file": (BytesIO(b"HEADER empty file\nEND\n"), "empty.pdb")
            },
            content_type="multipart/form-data"
        )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Uploaded file does not contain ATOM or HETATM records.", response.data)

    def test_analyze_missing_ligand_suggests_detected_ligands(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        with tempfile.TemporaryDirectory() as upload_dir:
            with patch.object(app_module, "UPLOAD_FOLDER", upload_dir):
                with open(pdb_path, "rb") as pdb_file:
                    response = self.client.post(
                        "/analyze",
                        data={
                            "ligand_name": "XXX",
                            "pdb_file": (BytesIO(pdb_file.read()), "1HSG.pdb")
                        },
                        content_type="multipart/form-data"
                    )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"No ligand named XXX found", response.data)
        self.assertIn(b"Available ligand candidates:", response.data)
        self.assertIn(b"MK1", response.data)

    def test_analyze_json_upload_uses_analyze_route(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        with tempfile.TemporaryDirectory() as upload_dir:
            with patch.object(app_module, "UPLOAD_FOLDER", upload_dir):
                with open(pdb_path, "rb") as pdb_file:
                    response = self.client.post(
                        "/analyze",
                        data={
                            "ligand_name": "MK1",
                            "pdb_file": (BytesIO(pdb_file.read()), "1HSG.pdb")
                        },
                        content_type="multipart/form-data",
                        headers={
                            "Accept": "application/json",
                            "X-Requested-With": "XMLHttpRequest"
                        }
                    )

        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertTrue(data["success"])
        self.assertIn("result_text", data)
        self.assertIn("pdb_url", data)

    def test_analyze_json_fetched_filename_uses_analyze_route(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        with tempfile.TemporaryDirectory() as upload_dir:
            fetched_name = "RCSB_1HSG_test.pdb"
            fetched_path = os.path.join(upload_dir, fetched_name)
            with open(pdb_path, "rb") as src, open(fetched_path, "wb") as dst:
                dst.write(src.read())

            with patch.object(app_module, "UPLOAD_FOLDER", upload_dir):
                response = self.client.post(
                    "/analyze",
                    data={
                        "ligand_name": "MK1",
                        "pdb_filename": fetched_name
                    },
                    headers={
                        "Accept": "application/json",
                        "X-Requested-With": "XMLHttpRequest"
                    }
                )

        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertTrue(data["success"])
        self.assertIn(f"/uploads/{fetched_name}", data["pdb_url"])

    def test_analyze_json_skip_ligand_runs_protein_only_analysis(self):
        pdb_text = (
            "ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N\n"
            "ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C\n"
            "ATOM      3  C   ALA A   2      15.685  12.755   5.133  1.00 10.00           C\n"
            "ATOM      4  CA  ASP B   1      12.000  11.000   7.000  1.00 10.00           C\n"
            "END\n"
        )

        with tempfile.TemporaryDirectory() as upload_dir:
            with patch.object(app_module, "UPLOAD_FOLDER", upload_dir):
                response = self.client.post(
                    "/analyze",
                    data={
                        "skipLigand": "true",
                        "ligand_name": "",
                        "pdb_file": (BytesIO(pdb_text.encode("utf-8")), "protein_only.pdb")
                    },
                    content_type="multipart/form-data",
                    headers={
                        "Accept": "application/json",
                        "X-Requested-With": "XMLHttpRequest"
                    }
                )

        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertTrue(data["success"])
        self.assertEqual(data["analysis_mode"], "protein_only")
        self.assertEqual(data["result_title"], "Protein-only structural overview")
        self.assertEqual(data["protein_summary"]["chain_count"], 2)
        self.assertEqual(data["protein_summary"]["residue_count"], 3)
        self.assertIn("pdb_url", data)

    def test_mutation_scan_route(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        with tempfile.TemporaryDirectory() as upload_dir:
            with patch.object(app_module, "UPLOAD_FOLDER", upload_dir):
                with open(pdb_path, "rb") as pdb_file:
                    response = self.client.post(
                        "/mutation_scan",
                        data={
                            "mutation_ligand_name": "MK1",
                            "mutation_text": "D25A",
                            "mutation_chain_id": "A",
                            "mutation_pdb_file": (BytesIO(pdb_file.read()), "1HSG.pdb")
                        },
                        content_type="multipart/form-data"
                    )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Mutation Impact Summary", response.data)
        self.assertIn(b"D25A", response.data)
        self.assertIn(b"Property Changes", response.data)

    def test_mutation_scan_accepts_loaded_pdb_filename(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        with tempfile.TemporaryDirectory() as upload_dir:
            loaded_name = "RCSB_1HSG_test.pdb"
            loaded_path = os.path.join(upload_dir, loaded_name)
            with open(pdb_path, "rb") as src, open(loaded_path, "wb") as dst:
                dst.write(src.read())

            with patch.object(app_module, "UPLOAD_FOLDER", upload_dir):
                response = self.client.post(
                    "/mutation_scan",
                    data={
                        "mutation_ligand_name": "MK1",
                        "mutation_text": "D25A",
                        "mutation_chain_id": "A",
                        "pdb_filename": loaded_name
                    }
                )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"Mutation Impact Summary", response.data)
        self.assertIn(b"D25A", response.data)

    def test_compare_accepts_loaded_structure_and_mutation(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        with tempfile.TemporaryDirectory() as upload_dir:
            loaded_name = "RCSB_1HSG_test.pdb"
            loaded_path = os.path.join(upload_dir, loaded_name)
            with open(pdb_path, "rb") as src, open(loaded_path, "wb") as dst:
                dst.write(src.read())

            with patch.object(app_module, "UPLOAD_FOLDER", upload_dir):
                response = self.client.post(
                    "/compare",
                    data={
                        "compare_ligand_name": "MK1",
                        "compare_mutation_text": "D25A",
                        "compare_chain_id": "A",
                        "pdb_filename": loaded_name
                    }
                )

        self.assertEqual(response.status_code, 200)
        self.assertIn(b"WT vs Mutant Heuristic Comparison", response.data)
        self.assertIn(b"D25A", response.data)


if __name__ == "__main__":
    unittest.main()
