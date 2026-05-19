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


if __name__ == "__main__":
    unittest.main()
