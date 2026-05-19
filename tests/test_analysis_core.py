import os
import unittest

from analysis_core import (
    analyze_ligand_pocket,
    classify_interaction,
    classify_residue,
    format_ligand_suggestions,
    get_hotspot_residues,
    infer_element,
    list_ligands,
    residue_keys_to_json,
)


ROOT_DIR = os.path.dirname(os.path.dirname(__file__))


class AnalysisCoreTest(unittest.TestCase):
    def test_classify_residue(self):
        self.assertEqual(classify_residue("VAL"), "hydrophobic")
        self.assertEqual(classify_residue("SER"), "polar")
        self.assertEqual(classify_residue("LYS"), "positive")
        self.assertEqual(classify_residue("ASP"), "negative")
        self.assertEqual(classify_residue("UNK"), "other")

    def test_infer_element_from_pdb_atom_name(self):
        self.assertEqual(infer_element(" CA "), "C")
        self.assertEqual(infer_element(" OD1"), "O")
        self.assertEqual(infer_element(" N1 "), "N")

    def test_classify_interaction(self):
        ligand_atom = {"atom_name": "C1"}
        protein_atom = {"atom_name": "CB", "res_name": "VAL"}

        interaction_type, color = classify_interaction(ligand_atom, protein_atom, 3.8)

        self.assertEqual(interaction_type, "hydrophobic contact")
        self.assertEqual(color, "#facc15")

    def test_analyze_ligand_pocket_finds_contacts(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        contacts, counts, interpretation, interactions = analyze_ligand_pocket(
            pdb_path,
            "MK1"
        )

        self.assertIsNotNone(contacts)
        self.assertGreater(len(contacts), 0)
        self.assertEqual(sum(counts.values()), len(contacts))
        self.assertIn("MK1 binding pocket", interpretation)
        self.assertLessEqual(len(interactions), 15)
        self.assertIn("chain_id", interactions[0])
        self.assertIn("res_name", interactions[0])
        self.assertIn("res_id", interactions[0])
        self.assertIn("atom_name", interactions[0])
        self.assertIn("element", interactions[0])
        self.assertIn("interaction_type", interactions[0])
        self.assertIn("color", interactions[0])
        self.assertIn("ligand_atom", interactions[0])

    def test_list_ligands_suggests_available_hetatm_ligands(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        ligands = list_ligands(pdb_path)
        suggestion_text = format_ligand_suggestions(ligands)

        self.assertTrue(any(ligand["res_name"] == "MK1" for ligand in ligands))
        self.assertIn("Available ligand candidates:", suggestion_text)
        self.assertIn("MK1", suggestion_text)

    def test_residue_keys_to_json_uses_standard_schema(self):
        residues = residue_keys_to_json([("A", "VAL", "82")])

        self.assertEqual(
            residues[0],
            {
                "chain_id": "A",
                "res_name": "VAL",
                "res_id": "82",
                "atom_name": "",
                "element": ""
            }
        )

    def test_hotspots_include_interaction_type(self):
        interactions = [
            {
                "chain_id": "A",
                "res_name": "VAL",
                "res_id": "82",
                "atom_name": "CB",
                "element": "C",
                "distance": 3.2,
                "interaction_type": "hydrophobic contact",
                "color": "#facc15"
            }
        ]

        hotspots = get_hotspot_residues(interactions)

        self.assertEqual(hotspots[0]["chain_id"], "A")
        self.assertEqual(hotspots[0]["res_name"], "VAL")
        self.assertEqual(hotspots[0]["res_id"], "82")
        self.assertEqual(hotspots[0]["atom_name"], "CB")
        self.assertEqual(hotspots[0]["element"], "C")
        self.assertEqual(hotspots[0]["interaction_type"], "hydrophobic contact")
        self.assertEqual(hotspots[0]["color"], "#facc15")

    def test_analyze_ligand_pocket_returns_none_for_missing_ligand(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        contacts, counts, interpretation, interactions = analyze_ligand_pocket(
            pdb_path,
            "XXX"
        )

        self.assertIsNone(contacts)
        self.assertIsNone(counts)
        self.assertIsNone(interpretation)
        self.assertIsNone(interactions)


if __name__ == "__main__":
    unittest.main()
