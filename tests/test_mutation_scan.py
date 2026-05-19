import json
import os
import subprocess
import sys
import unittest

from reports import build_mutation_scan_report
from services.mutation_scan import (
    analyze_mutation_scan,
    classify_interaction_impact,
    compare_properties,
    parse_mutation,
)


ROOT_DIR = os.path.dirname(os.path.dirname(__file__))


class MutationScanTest(unittest.TestCase):
    def test_parse_mutation(self):
        parsed = parse_mutation("R273H")

        self.assertEqual(parsed["mutation"], "R273H")
        self.assertEqual(parsed["wild_type_residue"], "ARG")
        self.assertEqual(parsed["residue_id"], "273")
        self.assertEqual(parsed["mutant_residue"], "HIS")

    def test_property_analysis(self):
        wild_type, mutant, changes = compare_properties("LYS", "ALA")

        self.assertEqual(wild_type["charge"], "positive")
        self.assertEqual(mutant["charge"], "neutral")
        self.assertEqual(changes["charge"], "positive -> neutral")
        self.assertEqual(changes["hydrophobicity"], "hydrophilic -> hydrophobic")

    def test_interaction_impact_predicts_electrostatic_loss(self):
        wild_type, mutant, changes = compare_properties("ASP", "ALA")
        interactions = [
            {
                "chain_id": "A",
                "res_name": "ASP",
                "res_id": "25",
                "atom_name": "OD1",
                "element": "O",
                "interaction_type": "charged / electrostatic",
                "distance": 2.8,
            }
        ]

        impact = classify_interaction_impact(interactions, changes, mutant)

        self.assertEqual(len(impact["possible_loss"]), 1)
        self.assertEqual(impact["possible_loss"][0]["interaction_type"], "charged / electrostatic")

    def test_analyze_mutation_scan(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")

        result = analyze_mutation_scan(pdb_path, "MK1", "D25A", chain_id="A")

        self.assertEqual(result["mutation"], "D25A")
        self.assertEqual(result["chain_id"], "A")
        self.assertEqual(result["original_residue"]["res_name"], "ASP")
        self.assertEqual(result["mutant_residue"]["res_name"], "ALA")
        self.assertIn("charge", result["property_changes"])
        self.assertIn("possible_loss", result["interaction_impact"])
        self.assertIn("Mutation Impact Summary", result["ai_interpretation"])

    def test_mutation_report(self):
        pdb_path = os.path.join(ROOT_DIR, "data", "1HSG.pdb")
        result = analyze_mutation_scan(pdb_path, "MK1", "D25A", chain_id="A")

        report = build_mutation_scan_report(result)

        self.assertIn("Mutation Impact Summary", report)
        self.assertIn("D25A", report)
        self.assertIn("Physicochemical Property Changes", report)

    def test_cli_mutation_mode_outputs_json(self):
        completed = subprocess.run(
            [
                sys.executable,
                os.path.join(ROOT_DIR, "scripts", "run_pipeline.py"),
                os.path.join(ROOT_DIR, "data", "1HSG.pdb"),
                "MK1",
                "--mutation",
                "D25A",
                "--chain-id",
                "A",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        result = json.loads(completed.stdout)

        self.assertEqual(result["mutation"], "D25A")
        self.assertEqual(result["chain_id"], "A")
        self.assertEqual(result["original_residue"]["res_name"], "ASP")


if __name__ == "__main__":
    unittest.main()

