import unittest
from unittest.mock import patch

from reports import build_report, count_interaction_types


class ReportsTest(unittest.TestCase):
    def test_count_interaction_types(self):
        interactions = [
            {"interaction_type": "hydrophobic contact"},
            {"interaction_type": "hydrophobic contact"},
            {"interaction_type": "polar / possible H-bond"},
        ]

        self.assertEqual(
            count_interaction_types(interactions),
            {
                "hydrophobic contact": 2,
                "polar / possible H-bond": 1
            }
        )

    @patch("reports.generate_ai_interpretation")
    def test_build_report_includes_interaction_types(self, mock_ai):
        mock_ai.return_value = "Mock AI interpretation"

        report_text, ai_text = build_report(
            "LIG",
            {("A", "VAL", "82"): {}},
            {
                "hydrophobic": 1,
                "polar": 0,
                "positive": 0,
                "negative": 0,
                "other": 0
            },
            "Primary interpretation",
            "highlight_test.pml",
            [
                {
                    "chain_id": "A",
                    "res_name": "VAL",
                    "res_id": "82",
                    "atom_name": "CB",
                    "element": "C",
                    "ligand_atom": "C1",
                    "distance": 3.8,
                    "interaction_type": "hydrophobic contact"
                }
            ]
        )

        self.assertEqual(ai_text, "Mock AI interpretation")
        self.assertIn("[2] Interaction Type Summary", report_text)
        self.assertIn("- hydrophobic contact: 1", report_text)
        self.assertIn("ligand atom C1 -> VAL82 atom CB", report_text)


if __name__ == "__main__":
    unittest.main()
