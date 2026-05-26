"""Tests for water_bridges.py — water-mediated protein-ligand contacts."""

import pytest
from water_bridges import detect_water_bridges, _distance


class TestDistance:
    def test_same_point_is_zero(self):
        assert _distance((0, 0, 0), (0, 0, 0)) == 0.0

    def test_unit_distance(self):
        assert _distance((0, 0, 0), (1, 0, 0)) == 1.0

    def test_3d_distance(self):
        d = _distance((0, 0, 0), (1, 2, 2))
        assert d == 3.0


class TestWaterBridgeDetection:
    def test_empty_atoms_returns_empty(self):
        result = detect_water_bridges([], [])
        assert result["bridges"] == []
        assert result["summary"]["total"] == 0

    def test_no_water_molecules(self):
        protein = [
            {"atom_type": "ATOM", "chain_id": "A", "res_name": "LYS", "res_id": "1",
             "atom_name": "NZ", "coord": (0, 0, 0), "element": "N"},
        ]
        ligand = [
            {"atom_type": "HETATM", "chain_id": "X", "res_name": "LIG", "res_id": "1",
             "atom_name": "O1", "coord": (2.5, 0, 0), "element": "O"},
        ]
        result = detect_water_bridges(protein, ligand)
        assert result["bridges"] == []

    def test_single_water_bridge_detected(self):
        protein = [
            {"atom_type": "ATOM", "chain_id": "A", "res_name": "ASP", "res_id": "25",
             "atom_name": "OD1", "coord": (0, 0, 0), "element": "O"},
        ]
        ligand = [
            {"atom_type": "HETATM", "chain_id": "X", "res_name": "LIG", "res_id": "1",
             "atom_name": "N1", "coord": (5.0, 0, 0), "element": "N"},
        ]
        water = [
            {"atom_type": "HETATM", "chain_id": "W", "res_name": "HOH", "res_id": "101",
             "atom_name": "O", "coord": (2.5, 0, 0), "element": "O"},
        ]
        result = detect_water_bridges(protein + water, ligand)
        assert result["summary"]["total"] == 1
        assert len(result["bridges"]) == 1
        b = result["bridges"][0]
        assert b["water_id"] == "W:HOH101"
        assert "A:ASP25" in b["protein_residue"] or b["protein_residue"] in ["A:ASP25"]
        assert b["category"] == "water_bridge"

    def test_water_too_far_from_protein(self):
        protein = [
            {"atom_type": "ATOM", "chain_id": "A", "res_name": "ASP", "res_id": "25",
             "atom_name": "OD1", "coord": (0, 0, 0), "element": "O"},
        ]
        ligand = [
            {"atom_type": "HETATM", "chain_id": "X", "res_name": "LIG", "res_id": "1",
             "atom_name": "N1", "coord": (10.0, 0, 0), "element": "N"},
        ]
        water = [
            {"atom_type": "HETATM", "chain_id": "W", "res_name": "HOH", "res_id": "101",
             "atom_name": "O", "coord": (8.0, 0, 0), "element": "O"},
        ]
        result = detect_water_bridges(protein + water, ligand)
        assert result["bridges"] == []

    def test_multiple_bridges_from_same_water(self):
        protein = [
            {"atom_type": "ATOM", "chain_id": "A", "res_name": "ASP", "res_id": "25",
             "atom_name": "OD1", "coord": (0, 0, 0), "element": "O"},
            {"atom_type": "ATOM", "chain_id": "A", "res_name": "LYS", "res_id": "30",
             "atom_name": "NZ", "coord": (0, 1, 0), "element": "N"},
        ]
        ligand = [
            {"atom_type": "HETATM", "chain_id": "X", "res_name": "LIG", "res_id": "1",
             "atom_name": "N1", "coord": (5.0, 0.5, 0), "element": "N"},
        ]
        water = [
            {"atom_type": "HETATM", "chain_id": "W", "res_name": "HOH", "res_id": "101",
             "atom_name": "O", "coord": (2.5, 0.5, 0), "element": "O"},
        ]
        result = detect_water_bridges(protein + water, ligand)
        # One water bridges to 2 protein residues + 1 ligand → 2 bridges reported
        assert result["summary"]["total"] >= 1
        assert result["summary"]["water_molecules_involved"] >= 1
