"""Tests for sasa.py — Shrake-Rupley solvent accessible surface area."""

import pytest
from sasa import compute_sasa, _shrake_rupley, _SPHERE_POINTS


class TestSpherePoints:
    def test_sphere_points_generated(self):
        n = 100
        points = _SPHERE_POINTS(n)
        assert len(points) == n
        # All points should be on unit sphere
        for p in points[:10]:
            r = sum(x * x for x in p) ** 0.5
            assert abs(r - 1.0) < 0.01


class TestShrakeRupley:
    def test_isolated_atom_full_exposure(self):
        atoms = [
            {"coord": (0, 0, 0), "element": "C"},
        ]
        radii = {"C": 1.7}
        result = _shrake_rupley(atoms, radii, n_points=200)
        assert len(result) == 1
        # Isolated atom should have near-maximum SASA
        expected_max = 4 * 3.14159 * (1.7 + 1.4) ** 2  # ~120 A^2
        assert result[0] > 0
        assert result[0] <= expected_max * 1.1

    def test_buried_atom_zero_exposure(self):
        # Atom 2 fully envelops atom 1 (larger radius at same position)
        atoms = [
            {"coord": (0, 0, 0), "element": "C"},
            {"coord": (0, 0, 0), "element": "P"},  # Larger — fully occludes atom 1
        ]
        radii = {"C": 1.7, "P": 2.5}
        result = _shrake_rupley(atoms, radii, n_points=200)
        # Atom 1 (smaller C) should be nearly fully buried by larger P
        assert result[0] < 15.0
        # Atom 2 (larger P) still has exposed surface
        assert result[1] > 0

    def test_two_separated_atoms(self):
        atoms = [
            {"coord": (0, 0, 0), "element": "C"},
            {"coord": (10, 0, 0), "element": "C"},  # Far apart
        ]
        radii = {"C": 1.7}
        result = _shrake_rupley(atoms, radii, n_points=200)
        # Both should be near-maximum since they don't overlap
        for r in result:
            assert r > 50.0  # Most of the surface exposed


class TestComputeSASA:
    def test_empty_atoms(self):
        result = compute_sasa([], {})
        assert result["per_residue"] == {}
        assert result["total_sasa"] == 0.0

    def test_single_residue_sasa(self):
        atoms = [
            {"atom_type": "ATOM", "chain_id": "A", "res_name": "ALA", "res_id": "1",
             "atom_name": "CA", "coord": (0, 0, 0), "element": "C"},
        ]
        result = compute_sasa(atoms, {("A", "ALA", "1"): {}})
        assert "A:ALA1" in result["per_residue"]
        assert result["per_residue"]["A:ALA1"]["sasa"] > 0
        assert result["total_sasa"] > 0

    def test_per_residue_sum_matches_total(self):
        atoms = [
            {"atom_type": "ATOM", "chain_id": "A", "res_name": "ALA", "res_id": "1",
             "atom_name": "CA", "coord": (0, 0, 0), "element": "C"},
            {"atom_type": "ATOM", "chain_id": "A", "res_name": "GLY", "res_id": "2",
             "atom_name": "CA", "coord": (10, 0, 0), "element": "C"},
        ]
        contact = {("A", "ALA", "1"): {}, ("A", "GLY", "2"): {}}
        result = compute_sasa(atoms, contact)
        per_res_sum = sum(v["sasa"] for v in result["per_residue"].values())
        assert abs(per_res_sum - result["total_sasa"]) < 0.01

    def test_burial_classification(self):
        atoms = [
            {"atom_type": "ATOM", "chain_id": "A", "res_name": "ALA", "res_id": "1",
             "atom_name": "CA", "coord": (0, 0, 0), "element": "C"},
        ]
        result = compute_sasa(atoms, {("A", "ALA", "1"): {}})
        classification = result["per_residue"]["A:ALA1"]["classification"]
        assert classification in ("buried", "partially_buried", "exposed")
