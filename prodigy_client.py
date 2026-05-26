"""
Prodigy binding affinity prediction client.

Prodigy (PROtein binDIng enerGY prediction) is a free web service from
Utrecht University that predicts protein-ligand binding affinity (ΔG, Kd)
from a PDB structure.

API: POST PDB file → returns ΔG (kcal/mol) and Kd (M) at 25°C.
"""

import requests

PRODIGY_URL = "https://wenmr.science.uu.nl/prodigy/run_single"
REQUEST_TIMEOUT = 20


def predict_binding_affinity(pdb_path, ligand_name):
    """Predict binding affinity using Prodigy API.

    Args:
        pdb_path: path to PDB file
        ligand_name: 3-letter residue name of ligand

    Returns:
        dict with {delta_g, kd, kd_units, temperature} or None on failure.
    """
    try:
        with open(pdb_path, "rb") as f:
            files = {"pdb_file": (f"structure.pdb", f, "application/octet-stream")}
            data = {
                "ligand_residue": ligand_name.upper(),
                "temperature": "298",
            }
            resp = requests.post(
                PRODIGY_URL,
                files=files,
                data=data,
                timeout=REQUEST_TIMEOUT,
            )
    except requests.RequestException:
        return None

    if resp.status_code != 200:
        return None

    try:
        result = resp.json()
    except ValueError:
        return None

    if not result:
        return None

    delta_g = result.get("delta_g") or result.get("binding_affinity") or result.get("dg")
    kd = result.get("kd") or result.get("dissociation_constant")

    if delta_g is None and kd is None:
        return None

    return {
        "delta_g": round(float(delta_g), 2) if delta_g else None,
        "kd": _format_kd(kd) if kd else None,
        "temperature": 298,
        "source": "Prodigy (Utrecht University)",
    }


def _format_kd(kd_value):
    """Format Kd value with appropriate units."""
    try:
        kd = float(kd_value)
    except (ValueError, TypeError):
        return str(kd_value)

    if kd < 1e-9:
        return f"{kd * 1e12:.1f} pM"
    elif kd < 1e-6:
        return f"{kd * 1e9:.1f} nM"
    elif kd < 1e-3:
        return f"{kd * 1e6:.1f} uM"
    elif kd < 1:
        return f"{kd * 1e3:.1f} mM"
    else:
        return f"{kd:.2f} M"
