import csv
import io
import json
import os
import re
import ssl
import urllib.error
import urllib.request
from uuid import uuid4

from dotenv import load_dotenv
from flask import Flask, jsonify, render_template, request, send_from_directory, url_for
from markupsafe import escape
from werkzeug.utils import secure_filename

from analysis_core import (
    analyze_ligand_pocket,
    classify_residue,
    format_ligand_suggestions,
    get_hotspot_residues,
    is_pdb_file,
    list_ligands,
    parse_pdb_atoms,
    residue_keys_to_json,
)
from reports import build_comparison_report, build_report, generate_pymol_script
from reports import build_mutation_scan_report
from services.mutation_scan import MutationScanError, analyze_mutation_scan
from services.stats import get_stats, increment_analysis, increment_comparison, increment_mutation_scan


load_dotenv()

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 32 * 1024 * 1024

UPLOAD_FOLDER = "uploads"
RESULT_FOLDER = "results"
ALLOWED_EXTENSIONS = {".pdb"}

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULT_FOLDER, exist_ok=True)


def make_urlopen_context():
    try:
        import certifi
    except ImportError:
        return ssl.create_default_context()

    return ssl.create_default_context(cafile=certifi.where())


def empty_page_context(**overrides):
    context = {
        "result_text": None,
        "ai_html": "",
        "pdb_url": None,
        "report_download_url": None,
        "json_download_url": None,
        "csv_download_url": None,
        "analysis_mode": None,
        "protein_summary": None,
        "interaction_data": [],
        "comparison_text": None,
        "mutation_scan_result": None,
        "mutation_scan_text": None,
        "lost_residues": [],
        "gained_residues": [],
        "hotspot_residues": []
    }
    context.update(overrides)
    return context


def render_index(**overrides):
    return render_template("index.html", **empty_page_context(**overrides))


def make_html(text):
    if not text:
        return ""
    escaped = str(escape(text))
    lines = escaped.split("\n")
    result = []
    for line in lines:
        if line.startswith("## "):
            result.append(
                '<div style="font-weight:700;font-size:15px;color:#60a5fa;'
                'margin-top:16px;margin-bottom:6px;">' + line[3:] + "</div>"
            )
        else:
            result.append(line)
    return "<br>".join(result)


def save_uploaded_pdb(file_storage, prefix=""):
    original_name = secure_filename(file_storage.filename or "")

    if not original_name:
        return None, None, "Please upload a PDB file."

    extension = os.path.splitext(original_name)[1].lower()

    if extension not in ALLOWED_EXTENSIONS:
        return None, None, "Only .pdb files are supported."

    filename = f"{prefix}{uuid4().hex}_{original_name}"
    path = os.path.join(UPLOAD_FOLDER, filename)
    file_storage.save(path)

    if not is_pdb_file(path):
        os.remove(path)
        return None, None, "Uploaded file does not contain ATOM or HETATM records."

    return filename, path, None


def resolve_loaded_structure(file_fields, prefix=""):
    pdb_filename = request.form.get("pdb_filename", "").strip()

    if pdb_filename:
        filename = secure_filename(pdb_filename)
        pdb_path = os.path.join(UPLOAD_FOLDER, filename)
        if not os.path.isfile(pdb_path):
            return None, None, "Fetched PDB file no longer available. Please re-fetch."
        if not is_pdb_file(pdb_path):
            return None, None, "Loaded PDB file does not contain ATOM or HETATM records."
        return filename, pdb_path, None

    for field in file_fields:
        pdb_file = request.files.get(field)
        if pdb_file and pdb_file.filename:
            return save_uploaded_pdb(pdb_file, prefix=prefix)

    return None, None, "Please load a structure from RCSB or upload a local PDB file first."


def resolve_structure_input(filename_field, file_field, prefix, missing_message):
    pdb_filename = request.form.get(filename_field, "").strip()

    if pdb_filename:
        filename = secure_filename(pdb_filename)
        pdb_path = os.path.join(UPLOAD_FOLDER, filename)
        if not os.path.isfile(pdb_path):
            return None, None, "Fetched PDB file no longer available. Please re-fetch."
        if not is_pdb_file(pdb_path):
            return None, None, "Loaded PDB file does not contain ATOM or HETATM records."
        return filename, pdb_path, None

    pdb_file = request.files.get(file_field)
    if pdb_file and pdb_file.filename:
        return save_uploaded_pdb(pdb_file, prefix=prefix)

    return None, None, missing_message


def write_result_file(filename, text):
    path = os.path.join(RESULT_FOLDER, filename)

    with open(path, "w", encoding="utf-8") as f:
        f.write(text)

    return filename


def mutation_residue_highlights(mutation_result):
    lost = []
    gained = []

    if mutation_result["interaction_impact"]["possible_loss"]:
        lost.append(mutation_result["original_residue"])

    if mutation_result["interaction_impact"]["possible_gain"]:
        gained.append(mutation_result["mutant_residue"])

    return lost, gained


@app.route("/", methods=["GET"])
def index():
    return render_index()


@app.route("/health", methods=["GET"])
def health():
    return jsonify({"status": "ok"})


@app.route("/stats", methods=["GET"])
def stats():
    return jsonify(get_stats())


@app.route("/analyze", methods=["POST"])
def analyze():
    wants_json = (
        request.headers.get("X-Requested-With") == "XMLHttpRequest" or
        request.accept_mimetypes.best == "application/json"
    )
    pdb_filename = request.form.get("pdb_filename", "").strip()
    pdb_file = request.files.get("pdb_file")
    ligand_name = request.form.get("ligand_name", "").strip().upper()
    skip_ligand = request.form.get("skipLigand", "").strip().lower() in {"1", "true", "yes", "on"}

    if pdb_filename:
        pdb_path = os.path.join(UPLOAD_FOLDER, secure_filename(pdb_filename))
        if not os.path.isfile(pdb_path):
            if wants_json:
                return jsonify({
                    "success": False,
                    "error_text": "Fetched PDB file no longer available. Please re-fetch.",
                    "ai_html": make_html("Fetched PDB file no longer available.")
                })
            return render_index(
                result_text="Fetched PDB file no longer available. Please re-fetch.",
                ai_html=make_html("Fetched PDB file no longer available.")
            )
    else:
        if not pdb_file or pdb_file.filename == "":
            error_text = "Please upload a PDB file."
            if wants_json:
                return jsonify({
                    "success": False,
                    "error_text": error_text,
                    "ai_html": make_html(error_text)
                })
            return render_index(result_text=error_text, ai_html=make_html(error_text))

        if not ligand_name and not skip_ligand:
            error_text = "Please enter ligand name, for example: MK1 / CLR."
            if wants_json:
                return jsonify({
                    "success": False,
                    "error_text": error_text,
                    "ai_html": make_html(error_text)
                })
            return render_index(result_text=error_text, ai_html=make_html(error_text))

        pdb_filename, pdb_path, error_text = save_uploaded_pdb(pdb_file)

        if error_text:
            if wants_json:
                return jsonify({
                    "success": False,
                    "error_text": error_text,
                    "ai_html": make_html(error_text)
                })
            return render_index(result_text=error_text, ai_html=make_html(error_text))

    if skip_ligand:
        result = _build_protein_only_result(pdb_path, pdb_filename)
        if wants_json:
            return jsonify(result)
        return render_index(**result)

    if not ligand_name:
        error_text = "Please enter ligand name, for example: MK1 / CLR."
        if wants_json:
            return jsonify({
                "success": False,
                "error_text": error_text,
                "ai_html": make_html(error_text)
            })
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    result = _build_analyze_result(pdb_path, pdb_filename, ligand_name)
    if wants_json:
        return jsonify(result)
    if not result["success"]:
        return render_index(
            result_text=result["error_text"],
            ai_html=result["ai_html"],
            pdb_url=result.get("pdb_url")
        )
    return render_index(**result)


@app.route("/analyze_async", methods=["POST"])
def analyze_async():
    pdb_filename = request.form.get("pdb_filename", "").strip()
    ligand_name = request.form.get("ligand_name", "").strip().upper()
    skip_ligand = request.form.get("skipLigand", "").strip().lower() in {"1", "true", "yes", "on"}

    if pdb_filename:
        pdb_path = os.path.join(UPLOAD_FOLDER, secure_filename(pdb_filename))
        if not os.path.isfile(pdb_path):
            return jsonify({
                "success": False,
                "error_text": "Fetched PDB file no longer available. Please re-fetch.",
                "ai_html": make_html("Fetched PDB file no longer available.")
            })
        if skip_ligand:
            return jsonify(_build_protein_only_result(pdb_path, pdb_filename))
        result = _build_analyze_result(pdb_path, pdb_filename, ligand_name)
        return jsonify(result)

    pdb_file = request.files.get("pdb_file")
    if not pdb_file or pdb_file.filename == "":
        return jsonify({
            "success": False,
            "error_text": "Please upload a PDB file or fetch from RCSB.",
            "ai_html": make_html("Please upload a PDB file or fetch from RCSB.")
        })

    if not ligand_name and not skip_ligand:
        return jsonify({
            "success": False,
            "error_text": "Please enter ligand name, for example: MK1 / CLR.",
            "ai_html": make_html("Please enter ligand name, for example: MK1 / CLR.")
        })

    filename, pdb_path, error_text = save_uploaded_pdb(pdb_file)

    if error_text:
        return jsonify({
            "success": False,
            "error_text": error_text,
            "ai_html": make_html(error_text)
        })

    if skip_ligand:
        result = _build_protein_only_result(pdb_path, filename)
    else:
        result = _build_analyze_result(pdb_path, filename, ligand_name)
    return jsonify(result)


def _build_protein_only_result(pdb_path, filename):
    atoms = parse_pdb_atoms(pdb_path)
    protein_atoms = [atom for atom in atoms if atom["atom_type"] == "ATOM"]
    residue_map = {}
    chain_map = {}

    for atom in protein_atoms:
        chain_id = atom["chain_id"] or "(blank)"
        residue_key = (chain_id, atom["res_id"], atom["res_name"])
        residue_map[residue_key] = atom["res_name"]

        chain_data = chain_map.setdefault(chain_id, {
            "chain_id": chain_id,
            "atom_count": 0,
            "residue_keys": set(),
            "hydrophobic": 0,
            "polar": 0,
            "positive": 0,
            "negative": 0,
            "other": 0
        })
        chain_data["atom_count"] += 1
        chain_data["residue_keys"].add(residue_key)

    residue_type_counts = {
        "hydrophobic": 0,
        "polar": 0,
        "positive": 0,
        "negative": 0,
        "other": 0
    }

    for residue_key, res_name in residue_map.items():
        residue_type = classify_residue(res_name)
        if residue_type not in residue_type_counts:
            residue_type = "other"
        residue_type_counts[residue_type] += 1

        chain_id = residue_key[0]
        chain_map[chain_id][residue_type] += 1

    chains = []
    for chain_id in sorted(chain_map):
        chain_data = chain_map[chain_id]
        chains.append({
            "chain_id": chain_id,
            "atom_count": chain_data["atom_count"],
            "residue_count": len(chain_data["residue_keys"]),
            "hydrophobic": chain_data["hydrophobic"],
            "polar": chain_data["polar"],
            "positive": chain_data["positive"],
            "negative": chain_data["negative"],
            "other": chain_data["other"]
        })

    ligand_candidates = sorted(set(ligand["res_name"] for ligand in list_ligands(pdb_path)))
    summary = {
        "atom_count": len(protein_atoms),
        "residue_count": len(residue_map),
        "chain_count": len(chains),
        "hetatm_count": len([atom for atom in atoms if atom["atom_type"] == "HETATM"]),
        "chains": chains,
        "residue_type_counts": residue_type_counts,
        "ligand_candidates": ligand_candidates,
        "secondary_structure": {
            "status": "not_calculated",
            "note": "Secondary structure statistics are reserved for a future DSSP-style analysis step."
        }
    }

    result_text = (
        "Protein-only structural overview\n\n"
        f"Structure: {filename}\n"
        f"Chains: {summary['chain_count']}\n"
        f"Protein residues: {summary['residue_count']}\n"
        f"Protein atoms: {summary['atom_count']}\n\n"
        "Residue composition:\n"
        f"- Hydrophobic: {residue_type_counts['hydrophobic']}\n"
        f"- Polar: {residue_type_counts['polar']}\n"
        f"- Positively charged: {residue_type_counts['positive']}\n"
        f"- Negatively charged: {residue_type_counts['negative']}\n"
        f"- Other: {residue_type_counts['other']}\n"
    )

    ai_text = (
        "Protein-only analysis mode. Ligand pocket contact analysis was skipped. "
        "This overview summarizes chain counts, residue counts, broad residue chemistry, "
        "and provides a 3D structure view for protein-only structures."
    )

    report_filename = write_result_file(
        f"protein_only_report_{uuid4().hex}.txt",
        result_text
    )

    json_filename = write_result_file(
        f"protein_only_data_{uuid4().hex}.json",
        json.dumps(summary, indent=2)
    )

    increment_analysis()

    return {
        "success": True,
        "analysis_mode": "protein_only",
        "result_title": "Protein-only structural overview",
        "result_text": result_text,
        "ai_html": make_html(ai_text),
        "pdb_url": url_for("uploaded_file", filename=filename),
        "report_download_url": url_for("download_report", filename=report_filename),
        "json_download_url": url_for("download_report_json", filename=json_filename),
        "protein_summary": summary,
        "interaction_data": [],
        "hotspot_residues": [],
        "lost_residues": [],
        "gained_residues": []
    }


def _build_analyze_result(pdb_path, filename, ligand_name):

    contact_residues, counts, primary_interpretation, interactions = analyze_ligand_pocket(
        pdb_path,
        ligand_name
    )

    if contact_residues is None:
        result_text = (
            f"No ligand named {ligand_name} found in this PDB file.\n"
            f"{format_ligand_suggestions(list_ligands(pdb_path))}"
        )
        return {
            "success": False,
            "error_text": result_text,
            "ai_html": make_html(result_text),
            "pdb_url": url_for("uploaded_file", filename=filename)
        }

    pymol_filename = generate_pymol_script(
        os.path.join(UPLOAD_FOLDER, filename),
        ligand_name,
        contact_residues,
        RESULT_FOLDER,
        output_prefix=os.path.splitext(filename)[0]
    )

    report_text, ai_text = build_report(
        ligand_name,
        contact_residues,
        counts,
        primary_interpretation,
        pymol_filename,
        interactions
    )

    report_filename = write_result_file(
        f"pocket_report_{uuid4().hex}.txt",
        report_text
    )

    json_filename = write_result_file(
        f"pocket_data_{uuid4().hex}.json",
        json.dumps(interactions, indent=2)
    )

    csv_buf = io.StringIO()
    csv_writer = csv.writer(csv_buf)
    csv_writer.writerow(["chain_id", "res_name", "res_id", "interaction_type",
                         "distance", "ligand_atom", "protein_atom", "element"])
    for item in interactions:
        csv_writer.writerow([
            item["chain_id"], item["res_name"], item["res_id"],
            item["interaction_type"], item["distance"],
            item["ligand_atom"], item["atom_name"], item["element"]
        ])
    csv_filename = write_result_file(
        f"pocket_data_{uuid4().hex}.csv",
        csv_buf.getvalue()
    )

    increment_analysis()

    return {
        "success": True,
        "analysis_mode": "ligand",
        "result_text": report_text,
        "ai_html": make_html(ai_text),
        "pdb_url": url_for("uploaded_file", filename=filename),
        "report_download_url": url_for("download_report", filename=report_filename),
        "json_download_url": url_for("download_report_json", filename=json_filename),
        "csv_download_url": url_for("download_report_csv", filename=csv_filename),
        "interaction_data": interactions,
        "hotspot_residues": get_hotspot_residues(interactions),
        "lost_residues": [],
        "gained_residues": []
    }


@app.route("/compare", methods=["POST"])
def compare():
    wants_json = (
        request.headers.get("X-Requested-With") == "XMLHttpRequest" or
        request.accept_mimetypes.best == "application/json"
    )
    ligand_name = request.form.get("compare_ligand_name", "").strip().upper()

    if not ligand_name:
        error_text = "Please enter ligand name for comparison, for example: MK1 / CLR."
        if wants_json:
            return jsonify({"success": False, "error_text": error_text, "ai_html": make_html(error_text)})
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    wt_filename, wt_path, error_text = resolve_structure_input(
        "wt_pdb_filename",
        "wt_file",
        "WT_",
        "Please load WT structure from RCSB or upload a local WT PDB file first."
    )

    if error_text:
        if wants_json:
            return jsonify({"success": False, "error_text": error_text, "ai_html": make_html(error_text)})
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    mut_filename, mut_path, error_text = resolve_structure_input(
        "mut_pdb_filename",
        "mut_file",
        "MUT_",
        "Please load Mutant structure from RCSB or upload a local Mutant PDB file first."
    )

    if error_text:
        if wants_json:
            return jsonify({"success": False, "error_text": error_text, "ai_html": make_html(error_text)})
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    wt_contacts, wt_counts, wt_primary_interpretation, wt_interactions = analyze_ligand_pocket(
        wt_path,
        ligand_name
    )

    mut_contacts, mut_counts, mut_primary_interpretation, mut_interactions = analyze_ligand_pocket(
        mut_path,
        ligand_name
    )

    if wt_contacts is None:
        result_text = (
            f"No ligand named {ligand_name} found in WT PDB file.\n"
            f"{format_ligand_suggestions(list_ligands(wt_path))}"
        )
        if wants_json:
            return jsonify({"success": False, "error_text": result_text, "ai_html": make_html(result_text)})
        return render_index(result_text=result_text, ai_html=make_html(result_text))

    if mut_contacts is None:
        result_text = (
            f"No ligand named {ligand_name} found in Mutant PDB file.\n"
            f"{format_ligand_suggestions(list_ligands(mut_path))}"
        )
        if wants_json:
            return jsonify({
                "success": False,
                "error_text": result_text,
                "ai_html": make_html(result_text),
                "pdb_url": url_for("uploaded_file", filename=wt_filename),
                "interaction_data": wt_interactions,
                "hotspot_residues": get_hotspot_residues(wt_interactions)
            })
        return render_index(
            result_text=result_text,
            ai_html=make_html(result_text),
            pdb_url=url_for("uploaded_file", filename=wt_filename),
            interaction_data=wt_interactions,
            hotspot_residues=get_hotspot_residues(wt_interactions)
        )

    wt_report_text, wt_ai_text = build_report(
        ligand_name,
        wt_contacts,
        wt_counts,
        wt_primary_interpretation,
        "WT_structure",
        wt_interactions
    )

    comparison_text, lost, gained = build_comparison_report(
        ligand_name,
        wt_contacts,
        mut_contacts,
        wt_counts,
        mut_counts
    )

    write_result_file(
        f"mutation_comparison_report_{uuid4().hex}.txt",
        comparison_text
    )

    increment_comparison()

    result = {
        "success": True,
        "analysis_mode": "comparison",
        "result_text": wt_report_text,
        "ai_html": make_html(wt_ai_text),
        "pdb_url": url_for("uploaded_file", filename=wt_filename),
        "interaction_data": wt_interactions,
        "comparison_text": comparison_text,
        "lost_residues": residue_keys_to_json(lost),
        "gained_residues": residue_keys_to_json(gained),
        "hotspot_residues": get_hotspot_residues(wt_interactions)
    }
    if wants_json:
        return jsonify(result)
    return render_index(**result)


@app.route("/mutation_scan", methods=["POST"])
def mutation_scan():
    wants_json = (
        request.headers.get("X-Requested-With") == "XMLHttpRequest" or
        request.accept_mimetypes.best == "application/json"
    )
    ligand_name = request.form.get("mutation_ligand_name", "").strip().upper()
    mutation_text = request.form.get("mutation_text", "").strip().upper()
    chain_id = request.form.get("mutation_chain_id", "").strip()

    filename, pdb_path, error_text = resolve_loaded_structure(
        ("pdb_file", "mutation_pdb_file"),
        prefix="MUTSCAN_"
    )

    if error_text:
        if wants_json:
            return jsonify({"success": False, "error_text": error_text, "ai_html": make_html(error_text)})
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    if not ligand_name:
        error_text = "Please enter ligand name for mutation scan."
        if wants_json:
            return jsonify({"success": False, "error_text": error_text, "ai_html": make_html(error_text)})
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    if not mutation_text:
        error_text = "Please enter mutation, for example: R273H."
        if wants_json:
            return jsonify({"success": False, "error_text": error_text, "ai_html": make_html(error_text)})
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    try:
        mutation_result = analyze_mutation_scan(
            pdb_path,
            ligand_name,
            mutation_text,
            chain_id=chain_id or None
        )
    except MutationScanError as error:
        error_text = str(error)
        if wants_json:
            return jsonify({
                "success": False,
                "error_text": error_text,
                "ai_html": make_html(error_text),
                "pdb_url": url_for("uploaded_file", filename=filename)
            })
        return render_index(
            result_text=error_text,
            ai_html=make_html(error_text),
            pdb_url=url_for("uploaded_file", filename=filename)
        )

    contact_residues, counts, primary_interpretation, interactions = analyze_ligand_pocket(
        pdb_path,
        ligand_name
    )

    mutation_scan_text = build_mutation_scan_report(mutation_result)

    increment_mutation_scan()

    result = {
        "success": True,
        "analysis_mode": "mutation",
        "result_text": mutation_scan_text,
        "ai_html": make_html(mutation_result["ai_interpretation"]),
        "pdb_url": url_for("uploaded_file", filename=filename),
        "interaction_data": interactions,
        "hotspot_residues": get_hotspot_residues(interactions),
        "mutation_scan_result": mutation_result,
        "mutation_scan_text": mutation_scan_text
    }
    if wants_json:
        return jsonify(result)
    return render_index(**result)


def _fetch_pdb_from_rcsb(pdb_id):
    pdb_id = pdb_id.strip().upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "ProteinStructureCopilot/1.0"})
        with urllib.request.urlopen(req, timeout=30, context=make_urlopen_context()) as resp:
            data = resp.read().decode("utf-8", errors="ignore")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None, None, f"PDB ID {pdb_id} not found in RCSB."
        return None, None, f"RCSB server error (HTTP {e.code})."
    except Exception as e:
        return None, None, f"Could not connect to RCSB: {str(e)}"

    if not data.strip():
        return None, None, f"Downloaded file for {pdb_id} is empty."

    filename = f"RCSB_{pdb_id}_{uuid4().hex[:8]}.pdb"
    path = os.path.join(UPLOAD_FOLDER, filename)
    with open(path, "w", encoding="utf-8") as f:
        f.write(data)

    if not is_pdb_file(path):
        os.remove(path)
        return None, None, f"Downloaded file for {pdb_id} does not contain valid PDB data."

    return filename, path, None


@app.route("/fetch_pdb", methods=["POST"])
def fetch_pdb():
    pdb_id = request.form.get("pdb_id", "").strip().upper()
    if not pdb_id:
        return jsonify({"success": False, "error_text": "Please enter a PDB ID."})

    if not re.fullmatch(r"[A-Za-z0-9]{4}", pdb_id):
        return jsonify({"success": False, "error_text": "PDB ID must be exactly 4 characters (e.g. 1HSG)."})

    filename, pdb_path, error_text = _fetch_pdb_from_rcsb(pdb_id)
    if error_text:
        return jsonify({"success": False, "error_text": error_text})

    ligands = list_ligands(pdb_path)
    ligand_names = sorted(set(lig["res_name"] for lig in ligands))

    return jsonify({
        "success": True,
        "pdb_id": pdb_id,
        "filename": filename,
        "pdb_url": url_for("uploaded_file", filename=filename),
        "ligands": ligand_names[:12],
        "ligand_count": len(ligand_names)
    })


@app.route("/uploads/<path:filename>")
def uploaded_file(filename):
    return send_from_directory(UPLOAD_FOLDER, filename)


@app.route("/download_report/<path:filename>")
def download_report(filename):
    return send_from_directory(RESULT_FOLDER, filename, as_attachment=True)


@app.route("/download_report_json/<path:filename>")
def download_report_json(filename):
    return send_from_directory(RESULT_FOLDER, filename, as_attachment=True)


@app.route("/download_report_csv/<path:filename>")
def download_report_csv(filename):
    return send_from_directory(RESULT_FOLDER, filename, as_attachment=True)


if __name__ == "__main__":
    app.run(debug=os.getenv("FLASK_DEBUG") == "1")
