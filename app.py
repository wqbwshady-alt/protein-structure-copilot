import os
from uuid import uuid4

from dotenv import load_dotenv
from flask import Flask, render_template, request, send_from_directory, url_for
from markupsafe import escape
from werkzeug.utils import secure_filename

from analysis_core import (
    analyze_ligand_pocket,
    format_ligand_suggestions,
    get_hotspot_residues,
    is_pdb_file,
    list_ligands,
    residue_keys_to_json,
)
from reports import build_comparison_report, build_report, generate_pymol_script
from reports import build_mutation_scan_report
from services.mutation_scan import MutationScanError, analyze_mutation_scan


load_dotenv()

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 32 * 1024 * 1024

UPLOAD_FOLDER = "uploads"
RESULT_FOLDER = "results"
ALLOWED_EXTENSIONS = {".pdb"}

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULT_FOLDER, exist_ok=True)


def empty_page_context(**overrides):
    context = {
        "result_text": None,
        "ai_html": "",
        "pdb_url": None,
        "report_download_url": None,
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
    return str(escape(text)).replace("\n", "<br>")


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


def write_result_file(filename, text):
    path = os.path.join(RESULT_FOLDER, filename)

    with open(path, "w", encoding="utf-8") as f:
        f.write(text)

    return filename


@app.route("/", methods=["GET"])
def index():
    return render_index()


@app.route("/analyze", methods=["POST"])
def analyze():
    pdb_file = request.files.get("pdb_file")
    ligand_name = request.form.get("ligand_name", "").strip().upper()

    if not pdb_file or pdb_file.filename == "":
        error_text = "Please upload a PDB file."
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    if not ligand_name:
        error_text = "Please enter ligand name, for example: MK1 / CLR."
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    filename, pdb_path, error_text = save_uploaded_pdb(pdb_file)

    if error_text:
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    contact_residues, counts, primary_interpretation, interactions = analyze_ligand_pocket(
        pdb_path,
        ligand_name
    )

    if contact_residues is None:
        result_text = (
            f"No ligand named {ligand_name} found in this PDB file.\n"
            f"{format_ligand_suggestions(list_ligands(pdb_path))}"
        )
        return render_index(
            result_text=result_text,
            ai_html=make_html(result_text),
            pdb_url=url_for("uploaded_file", filename=filename)
        )

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

    return render_index(
        result_text=report_text,
        ai_html=make_html(ai_text),
        pdb_url=url_for("uploaded_file", filename=filename),
        report_download_url=url_for("download_report", filename=report_filename),
        interaction_data=interactions,
        hotspot_residues=get_hotspot_residues(interactions)
    )


@app.route("/compare", methods=["POST"])
def compare():
    wt_file = request.files.get("wt_file")
    mut_file = request.files.get("mut_file")
    ligand_name = request.form.get("compare_ligand_name", "").strip().upper()

    if not wt_file or wt_file.filename == "":
        error_text = "Please upload WT PDB file."
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    if not mut_file or mut_file.filename == "":
        error_text = "Please upload Mutant PDB file."
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    if not ligand_name:
        error_text = "Please enter ligand name for comparison, for example: MK1 / CLR."
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    wt_filename, wt_path, error_text = save_uploaded_pdb(wt_file, prefix="WT_")

    if error_text:
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    mut_filename, mut_path, error_text = save_uploaded_pdb(mut_file, prefix="MUT_")

    if error_text:
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
        return render_index(result_text=result_text, ai_html=make_html(result_text))

    if mut_contacts is None:
        result_text = (
            f"No ligand named {ligand_name} found in Mutant PDB file.\n"
            f"{format_ligand_suggestions(list_ligands(mut_path))}"
        )
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

    return render_index(
        result_text=wt_report_text,
        ai_html=make_html(wt_ai_text),
        pdb_url=url_for("uploaded_file", filename=wt_filename),
        interaction_data=wt_interactions,
        comparison_text=comparison_text,
        lost_residues=residue_keys_to_json(lost),
        gained_residues=residue_keys_to_json(gained),
        hotspot_residues=get_hotspot_residues(wt_interactions)
    )


@app.route("/mutation_scan", methods=["POST"])
def mutation_scan():
    pdb_file = request.files.get("mutation_pdb_file")
    ligand_name = request.form.get("mutation_ligand_name", "").strip().upper()
    mutation_text = request.form.get("mutation_text", "").strip().upper()
    chain_id = request.form.get("mutation_chain_id", "").strip()

    if not pdb_file or pdb_file.filename == "":
        error_text = "Please upload a PDB file for mutation scan."
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    if not ligand_name:
        error_text = "Please enter ligand name for mutation scan."
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    if not mutation_text:
        error_text = "Please enter mutation, for example: R273H."
        return render_index(result_text=error_text, ai_html=make_html(error_text))

    filename, pdb_path, error_text = save_uploaded_pdb(pdb_file, prefix="MUTSCAN_")

    if error_text:
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

    return render_index(
        result_text=mutation_scan_text,
        ai_html=make_html(mutation_result["ai_interpretation"]),
        pdb_url=url_for("uploaded_file", filename=filename),
        interaction_data=interactions,
        hotspot_residues=get_hotspot_residues(interactions),
        mutation_scan_result=mutation_result,
        mutation_scan_text=mutation_scan_text
    )


@app.route("/uploads/<path:filename>")
def uploaded_file(filename):
    return send_from_directory(UPLOAD_FOLDER, filename)


@app.route("/download_report/<path:filename>")
def download_report(filename):
    return send_from_directory(RESULT_FOLDER, filename, as_attachment=True)


if __name__ == "__main__":
    app.run(debug=os.getenv("FLASK_DEBUG") == "1")
