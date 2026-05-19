# Protein Structure Copilot

Protein Structure Copilot is a small Flask web app for inspecting protein-ligand
binding pockets from PDB files. It supports single-structure pocket analysis,
WT/mutant pocket comparison, browser-based 3D visualization with 3Dmol.js, and
optional DeepSeek-powered structural interpretation.

## Features

- Upload a PDB file and analyze residues within 5 A of a target ligand.
- Classify contact residues as hydrophobic, polar, positive, negative, or other.
- Label closest ligand-residue interaction lines as hydrophobic, polar / possible
  H-bond, charged / electrostatic, or van der Waals contacts.
- Show color-coded interaction lines and hotspot residues in the 3D viewer.
- Run heuristic mutation scans such as `R273H`, `K120A`, or `Y220C` without
  sidechain remodeling.
- Compare WT and mutant structures to identify lost and gained pocket contacts.
- Suggest detected ligand names when the requested ligand is missing.
- Generate text reports and PyMOL scripts in `results/`.
- Fall back to local rule-based interpretation when DeepSeek is unavailable.

## Project Layout

- `app.py`: Flask routes, upload handling, response rendering, downloads.
- `analysis_core.py`: PDB parsing, distance calculation, residue classification,
  pocket analysis, hotspot helpers.
- `ai_client.py`: DeepSeek API call and local fallback interpretation.
- `reports.py`: Text report, comparison report, and PyMOL script generation.
- `services/mutation_scan.py`: Rule-based mutation parsing, residue lookup,
  property comparison, and interaction impact estimation.
- `templates/index.html`: Single-page UI and 3Dmol.js viewer logic.
- `scripts/run_pipeline.py`: Command-line pocket analysis using the same core
  logic as the web app.
- `data/`: Example PDB files.
- `uploads/`: Uploaded PDB files created by the web app.
- `results/`: Generated reports and PyMOL scripts.

## Setup

Create and activate a Python environment, then install dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

To enable DeepSeek interpretation, create a `.env` file:

```bash
DEEPSEEK_API_KEY=your_api_key_here
```

The app still works without a valid key; it will show a local rule-based
interpretation instead.

## Run

Local Flask development server:

```bash
python app.py
```

Open:

```text
http://127.0.0.1:5000
```

Debug mode is disabled by default. For local development, enable it explicitly:

```bash
FLASK_DEBUG=1 python app.py
```

Production-style local server:

```bash
gunicorn app:app
```

`gunicorn.conf.py` binds to `0.0.0.0:$PORT`, defaulting to port `8000` when
`PORT` is not set.

The front-end viewer loads 3Dmol.js from a CDN, so the browser needs network
access for 3D rendering.

## Mutation Scan

Use the Mutation Scan form to upload a PDB file, enter a ligand name, and provide
a mutation such as:

```text
R273H
K120A
Y220C
```

The first version is heuristic: it does not perform sidechain remodeling,
Rosetta/FoldX scoring, or mutant structure generation. It compares amino-acid
properties and the current ligand interaction set to estimate possible changes
in charge, polarity, hydrophobicity, aromaticity, sidechain size, interaction
gain/loss, and binding stability.

## CLI

Pocket analysis:

```bash
python scripts/run_pipeline.py data/1HSG.pdb MK1
```

Mutation scan JSON output:

```bash
python scripts/run_pipeline.py data/1HSG.pdb MK1 --mutation D25A --chain-id A
```

## Test

```bash
python -m pytest
```

## Deploy To Render

1. Push this repository to GitHub.
2. In Render, create a new Web Service from the GitHub repository.
3. Use these settings:

```text
Build Command: pip install -r requirements.txt
Start Command: gunicorn app:app
```

4. Add environment variables as needed:

```text
DEEPSEEK_API_KEY=your_api_key_here
```

`uploads/` and `results/` are created automatically at app startup. They are
ignored by Git and should be treated as ephemeral runtime storage on Render.

## Suggested Next Development Steps

- Add richer structure chemistry detection, such as hydrogen bonds, salt bridges,
  aromatic contacts, and ligand atom typing.
- Add result history instead of writing only downloadable text artifacts.
- Move long-running analysis and AI calls into background jobs for larger PDBs.
- Add front-end states for upload progress, invalid ligand names, and empty
  interaction sets.
