# CLAUDE.md

This is the Protein Structure Copilot project — a Flask web app for protein-ligand binding pocket analysis with 3D visualization.

## Guidelines

- Do not rewrite the entire project without explicit approval.
- Before making changes, analyze the problem first and explain the plan.
- Preserve existing functionality: ligand detection, pocket analysis, 3D viewer (3Dmol.js), PyMOL script generation, mutation scan, WT/mutant comparison.
- Do not remove existing API routes (`/analyze`, `/compare`, `/mutation_scan`, `/health`, `/uploads/<path>`, `/download_report/<path>`).
- Prefer small, incremental changes over large refactors.
- After each change, list the files that were modified.
- Keep the Flask backend and frontend structure stable. Avoid introducing new frameworks or build tools unless the user explicitly requests them.
