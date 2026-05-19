import argparse
import json
import os
import sys

ROOT_DIR = os.path.dirname(os.path.dirname(__file__))

if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from analysis_core import analyze_ligand_pocket, format_ligand_suggestions, list_ligands, sort_residue_keys
from reports import count_interaction_types, generate_pymol_script
from services.mutation_scan import analyze_mutation_scan


def main():
    parser = argparse.ArgumentParser(
        description="Analyze ligand-binding pocket contacts in a PDB file."
    )
    parser.add_argument("pdb_file", help="Path to the PDB file.")
    parser.add_argument("ligand_name", help="Three-letter ligand residue name, for example MK1 or CLR.")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Contact cutoff distance in Angstroms.")
    parser.add_argument("--mutation", help="Optional mutation scan, for example R273H.")
    parser.add_argument("--chain-id", help="Optional chain ID for mutation scan.")

    args = parser.parse_args()
    pdb_path = os.path.abspath(args.pdb_file)
    ligand_name = args.ligand_name.strip().upper()

    if args.mutation:
        mutation_result = analyze_mutation_scan(
            pdb_path,
            ligand_name,
            args.mutation,
            chain_id=args.chain_id
        )
        print(json.dumps(mutation_result, indent=2))
        return

    contacts, counts, interpretation, interactions = analyze_ligand_pocket(
        pdb_path,
        ligand_name,
        cutoff=args.cutoff
    )

    os.makedirs("results", exist_ok=True)

    if contacts is None:
        print("========== Structure Analysis Report ==========\n")
        print(f"Structure file: {pdb_path}")
        print(f"Target ligand: {ligand_name}")
        print(f"\nNo ligand named {ligand_name} found in this PDB file.")
        print(format_ligand_suggestions(list_ligands(pdb_path)))
        sys.exit(1)

    contact_data = [
        {
            "chain_id": chain_id,
            "res_name": res_name,
            "res_id": res_id,
            "atom_name": "",
            "element": ""
        }
        for chain_id, res_name, res_id in sort_residue_keys(contacts.keys())
    ]

    with open("results/latest_contacts.json", "w", encoding="utf-8") as f:
        json.dump(contact_data, f, indent=2)

    pymol_filename = generate_pymol_script(
        pdb_path,
        ligand_name,
        contacts,
        "results",
        output_prefix=os.path.splitext(os.path.basename(pdb_path))[0]
    )

    print("========== Structure Analysis Report ==========\n")
    print(f"Structure file: {pdb_path}")
    print(f"Target ligand: {ligand_name}")

    print("\n[1] Ligand Pocket Summary")
    print(f"Contact residues within {args.cutoff:g}A: {len(contacts)}")
    print(f"Hydrophobic residues: {counts['hydrophobic']}")
    print(f"Polar residues: {counts['polar']}")
    print(f"Positive residues: {counts['positive']}")
    print(f"Negative residues: {counts['negative']}")

    print("\n[2] Interaction Type Summary")
    for interaction_type, count in sorted(count_interaction_types(interactions).items()):
        print(f"- {interaction_type}: {count}")

    print("\n[3] Contact Residues")
    for chain_id, res_name, res_id in sort_residue_keys(contacts.keys()):
        print(f"- Chain {chain_id}: {res_name}{res_id}")

    print("\n[4] Closest Interaction Lines")
    for interaction in interactions:
        print(
            f"- Chain {interaction['chain_id']}: ligand atom {interaction.get('ligand_atom', '?')} -> "
            f"{interaction['res_name']}{interaction['res_id']} atom {interaction['atom_name']} | "
            f"{interaction.get('interaction_type', 'contact')} | {interaction['distance']} A"
        )

    print("\n[5] Interpretation")
    print(interpretation)

    print("\n[6] Visualization")
    print(f"PyMOL script: results/{pymol_filename}")


if __name__ == "__main__":
    main()
