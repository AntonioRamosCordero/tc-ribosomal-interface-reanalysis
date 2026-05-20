from Bio.PDB import PDBParser, SASA
import sys
from pathlib import Path

if len(sys.argv) != 3:
    print("Usage: python calculate_sasa_by_residue.py input.pdb output.tsv")
    sys.exit(1)

input_pdb = sys.argv[1]
output_file = sys.argv[2]

parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", input_pdb)

sr = SASA.ShrakeRupley()
sr.compute(structure, level="R")

with open(output_file, "w") as f:
    f.write("Residue\tSASA\n")

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":
                    resname = residue.resname
                    resid = residue.id[1]
                    sasa = getattr(residue, "sasa", 0.0)

                    f.write(f"{chain.id}:{resname}{resid}\t{sasa:.2f}\n")

print(f"SASA table saved to: {output_file}")
