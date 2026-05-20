from Bio.PDB import PDBParser
import numpy as np
import pandas as pd
from pathlib import Path

pdb = "01_PDB/5OPT.pdb"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("5OPT", pdb)

hotspot_chain = "P"
hotspot_resids = [149,152,153,161,162,163]

hotspot_atoms = []
other_atoms = []

for model in structure:
    for chain in model:
        for res in chain:
            for atom in res:
                item = {
                    "chain": chain.id,
                    "resid": res.id[1],
                    "resname": res.resname,
                    "atom": atom.name,
                    "coord": atom.coord,
                    "hetflag": res.id[0]
                }

                if chain.id == hotspot_chain and res.id[1] in hotspot_resids:
                    hotspot_atoms.append(item)
                else:
                    other_atoms.append(item)
    break

rows = []

for ha in hotspot_atoms:
    for oa in other_atoms:
        d = float(np.linalg.norm(ha["coord"] - oa["coord"]))
        if d <= 4.0:
            rows.append({
                "hotspot_residue": f'{ha["resname"]}{ha["resid"]}',
                "hotspot_atom": ha["atom"],
                "contact_chain": oa["chain"],
                "contact_residue": f'{oa["resname"]}{oa["resid"]}',
                "contact_atom": oa["atom"],
                "distance_A": round(d,3),
                "contact_type": "rRNA_or_nucleotide" if oa["resname"] in ["A","U","G","C","DA","DT","DG","DC"] else "protein_or_other"
            })

df = pd.DataFrame(rows)
out = Path("05_RESULTS/5OPT_hotspot_contacts_4A.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print(df.head(80).to_string(index=False))
print("\nN contacts:", len(df))
print("Saved:", out)

print("\nContact chains summary:")
print(df.groupby(["contact_chain","contact_type"]).size().reset_index(name="n").to_string(index=False))
