from Bio.PDB import PDBParser, Superimposer
import numpy as np
import pandas as pd
from pathlib import Path

af_pdb = "13_VALIDATED_S6_MINI_REANALYSIS/01_INPUT/AF_Q4DSU0_validated_S6.pdb"
cryo_pdb = "13_VALIDATED_S6_MINI_REANALYSIS/01_INPUT/5OPT_chainP_validated_S6.pdb"

parser = PDBParser(QUIET=True)
af = parser.get_structure("AF", af_pdb)
cryo = parser.get_structure("5OPT", cryo_pdb)

def ca_dict(structure):
    d = {}
    for r in structure.get_residues():
        if r.id[0] == " " and "CA" in r:
            d[r.id[1]] = r["CA"]
    return d

af_ca = ca_dict(af)
cryo_ca = ca_dict(cryo)

common = sorted(set(af_ca) & set(cryo_ca))
fixed = [cryo_ca[i] for i in common]
moving = [af_ca[i] for i in common]

sup = Superimposer()
sup.set_atoms(fixed, moving)
global_rmsd = sup.rms

hotspot = [149,152,153,161,162,163]
fixed_h = [cryo_ca[i] for i in hotspot if i in cryo_ca and i in af_ca]
moving_h = [af_ca[i] for i in hotspot if i in cryo_ca and i in af_ca]

sup_h = Superimposer()
sup_h.set_atoms(fixed_h, moving_h)
hotspot_rmsd = sup_h.rms

rows = []
for i in common:
    d = np.linalg.norm(af_ca[i].coord - cryo_ca[i].coord)
    rows.append({
        "resid": i,
        "AF_residue": af_ca[i].get_parent().resname,
        "5OPT_residue": cryo_ca[i].get_parent().resname,
        "CA_distance_unaligned_A": d,
        "is_hotspot": i in hotspot
    })

df = pd.DataFrame(rows)

outdir = Path("13_VALIDATED_S6_MINI_REANALYSIS/06_RESULTS")
outdir.mkdir(parents=True, exist_ok=True)

df.to_csv(outdir / "AF_Q4DSU0_vs_5OPT_CA_distance_unaligned.csv", index=False)

with open(outdir / "AF_Q4DSU0_vs_5OPT_structural_comparison_summary.txt", "w") as f:
    f.write(f"Common CA residues: {len(common)}\n")
    f.write(f"Global CA RMSD after superposition: {global_rmsd:.3f} A\n")
    f.write(f"Hotspot CA RMSD after local superposition: {hotspot_rmsd:.3f} A\n")
    f.write(f"Hotspot residues used: {hotspot}\n")

print("Common CA residues:", len(common))
print("Global CA RMSD after superposition:", round(global_rmsd,3), "A")
print("Hotspot CA RMSD after local superposition:", round(hotspot_rmsd,3), "A")
print("\nSaved results in 13_VALIDATED_S6_MINI_REANALYSIS/06_RESULTS/")
