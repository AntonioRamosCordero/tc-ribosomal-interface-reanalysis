from Bio.PDB import PDBParser, Superimposer, PDBIO
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio import pairwise2
import pandas as pd
import numpy as np
from pathlib import Path

tc_pdb = "01_INPUT/AF_Q4DSU0_validated_S6.pdb"
hu_pdb = "05_COMPARISON/HUMAN_RPS6/01_INPUT/AF_P62753_HUMAN_RPS6.pdb"

parser = PDBParser(QUIET=True)
tc = parser.get_structure("TcS6", tc_pdb)
hu = parser.get_structure("HumanRPS6", hu_pdb)

def get_seq_map(structure):
    seq = ""
    mapping = []
    for model in structure:
        for chain in model:
            for r in chain:
                if r.id[0] != " ":
                    continue
                aa = protein_letters_3to1.get(r.resname, "X")
                seq += aa
                mapping.append({
                    "seq_index": len(seq),
                    "chain": chain.id,
                    "resid": r.id[1],
                    "resname": r.resname,
                    "aa": aa,
                    "residue": r
                })
            break
        break
    return seq, mapping

tc_seq, tc_map = get_seq_map(tc)
hu_seq, hu_map = get_seq_map(hu)

aln = pairwise2.align.globalxx(tc_seq, hu_seq, one_alignment_only=True)[0]
a_tc, a_hu, score, start, end = aln

i_tc = 0
i_hu = 0
pairs = []

for pos, (ct, ch) in enumerate(zip(a_tc, a_hu), start=1):
    tc_info = None
    hu_info = None

    if ct != "-":
        tc_info = tc_map[i_tc]
        i_tc += 1

    if ch != "-":
        hu_info = hu_map[i_hu]
        i_hu += 1

    if tc_info and hu_info:
        if "CA" in tc_info["residue"] and "CA" in hu_info["residue"]:
            pairs.append({
                "alignment_pos": pos,
                "tc_resid": tc_info["resid"],
                "tc_resname": tc_info["resname"],
                "tc_aa": tc_info["aa"],
                "human_resid": hu_info["resid"],
                "human_resname": hu_info["resname"],
                "human_aa": hu_info["aa"],
                "match": tc_info["aa"] == hu_info["aa"],
                "tc_CA": tc_info["residue"]["CA"],
                "hu_CA": hu_info["residue"]["CA"],
            })

fixed = [p["tc_CA"] for p in pairs]
moving = [p["hu_CA"] for p in pairs]

sup = Superimposer()
sup.set_atoms(fixed, moving)
rmsd = sup.rms

# apply human -> Tc coordinate frame
sup.apply(hu.get_atoms())

outdir = Path("05_COMPARISON/HUMAN_RPS6/05_RESULTS")
outdir.mkdir(parents=True, exist_ok=True)

# save aligned human
io = PDBIO()
io.set_structure(hu)
io.save(str(outdir / "AF_P62753_HUMAN_RPS6_aligned_to_TcS6_Q4DSU0.pdb"))

rows = []
for p in pairs:
    dist = np.linalg.norm(p["tc_CA"].coord - p["hu_CA"].coord)
    rows.append({
        "tc_resid": p["tc_resid"],
        "tc_resname": p["tc_resname"],
        "tc_aa": p["tc_aa"],
        "human_resid": p["human_resid"],
        "human_resname": p["human_resname"],
        "human_aa": p["human_aa"],
        "match": p["match"],
        "CA_distance_after_alignment_A": dist
    })

df = pd.DataFrame(rows)
df.to_csv(outdir / "TcS6_Q4DSU0_vs_HUMAN_RPS6_residue_mapping.csv", index=False)

summary = {
    "TcS6_length": len(tc_seq),
    "Human_RPS6_length": len(hu_seq),
    "aligned_CA_pairs": len(pairs),
    "sequence_identity_on_aligned_pairs": df["match"].mean(),
    "global_CA_RMSD_after_alignment_A": rmsd
}

with open(outdir / "TcS6_Q4DSU0_vs_HUMAN_RPS6_alignment_summary.txt", "w") as f:
    for k, v in summary.items():
        f.write(f"{k}: {v}\n")

print("TcS6 length:", len(tc_seq))
print("Human length:", len(hu_seq))
print("Aligned CA pairs:", len(pairs))
print("Sequence identity:", round(df["match"].mean(), 3))
print("Global CA RMSD:", round(rmsd, 3), "A")
print("\nSaved:")
print(outdir / "AF_P62753_HUMAN_RPS6_aligned_to_TcS6_Q4DSU0.pdb")
print(outdir / "TcS6_Q4DSU0_vs_HUMAN_RPS6_residue_mapping.csv")
print(outdir / "TcS6_Q4DSU0_vs_HUMAN_RPS6_alignment_summary.txt")

# print mapping for Tc electropositive pockets residues of interest
interest = [220,221,224,227,230,231,235,152,162,163]
print("\nResidue mapping for Tc residues of interest:")
print(df[df["tc_resid"].isin(interest)][[
    "tc_resid","tc_resname","human_resid","human_resname","match","CA_distance_after_alignment_A"
]].to_string(index=False))
