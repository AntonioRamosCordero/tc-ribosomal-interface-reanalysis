from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

files = {
    "AF_TcS6_like": "/mnt/c/Users/panco/Desktop/PPS6_PROJECT_CLEAN/LIGAND_DESIGN_STAGE1/03_RECEPTOR_PREPARATION/TC_PPS6/PPS6_apbs_from_pqr.pdb",
    "5OPT_chainP_S6": "03_EXTRACTED_S6/5OPT_chainP_S6.pdb"
}

def get_seq_and_map(pdb_path):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", pdb_path)

    seq = ""
    mapping = []

    for model in s:
        for chain in model:
            for r in chain:
                if r.id[0] != " ":
                    continue
                resn = r.resname
                aa = protein_letters_3to1.get(resn, "X")
                seq += aa
                mapping.append((len(seq), chain.id, r.id[1], resn, aa))
            break
        break

    return seq, mapping

seq1, map1 = get_seq_and_map(files["AF_TcS6_like"])
seq2, map2 = get_seq_and_map(files["5OPT_chainP_S6"])

aln = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
a1, a2, score, start, end = aln

print("AF length:", len(seq1))
print("5OPT length:", len(seq2))
print("Alignment score:", score)
print(format_alignment(*aln))

# Build alignment position map
i1 = 0
i2 = 0
rows = []

hotspot_AF = {96,97,100,101,102,103,104,105,106}

for pos, (c1, c2) in enumerate(zip(a1, a2), start=1):
    af_info = None
    opt_info = None

    if c1 != "-":
        af_info = map1[i1]
        i1 += 1

    if c2 != "-":
        opt_info = map2[i2]
        i2 += 1

    if af_info and opt_info:
        af_resid = af_info[2]
        opt_resid = opt_info[2]

        if af_resid in hotspot_AF:
            rows.append({
                "AF_resid": af_resid,
                "AF_resn": af_info[3],
                "AF_aa": af_info[4],
                "5OPT_resid": opt_resid,
                "5OPT_resn": opt_info[3],
                "5OPT_aa": opt_info[4],
                "match": af_info[4] == opt_info[4]
            })

print("\nHotspot mapping AF -> 5OPT:")
for r in rows:
    print(r)
