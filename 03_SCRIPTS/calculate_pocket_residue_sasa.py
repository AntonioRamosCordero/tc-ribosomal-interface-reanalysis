from Bio.PDB import PDBParser, ShrakeRupley
import pandas as pd
from pathlib import Path
import re

systems = {
    "TcS6_Q4DSU0": {
        "pdb": "01_INPUT/AF_Q4DSU0_validated_S6.pdb",
        "ranking": "06_RESULTS/TcS6_vs_Human_pocket_candidate_ranking.csv",
        "pockets_dir": "04_POCKETS/FPOCKET_Q4DSU0/AF_Q4DSU0_validated_S6_out/pockets"
    },
    "Human_RPS6": {
        "pdb": "05_COMPARISON/HUMAN_RPS6/01_INPUT/AF_P62753_HUMAN_RPS6.pdb",
        "ranking": "06_RESULTS/TcS6_vs_Human_pocket_candidate_ranking.csv",
        "pockets_dir": "05_COMPARISON/HUMAN_RPS6/04_POCKETS/AF_P62753_HUMAN_RPS6_out/pockets"
    }
}

def get_pocket_residues(pocket_file):
    residues = set()
    for line in Path(pocket_file).read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            chain = line[21].strip()
            resid = int(line[22:26].strip())
            residues.add((chain, resid, resname))
    return residues

rows = []

for species, cfg in systems.items():
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(species, cfg["pdb"])

    sr = ShrakeRupley()
    sr.compute(structure, level="R")

    residue_sasa = {}
    for res in structure.get_residues():
        if res.id[0] != " ":
            continue
        chain = res.get_parent().id
        resid = res.id[1]
        resname = res.resname
        residue_sasa[(chain, resid, resname)] = getattr(res, "sasa", 0.0)

    pockets_dir = Path(cfg["pockets_dir"])

    for pocket_file in sorted(pockets_dir.glob("pocket*_atm.pdb")):
        m = re.search(r"pocket(\d+)_atm", pocket_file.name)
        if not m:
            continue
        pocket_id = int(m.group(1))

        residues = get_pocket_residues(pocket_file)
        sasas = [residue_sasa.get(r, 0.0) for r in residues]

        basic = [r for r in residues if r[2] in ["ARG", "LYS", "HIS"]]
        basic_sasas = [residue_sasa.get(r, 0.0) for r in basic]

        rows.append({
            "species": species,
            "pocket": pocket_id,
            "n_pocket_residues": len(residues),
            "mean_residue_SASA": sum(sasas)/len(sasas) if sasas else 0,
            "median_residue_SASA": sorted(sasas)[len(sasas)//2] if sasas else 0,
            "total_pocket_residue_SASA": sum(sasas),
            "n_basic_residues": len(basic),
            "mean_basic_residue_SASA": sum(basic_sasas)/len(basic_sasas) if basic_sasas else 0,
            "total_basic_residue_SASA": sum(basic_sasas),
            "residue_SASA_list": "; ".join(
                f"{resname}{resid}{chain}:{residue_sasa.get((chain,resid,resname),0.0):.1f}"
                for chain, resid, resname in sorted(residues, key=lambda x: x[1])
            )
        })

df = pd.DataFrame(rows)

profiles = pd.read_csv("06_RESULTS/TcS6_vs_Human_pocket_candidate_ranking.csv")

out = profiles.merge(df, on=["species", "pocket"], how="left")

out_path = "06_RESULTS/TcS6_vs_Human_pocket_profiles_with_SASA.csv"
out.to_csv(out_path, index=False)

cols = [
    "species", "pocket",
    "overall_candidate_score",
    "electropositive_patch_overlap_score",
    "druggability_score",
    "volume",
    "total_SASA",
    "mean_residue_SASA",
    "total_pocket_residue_SASA",
    "n_basic_residues_x",
    "mean_basic_residue_SASA",
    "total_basic_residue_SASA",
    "residue_list"
]

print(out[cols].head(20).to_string(index=False))
print("\nSaved:", out_path)
