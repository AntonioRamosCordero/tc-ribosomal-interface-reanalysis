from pathlib import Path
import subprocess
import re
import pandas as pd

base = Path("09_TC_POCKET2_COUNTERDOCKING")

ligand_dir = Path(
    "/mnt/c/Users/panco/Desktop/PPS6_PROJECT_CLEAN/09_SBVS_FDA_APPROVED/05_LIGANDS_PDBQT/FDA_ALL"
)

ligands = sorted(ligand_dir.glob("*.pdbqt"))

configs = {
    "TcS6_pocket2": base / "02_BOXES/vina_config_TcS6_pocket2.txt",

    "TcS6_pocket11": base / "02_BOXES/vina_config_TcS6_pocket11.txt",

    "Human_pocket9": base / "02_BOXES/vina_config_Human_pocket9.txt",

    "Human_pocket12": base / "02_BOXES/vina_config_Human_pocket12.txt",

    "Human_pocket5": base / "02_BOXES/vina_config_Human_pocket5.txt",
}

outdir = base / "04_DOCKING/FDA_ALL_TcPocket2_11_vs_Human"

outdir.mkdir(parents=True, exist_ok=True)

rows = []

def parse_best_score(log_text):

    for line in log_text.splitlines():

        if re.match(r"\s*1\s+", line):

            parts = line.split()

            try:
                return float(parts[1])

            except:
                return None

    return None

for i, lig in enumerate(ligands, start=1):

    lig_name = lig.stem

    if i % 25 == 0:

        print(f"Progress: {i}/{len(ligands)} ligands")

    for pocket_name, cfg in configs.items():

        out_pdbqt = outdir / f"{lig_name}__{pocket_name}_out.pdbqt"

        log_file = outdir / f"{lig_name}__{pocket_name}_log.txt"

        if log_file.exists():

            log_text = log_file.read_text(errors="ignore")

            score = parse_best_score(log_text)

        else:

            cmd = [
                "vina",
                "--config", str(cfg),
                "--ligand", str(lig),
                "--out", str(out_pdbqt),
            ]

            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )

            log_file.write_text(result.stdout)

            score = parse_best_score(result.stdout)

        rows.append({
            "ligand": lig_name,
            "pocket": pocket_name,
            "best_score_kcal_mol": score,
            "out_pdbqt": str(out_pdbqt),
            "log_file": str(log_file)
        })

raw = pd.DataFrame(rows)

raw_out = base / "05_RESULTS/FDA_ALL_counterdocking_Tc2_11_vs_Human_raw.csv"

raw.to_csv(raw_out, index=False)

pivot = raw.pivot(
    index="ligand",
    columns="pocket",
    values="best_score_kcal_mol"
).reset_index()

pivot["best_Tc_score"] = pivot[
    ["TcS6_pocket2", "TcS6_pocket11"]
].min(axis=1)

pivot["best_Tc_pocket"] = pivot[
    ["TcS6_pocket2", "TcS6_pocket11"]
].idxmin(axis=1)

pivot["best_human_score"] = pivot[
    ["Human_pocket9", "Human_pocket12", "Human_pocket5"]
].min(axis=1)

pivot["best_human_pocket"] = pivot[
    ["Human_pocket9", "Human_pocket12", "Human_pocket5"]
].idxmin(axis=1)

pivot["delta_best_human_minus_best_Tc"] = (
    pivot["best_human_score"] - pivot["best_Tc_score"]
)

pivot["Tc_selective_candidate"] = (
    pivot["delta_best_human_minus_best_Tc"] > 0
)

out = base / "05_RESULTS/FDA_ALL_TcPocket2_11_vs_Human_selectivity.csv"

pivot.sort_values(
    ["delta_best_human_minus_best_Tc", "best_Tc_score"],
    ascending=[False, True]
).to_csv(
    out,
    index=False
)

print("\nTop 40 Tc-selective / least human-favored candidates:")

print(
    pivot.sort_values(
        ["delta_best_human_minus_best_Tc", "best_Tc_score"],
        ascending=[False, True]
    ).head(40).to_string(index=False)
)

print("\nSaved:")
print(raw_out)
print(out)
