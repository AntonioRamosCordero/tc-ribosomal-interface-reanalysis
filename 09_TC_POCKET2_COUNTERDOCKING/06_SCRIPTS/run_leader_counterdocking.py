from pathlib import Path
import subprocess
import re
import pandas as pd

base = Path("09_TC_POCKET2_COUNTERDOCKING")

ligands = sorted((base / "03_LIGANDS").glob("*.pdbqt"))

configs = {
    "TcS6_pocket2": base / "02_BOXES/vina_config_TcS6_pocket2.txt",

    "TcS6_pocket11": base / "02_BOXES/vina_config_TcS6_pocket11.txt",

    "Human_pocket9": base / "02_BOXES/vina_config_Human_pocket9.txt",

    "Human_pocket12": base / "02_BOXES/vina_config_Human_pocket12.txt",

    "Human_pocket5": base / "02_BOXES/vina_config_Human_pocket5.txt",
}

outdir = base / "04_DOCKING/leader_test"
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

for lig in ligands:

    lig_name = lig.stem

    for pocket_name, cfg in configs.items():

        out_pdbqt = outdir / f"{lig_name}__{pocket_name}_out.pdbqt"

        log_file = outdir / f"{lig_name}__{pocket_name}_log.txt"

        cmd = [
            "vina",
            "--config", str(cfg),
            "--ligand", str(lig),
            "--out", str(out_pdbqt),
        ]

        print("Running:", lig_name, "->", pocket_name)

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

df = pd.DataFrame(rows)

pivot = df.pivot(
    index="ligand",
    columns="pocket",
    values="best_score_kcal_mol"
).reset_index()

pivot["best_human_score"] = pivot[
    ["Human_pocket9","Human_pocket12","Human_pocket5"]
].min(axis=1)

pivot["best_human_pocket"] = pivot[
    ["Human_pocket9","Human_pocket12","Human_pocket5"]
].idxmin(axis=1)

pivot["best_Tc_score"] = pivot[
    ["TcS6_pocket2","TcS6_pocket11"]
].min(axis=1)

pivot["best_Tc_pocket"] = pivot[
    ["TcS6_pocket2","TcS6_pocket11"]
].idxmin(axis=1)

pivot["delta_best_human_minus_best_Tc"] = (
    pivot["best_human_score"] - pivot["best_Tc_score"]
)

out_csv = base / "05_RESULTS/leader_counterdocking_TcPocket2_11_vs_HumanPockets.csv"

df.to_csv(
    base / "05_RESULTS/leader_counterdocking_raw.csv",
    index=False
)

pivot.sort_values(
    "delta_best_human_minus_best_Tc",
    ascending=False
).to_csv(
    out_csv,
    index=False
)

print("\nSUMMARY:")

print(
    pivot.sort_values(
        "delta_best_human_minus_best_Tc",
        ascending=False
    ).to_string(index=False)
)

print("\nSaved:", out_csv)
