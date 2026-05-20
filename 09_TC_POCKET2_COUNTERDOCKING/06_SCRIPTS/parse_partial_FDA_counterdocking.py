from pathlib import Path
import re
import pandas as pd

base = Path("09_TC_POCKET2_COUNTERDOCKING")

dockdir = base / "04_DOCKING/FDA_ALL_TcPocket2_11_vs_Human"

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

for log in sorted(dockdir.glob("*_log.txt")):
    name = log.name.replace("_log.txt", "")
    if "__" not in name:
        continue

    ligand, pocket = name.split("__", 1)

    score = parse_best_score(log.read_text(errors="ignore"))

    rows.append({
        "ligand": ligand,
        "pocket": pocket,
        "best_score_kcal_mol": score,
        "log_file": str(log)
    })

raw = pd.DataFrame(rows)

raw_out = base / "05_RESULTS/PARTIAL_FDA_counterdocking_Tc2_11_vs_Human_raw.csv"
raw.to_csv(raw_out, index=False)

pivot = raw.pivot_table(
    index="ligand",
    columns="pocket",
    values="best_score_kcal_mol",
    aggfunc="first"
).reset_index()

required = [
    "TcS6_pocket2",
    "TcS6_pocket11",
    "Human_pocket9",
    "Human_pocket12",
    "Human_pocket5"
]

for col in required:
    if col not in pivot.columns:
        pivot[col] = None

complete = pivot.dropna(subset=required).copy()

complete["best_Tc_score"] = complete[
    ["TcS6_pocket2", "TcS6_pocket11"]
].min(axis=1)

complete["best_Tc_pocket"] = complete[
    ["TcS6_pocket2", "TcS6_pocket11"]
].idxmin(axis=1)

complete["best_human_score"] = complete[
    ["Human_pocket9", "Human_pocket12", "Human_pocket5"]
].min(axis=1)

complete["best_human_pocket"] = complete[
    ["Human_pocket9", "Human_pocket12", "Human_pocket5"]
].idxmin(axis=1)

complete["delta_best_human_minus_best_Tc"] = (
    complete["best_human_score"] - complete["best_Tc_score"]
)

complete["Tc_selective_candidate"] = (
    complete["delta_best_human_minus_best_Tc"] > 0
)

out = base / "05_RESULTS/PARTIAL_FDA_TcPocket2_11_vs_Human_selectivity.csv"

complete.sort_values(
    ["delta_best_human_minus_best_Tc", "best_Tc_score"],
    ascending=[False, True]
).to_csv(out, index=False)

print("Total logs:", len(raw))
print("Ligands with complete 5-pocket docking:", len(complete))

print("\nTop partial candidates:")
print(
    complete.sort_values(
        ["delta_best_human_minus_best_Tc", "best_Tc_score"],
        ascending=[False, True]
    ).head(40).to_string(index=False)
)

print("\nSaved:")
print(raw_out)
print(out)
