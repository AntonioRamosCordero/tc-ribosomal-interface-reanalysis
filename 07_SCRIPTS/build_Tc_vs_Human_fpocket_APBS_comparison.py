import pandas as pd
from pathlib import Path

tc = pd.read_csv("06_RESULTS/Q4DSU0_fpocket_APBS_potential_ranking.csv")
hu = pd.read_csv("05_COMPARISON/HUMAN_RPS6/05_RESULTS/HUMAN_fpocket_APBS_potential_ranking.csv")

tc["species"] = "TcS6_Q4DSU0"
hu["species"] = "Human_RPS6"

# Keep comparable electrostatic/druggability columns currently available from APBS ranking
cols = [
    "species",
    "pocket",
    "n_atoms",
    "n_residues",
    "mean_APBS_0p15M",
    "median_APBS_0p15M",
    "max_APBS_0p15M",
    "fraction_positive",
    "fraction_gt_50",
    "fraction_gt_100",
    "residue_list"
]

merged = pd.concat([tc[cols], hu[cols]], ignore_index=True)

out_all = Path("06_RESULTS/TcS6_vs_Human_fpocket_APBS_all_ranked.csv")
merged.sort_values(
    ["fraction_gt_50", "median_APBS_0p15M", "mean_APBS_0p15M"],
    ascending=[False, False, False]
).to_csv(out_all, index=False)

# top pockets by electropositive enrichment
top = merged.sort_values(
    ["fraction_gt_50", "median_APBS_0p15M", "mean_APBS_0p15M"],
    ascending=[False, False, False]
).groupby("species").head(6)

out_top = Path("06_RESULTS/TcS6_vs_Human_top_electropositive_pockets.csv")
top.to_csv(out_top, index=False)

print(top[[
    "species", "pocket", "n_residues",
    "median_APBS_0p15M", "fraction_gt_50", "fraction_gt_100",
    "residue_list"
]].to_string(index=False))

print("\nSaved:")
print(out_all)
print(out_top)
