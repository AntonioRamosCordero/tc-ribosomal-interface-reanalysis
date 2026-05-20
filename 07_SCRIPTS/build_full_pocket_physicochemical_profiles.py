import re
from pathlib import Path
import pandas as pd

def parse_fpocket_info(path, species):
    lines = Path(path).read_text(errors="ignore").splitlines()

    rows = []
    current = None

    keymap = {
        "Score": "fpocket_score",
        "Druggability Score": "druggability_score",
        "Number of Alpha Spheres": "n_alpha_spheres",
        "Total SASA": "total_SASA",
        "Polar SASA": "polar_SASA",
        "Apolar SASA": "apolar_SASA",
        "Volume": "volume",
        "Mean local hydrophobic density": "mean_local_hydrophobic_density",
        "Mean alpha sphere radius": "mean_alpha_sphere_radius",
        "Mean alp. sph. solvent access": "mean_alpha_sphere_solvent_access",
        "Apolar alpha sphere proportion": "apolar_alpha_sphere_proportion",
        "Hydrophobicity score": "hydrophobicity_score",
        "Volume score": "volume_score",
        "Polarity score": "polarity_score",
        "Charge score": "charge_score",
        "Proportion of polar atoms": "proportion_polar_atoms",
        "Alpha sphere density": "alpha_sphere_density",
        "Flexibility": "flexibility",
    }

    for line in lines:

        m = re.match(r"Pocket\s+(\d+)\s+:", line)

        if m:
            if current:
                rows.append(current)

            current = {
                "species": species,
                "pocket": int(m.group(1))
            }

            continue

        if current is None:
            continue

        if ":" not in line:
            continue

        k, v = line.split(":", 1)

        k = k.strip()
        v = v.strip()

        if k in keymap:

            try:
                current[keymap[k]] = float(v)

            except:
                current[keymap[k]] = v

    if current:
        rows.append(current)

    return pd.DataFrame(rows)


def classify_residues(residue_list):

    basics = ["ARG", "LYS", "HIS"]
    acids = ["ASP", "GLU"]
    polars = ["ASN", "GLN", "SER", "THR", "CYS", "TYR", "HIS"]
    aromatics = ["PHE", "TYR", "TRP", "HIS"]
    hydrophobics = ["ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP"]

    residues = [
        x.strip()
        for x in str(residue_list).split(";")
        if x.strip()
    ]

    names = [r[:3] for r in residues]

    n = len(names)

    if n == 0:
        n = 1

    return pd.Series({

        "n_basic_residues":
            sum(x in basics for x in names),

        "n_acidic_residues":
            sum(x in acids for x in names),

        "n_polar_residues":
            sum(x in polars for x in names),

        "n_aromatic_residues":
            sum(x in aromatics for x in names),

        "n_hydrophobic_residues":
            sum(x in hydrophobics for x in names),

        "frac_basic_residues":
            sum(x in basics for x in names) / n,

        "frac_acidic_residues":
            sum(x in acids for x in names) / n,

        "frac_polar_residues":
            sum(x in polars for x in names) / n,

        "frac_aromatic_residues":
            sum(x in aromatics for x in names) / n,

        "frac_hydrophobic_residues":
            sum(x in hydrophobics for x in names) / n,
    })


tc_info = parse_fpocket_info(
    "04_POCKETS/FPOCKET_Q4DSU0/AF_Q4DSU0_validated_S6_out/AF_Q4DSU0_validated_S6_info.txt",
    "TcS6_Q4DSU0"
)

hu_info = parse_fpocket_info(
    "05_COMPARISON/HUMAN_RPS6/04_POCKETS/AF_P62753_HUMAN_RPS6_out/AF_P62753_HUMAN_RPS6_info.txt",
    "Human_RPS6"
)

tc_apbs = pd.read_csv(
    "06_RESULTS/Q4DSU0_fpocket_APBS_potential_ranking.csv"
)
tc_apbs["species"] = "TcS6_Q4DSU0"

hu_apbs = pd.read_csv(
    "05_COMPARISON/HUMAN_RPS6/05_RESULTS/HUMAN_fpocket_APBS_potential_ranking.csv"
)
hu_apbs["species"] = "Human_RPS6"

tc = tc_info.merge(
    tc_apbs,
    on=["species", "pocket"],
    how="left"
)

hu = hu_info.merge(
    hu_apbs,
    on=["species", "pocket"],
    how="left"
)

df = pd.concat([tc, hu], ignore_index=True)

chem = df["residue_list"].apply(classify_residues)

df = pd.concat([df, chem], axis=1)

df["electropositive_patch_overlap_score"] = (
    df["fraction_gt_50"].fillna(0) * 0.5
    + df["fraction_gt_100"].fillna(0) * 0.3
    + df["fraction_positive"].fillna(0) * 0.2
)

df["geometry_exposure_score"] = (
    df["druggability_score"].fillna(0) * 0.35
    + (
        df["volume"].fillna(0)
        / df["volume"].fillna(0).max()
    ) * 0.25
    + df["mean_alpha_sphere_solvent_access"].fillna(0) * 0.25
    + (
        df["total_SASA"].fillna(0)
        / df["total_SASA"].fillna(0).max()
    ) * 0.15
)

df["overall_candidate_score"] = (
    df["electropositive_patch_overlap_score"] * 0.45
    + df["geometry_exposure_score"] * 0.35
    + df["frac_basic_residues"] * 0.20
)

out_all = (
    "06_RESULTS/"
    "TcS6_vs_Human_full_pocket_physicochemical_profiles.csv"
)

df.to_csv(out_all, index=False)

ranked = df.sort_values(
    [
        "overall_candidate_score",
        "electropositive_patch_overlap_score",
        "druggability_score"
    ],
    ascending=[False, False, False]
)

out_ranked = (
    "06_RESULTS/"
    "TcS6_vs_Human_pocket_candidate_ranking.csv"
)

ranked.to_csv(out_ranked, index=False)

cols = [
    "species",
    "pocket",
    "overall_candidate_score",
    "electropositive_patch_overlap_score",
    "geometry_exposure_score",
    "druggability_score",
    "volume",
    "total_SASA",
    "polar_SASA",
    "apolar_SASA",
    "fraction_positive",
    "fraction_gt_50",
    "fraction_gt_100",
    "n_basic_residues",
    "n_acidic_residues",
    "n_polar_residues",
    "n_aromatic_residues",
    "n_hydrophobic_residues",
    "residue_list"
]

print(
    ranked[cols]
    .head(20)
    .to_string(index=False)
)

print("\nSaved:")
print(out_all)
print(out_ranked)
