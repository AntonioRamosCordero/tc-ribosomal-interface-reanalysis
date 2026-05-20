from pathlib import Path
import numpy as np
import pandas as pd
import re

dx_path = Path("05_COMPARISON/HUMAN_RPS6/03_APBS/HUMAN_RPS6_0.15M.dx")
pocket_dir = Path("05_COMPARISON/HUMAN_RPS6/04_POCKETS/AF_P62753_HUMAN_RPS6_out/pockets")

def read_dx(path):
    lines = path.read_text(errors="ignore").splitlines()
    counts = None
    origin = None
    deltas = []
    data_start = None

    for i, line in enumerate(lines):
        if line.startswith("object 1 class gridpositions counts"):
            counts = tuple(map(int, line.split()[-3:]))
        elif line.startswith("origin"):
            origin = np.array(list(map(float, line.split()[1:4])))
        elif line.startswith("delta"):
            deltas.append(np.array(list(map(float, line.split()[1:4]))))
        elif line.startswith("object 3 class array"):
            data_start = i + 1
            break

    data = []
    for line in lines[data_start:]:
        if line.startswith("attribute") or line.startswith("object"):
            break
        data.extend([float(x) for x in line.split()])

    grid = np.array(data).reshape(counts)
    spacing = np.array([d[np.nonzero(d)][0] for d in deltas])
    return grid, origin, spacing, counts

grid, origin, spacing, counts = read_dx(dx_path)

def sample(coord):
    idx = np.round((coord - origin) / spacing).astype(int)
    idx = np.clip(idx, [0,0,0], np.array(counts)-1)
    return float(grid[tuple(idx)])

rows = []

for pdb in sorted(pocket_dir.glob("pocket*_atm.pdb")):
    m = re.search(r"pocket(\d+)_atm", pdb.name)
    pocket_id = int(m.group(1)) if m else -1

    coords = []
    residues = set()

    for line in pdb.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM","HETATM")):
            try:
                coords.append(np.array([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54])
                ]))
                resname = line[17:20].strip()
                chain = line[21].strip()
                resid = line[22:26].strip()
                residues.add(f"{resname}{resid}{chain}")
            except:
                pass

    if not coords:
        continue

    values = np.array([sample(c) for c in coords])

    rows.append({
        "species": "human",
        "pocket": pocket_id,
        "n_atoms": len(coords),
        "n_residues": len(residues),
        "mean_APBS_0p15M": values.mean(),
        "median_APBS_0p15M": np.median(values),
        "min_APBS_0p15M": values.min(),
        "max_APBS_0p15M": values.max(),
        "fraction_positive": float((values > 0).mean()),
        "fraction_gt_25": float((values > 25).mean()),
        "fraction_gt_50": float((values > 50).mean()),
        "fraction_gt_100": float((values > 100).mean()),
        "residue_list": "; ".join(sorted(residues))
    })

df = pd.DataFrame(rows)

out = Path("05_COMPARISON/HUMAN_RPS6/05_RESULTS/HUMAN_fpocket_APBS_potential_ranking.csv")
df.sort_values(
    ["fraction_gt_50","median_APBS_0p15M","mean_APBS_0p15M"],
    ascending=[False,False,False]
).to_csv(out, index=False)

print(
    df.sort_values(
        ["fraction_gt_50","median_APBS_0p15M","mean_APBS_0p15M"],
        ascending=[False,False,False]
    )[
        [
            "pocket","n_atoms","n_residues",
            "mean_APBS_0p15M","median_APBS_0p15M",
            "max_APBS_0p15M","fraction_positive",
            "fraction_gt_50","fraction_gt_100",
            "residue_list"
        ]
    ].head(20).to_string(index=False)
)

print("\nSaved:", out)
