from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

pdb_path = Path("01_INPUT/AF_Q4DSU0_validated_S6.pdb")
dx_path = Path("03_APBS/validated_S6_0.15M.dx")
out_csv = Path("06_RESULTS/validated_S6_hotspot_APBS_atom_samples.csv")

hotspot_resids = [149,152,153,161,162,163]

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
    return grid[tuple(idx)], idx

parser = PDBParser(QUIET=True)
s = parser.get_structure("s6", pdb_path)

rows = []
for r in s.get_residues():
    if r.id[0] != " " or r.id[1] not in hotspot_resids:
        continue
    for atom in r:
        val, idx = sample(atom.coord)
        rows.append({
            "resid": r.id[1],
            "resname": r.resname,
            "atom": atom.name,
            "x": atom.coord[0],
            "y": atom.coord[1],
            "z": atom.coord[2],
            "apbs_value_0p15M": val,
            "grid_i": idx[0],
            "grid_j": idx[1],
            "grid_k": idx[2]
        })

df = pd.DataFrame(rows)
df.to_csv(out_csv, index=False)

summary = (
    df.groupby(["resid","resname"])
      .agg(
          n_atoms=("apbs_value_0p15M","count"),
          mean_APBS=("apbs_value_0p15M","mean"),
          median_APBS=("apbs_value_0p15M","median"),
          min_APBS=("apbs_value_0p15M","min"),
          max_APBS=("apbs_value_0p15M","max")
      )
      .reset_index()
)

summary_out = Path("06_RESULTS/validated_S6_hotspot_APBS_residue_summary.csv")
summary.to_csv(summary_out, index=False)

print(summary.to_string(index=False))
print("\nSaved:")
print(out_csv)
print(summary_out)
