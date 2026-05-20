import pandas as pd
from pathlib import Path

boxes = pd.read_csv("09_TC_POCKET2_COUNTERDOCKING/02_BOXES/pocket_docking_boxes.csv")

outdir = Path("09_TC_POCKET2_COUNTERDOCKING/02_BOXES")
outdir.mkdir(parents=True, exist_ok=True)

for _, row in boxes.iterrows():
    pocket = row["pocket"]

    if pocket.startswith("TcS6"):
        receptor = "09_TC_POCKET2_COUNTERDOCKING/00_RECEPTORS/TcS6_Q4DSU0_receptor.pdbqt"
    else:
        receptor = "09_TC_POCKET2_COUNTERDOCKING/00_RECEPTORS/Human_RPS6_receptor.pdbqt"

    txt = f"""receptor = {receptor}

center_x = {row['center_x']:.3f}
center_y = {row['center_y']:.3f}
center_z = {row['center_z']:.3f}

size_x = {row['size_x']:.3f}
size_y = {row['size_y']:.3f}
size_z = {row['size_z']:.3f}

exhaustiveness = 8
num_modes = 9
energy_range = 3
"""

    path = outdir / f"vina_config_{pocket}.txt"
    path.write_text(txt)

    print("Saved:", path)
