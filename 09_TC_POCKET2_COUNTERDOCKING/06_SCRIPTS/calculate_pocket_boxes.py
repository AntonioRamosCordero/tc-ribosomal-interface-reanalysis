from pathlib import Path
import numpy as np
import pandas as pd

pockets = {
    "TcS6_pocket2": "09_TC_POCKET2_COUNTERDOCKING/01_POCKETS/TcS6_pocket2_atm.pdb",

    "TcS6_pocket11": "09_TC_POCKET2_COUNTERDOCKING/01_POCKETS/TcS6_pocket11_atm.pdb",

    "Human_pocket9": "09_TC_POCKET2_COUNTERDOCKING/01_POCKETS/Human_pocket9_atm.pdb",

    "Human_pocket12": "09_TC_POCKET2_COUNTERDOCKING/01_POCKETS/Human_pocket12_atm.pdb",

    "Human_pocket5": "09_TC_POCKET2_COUNTERDOCKING/01_POCKETS/Human_pocket5_atm.pdb",
}

rows = []

for name, path in pockets.items():

    coords = []

    for line in Path(path).read_text(errors="ignore").splitlines():

        if line.startswith(("ATOM","HETATM")):

            coords.append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])

    coords = np.array(coords)

    center = coords.mean(axis=0)

    mins = coords.min(axis=0)

    maxs = coords.max(axis=0)

    span = maxs - mins

    # margen razonable para ligandos FDA medianos
    size = np.maximum(span + 8.0, 16.0)

    rows.append({
        "pocket": name,
        "center_x": center[0],
        "center_y": center[1],
        "center_z": center[2],
        "size_x": size[0],
        "size_y": size[1],
        "size_z": size[2],
        "span_x": span[0],
        "span_y": span[1],
        "span_z": span[2],
        "n_atoms": len(coords)
    })

df = pd.DataFrame(rows)

out = "09_TC_POCKET2_COUNTERDOCKING/02_BOXES/pocket_docking_boxes.csv"

df.to_csv(out, index=False)

print(df.to_string(index=False))

print("\nSaved:", out)
