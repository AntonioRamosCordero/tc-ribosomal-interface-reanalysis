from Bio.PDB import PDBParser, Superimposer
import numpy as np
import pandas as pd
from pathlib import Path

af_pdb = "05_COMPARISON/HUMAN_RPS6/01_INPUT/AF_P62753_HUMAN_RPS6.pdb"
cryo_pdb = "05_COMPARISON/HUMAN_RPS6/07_CRYOEM_7R4X/7R4X_chainG_RPS6.pdb"
full_pdb = "05_COMPARISON/HUMAN_RPS6/07_CRYOEM_7R4X/7R4X.pdb"

pocket_files = {
    "Human_pocket5": "05_COMPARISON/HUMAN_RPS6/04_POCKETS/AF_P62753_HUMAN_RPS6_out/pockets/pocket5_atm.pdb",
    "Human_pocket9": "05_COMPARISON/HUMAN_RPS6/04_POCKETS/AF_P62753_HUMAN_RPS6_out/pockets/pocket9_atm.pdb",
    "Human_pocket12": "05_COMPARISON/HUMAN_RPS6/04_POCKETS/AF_P62753_HUMAN_RPS6_out/pockets/pocket12_atm.pdb",
}

parser = PDBParser(QUIET=True)
af = parser.get_structure("AF_HUMAN", af_pdb)
cryo = parser.get_structure("7R4X_G", cryo_pdb)
full = parser.get_structure("7R4X_full", full_pdb)

def ca_dict(structure):
    d = {}
    for r in structure.get_residues():
        if r.id[0] == " " and "CA" in r:
            d[r.id[1]] = r["CA"]
    return d

af_ca = ca_dict(af)
cryo_ca = ca_dict(cryo)
common = sorted(set(af_ca) & set(cryo_ca))

fixed = [cryo_ca[i] for i in common]
moving = [af_ca[i] for i in common]

sup = Superimposer()
sup.set_atoms(fixed, moving)
rot, tran = sup.rotran

def transform(coords):
    return np.dot(coords, rot) + tran

def read_pocket_coords(path):
    coords = []
    residues = set()
    for line in Path(path).read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM","HETATM")):
            coords.append(np.array([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ]))
            resname = line[17:20].strip()
            resid = int(line[22:26])
            chain = line[21].strip()
            residues.add(f"{resname}{resid}{chain}")
    return np.array(coords), sorted(residues)

rna_resnames = {"A","U","G","C","DA","DT","DG","DC"}
rna_atoms = []

for atom in full.get_atoms():
    res = atom.get_parent()
    chain = res.get_parent()
    if res.resname.strip() in rna_resnames:
        rna_atoms.append({
            "chain": chain.id,
            "resname": res.resname.strip(),
            "resid": res.id[1],
            "atom": atom.name,
            "coord": atom.coord
        })

rows = []
contact_rows = []

for pocket_name, path in pocket_files.items():
    coords, residues = read_pocket_coords(path)
    coords_7r4x = transform(coords)
    center = coords_7r4x.mean(axis=0)

    min_d = 999
    closest = None
    contacts = []

    for ra in rna_atoms:
        dists = np.linalg.norm(coords_7r4x - ra["coord"], axis=1)
        d = float(dists.min())

        if d < min_d:
            min_d = d
            closest = ra

        if d <= 5.0:
            contacts.append((d, ra))

    rows.append({
        "pocket": pocket_name,
        "n_pocket_atoms": len(coords),
        "pocket_residues": "; ".join(residues),
        "AF_to_7R4X_CA_RMSD_A": sup.rms,
        "n_common_CA": len(common),
        "center_x_7R4X": center[0],
        "center_y_7R4X": center[1],
        "center_z_7R4X": center[2],
        "min_distance_to_rRNA_A": min_d,
        "closest_rRNA": f"{closest['resname']}{closest['resid']}{closest['chain']}:{closest['atom']}" if closest else None,
        "n_rRNA_atoms_within_5A": len(contacts)
    })

    for d, ra in sorted(contacts, key=lambda x: x[0])[:100]:
        contact_rows.append({
            "pocket": pocket_name,
            "distance_A": d,
            "rRNA_chain": ra["chain"],
            "rRNA_residue": f"{ra['resname']}{ra['resid']}",
            "rRNA_atom": ra["atom"]
        })

summary = pd.DataFrame(rows)
contacts = pd.DataFrame(contact_rows)

summary_out = "05_COMPARISON/HUMAN_RPS6/05_RESULTS/Human_pockets5_9_12_distance_to_7R4X_rRNA.csv"
contacts_out = "05_COMPARISON/HUMAN_RPS6/05_RESULTS/Human_pockets5_9_12_rRNA_contacts_5A.csv"

summary.to_csv(summary_out, index=False)
contacts.to_csv(contacts_out, index=False)

print(summary.to_string(index=False))
print("\nTop contacts:")
print(contacts.head(60).to_string(index=False) if len(contacts) else "No contacts within 5A")

print("\nSaved:")
print(summary_out)
print(contacts_out)
