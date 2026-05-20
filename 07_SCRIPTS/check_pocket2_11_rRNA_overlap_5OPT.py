from Bio.PDB import PDBParser, Superimposer
import numpy as np
import pandas as pd
from pathlib import Path

tc_af_pdb = "01_INPUT/AF_Q4DSU0_validated_S6.pdb"
s6_5opt_pdb = "5OPT_chainP_validated_S6.pdb"
full_5opt_pdb = "../12_CRYOEM_5OPT_CONTEXT/01_PDB/5OPT.pdb"

# fallback if running from project root
if not Path(full_5opt_pdb).exists():
    full_5opt_pdb = "/mnt/c/Users/panco/Desktop/PPS6_PROJECT_CLEAN/12_CRYOEM_5OPT_CONTEXT/01_PDB/5OPT.pdb"

if not Path(s6_5opt_pdb).exists():
    s6_5opt_pdb = "01_INPUT/5OPT_chainP_validated_S6.pdb"

pocket_files = {
    "TcS6_pocket2": "04_POCKETS/FPOCKET_Q4DSU0/AF_Q4DSU0_validated_S6_out/pockets/pocket2_atm.pdb",
    "TcS6_pocket11": "04_POCKETS/FPOCKET_Q4DSU0/AF_Q4DSU0_validated_S6_out/pockets/pocket11_atm.pdb",
}

parser = PDBParser(QUIET=True)

af = parser.get_structure("AF_Q4DSU0", tc_af_pdb)
cryo_s6 = parser.get_structure("5OPT_chainP", s6_5opt_pdb)
full = parser.get_structure("5OPT_full", full_5opt_pdb)

def ca_dict(structure):
    d = {}
    for r in structure.get_residues():
        if r.id[0] == " " and "CA" in r:
            d[r.id[1]] = r["CA"]
    return d

af_ca = ca_dict(af)
cryo_ca = ca_dict(cryo_s6)
common = sorted(set(af_ca) & set(cryo_ca))

fixed = [cryo_ca[i] for i in common]
moving = [af_ca[i] for i in common]

sup = Superimposer()
sup.set_atoms(fixed, moving)

# Load pocket coordinates from AF coordinate frame and transform into 5OPT frame
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

rot, tran = sup.rotran

# Biopython applies coord = coord @ rot + tran
def transform(coords):
    return np.dot(coords, rot) + tran

# rRNA atoms in full 5OPT
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

for pocket_name, pocket_path in pocket_files.items():
    coords, residues = read_pocket_coords(pocket_path)
    coords_5opt = transform(coords)
    center = coords_5opt.mean(axis=0)

    min_d = 999
    closest = None
    contacts = []

    for ra in rna_atoms:
        dists = np.linalg.norm(coords_5opt - ra["coord"], axis=1)
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
        "center_x_5OPT": center[0],
        "center_y_5OPT": center[1],
        "center_z_5OPT": center[2],
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

summary_out = "06_RESULTS/Q4DSU0_pocket2_11_distance_to_5OPT_rRNA.csv"
contacts_out = "06_RESULTS/Q4DSU0_pocket2_11_rRNA_contacts_5A.csv"

summary.to_csv(summary_out, index=False)
contacts.to_csv(contacts_out, index=False)

print(summary.to_string(index=False))
print("\nTop contacts:")
print(contacts.head(40).to_string(index=False) if len(contacts) else "No contacts within 5A")

print("\nSaved:")
print(summary_out)
print(contacts_out)
