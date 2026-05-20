from Bio.PDB import PDBParser, PDBIO, Select
from pathlib import Path

inp = Path("05_COMPARISON/HUMAN_RPS6/07_CRYOEM_7R4X/7R4X.pdb")
out = Path("05_COMPARISON/HUMAN_RPS6/07_CRYOEM_7R4X/7R4X_chainG_RPS6.pdb")

class ChainSelect(Select):
    def accept_chain(self, chain):
        return chain.id == "G"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("7R4X", inp)

io = PDBIO()
io.set_structure(structure)
io.save(str(out), ChainSelect())

print("Saved:", out)
