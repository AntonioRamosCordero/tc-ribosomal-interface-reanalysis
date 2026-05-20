from Bio.PDB import PDBParser, PDBIO, Select
from pathlib import Path

inp = Path("01_PDB/5OPT.pdb")
out = Path("03_EXTRACTED_S6/5OPT_chainP_S6.pdb")
out.parent.mkdir(parents=True, exist_ok=True)

class ChainSelect(Select):
    def accept_chain(self, chain):
        return chain.id == "P"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("5OPT", inp)

io = PDBIO()
io.set_structure(structure)
io.save(str(out), ChainSelect())

print("Saved:", out)
