# Methods Overview

## Structural models

The project initially used the AlphaFold-predicted structure associated with the *Trypanosoma cruzi* UniProt/AlphaFold entry Q4DSU0, initially explored as a putative TcS6-like ribosomal protein candidate.

Human RPS6 AlphaFold structures were used for comparative analyses.

Later stages incorporated experimentally resolved cryo-EM ribosome structures for structural reevaluation and RNA-interface validation.

Cryo-EM analyses suggested that several initially proposed structural interpretations required reevaluation, particularly regarding functional accessibility and RNA-interface overlap.

---

## Electrostatic calculations

Electrostatic potentials were calculated using:

- PDB2PQR v3.6.1
- APBS v3.4.1
- PARSE force field
- pH 7.0 conditions

Both structural preprocessing and electrostatic grids were generated under reproducible command-line workflows.

---

## Pocket prediction

Predicted pockets were identified using:

- fpocket
- P2Rank

Pocket rankings were compared using:
- electrostatic enrichment
- residue composition
- SASA
- overlap with electropositive regions
- proximity to functional Arg-rich regions

---

## Structural comparison

Comparisons between parasite and human structures included:

- residue mapping
- physicochemical profiling
- solvent exposure
- electrostatic characterization
- pocket ranking
- RNA-interface overlap

---

## Cryo-EM reevaluation

Cryo-EM ribosomal structures were incorporated during later stages of the project to evaluate whether predicted pockets overlapped experimentally resolved rRNA interfaces.

Structures analyzed included:
- 5OPT (*Trypanosoma cruzi* ribosome)
- 7R4X (human ribosome)

---

## Docking and counter-docking

Exploratory docking experiments were performed using AutoDock Vina.

Docking analyses included:
- parasite pockets
- human comparative pockets
- counter-docking comparisons
- preliminary selectivity estimation

Docking results should be interpreted as exploratory structural observations rather than validated binding evidence.

---

## Reproducibility

All analyses were performed through reproducible command-line workflows under Linux/WSL environments using Python-based scripts and open-source structural bioinformatics software.
