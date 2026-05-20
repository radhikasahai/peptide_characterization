# Peptide characterization

Small Python toolkit that connects **one-letter peptide sequences** to **RDKit molecules**: composition summaries, canonical SMILES from **iterative amide coupling**, and cheminformatics descriptors / visualization helpers.

## Chemistry model (what the code actually does)

- **Standard 20 L-amino acids** only (validated single-letter alphabet).
- **Linear peptides**: each extension couples the chain **C-terminal carboxyl** to the incoming residue **α-amine** using RDKit reaction SMARTS (handles **Gly** and **Pro** as distinct cases).
- **Asp / Glu**: side-chain carboxylates are represented as **methyl esters in the building blocks** so backbone coupling stays **regioselective**. After assembly, methyl esters are **hydrolyzed back to acids by default** (`hydrolyze_sidechains=True` in `sequence_to_rdkit_mol` / `sequence_to_smiles`).


## Setup

Create the conda environment (see `environment.yml`):

```bash
conda env create -f environment.yml
conda activate peptide
```

## Command-line demos

```bash
python peptide_parser.py
python peptide_visual.py
```

## Web UI (Streamlit)

```bash
streamlit run app/streamlit_app.py
```

## Tests & CI

```bash
pytest tests/
```

GitHub Actions runs the same tests using Micromamba (`.github/workflows/ci.yml`).

## Layout

| Module | Role |
|--------|------|
| `chem_utils.py` | Minimal SMILES helpers shared across scripts |
| `peptide_parser.py` | Sequence validation, composition, **RDKit peptide assembly** |
| `peptide_visual.py` | Descriptors, fingerprints, 2D/3D export helpers |
| `app/streamlit_app.py` | Browser UI combining parser + visualization |
