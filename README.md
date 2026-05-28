# Peptide characterization

Small Python toolkit that connects **one-letter peptide sequences** to **RDKit molecules**: composition summaries, canonical SMILES from **iterative amide coupling**, cheminformatics descriptors, and batch export of 2D/3D visuals.

**Live app:** https://radhikasahai.github.io/peptide_characterization/  
**Presentation:** https://radhikasahai.github.io/peptide_characterization/presentation.html

## Chemistry model

- **Standard 20 L-amino acids** only (validated single-letter alphabet).
- **Linear peptides**: each extension couples the chain **C-terminal carboxyl** to the incoming residue **α-amine** using RDKit reaction SMARTS (handles **Gly** and **Pro** as distinct cases).
- **Asp / Glu**: side-chain carboxylates are represented as **methyl esters in the building blocks** so backbone coupling stays **regioselective**. After assembly, methyl esters are **hydrolyzed back to acids by default** (`hydrolyze_sidechains=True` in `sequence_to_rdkit_mol` / `sequence_to_smiles`).

## Setup

From the repository root:

```bash
conda env create -f environment.yml
conda activate peptide
```

## Command-line usage

Run all commands from the **repository root** after activating the `peptide` environment.

### Parse sequences and build SMILES

```bash
python peptide_parser.py
```

### Generate sequences and build structures

```bash
python peptide_synthesis.py
```

### Descriptors, fingerprints, and example exports (aspirin demo)

```bash
python peptide_visual.py
```

### Batch-export peptide visuals from CSV

Input CSV: header row required; columns `id` or `name`, plus `smiles` and/or `sequence` (see `data/peptides_example.csv`).

```bash
python scripts/export_peptide_visuals.py --csv data/peptides_example.csv --out output/peptides
python scripts/export_peptide_visuals.py --csv data/peptides_example.csv --out output/peptides --skip-3d --grid
```

Per row in `--out`:

| File | Description |
|------|-------------|
| `{id}.png` | 2D structure image |
| `{id}_3d.sdf` | 3D conformer (UFF), unless `--skip-3d` |
| `{id}_ligand.sdf` | Docking-style 3D (MMFF), unless `--skip-3d` |
| `export_summary.csv` | Per-row success/failure log |
| `descriptors_summary.csv` | RDKit descriptors for rows with valid SMILES |
| `library_grid.png` | Optional grid (`--grid`) |

### Benchmark and tests

Frozen regression suite: `data/benchmark_sequences.csv` (47 rows: 39 valid assemblies, 8 invalid inputs).

```bash
python scripts/generate_benchmark.py   # regenerate after changing assembly logic
pytest tests/
```

GitHub Actions runs the same tests using Micromamba (`.github/workflows/ci.yml`).

## Web UI (GitHub Pages + FastAPI)

The interactive UI lives in `docs/` (static frontend). Structure building runs on the FastAPI backend in `api/`.

### Run locally

Terminal 1 — API:

```bash
conda activate peptide
pip install -r requirements-api-local.txt
uvicorn api.main:app --reload --port 8000
```

Terminal 2 — frontend:

```bash
cd docs
python -m http.server 5500
```

Open http://localhost:5500 (the frontend defaults to `http://localhost:8000`).

### Deploy

1. **GitHub Pages** — push to `main`; the **Deploy GitHub Pages** workflow publishes `docs/`. Enable Pages once under **Settings → Pages → Build and deployment → Source: GitHub Actions**.
2. **Render API** — create a Web Service from this repo using `render.yaml` / `Dockerfile`. Set `ALLOWED_ORIGINS=https://radhikasahai.github.io`.
3. **Connect frontend** — in the repo, go to **Settings → Secrets and variables → Actions → Variables** and set **`PEPTIDE_API_URL`** to your Render service URL (e.g. `https://peptide-characterization-api.onrender.com`). Re-run the Pages deploy workflow so `docs/js/config.js` is updated.

Verify: open the live site, enter `GG`, and confirm the 2D structure and descriptors load.

### Streamlit (legacy demo)

Run from the **repository root** (so imports resolve):

```bash
conda activate peptide
streamlit run app/streamlit_app.py
```

## Project layout

| Module | Role |
|--------|------|
| `chem_utils.py` | Minimal SMILES helpers |
| `peptide_parser.py` | Sequence validation, composition, RDKit peptide assembly |
| `peptide_synthesis.py` | Random/motif/combinatorial sequences and structure builds |
| `peptide_visual.py` | Descriptors, fingerprints, 2D/3D export |
| `api/main.py` | FastAPI backend for the web UI |
| `docs/` | Static frontend (GitHub Pages) |
| `app/streamlit_app.py` | Legacy Streamlit demo UI |
| `scripts/export_peptide_visuals.py` | Batch PNG/SDF export from CSV |
| `scripts/generate_benchmark.py` | Regenerate benchmark CSV |
| `data/benchmark_sequences.csv` | Committed regression benchmark |
| `data/peptides_example.csv` | Sample CSV for visual export |

## References

Published methods and software behind this codebase. Amino-acid monomer SMILES and amide-coupling reaction SMARTS in `peptide_parser.py` are **project-defined** building blocks (regression-checked in `data/benchmark_sequences.csv`); they implement standard peptide-bond chemistry rather than reproducing a single published structure database.

### Cheminformatics platform and notation

| Topic | Used in | Reference |
|-------|---------|-----------|
| RDKit (parsing, reactions, descriptors, depictions, fingerprints) | `peptide_parser.py`, `peptide_visual.py`, `api/main.py`, `chem_utils.py` | Landrum, G. *RDKit: Open-Source Cheminformatics*; [https://www.rdkit.org](https://www.rdkit.org) — see [RDKit citation info](https://www.rdkit.org/docs/Overview.html#citation-information) (e.g. Zenodo [10.5281/zenodo.591637](https://doi.org/10.5281/zenodo.591637)). |
| SMILES line notation | All SMILES I/O | Weininger, D. *J. Chem. Inf. Comput. Sci.* **1988**, *28*, 31–36. Weininger, D.; Weininger, A.; Weininger, J. L. *J. Chem. Inf. Comput. Sci.* **1989**, *29*, 97–101. |
| SMARTS reaction transforms | `peptide_parser.py` (`ReactionFromSmarts`, amide coupling, methyl-ester hydrolysis) | Daylight SMARTS language; RDKit reaction machinery (Landrum, G., et al., as above). |

### In silico linear peptide assembly

| Topic | Used in | Reference |
|-------|---------|-----------|
| Iterative peptide-bond (amide) formation | `sequence_to_rdkit_mol`, `build_peptide` | Standard solution-phase peptide coupling concept; solid-phase foundation: Merrifield, R. B. *J. Am. Chem. Soc.* **1963**, *85*, 2149–2154. |
| Asp/Glu side-chain masking as methyl esters, then deprotection | `AA_BUILDING_BLOCK_SMILES`, `hydrolyze_sidechains` | Common protecting-group strategy in peptide synthesis; overview: Bodanszky, M. *Principles of Peptide Synthesis* (Springer); Fields, G. B.; Noble, R. L. *Int. J. Pept. Protein Res.* **1990**, *35*, 161–214. |
| One-letter amino acid codes | `validate_sequence`, composition | IUPAC-IUBMB one-letter symbols for the standard amino acids. |

### Molecular descriptors (`calculate_basic_descriptors`)

| Descriptor (code) | Reference |
|-------------------|-----------|
| `logp` (Crippen `MolLogP`) | Wildman, M. D.; Crippen, G. M. *J. Chem. Inf. Comput. Sci.* **1999**, *39*, 868–873. |
| `tpsa` | Ertl, P.; Rohde, B.; Selzer, P. *J. Med. Chem.* **2000**, *43*, 3714–3717. |
| `hbond_donors`, `hbond_acceptors`, `rotatable_bonds` (Lipinski) | Lipinski, C. A.; Lombardo, F.; Dominy, B. W.; Feeney, P. J. *Adv. Drug Delivery Rev.* **1997**, *23*, 3–25. Lipinski, C. A. *J. Pharmacol. Toxicol. Methods* **2001**, *44*, 235–249. |
| `fraction_csp3` | Lovering, F.; Bikker, J.; Humblet, C. *J. Med. Chem.* **2009**, *52*, 6752–6756 (Fsp³ concept; RDKit implementation `CalcFractionCSP3`). |
| Molecular weight / heavy atoms | Standard RDKit `Descriptors` (see RDKit citation). |

### 3D structures and similarity (`peptide_visual.py`)

| Topic | Used in | Reference |
|-------|---------|-----------|
| 3D embedding (`EmbedMolecule`) | `_embed_3d`, `_embed_docking_3d` | Riniker, S.; Landrum, G. A. *J. Chem. Inf. Model.* **2015**, *55*, 2562–2574 (ETKDG family; RDKit default embed protocol). |
| UFF relaxation | `{id}_3d.sdf` | Rappe, A. K.; Casewit, C. J.; Colwell, K. S.; et al. *J. Am. Chem. Soc.* **1992**, *114*, 10024–10035. |
| MMFF relaxation | `{id}_ligand.sdf` | Halgren, T. A. *J. Comput. Chem.* **1996**, *17*, 490–519. |
| Morgan fingerprint (`GetMorganGenerator`) | `generate_morgan_fingerprint` | Rogers, D.; Hahn, M. *J. Chem. Inf. Model.* **2010**, *50*, 742–754. |
| MACCS keys | `generate_maccs_fingerprint` | Durant, J. L.; Leland, B. A.; Henry, D. R.; Hounshell, J. D. *J. Chem. Inf. Comput. Sci.* **2002**, *42*, 1273–1280. |
| Tanimoto similarity | `calculate_similarity` | Tanimoto, T. T. *IBM Technical Report Series* **1957** (similarity coefficient); widely applied to molecular fingerprints. |
