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

## Synthetic peptides (Phase 1)

Generate sequences and build structures with `peptide_synthesis.py`:

```python
from peptide_synthesis import (
    random_sequence,
    from_motif,
    combinatorial_library,
    build_peptide,
    build_peptide_library,
)

seq = random_sequence(10)  # optional: weights={"K": 2.0, "R": 2.0}
motif = from_motif("RGD", repeats=2, n_term_flank="A")
variants = combinatorial_library("AAAAA", {3: "DE"})  # 1-based positions

result = build_peptide(motif)  # -> valid, sequence, smiles, mol, errors
library = build_peptide_library(variants)
```

```bash
python peptide_synthesis.py
```

## Export peptide visuals from CSV

Batch-export the same assets as the aspirin demo (2D PNG, 3D SDF, docking SDF) for
each row in a peptide library CSV.

### Input CSV format

Header row required. Column names are case-insensitive.

| Column | Required | Description |
|--------|----------|-------------|
| `id` or `name` | Yes | Base filename for outputs (e.g. `gly_dipeptide` → `gly_dipeptide.png`) |
| `smiles` | Recommended | Valid SMILES; used as-is when present |
| `sequence` | Optional | One-letter sequence; used when `smiles` is empty (built via `sequence_to_smiles`) |

Example (`data/peptides_example.csv`):

```csv
id,sequence,smiles,notes
gly_dipeptide,GG,NCC(=O)NCC(=O)O,short test peptide
motif_rgd,RGD,,build SMILES from sequence
```

### Outputs (per row in `--out` directory)

| File | Description |
|------|-------------|
| `{id}.png` | 2D structure image |
| `{id}_3d.sdf` | 3D conformer (UFF), unless `--skip-3d` |
| `{id}_ligand.sdf` | Docking-style 3D (MMFF), unless `--skip-3d` |
| `export_summary.csv` | Per-row success/failure log |
| `descriptors_summary.csv` | RDKit descriptors for all rows with valid SMILES |
| `library_grid.png` | Optional grid (`--grid`) |

### Commands

```bash
python scripts/export_peptide_visuals.py --csv data/peptides_example.csv --out output/peptides
python scripts/export_peptide_visuals.py --csv data/peptides_example.csv --out output/peptides --skip-3d --grid
```

From Python:

```python
from peptide_visual import export_peptides_from_csv

report = export_peptides_from_csv("data/peptides_example.csv", "output/peptides")
print(report["summary"])
```

Build a CSV from a synthesized library:

```python
import csv
from peptide_synthesis import build_peptide_library, combinatorial_library

variants = combinatorial_library("AAAAA", {3: "DE"})
lib = build_peptide_library(variants, names=[f"v{i}" for i in range(len(variants))])
with open("data/my_library.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["id", "sequence", "smiles", "valid"])
    w.writeheader()
    for entry in lib:
        w.writerow({
            "id": entry.get("name") or entry["sequence"],
            "sequence": entry["sequence"],
            "smiles": entry.get("smiles") or "",
            "valid": entry["valid"],
        })
```

## Command-line demos

```bash
python peptide_parser.py
python peptide_visual.py
```

## Web UI

Interactive app: **GitHub Pages frontend** (`docs/`) + **FastAPI backend** (`api/`).

Live site: **https://radhikasahai.github.io/peptide_characterization/**

### Local development

Terminal 1 — API (from repo root):

```bash
conda activate peptide
pip install -r requirements-api-local.txt
uvicorn api.main:app --reload --port 8000
```

Do **not** run `pip install -r requirements-api.txt` locally if you already have conda RDKit — it installs a conflicting `rdkit-pypi` copy.

Terminal 2 — frontend:

```bash
cd docs
python -m http.server 5500
```

Open **http://localhost:5500**.

### Deploy to GitHub Pages + Render

**Step 1 — Enable Pages (one-time)**

1. Repo → **Settings** → **Pages**
2. **Build and deployment → Source:** **GitHub Actions**

**Step 2 — Deploy the API on Render**

1. [Render](https://render.com) → **New** → **Blueprint** (or Web Service from this repo)
2. Use the repo `render.yaml` / `Dockerfile`
3. Confirm env var `ALLOWED_ORIGINS=https://radhikasahai.github.io`
4. Wait for deploy; copy the service URL (e.g. `https://peptide-characterization-api.onrender.com`)

**Step 3 — Connect frontend to API**

Either:

- **Quick (no redeploy):** open the live site → **Connect API backend** panel → paste your Render URL → **Save & connect**
- **Permanent:** repo → **Settings → Secrets and variables → Actions → Variables** → add **`PEPTIDE_API_URL`** → re-run **Deploy GitHub Pages**

**Step 4 — Verify**

- Pages site shows the dark interactive UI (not README)
- Enter `GG` → 2D structure and descriptors appear

**Fallback (if Actions deploy fails):** Settings → Pages → deploy from branch **`main`** / folder **`/docs`**.

### Streamlit (legacy demo)

Run from the **repository root** (so imports resolve):

```bash
cd peptide_characterization
conda activate peptide
streamlit run app/streamlit_app.py
```

## Benchmark dataset

`data/benchmark_sequences.csv` is a **frozen regression suite** (47 rows: 39 valid assemblies, 8 invalid inputs). Each valid row stores golden **canonical SMILES** and **composition counts** generated once by RDKit.

Regenerate after changing assembly logic :

```bash
python scripts/generate_benchmark.py
pytest tests/
```

## Tests & CI

```bash
pytest tests/
```

GitHub Actions runs the same tests using Micromamba (`.github/workflows/ci.yml`).

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



## Layout

| Module | Role |
|--------|------|
| `chem_utils.py` | Minimal SMILES helpers shared across scripts |
| `peptide_parser.py` | Sequence validation, composition, **RDKit peptide assembly** |
| `peptide_synthesis.py` | **Generate** sequences (random/motif/library) and **build** structures |
| `peptide_visual.py` | Descriptors, fingerprints, 2D/3D export helpers |
| `api/main.py` | FastAPI backend for the GitHub Pages web UI |
| `docs/` | Static frontend (GitHub Pages) |
| `app/streamlit_app.py` | Legacy Streamlit demo UI |
| `data/benchmark_sequences.csv` | Committed regression benchmark (golden SMILES + composition) |
| `scripts/generate_benchmark.py` | Regenerate benchmark CSV from manifest |
| `scripts/export_peptide_visuals.py` | Batch PNG/SDF export from peptide CSV |
| `data/peptides_example.csv` | Sample input for visual export |
| `benchmark_data.py` | Load benchmark rows (Streamlit examples, summaries) |
