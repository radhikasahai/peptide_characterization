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

Regenerate after changing assembly logic (review diffs before committing):

```bash
python scripts/generate_benchmark.py
pytest tests/
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
