"""FastAPI backend for the peptide characterization web UI."""

from __future__ import annotations

import base64
import io
import os
import sys
from pathlib import Path
from typing import Any, Literal, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from rdkit.Chem.Draw import MolToImage

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from benchmark_data import benchmark_summary, load_valid_examples
from peptide_parser import classify_sequence, validate_sequence
from peptide_synthesis import (
    build_peptide,
    build_peptide_library,
    combinatorial_library,
    from_motif,
    library_summary,
    random_sequence,
)
from peptide_visual import calculate_basic_descriptors

_DEFAULT_ORIGINS = [
    "http://localhost:8000",
    "http://127.0.0.1:8000",
    "http://localhost:5500",
    "http://127.0.0.1:5500",
    "https://radhikasahai.github.io",
]


def _allowed_origins() -> list[str]:
    extra = os.environ.get("ALLOWED_ORIGINS", "")
    origins = list(_DEFAULT_ORIGINS)
    for item in extra.split(","):
        item = item.strip()
        if item and item not in origins:
            origins.append(item)
    return origins


app = FastAPI(title="Peptide characterization API", version="1.0.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=_allowed_origins(),
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class CharacterizeRequest(BaseModel):
    sequence: str = ""
    hydrolyze_sidechains: bool = True


class GenerateRequest(BaseModel):
    mode: Literal["random", "motif", "combinatorial"]
    length: int = Field(default=8, ge=1, le=100)
    motif: str = "RGD"
    repeats: int = Field(default=2, ge=1, le=20)
    n_term_flank: str = "A"
    base: str = "AAAAA"
    position: int = Field(default=3, ge=1, le=50)
    allowed: str = "DE"
    max_variants: int = Field(default=50, ge=1, le=500)
    include_library_summary: bool = False
    hydrolyze_sidechains: bool = True


def _mol_png_base64(mol, size: tuple[int, int] = (520, 420)) -> str:
    img = MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _composition_payload(classification: dict[str, Any]) -> dict[str, Any]:
    return {
        "length": classification["length"],
        "hydrophobic": classification["hydrophobic_count"],
        "charged": classification["charged_count"],
        "positive": classification["positive_count"],
        "negative": classification["negative_count"],
        "polar": classification["polar_count"],
        "aromatic": classification["aromatic_count"],
        "counts": classification["composition"],
    }


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/api/benchmark")
def get_benchmark() -> dict[str, Any]:
    summary = benchmark_summary()
    examples = [
        {"id": row["id"], "sequence": row["sequence"], "notes": row.get("notes", "")}
        for row in load_valid_examples()
    ]
    examples.sort(key=lambda row: row["id"])
    return {"summary": summary, "examples": examples}


@app.post("/api/characterize")
def characterize(body: CharacterizeRequest) -> dict[str, Any]:
    sequence = body.sequence.strip()
    validation = validate_sequence(sequence)
    if not validation["valid"]:
        raise HTTPException(
            status_code=400,
            detail={
                "message": "Invalid sequence",
                "invalid_residues": validation["invalid_residues"],
            },
        )

    built = build_peptide(sequence, hydrolyze_sidechains=body.hydrolyze_sidechains)
    if not built["valid"]:
        raise HTTPException(
            status_code=400,
            detail={"message": "Build failed", "errors": built["errors"]},
        )

    classification = classify_sequence(sequence)
    smiles = built["smiles"]
    mol = built["mol"]
    descriptors = calculate_basic_descriptors(smiles) if smiles else None

    return {
        "valid": True,
        "sequence": validation["sequence"],
        "hydrolyze_sidechains": body.hydrolyze_sidechains,
        "composition": _composition_payload(classification),
        "smiles": smiles,
        "structure_png_base64": _mol_png_base64(mol) if mol is not None else None,
        "descriptors": descriptors,
    }


@app.post("/api/generate")
def generate(body: GenerateRequest) -> dict[str, Any]:
    try:
        if body.mode == "random":
            sequence = random_sequence(body.length)
            variants: Optional[list[str]] = None
        elif body.mode == "motif":
            sequence = from_motif(
                body.motif,
                repeats=body.repeats,
                n_term_flank=body.n_term_flank,
            )
            variants = None
        else:
            variants = combinatorial_library(
                body.base,
                {body.position: body.allowed},
                max_variants=body.max_variants,
            )
            sequence = variants[0]
    except ValueError as exc:
        raise HTTPException(status_code=400, detail={"message": str(exc)}) from exc

    result: dict[str, Any] = {"sequence": sequence, "variants": variants}
    if body.mode == "combinatorial" and body.include_library_summary and variants:
        preview = variants[:10]
        result["library_summary"] = library_summary(
            build_peptide_library(preview, hydrolyze_sidechains=body.hydrolyze_sidechains)
        )
    return result
