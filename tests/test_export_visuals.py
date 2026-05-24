"""Tests for CSV peptide export (peptide_visual.export_peptides_from_csv)."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from peptide_visual import export_peptides_from_csv


@pytest.fixture
def tiny_csv(tmp_path: Path) -> Path:
    path = tmp_path / "peptides.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["id", "sequence", "smiles"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "id": "gg",
                "sequence": "GG",
                "smiles": "NCC(=O)NCC(=O)O",
            }
        )
        writer.writerow(
            {
                "id": "bad",
                "sequence": "",
                "smiles": "not_smiles",
            }
        )
    return path


def test_export_peptides_from_csv_png_and_summary(tiny_csv: Path, tmp_path: Path):
    out = tmp_path / "export_out"
    report = export_peptides_from_csv(
        tiny_csv,
        out,
        skip_3d=True,
        write_grid=False,
        write_descriptor_summary=True,
    )

    assert report["summary"]["total"] == 2
    assert (out / "gg.png").is_file()
    assert (out / "export_summary.csv").is_file()
    assert (out / "descriptors_summary.csv").is_file()

    gg_row = next(r for r in report["results"] if r["id"] == "gg")
    assert gg_row["png"] is True
    assert gg_row["smiles"] == "NCC(=O)NCC(=O)O"

    bad_row = next(r for r in report["results"] if r["id"] == "bad")
    assert bad_row["png"] is False
    assert bad_row["errors"]


def test_export_builds_smiles_from_sequence_only(tmp_path: Path):
    path = tmp_path / "seq_only.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "sequence"])
        writer.writeheader()
        writer.writerow({"name": "dipeptide", "sequence": "GG"})

    out = tmp_path / "out"
    report = export_peptides_from_csv(path, out, skip_3d=True, write_descriptor_summary=False)
    assert report["results"][0]["smiles"] == "NCC(=O)NCC(=O)O"
    assert (out / "dipeptide.png").is_file()
