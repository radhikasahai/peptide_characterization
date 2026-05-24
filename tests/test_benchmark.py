"""Regression tests driven by data/benchmark_sequences.csv."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from chem_utils import validate_smiles
from peptide_parser import classify_sequence, sequence_to_smiles, validate_sequence

BENCHMARK_CSV = Path(__file__).resolve().parent.parent / "data" / "benchmark_sequences.csv"


def _load_benchmark_rows() -> list[dict[str, str]]:
    with BENCHMARK_CSV.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _expect_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def _optional_int(value: str) -> int | None:
    value = (value or "").strip()
    if not value:
        return None
    return int(value)


BENCHMARK_ROWS = _load_benchmark_rows()


@pytest.mark.parametrize("row", BENCHMARK_ROWS, ids=lambda r: r["id"])
def test_benchmark_validation(row: dict[str, str]) -> None:
    expect_valid = _expect_bool(row["expect_valid"])
    result = validate_sequence(row["sequence"])
    assert result["valid"] == expect_valid


@pytest.mark.parametrize(
    "row",
    [r for r in BENCHMARK_ROWS if _expect_bool(r["expect_valid"])],
    ids=lambda r: r["id"],
)
def test_benchmark_smiles(row: dict[str, str]) -> None:
    hydrolyze = _expect_bool(row["hydrolyze_sidechains"])
    expected = row["expect_smiles"].strip()
    assert expected, f"{row['id']}: missing expect_smiles in benchmark CSV"
    actual = sequence_to_smiles(row["sequence"], hydrolyze_sidechains=hydrolyze)
    assert actual == expected
    assert validate_smiles(actual)


@pytest.mark.parametrize(
    "row",
    [r for r in BENCHMARK_ROWS if _expect_bool(r["expect_valid"])],
    ids=lambda r: r["id"],
)
def test_benchmark_composition(row: dict[str, str]) -> None:
    classification = classify_sequence(row["sequence"])
    assert classification["valid"] is True

    for field, key in (
        ("expect_length", "length"),
        ("expect_hydrophobic_count", "hydrophobic_count"),
        ("expect_charged_count", "charged_count"),
        ("expect_positive_count", "positive_count"),
        ("expect_negative_count", "negative_count"),
        ("expect_polar_count", "polar_count"),
        ("expect_aromatic_count", "aromatic_count"),
    ):
        expected = _optional_int(row[field])
        assert expected is not None, f"{row['id']}: missing {field}"
        assert classification[key] == expected


def test_benchmark_file_has_expected_scale() -> None:
    """Guard against accidentally shrinking the benchmark."""
    assert len(BENCHMARK_ROWS) >= 40
    invalid = sum(1 for r in BENCHMARK_ROWS if not _expect_bool(r["expect_valid"]))
    valid = len(BENCHMARK_ROWS) - invalid
    assert invalid >= 8
    assert valid >= 30
