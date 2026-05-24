"""Tests for peptide_synthesis (sequence generation + build wrappers)."""

from __future__ import annotations

import pytest

from peptide_parser import validate_sequence
from peptide_synthesis import (
    build_peptide,
    build_peptide_library,
    combinatorial_library,
    from_motif,
    library_summary,
    random_sequence,
    sequences_from_benchmark,
)

pytest.importorskip("rdkit")


def test_random_sequence_length_and_alphabet():
    rng = __import__("random").Random(42)
    seq = random_sequence(12, rng=rng)
    assert len(seq) == 12
    assert validate_sequence(seq)["valid"]


def test_random_sequence_weighted():
    rng = __import__("random").Random(0)
    seq = random_sequence(20, weights={"A": 1.0}, rng=rng)
    assert seq == "A" * 20


def test_random_sequence_invalid_length():
    with pytest.raises(ValueError, match="length"):
        random_sequence(0)


def test_from_motif_with_flanks():
    assert from_motif("RGD", repeats=2, n_term_flank="A") == "ARGDRGD"


def test_from_motif_invalid_residue():
    with pytest.raises(ValueError, match="unsupported"):
        from_motif("AX")


def test_combinatorial_library_size():
    variants = combinatorial_library("AAAAA", {3: "DE", 5: "KR"})
    assert len(variants) == 4
    assert set(v[2] for v in variants) == {"D", "E"}
    assert set(v[4] for v in variants) == {"K", "R"}


def test_combinatorial_library_max_variants_guard():
    with pytest.raises(ValueError, match="max_variants"):
        combinatorial_library(
            "A" * 5,
            {i: "ACDEFGHIKLMNPQRSTWYV" for i in range(1, 6)},
            max_variants=10,
        )


def test_build_peptide_valid():
    result = build_peptide("GG", name="dipeptide")
    assert result["valid"] is True
    assert result["name"] == "dipeptide"
    assert result["smiles"] == "NCC(=O)NCC(=O)O"
    assert result["mol"] is not None
    assert result["errors"] == []


def test_build_peptide_invalid():
    result = build_peptide("AX")
    assert result["valid"] is False
    assert result["smiles"] is None
    assert result["mol"] is None
    assert result["errors"]


def test_build_peptide_library():
    lib = build_peptide_library(["GG", "AX"], names=["ok", "bad"])
    assert lib[0]["valid"]
    assert lib[1]["valid"] is False
    summary = library_summary(lib)
    assert summary["total"] == 2
    assert summary["built"] == 1
    assert summary["failed"] == 1


def test_sequences_from_benchmark_tag():
    gly = sequences_from_benchmark(tag="gly")
    assert len(gly) >= 2
    assert all("gly" in row["tags"].lower() for row in gly)


def test_sequences_from_benchmark_limit():
    rows = sequences_from_benchmark(limit=2)
    assert len(rows) == 2
