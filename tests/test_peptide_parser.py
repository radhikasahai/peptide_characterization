"""Tests for peptide_parser (requires RDKit)."""

from __future__ import annotations

import pytest

pytest.importorskip("rdkit")

from peptide_parser import (
    classify_sequence,
    sequence_to_rdkit_mol,
    sequence_to_smiles,
    validate_sequence,
)


def test_validate_sequence_invalid_residue():
    out = validate_sequence("AZX")
    assert out["valid"] is False
    assert set(out["invalid_residues"]) == {"X", "Z"}


def test_validate_sequence_normalizes():
    out = validate_sequence("  acy  ")
    assert out["valid"] is True
    assert out["sequence"] == "ACY"


def test_validate_sequence_whitespace_only_invalid():
    out = validate_sequence("   ")
    assert out["valid"] is False
    assert out["sequence"] == ""


def test_classify_sequence_invalid():
    assert classify_sequence("B")["valid"] is False


def test_sequence_to_smiles_glycine_dipeptide():
    assert sequence_to_smiles("GG") == "NCC(=O)NCC(=O)O"


def test_sequence_to_smiles_retains_methyl_esters_when_requested():
    assert (
        sequence_to_smiles("DE", hydrolyze_sidechains=False)
        == "COC(=O)CC[C@H](NC(=O)[C@@H](N)CC(=O)OC)C(=O)O"
    )


def test_sequence_to_smiles_hydrolyzes_sidechains():
    assert (
        sequence_to_smiles("DE", hydrolyze_sidechains=True)
        == "N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)O"
    )


def test_sequence_to_rdkit_mol_roundtrip():
    mol = sequence_to_rdkit_mol("ACYDEKGP")
    assert mol is not None
    smi = sequence_to_smiles("ACYDEKGP")
    assert smi is not None
