#!/usr/bin/env python3
"""
Regenerate data/benchmark_sequences.csv from the manifest below.

Run from repo root (requires RDKit):
  python scripts/generate_benchmark.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from peptide_parser import (  # noqa: E402
    classify_sequence,
    sequence_to_smiles,
    validate_sequence,
)

OUTPUT = REPO_ROOT / "data" / "benchmark_sequences.csv"

# (id, sequence, expect_valid, hydrolyze_sidechains, tags, notes)
MANIFEST: list[tuple[str, str, bool, bool, str, str]] = [
    # --- validation: invalid ---
    ("empty", "", False, True, "validation,invalid", "empty input"),
    ("whitespace_only", "   ", False, True, "validation,invalid", "whitespace only"),
    ("invalid_B", "AB", False, True, "validation,invalid", "unsupported B"),
    ("invalid_X", "AX", False, True, "validation,invalid", "unsupported X"),
    ("invalid_Z", "AZ", False, True, "validation,invalid", "unsupported Z"),
    ("invalid_U", "AU", False, True, "validation,invalid", "unsupported U"),
    ("invalid_mixed", "ACXB", False, True, "validation,invalid", "mixed valid and invalid"),
    ("invalid_only_xxx", "XXX", False, True, "validation,invalid", "all invalid"),
    # --- validation: normalize ---
    ("normalize_lowercase", "  acy  ", True, True, "validation,normalize", "strip and uppercase"),
    # --- gly ---
    ("gly_single", "G", True, True, "gly", "single glycine"),
    ("gly_dipeptide", "GG", True, True, "gly", "gly-gly"),
    ("gly_tri", "GGG", True, True, "gly", "gly trimer"),
    # --- pro ---
    ("pro_single", "P", True, True, "pro", "single proline"),
    ("ala_pro", "AP", True, True, "pro", "ala-pro"),
    ("pro_gly", "PG", True, True, "pro", "pro-gly"),
    ("gly_pro", "GP", True, True, "pro", "gly-pro"),
    ("pro_pro", "PP", True, True, "pro", "pro-pro"),
    # --- asp / glu ---
    ("asp_glu_de", "DE", True, True, "charged,asp_glu", "asp-glu dipeptide hydrolyzed"),
    (
        "asp_glu_de_protected",
        "DE",
        True,
        False,
        "charged,asp_glu,protected",
        "asp-glu with methyl esters retained",
    ),
    ("asp_ala_da", "DA", True, True, "charged,asp_glu", "asp-ala"),
    ("ala_asp_ad", "AD", True, True, "charged,asp_glu", "ala-asp"),
    ("asp_asp_dd", "DD", True, True, "charged,asp_glu", "asp-asp"),
    ("glu_glu_ee", "EE", True, True, "charged,asp_glu", "glu-glu"),
    ("glu_lys_ek", "EK", True, True, "charged,asp_glu", "glu-lys"),
    ("lys_glu_ke", "KE", True, True, "charged,asp_glu", "lys-glu"),
    ("asp_in_motif", "RGD", True, True, "charged,asp_glu,motif", "RGD integrin motif"),
    # --- charged / aromatic ---
    ("lys_tetramer", "KKKK", True, True, "charged", "poly-lys"),
    ("arg_dipeptide", "RR", True, True, "charged", "arg-arg"),
    ("aromatic_tri", "FWY", True, True, "aromatic", "phe-trp-tyr"),
    ("charged_cluster", "HKRDE", True, True, "charged", "mixed charged"),
    # --- realistic motifs ---
    ("demo_acydekgp", "ACYDEKGP", True, True, "motif", "streamlit default demo"),
    ("all_twenty_once", "ARNDCQEGHILKMFPSTWYV", True, True, "motif", "each standard AA once"),
    ("ala_core", "AAAA", True, True, "hydrophobic", "poly-ala"),
    ("val_leu_ile", "VLI", True, True, "hydrophobic", "branched aliphatic"),
    ("cys_pair", "CC", True, True, "polar", "cys-cys (no disulfide)"),
    ("asn_gln", "NQ", True, True, "polar", "asn-gln"),
    ("ser_thr", "ST", True, True, "polar", "ser-thr"),
    ("met_phe", "MF", True, True, "hydrophobic,aromatic", "met-phe"),
    ("trp_tyr", "WY", True, True, "aromatic", "trp-tyr"),
    ("his_middle", "AHA", True, True, "charged", "his between ala"),
    # --- longer / stress ---
    (
        "stress_30_ala",
        "A" * 30,
        True,
        True,
        "stress",
        "30-residue poly-alanine",
    ),
    (
        "stress_50_mixed",
        "ACDEFGHIKLMNPQRSTVWY" * 2 + "ACDEFGHIKLMNPQRSTVWY"[:10],
        True,
        True,
        "stress",
        "50-residue repeating alphabet",
    ),
    # --- extra dipeptide coverage (charged neighbors) ---
    ("arg_lys_rk", "RK", True, True, "charged", "arg-lys"),
    ("asp_arg_dr", "DR", True, True, "charged,asp_glu", "asp-arg"),
    ("glu_arg_er", "ER", True, True, "charged,asp_glu", "glu-arg"),
    ("asn_asp_nd", "ND", True, True, "polar,asp_glu", "asn-asp"),
    ("gln_glu_qe", "QE", True, True, "polar,asp_glu", "gln-glu"),
]


FIELDNAMES = [
    "id",
    "sequence",
    "expect_valid",
    "expect_smiles",
    "hydrolyze_sidechains",
    "expect_length",
    "expect_hydrophobic_count",
    "expect_charged_count",
    "expect_positive_count",
    "expect_negative_count",
    "expect_polar_count",
    "expect_aromatic_count",
    "tags",
    "notes",
]


def _bool_str(value: bool) -> str:
    return "true" if value else "false"


def build_row(
    row_id: str,
    sequence: str,
    expect_valid: bool,
    hydrolyze: bool,
    tags: str,
    notes: str,
) -> dict[str, str]:
    row: dict[str, str] = {
        "id": row_id,
        "sequence": sequence,
        "expect_valid": _bool_str(expect_valid),
        "expect_smiles": "",
        "hydrolyze_sidechains": _bool_str(hydrolyze),
        "expect_length": "",
        "expect_hydrophobic_count": "",
        "expect_charged_count": "",
        "expect_positive_count": "",
        "expect_negative_count": "",
        "expect_polar_count": "",
        "expect_aromatic_count": "",
        "tags": tags,
        "notes": notes,
    }

    validation = validate_sequence(sequence)
    if expect_valid:
        if not validation["valid"]:
            raise RuntimeError(f"{row_id}: expected valid but validation failed")
        classification = classify_sequence(sequence)
        if not classification["valid"]:
            raise RuntimeError(f"{row_id}: classify_sequence failed")
        smi = sequence_to_smiles(sequence, hydrolyze_sidechains=hydrolyze)
        if smi is None:
            raise RuntimeError(f"{row_id}: sequence_to_smiles returned None")
        row["expect_smiles"] = smi
        row["expect_length"] = str(classification["length"])
        row["expect_hydrophobic_count"] = str(classification["hydrophobic_count"])
        row["expect_charged_count"] = str(classification["charged_count"])
        row["expect_positive_count"] = str(classification["positive_count"])
        row["expect_negative_count"] = str(classification["negative_count"])
        row["expect_polar_count"] = str(classification["polar_count"])
        row["expect_aromatic_count"] = str(classification["aromatic_count"])
    else:
        if validation["valid"]:
            raise RuntimeError(f"{row_id}: expected invalid but validation passed")

    return row


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        build_row(row_id, seq, valid, hydro, tags, notes)
        for row_id, seq, valid, hydro, tags, notes in MANIFEST
    ]
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    valid_n = sum(1 for r in rows if r["expect_valid"] == "true")
    invalid_n = len(rows) - valid_n
    print(f"Wrote {len(rows)} rows to {OUTPUT}")
    print(f"  valid: {valid_n}, invalid: {invalid_n}")
    print(json.dumps({"output": str(OUTPUT), "rows": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
