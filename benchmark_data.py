"""Load the committed peptide benchmark CSV."""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
from typing import Any

BENCHMARK_PATH = Path(__file__).resolve().parent / "data" / "benchmark_sequences.csv"


@lru_cache(maxsize=1)
def load_benchmark_rows() -> tuple[dict[str, str], ...]:
    with BENCHMARK_PATH.open(encoding="utf-8") as handle:
        return tuple(csv.DictReader(handle))


def load_valid_examples() -> list[dict[str, str]]:
    """Valid benchmark rows (for demos and documentation)."""
    return [
        row
        for row in load_benchmark_rows()
        if row.get("expect_valid", "").strip().lower() == "true"
    ]


def benchmark_summary() -> dict[str, Any]:
    rows = load_benchmark_rows()
    valid = sum(1 for r in rows if r.get("expect_valid", "").lower() == "true")
    return {
        "total": len(rows),
        "valid": valid,
        "invalid": len(rows) - valid,
    }
