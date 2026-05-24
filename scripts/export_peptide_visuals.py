#!/usr/bin/env python3
"""
Export 2D/3D assets for peptides listed in a CSV.

Example:
  python scripts/export_peptide_visuals.py \\
    --csv data/peptides_example.csv \\
    --out output/peptides
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from peptide_visual import export_peptides_from_csv  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export PNG and SDF files for peptides in a CSV."
    )
    parser.add_argument(
        "--csv",
        required=True,
        type=Path,
        help="Input CSV (columns: id or name; smiles and/or sequence)",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output directory for PNG/SDF files",
    )
    parser.add_argument(
        "--no-hydrolyze",
        action="store_true",
        help="When building from sequence, keep Asp/Glu methyl esters",
    )
    parser.add_argument(
        "--skip-3d",
        action="store_true",
        help="Only write PNG files (skip _3d.sdf and _ligand.sdf)",
    )
    parser.add_argument(
        "--grid",
        action="store_true",
        help="Also write library_grid.png for all successful PNGs",
    )
    parser.add_argument(
        "--no-descriptors",
        action="store_true",
        help="Skip descriptors_summary.csv",
    )

    args = parser.parse_args()

    report = export_peptides_from_csv(
        args.csv,
        args.out,
        hydrolyze_sidechains=not args.no_hydrolyze,
        skip_3d=args.skip_3d,
        write_grid=args.grid,
        write_descriptor_summary=not args.no_descriptors,
    )

    print(json.dumps(report["summary"], indent=2))
    for row in report["results"]:
        status = "ok" if not row["errors"] else "partial" if row["png"] else "fail"
        print(f"  [{status}] {row['id']}: {', '.join(row['errors']) or 'exported'}")


if __name__ == "__main__":
    main()
