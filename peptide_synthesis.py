"""
Synthetic linear peptide generation and structure building.

Phase 1 pipeline: generate sequence(s) -> build_peptide / build_peptide_library
(uses ``peptide_parser`` for RDKit assembly).
"""

from __future__ import annotations

import itertools
import random
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union

from rdkit import Chem

from benchmark_data import load_valid_examples
from peptide_parser import (
    STANDARD_AAS,
    list_standard_aas,
    sequence_to_rdkit_mol,
    sequence_to_smiles,
    validate_sequence,
)

# Re-export for callers configuring generators
DEFAULT_ALPHABET: tuple[str, ...] = tuple(sorted(STANDARD_AAS))


def _normalize_alphabet(alphabet: Optional[Sequence[str]]) -> List[str]:
    letters = [aa.upper() for aa in (alphabet or list_standard_aas())]
    invalid = sorted({aa for aa in letters if aa not in STANDARD_AAS})
    if invalid:
        raise ValueError(f"Unknown amino acid codes in alphabet: {invalid}")
    if not letters:
        raise ValueError("Alphabet must not be empty")
    return letters


def _validate_motif(motif: str) -> str:
    motif = motif.upper().strip()
    invalid = sorted({aa for aa in motif if aa not in STANDARD_AAS})
    if invalid:
        raise ValueError(f"Motif contains unsupported residues: {invalid}")
    if not motif:
        raise ValueError("Motif must not be empty")
    return motif


def random_sequence(
    length: int,
    *,
    alphabet: Optional[Sequence[str]] = None,
    weights: Optional[Mapping[str, float]] = None,
    rng: Optional[random.Random] = None,
) -> str:
    """
    Sample a random one-letter sequence of the given length.

    ``weights`` maps amino-acid code -> non-negative weight (need not sum to 1).
    """
    if length < 1:
        raise ValueError("length must be >= 1")

    letters = _normalize_alphabet(alphabet)
    rng = rng or random.Random()

    if weights:
        pool: List[str] = []
        weight_values: List[float] = []
        for aa in letters:
            w = weights.get(aa, 0.0)
            if w > 0:
                pool.append(aa)
                weight_values.append(float(w))
        if not pool:
            raise ValueError("weights must assign positive weight to at least one letter")
        return "".join(rng.choices(pool, weights=weight_values, k=length))

    return "".join(rng.choice(letters) for _ in range(length))


def from_motif(
    motif: str,
    *,
    repeats: int = 1,
    n_term_flank: str = "",
    c_term_flank: str = "",
) -> str:
    """
    Build a sequence from a repeated motif and optional N/C-terminal flanks.

    Example: ``from_motif("RGD", repeats=2, n_term_flank="A")`` -> ``ARGDRGD``.
    """
    if repeats < 1:
        raise ValueError("repeats must be >= 1")

    motif = _validate_motif(motif)
    n_term_flank = _validate_motif(n_term_flank) if n_term_flank else ""
    c_term_flank = _validate_motif(c_term_flank) if c_term_flank else ""

    return n_term_flank + motif * repeats + c_term_flank


def combinatorial_library(
    base: str,
    positions: Mapping[int, str],
    *,
    max_variants: int = 100,
) -> List[str]:
    """
    Enumerate sequences by varying ``base`` at selected **1-based** positions.

    ``positions`` maps position -> allowed amino acids, e.g. ``{3: "DE", 5: "KR"}``.
    Stops with ``ValueError`` if the Cartesian product exceeds ``max_variants``.
    """
    validation = validate_sequence(base)
    if not validation["valid"]:
        raise ValueError(
            f"Invalid base sequence: {validation['invalid_residues']}"
        )

    base_seq = validation["sequence"]
    n = len(base_seq)

    if not positions:
        return [base_seq]

    sorted_positions = sorted(positions.keys())
    for pos in sorted_positions:
        if pos < 1 or pos > n:
            raise ValueError(f"Position {pos} out of range for base length {n}")
        allowed = _validate_motif(positions[pos])
        if not allowed:
            raise ValueError(f"No allowed amino acids at position {pos}")

    option_lists = [_validate_motif(positions[p]) for p in sorted_positions]
    total = 1
    for opts in option_lists:
        total *= len(opts)
    if total > max_variants:
        raise ValueError(
            f"Library would have {total} variants (max_variants={max_variants})"
        )

    sequences: List[str] = []
    chars = list(base_seq)
    for combo in itertools.product(*option_lists):
        variant = chars.copy()
        for pos, aa in zip(sorted_positions, combo):
            variant[pos - 1] = aa
        sequences.append("".join(variant))

    return sequences


def sequences_from_benchmark(
    *,
    tag: Optional[str] = None,
    limit: Optional[int] = None,
) -> List[Dict[str, str]]:
    """
    Return valid benchmark entries as ``{id, sequence, notes, tags}``.

    If ``tag`` is set, only rows whose ``tags`` column contains that tag
    (comma-separated) are included.
    """
    rows = load_valid_examples()
    if tag:
        tag = tag.strip().lower()
        rows = [
            r
            for r in rows
            if tag in {t.strip().lower() for t in r.get("tags", "").split(",")}
        ]
    if limit is not None:
        rows = rows[:limit]
    return [
        {
            "id": r["id"],
            "sequence": r["sequence"],
            "notes": r.get("notes", ""),
            "tags": r.get("tags", ""),
        }
        for r in rows
    ]


def build_peptide(
    sequence: str,
    *,
    hydrolyze_sidechains: bool = True,
    name: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Validate a sequence and build a linear peptide structure in RDKit.

    Returns a dict with ``valid``, ``sequence``, ``smiles``, ``mol``, and ``errors``.
    """
    validation = validate_sequence(sequence)
    result: Dict[str, Any] = {
        "name": name,
        "sequence": validation["sequence"],
        "valid": False,
        "invalid_residues": validation["invalid_residues"],
        "smiles": None,
        "mol": None,
        "errors": [],
        "hydrolyze_sidechains": hydrolyze_sidechains,
    }

    if not validation["valid"]:
        result["errors"].append(
            f"Invalid residues: {validation['invalid_residues']}"
        )
        return result

    mol = sequence_to_rdkit_mol(
        validation["sequence"],
        hydrolyze_sidechains=hydrolyze_sidechains,
    )
    if mol is None:
        result["errors"].append("RDKit assembly failed")
        return result

    smiles = sequence_to_smiles(
        validation["sequence"],
        hydrolyze_sidechains=hydrolyze_sidechains,
    )
    if smiles is None:
        result["errors"].append("SMILES generation failed")
        return result

    result["smiles"] = smiles
    result["mol"] = mol
    result["valid"] = True
    return result


def build_peptide_library(
    sequences: Sequence[str],
    *,
    names: Optional[Sequence[Optional[str]]] = None,
    hydrolyze_sidechains: bool = True,
) -> List[Dict[str, Any]]:
    """
    Build many peptides. ``names`` is optional and aligned with ``sequences``.
    """
    seq_list = list(sequences)
    if names is not None:
        name_list = list(names)
        if len(name_list) != len(seq_list):
            raise ValueError("names must have the same length as sequences")
    else:
        name_list = [None] * len(seq_list)

    return [
        build_peptide(
            seq,
            hydrolyze_sidechains=hydrolyze_sidechains,
            name=name,
        )
        for seq, name in zip(seq_list, name_list)
    ]


def library_summary(library: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """Aggregate success/failure counts for a build_peptide_library result."""
    built = [entry for entry in library if entry.get("valid")]
    failed = [entry for entry in library if not entry.get("valid")]
    return {
        "total": len(library),
        "built": len(built),
        "failed": len(failed),
        "sequences": [e["sequence"] for e in built],
    }


if __name__ == "__main__":
    rng = random.Random(0)

    print("=" * 60)
    print("Random sequence")
    print("=" * 60)
    seq = random_sequence(8, rng=rng)
    print(seq)
    print(build_peptide(seq, name="random_8"))

    print("\n" + "=" * 60)
    print("Motif (RGD x2 + N-term Ala)")
    print("=" * 60)
    motif_seq = from_motif("RGD", repeats=2, n_term_flank="A")
    print(motif_seq)
    print(build_peptide(motif_seq, name="motif_rgd"))

    print("\n" + "=" * 60)
    print("Combinatorial library (base=AAAAA, pos3=DE)")
    print("=" * 60)
    variants = combinatorial_library("AAAAA", {3: "DE"})
    lib = build_peptide_library(variants, names=[f"v{i}" for i in range(len(variants))])
    print(library_summary(lib))

    print("\n" + "=" * 60)
    print("Benchmark tag=gly (first 3)")
    print("=" * 60)
    for row in sequences_from_benchmark(tag="gly", limit=3):
        built = build_peptide(row["sequence"], name=row["id"])
        print(row["id"], built["valid"], built.get("smiles", "")[:40], "...")
