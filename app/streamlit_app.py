"""Streamlit UI for peptide sequence characterization (demo site)."""

from __future__ import annotations

import io
import sys
from pathlib import Path
from typing import Any, Dict

# Repo root on sys.path (Streamlit runs this file from app/, not project root).
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import streamlit as st
from rdkit.Chem.Draw import MolToImage

from benchmark_data import benchmark_summary, load_valid_examples
from peptide_parser import classify_sequence, validate_sequence
from peptide_synthesis import (
    build_peptide,
    combinatorial_library,
    from_motif,
    library_summary,
    random_sequence,
)
from peptide_visual import calculate_basic_descriptors


def _mol_png(mol, size: tuple[int, int] = (520, 420)) -> bytes:
    img = MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def main() -> None:
    st.set_page_config(
        page_title="Peptide characterization",
        page_icon="🧬",
        layout="wide",
    )
    st.title("Peptide characterization")
    summary = benchmark_summary()
    st.markdown(
        "Built with **RDKit**: validation and composition from sequence, "
        "linear peptide assembly via **amide coupling**, optional Asp/Glu "
        "**side-chain methyl protection** during assembly (documented in the README). "
        f"Regression benchmark: **{summary['valid']}** valid and "
        f"**{summary['invalid']}** invalid sequences (`data/benchmark_sequences.csv`)."
    )

    examples = load_valid_examples()
    example_by_id = {row["id"]: row for row in examples}
    example_ids = sorted(example_by_id.keys())

    if "seq_input" not in st.session_state:
        st.session_state.seq_input = "ACYDEKGP"

    picked = st.selectbox(
        "Load benchmark example",
        options=["— custom —", *example_ids],
        format_func=lambda eid: (
            "— custom —"
            if eid == "— custom —"
            else f"{eid} — {example_by_id[eid]['notes']}"
        ),
    )
    if picked != "— custom —":
        st.session_state.seq_input = example_by_id[picked]["sequence"]

    with st.expander("Generate synthetic sequence (Phase 1)"):
        gen_mode = st.radio(
            "Mode",
            ["Random", "Motif repeat", "Combinatorial"],
            horizontal=True,
        )

        if gen_mode == "Random":
            syn_length = st.number_input(
                "Length", min_value=1, max_value=100, value=8, key="syn_random_len"
            )
        elif gen_mode == "Motif repeat":
            st.text_input("Motif", value="RGD", key="syn_motif")
            st.number_input("Repeats", min_value=1, max_value=20, value=2, key="syn_repeats")
            st.text_input("N-term flank", value="A", key="syn_flank")
        else:
            st.text_input("Base sequence", value="AAAAA", key="syn_base")
            st.number_input(
                "Vary position (1-based)",
                min_value=1,
                max_value=50,
                value=3,
                key="syn_pos",
            )
            st.text_input("Allowed AAs at position", value="DE", key="syn_allowed")

        if st.button("Generate and use sequence"):
            try:
                if gen_mode == "Random":
                    st.session_state.seq_input = random_sequence(
                        int(st.session_state.syn_random_len)
                    )
                elif gen_mode == "Motif repeat":
                    st.session_state.seq_input = from_motif(
                        st.session_state.syn_motif,
                        repeats=int(st.session_state.syn_repeats),
                        n_term_flank=st.session_state.syn_flank,
                    )
                else:
                    variants = combinatorial_library(
                        st.session_state.syn_base,
                        {int(st.session_state.syn_pos): st.session_state.syn_allowed},
                        max_variants=50,
                    )
                    st.session_state.seq_input = variants[0]
                    st.session_state.syn_variants = variants
                st.rerun()
            except ValueError as exc:
                st.error(str(exc))

        if gen_mode == "Combinatorial" and st.session_state.get("syn_variants"):
            st.caption(
                f"Library: {len(st.session_state.syn_variants)} sequences "
                "(first loaded in input)."
            )
            st.json(library_summary(build_peptide_library(st.session_state.syn_variants[:10])))

    seq = st.text_input(
        "One-letter amino acid sequence",
        key="seq_input",
        help="Standard 20 L-amino acids. Spaces are stripped; sequence is uppercased.",
    ).strip()

    col_a, col_b = st.columns(2)
    with col_a:
        hydrolyze = st.checkbox(
            "Hydrolyze Asp/Glu side-chain methyl esters after assembly",
            value=True,
            help=(
                "If disabled, Asp/Glu retain methyl esters used as protecting groups "
                "during coupling."
            ),
        )

    validation: Dict[str, Any] = validate_sequence(seq)
    if not validation["valid"]:
        st.error(
            "Invalid sequence: unsupported residues "
            f"`{validation['invalid_residues']}`."
        )
        return

    built = build_peptide(seq, hydrolyze_sidechains=hydrolyze)
    if not built["valid"]:
        st.error("Build failed: " + "; ".join(built["errors"]))
        return

    classification = classify_sequence(seq)
    smiles = built["smiles"]
    mol = built["mol"]

    with col_a:
        st.subheader("Composition")
        st.json(
            {
                "length": classification["length"],
                "hydrophobic": classification["hydrophobic_count"],
                "charged": classification["charged_count"],
                "positive": classification["positive_count"],
                "negative": classification["negative_count"],
                "polar": classification["polar_count"],
                "aromatic": classification["aromatic_count"],
                "counts": classification["composition"],
            }
        )

    with col_b:
        st.subheader("Structure")
        if smiles:
            st.code(smiles, language=None)
        if mol is not None:
            st.image(_mol_png(mol), caption="2D structure (RDKit depiction)")
        else:
            st.warning("Could not build an RDKit molecule from this sequence.")

    st.subheader("Whole-molecule descriptors")
    if smiles:
        desc = calculate_basic_descriptors(smiles)
        if desc:
            # Present as a compact table-friendly dict
            st.dataframe(
                [{"descriptor": k, "value": v} for k, v in sorted(desc.items())],
                hide_index=True,
                use_container_width=True,
            )
        else:
            st.warning("Descriptor calculation failed.")


if __name__ == "__main__":
    main()
