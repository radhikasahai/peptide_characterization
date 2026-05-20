"""Streamlit UI for peptide sequence characterization (portfolio / demo site)."""

from __future__ import annotations

import io
from typing import Any, Dict

import streamlit as st
from rdkit.Chem.Draw import MolToImage

from peptide_parser import (
    classify_sequence,
    sequence_to_rdkit_mol,
    sequence_to_smiles,
    validate_sequence,
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
    st.markdown(
        "Built with **RDKit**: validation and composition from sequence, "
        "linear peptide assembly via **amide coupling**, optional Asp/Glu "
        "**side-chain methyl protection** during assembly (documented in the README)."
    )

    seq = st.text_input(
        "One-letter amino acid sequence",
        value="ACYDEKGP",
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

    classification = classify_sequence(seq)
    smiles = sequence_to_smiles(seq, hydrolyze_sidechains=hydrolyze)
    mol = sequence_to_rdkit_mol(seq, hydrolyze_sidechains=hydrolyze)

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
