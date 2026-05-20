"""Shared RDKit helpers for SMILES validation and canonicalization."""

from __future__ import annotations

from typing import Optional

from rdkit import Chem


def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    """Parse SMILES into an RDKit molecule, or ``None`` if invalid."""
    if not smiles or not isinstance(smiles, str):
        return None
    return Chem.MolFromSmiles(smiles.strip())


def validate_smiles(smiles: str) -> bool:
    """Return True if RDKit can parse ``smiles``."""
    return smiles_to_mol(smiles) is not None


def canonicalize_smiles(smiles: str) -> Optional[str]:
    """Return canonical SMILES, or ``None`` if parsing fails."""
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)
