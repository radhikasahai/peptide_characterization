# peptide_parser.py

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from chem_utils import canonicalize_smiles, validate_smiles


# ============================================================
# Amino Acid Data
# ============================================================

AMINO_ACID_DATA = {
    "A": {
        "name": "Alanine",
        "charge": "neutral",
        "molar_mass": 89.09,
        "property": "hydrophobic",
    },
    "R": {
        "name": "Arginine",
        "charge": "positive",
        "molar_mass": 174.20,
        "property": "charged",
    },
    "N": {
        "name": "Asparagine",
        "charge": "neutral",
        "molar_mass": 132.12,
        "property": "polar",
    },
    "D": {
        "name": "Aspartic Acid",
        "charge": "negative",
        "molar_mass": 133.10,
        "property": "charged",
    },
    "C": {
        "name": "Cysteine",
        "charge": "neutral",
        "molar_mass": 121.15,
        "property": "polar",
    },
    "E": {
        "name": "Glutamic Acid",
        "charge": "negative",
        "molar_mass": 147.13,
        "property": "charged",
    },
    "Q": {
        "name": "Glutamine",
        "charge": "neutral",
        "molar_mass": 146.15,
        "property": "polar",
    },
    "G": {
        "name": "Glycine",
        "charge": "neutral",
        "molar_mass": 75.07,
        "property": "hydrophobic",
    },
    "H": {
        "name": "Histidine",
        "charge": "positive",
        "molar_mass": 155.16,
        "property": "charged",
    },
    "I": {
        "name": "Isoleucine",
        "charge": "neutral",
        "molar_mass": 131.17,
        "property": "hydrophobic",
    },
    "L": {
        "name": "Leucine",
        "charge": "neutral",
        "molar_mass": 131.17,
        "property": "hydrophobic",
    },
    "K": {
        "name": "Lysine",
        "charge": "positive",
        "molar_mass": 146.19,
        "property": "charged",
    },
    "M": {
        "name": "Methionine",
        "charge": "neutral",
        "molar_mass": 149.21,
        "property": "hydrophobic",
    },
    "F": {
        "name": "Phenylalanine",
        "charge": "neutral",
        "molar_mass": 165.19,
        "property": "aromatic",
    },
    "P": {
        "name": "Proline",
        "charge": "neutral",
        "molar_mass": 115.13,
        "property": "hydrophobic",
    },
    "S": {
        "name": "Serine",
        "charge": "neutral",
        "molar_mass": 105.09,
        "property": "polar",
    },
    "T": {
        "name": "Threonine",
        "charge": "neutral",
        "molar_mass": 119.12,
        "property": "polar",
    },
    "W": {
        "name": "Tryptophan",
        "charge": "neutral",
        "molar_mass": 204.23,
        "property": "aromatic",
    },
    "Y": {
        "name": "Tyrosine",
        "charge": "neutral",
        "molar_mass": 181.19,
        "property": "aromatic",
    },
    "V": {
        "name": "Valine",
        "charge": "neutral",
        "molar_mass": 117.15,
        "property": "hydrophobic",
    },
}


# ============================================================
# Amino Acid Classifications
# ============================================================

STANDARD_AAS = set(AMINO_ACID_DATA.keys())

HYDROPHOBIC_AAS = {"A", "G", "I", "L", "M", "P", "V"}
CHARGED_AAS = {"R", "D", "E", "H", "K"}
POSITIVE_AAS = {"R", "H", "K"}
NEGATIVE_AAS = {"D", "E"}
POLAR_AAS = {"N", "C", "Q", "S", "T", "Y"}
AROMATIC_AAS = {"F", "W", "Y"}


# ============================================================
# Amino Acid to SMILES
# ============================================================

AA_TO_SMILES = {
    "A": "N[C@@H](C)C(=O)O",
    "R": "N[C@@H](CCCNC(N)=N)C(=O)O",
    "N": "N[C@@H](CC(=O)N)C(=O)O",
    "D": "N[C@@H](CC(=O)O)C(=O)O",
    "C": "N[C@@H](CS)C(=O)O",
    "E": "N[C@@H](CCC(=O)O)C(=O)O",
    "Q": "N[C@@H](CCC(=O)N)C(=O)O",
    "G": "NCC(=O)O",
    "H": "N[C@@H](CC1=CN=CN1)C(=O)O",
    "I": "N[C@@H](C(C)CC)C(=O)O",
    "L": "N[C@@H](CC(C)C)C(=O)O",
    "K": "N[C@@H](CCCCN)C(=O)O",
    "M": "N[C@@H](CCSC)C(=O)O",
    "F": "N[C@@H](CC1=CC=CC=C1)C(=O)O",
    "P": "N1CCCC1C(=O)O",
    "S": "N[C@@H](CO)C(=O)O",
    "T": "N[C@@H](C(O)C)C(=O)O",
    "W": "N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)O",
    "Y": "N[C@@H](CC1=CC=C(O)C=C1)C(=O)O",
    "V": "N[C@@H](C(C)C)C(=O)O",
}


# Building blocks for linear assembly (amide coupling).
# Asp/Glu use side-chain methyl esters so backbone coupling is regioselective;
# methyl groups are optionally hydrolyzed after assembly (see ``sequence_to_rdkit_mol``).
AA_BUILDING_BLOCK_SMILES = dict(AA_TO_SMILES)
AA_BUILDING_BLOCK_SMILES["D"] = "N[C@@H](CC(=O)OC)C(=O)O"
AA_BUILDING_BLOCK_SMILES["E"] = "N[C@@H](CCC(=O)OC)C(=O)O"


_RXN_AMIDE_CHIRAL = rdChemReactions.ReactionFromSmarts(
    "[C:1](=[O:2])[OH1].[N;H2:3][C@@H:4]([*:5])[C:6](=[O:7])[OH1]>>"
    "[C:1](=[O:2])[N;H1:3][C@@H:4]([*:5])[C:6](=[O:7])[OH1]"
)
_RXN_AMIDE_GLY = rdChemReactions.ReactionFromSmarts(
    "[C:1](=[O:2])[OH1].[N;H2:3][CH2:4][C:5](=[O:6])[OH1]>>"
    "[C:1](=[O:2])[N;H1:3][CH2:4][C:5](=[O:6])[OH1]"
)
_RXN_AMIDE_PRO = rdChemReactions.ReactionFromSmarts(
    "[C:1](=[O:2])[OH1].[N;H1;R:3][C:4]>>[C:1](=[O:2])[N;H0;R:3][C:4]"
)
# Methyl carbon is mapped to a methane by-product so RDKit does not warn about
# unmapped reactant atoms (the previous single-product SMARTS left [CH3] unmapped).
_RXN_ME_ESTER_HYDROLYSIS = rdChemReactions.ReactionFromSmarts(
    "[C:1](=[O:2])[O:3][CH3:4]>>[C:1](=[O:2])[OH1:3].[CH4:4]"
)


def _amide_coupling_rxn_order(incoming_aa: str):
    """Prefer reactions that uniquely describe the incoming residue."""
    if incoming_aa == "G":
        return (_RXN_AMIDE_GLY, _RXN_AMIDE_CHIRAL, _RXN_AMIDE_PRO)
    if incoming_aa == "P":
        return (_RXN_AMIDE_PRO, _RXN_AMIDE_CHIRAL, _RXN_AMIDE_GLY)
    return (_RXN_AMIDE_CHIRAL, _RXN_AMIDE_GLY, _RXN_AMIDE_PRO)


def _couple_residue(chain_mol: Chem.Mol, aa_code: str) -> Optional[Chem.Mol]:
    """Form one peptide bond: chain C-terminus + incoming N-terminus."""
    bb_smiles = AA_BUILDING_BLOCK_SMILES.get(aa_code)
    if bb_smiles is None:
        return None
    aa_mol = Chem.MolFromSmiles(bb_smiles)
    if aa_mol is None:
        return None

    for rxn in _amide_coupling_rxn_order(aa_code):
        products = rxn.RunReactants((chain_mol, aa_mol))
        smiles_seen = set()
        last_good: Optional[Chem.Mol] = None
        for prod in products:
            mol = prod[0]
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                continue
            smi = Chem.MolToSmiles(mol, canonical=True)
            if smi not in smiles_seen:
                smiles_seen.add(smi)
                last_good = mol
        if len(smiles_seen) == 1 and last_good is not None:
            return last_good

    return None


def _hydrolyze_building_block_methyl_esters(mol: Chem.Mol) -> Chem.Mol:
    """Convert Asp/Glu side-chain methyl esters used during assembly into acids."""
    current = Chem.Mol(mol)
    while True:
        products = _RXN_ME_ESTER_HYDROLYSIS.RunReactants((current,))
        if not products:
            break
        nxt = products[0][0]
        Chem.SanitizeMol(nxt)
        current = nxt
    return current


# ============================================================
# Sequence Validation
# ============================================================

def validate_sequence(sequence: str) -> Dict[str, Any]:
    """
    Validate amino acid sequence.
    Returns:
        {
            "valid": bool,
            "sequence": str,
            "invalid_residues": list
        }
    """
    if not sequence:
        return {
            "valid": False,
            "sequence": "",
            "invalid_residues": [],
        }

    sequence = sequence.upper().strip()

    if not sequence:
        return {
            "valid": False,
            "sequence": "",
            "invalid_residues": [],
        }

    invalid = sorted(set([aa for aa in sequence if aa not in STANDARD_AAS]))

    return {
        "valid": len(invalid) == 0,
        "sequence": sequence,
        "invalid_residues": invalid,
    }


def get_sequence_length(sequence: str) -> int:
    """
    Get length of valid amino acid sequence.
    """
    result = validate_sequence(sequence)

    if not result["valid"]:
        return 0

    return len(result["sequence"])


# ============================================================
# Amino Acid Operations
# ============================================================

def get_amino_acid_properties(aa_code: str) -> Optional[Dict]:
    """
    Retrieve properties for a single amino acid.
    """
    aa_code = aa_code.upper()

    return AMINO_ACID_DATA.get(aa_code)


def get_all_amino_acids_info() -> Dict:
    """
    Return all amino acid information.
    """
    return AMINO_ACID_DATA


def classify_sequence(sequence: str) -> Dict[str, Any]:
    """
    Analyze amino acid composition.
    """
    validation = validate_sequence(sequence)

    if not validation["valid"]:
        return {
            "valid": False,
            "invalid_residues": validation["invalid_residues"],
        }

    sequence = validation["sequence"]
    counts = Counter(sequence)

    hydrophobic = sum(counts[aa] for aa in HYDROPHOBIC_AAS)
    charged = sum(counts[aa] for aa in CHARGED_AAS)
    positive = sum(counts[aa] for aa in POSITIVE_AAS)
    negative = sum(counts[aa] for aa in NEGATIVE_AAS)
    polar = sum(counts[aa] for aa in POLAR_AAS)
    aromatic = sum(counts[aa] for aa in AROMATIC_AAS)

    return {
        "valid": True,
        "length": len(sequence),
        "composition": dict(counts),
        "hydrophobic_count": hydrophobic,
        "charged_count": charged,
        "positive_count": positive,
        "negative_count": negative,
        "polar_count": polar,
        "aromatic_count": aromatic,
    }


# ============================================================
# Sequence Conversion
# ============================================================

def sequence_to_rdkit_mol(
    sequence: str,
    *,
    hydrolyze_sidechains: bool = True,
) -> Optional[Chem.Mol]:
    """
    Build a single linear peptide ``Mol`` via iterative amide coupling in RDKit.

    Standard L-amino acids only (same alphabet as ``validate_sequence``).

    Asp/Glu side chains are methyl-protected during assembly to avoid coupling at
    side-chain carboxylates; by default those esters are hydrolyzed back to acids
    in the returned structure.

    Limitations: no cyclization, no non-natural residues, no PTMs, no explicit
    stereochemistry beyond building-block templates.
    """
    validation = validate_sequence(sequence)

    if not validation["valid"]:
        return None

    seq = validation["sequence"]

    chain = Chem.MolFromSmiles(AA_BUILDING_BLOCK_SMILES[seq[0]])
    if chain is None:
        return None

    for aa in seq[1:]:
        chain = _couple_residue(chain, aa)
        if chain is None:
            return None

    if hydrolyze_sidechains:
        chain = _hydrolyze_building_block_methyl_esters(chain)

    return chain


def sequence_to_smiles(
    sequence: str,
    *,
    hydrolyze_sidechains: bool = True,
) -> Optional[str]:
    """
    Canonical SMILES for a linear peptide built from the sequence.

    See ``sequence_to_rdkit_mol`` for chemistry assumptions.
    """
    mol = sequence_to_rdkit_mol(sequence, hydrolyze_sidechains=hydrolyze_sidechains)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


# ============================================================
# Batch Processing
# ============================================================

def validate_batch_input(data: List[Dict]) -> Dict[str, Any]:
    """
    Validate batch peptide input.

    Expected format:
    [
        {
            "name": "pep1",
            "sequence": "ACDE"
        },
        {
            "name": "pep2",
            "smiles": "CCO"
        }
    ]
    """

    results = []
    seen_names = set()
    duplicate_names = set()

    for item in data:
        name = item.get("name")

        if name in seen_names:
            duplicate_names.add(name)

        seen_names.add(name)

        result = {
            "name": name,
            "valid": False,
            "errors": [],
        }

        if "sequence" in item:
            seq_result = validate_sequence(item["sequence"])

            if seq_result["valid"]:
                result["valid"] = True
            else:
                result["errors"].append(
                    f"Invalid residues: {seq_result['invalid_residues']}"
                )

        elif "smiles" in item:
            if validate_smiles(item["smiles"]):
                result["valid"] = True
            else:
                result["errors"].append("Invalid SMILES")

        else:
            result["errors"].append(
                "Missing sequence or smiles field"
            )

        results.append(result)

    return {
        "results": results,
        "duplicate_names": sorted(list(duplicate_names)),
    }


# ============================================================
# Utility Functions
# ============================================================

def count_standard_aas() -> int:
    """
    Return number of standard amino acids.
    """
    return len(STANDARD_AAS)


def list_standard_aas() -> List[str]:
    """
    Return sorted list of amino acid codes.
    """
    return sorted(list(STANDARD_AAS))


def get_aa_summary() -> str:
    """
    Generate formatted amino acid summary.
    """
    lines = []

    for aa in sorted(STANDARD_AAS):
        info = AMINO_ACID_DATA[aa]

        lines.append(
            f"{aa}: {info['name']} | "
            f"Charge={info['charge']} | "
            f"Mass={info['molar_mass']}"
        )

    return "\n".join(lines)


# ============================================================
# Example Usage
# ============================================================

if __name__ == "__main__":

    sequence = "ACDEFGHIK"

    print("Sequence Validation")
    print(validate_sequence(sequence))

    print("\nSequence Classification")
    print(classify_sequence(sequence))

    print("\nSequence Length")
    print(get_sequence_length(sequence))

    print("\nCanonical SMILES (small-molecule demo, not the peptide)")
    demo_smiles = "CC(O)C"
    print(f"  input:  {demo_smiles}")
    print(f"  output: {canonicalize_smiles(demo_smiles)}")

    print("\nPeptide SMILES (from sequence)")
    print(sequence_to_smiles(sequence))

    print("\nAA Summary")
    print(get_aa_summary())