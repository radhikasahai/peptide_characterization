from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Crippen,
    Draw,
    Lipinski,
    MACCSkeys,
    Descriptors,
    rdMolDescriptors,
)
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

from typing import Dict, Optional, List

from chem_utils import canonicalize_smiles, smiles_to_mol, validate_smiles


# ============================================================
# Basic Molecule Utilities
# ============================================================


def calculate_basic_descriptors(smiles: str) -> Optional[Dict]:
    """
    Calculate common molecular descriptors.
    """

    mol = smiles_to_mol(smiles)

    if mol is None:
        return None

    descriptors = {
        "molecular_weight": Descriptors.MolWt(mol),
        "exact_molecular_weight": Descriptors.ExactMolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "hbond_donors": Lipinski.NumHDonors(mol),
        "hbond_acceptors": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "ring_count": Lipinski.RingCount(mol),
        "heavy_atom_count": Lipinski.HeavyAtomCount(mol),
        "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
    }

    return descriptors



def calculate_full_rdkit_descriptors(smiles: str) -> Optional[Dict]:
    """
    Calculate all available RDKit descriptors.
    """

    mol = smiles_to_mol(smiles)

    if mol is None:
        return None

    descriptor_dict = {}

    for descriptor_name, function in Descriptors.descList:
        try:
            descriptor_dict[descriptor_name] = function(mol)
        except Exception:
            descriptor_dict[descriptor_name] = None

    return descriptor_dict


# ============================================================
# Fingerprint Generation
# ============================================================


def generate_morgan_fingerprint(
    smiles: str,
    radius: int = 2,
    n_bits: int = 2048,
):
    """
    Generate Morgan fingerprint.
    """

    mol = smiles_to_mol(smiles)

    if mol is None:
        return None

    generator = GetMorganGenerator(radius=radius, fpSize=n_bits)

    fingerprint = generator.GetFingerprint(mol)

    return fingerprint



def generate_maccs_fingerprint(smiles: str):
    """
    Generate MACCS fingerprint.
    """

    mol = smiles_to_mol(smiles)

    if mol is None:
        return None

    return MACCSkeys.GenMACCSKeys(mol)



def fingerprint_to_bitstring(fingerprint) -> Optional[str]:
    """
    Convert RDKit fingerprint to bitstring.
    """

    if fingerprint is None:
        return None

    return fingerprint.ToBitString()



def calculate_tanimoto_similarity(smiles1: str, smiles2: str) -> Optional[float]:
    """
    Calculate Tanimoto similarity between two molecules.
    """

    fp1 = generate_morgan_fingerprint(smiles1)
    fp2 = generate_morgan_fingerprint(smiles2)

    if fp1 is None or fp2 is None:
        return None

    return TanimotoSimilarity(fp1, fp2)


# ============================================================
# Molecular Visualization
# ============================================================


def save_molecule_image(
    smiles: str,
    output_path: str = "molecule.png",
    size: tuple = (500, 500),
):
    """
    Save 2D molecular structure image.
    """

    mol = smiles_to_mol(smiles)

    if mol is None:
        return False

    img = MolToImage(mol, size=size)

    img.save(output_path)

    return True



def save_grid_image(
    smiles_list: List[str],
    output_path: str = "molecule_grid.png",
    mols_per_row: int = 4,
):
    """
    Save grid image of molecules.
    """

    mols = []

    for smiles in smiles_list:
        mol = smiles_to_mol(smiles)

        if mol is not None:
            mols.append(mol)

    if len(mols) == 0:
        return False

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=(300, 300),
    )

    img.save(output_path)

    return True


# ============================================================
# 3D Conformer Generation
# ============================================================


def generate_3d_conformer(
    smiles: str,
    output_sdf: str = "molecule_3d.sdf",
):
    """
    Generate 3D conformer and save as SDF.
    """

    mol = smiles_to_mol(smiles)

    if mol is None:
        return False

    mol = Chem.AddHs(mol)

    status = AllChem.EmbedMolecule(mol)

    if status != 0:
        return False

    AllChem.UFFOptimizeMolecule(mol)

    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()

    return True


# ============================================================
# Docking Preparation
# ============================================================


def prepare_ligand_for_docking(
    smiles: str,
    output_sdf: str = "ligand.sdf",
):
    """
    Prepare ligand structure for docking workflows.
    """

    mol = smiles_to_mol(smiles)

    if mol is None:
        return False

    mol = Chem.AddHs(mol)

    status = AllChem.EmbedMolecule(mol)

    if status != 0:
        return False

    AllChem.MMFFOptimizeMolecule(mol)

    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()

    return True


# ============================================================
# Batch Processing
# ============================================================


def batch_calculate_descriptors(smiles_list: List[str]):
    """
    Calculate descriptors for multiple molecules.
    """

    results = []

    for smiles in smiles_list:
        result = {
            "smiles": smiles,
            "valid": validate_smiles(smiles),
        }

        if result["valid"]:
            result.update(calculate_basic_descriptors(smiles))

        results.append(result)

    return results


# ============================================================
# Example Usage
# ============================================================

if __name__ == "__main__":
    # Example small molecule (classic aspirin demo).
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"

    print("=" * 60)
    print("SMILES Validation (aspirin)")
    print("=" * 60)
    print(validate_smiles(aspirin))

    print("\n" + "=" * 60)
    print("Canonical SMILES")
    print("=" * 60)
    print(canonicalize_smiles(aspirin))

    print("\n" + "=" * 60)
    print("Basic Descriptors")
    print("=" * 60)
    for key, value in calculate_basic_descriptors(aspirin).items():
        print(f"{key}: {value}")

    print("\n" + "=" * 60)
    print("Peptide built from sequence (for RDKit + peptide demo)")
    print("=" * 60)

    from peptide_parser import sequence_to_smiles as seq2smi

    pep = "ACYDEKGP"
    pep_smi = seq2smi(pep)
    print(f"Sequence: {pep}")
    print(f"SMILES:   {pep_smi}")

    fp = generate_morgan_fingerprint(pep_smi)
    print(type(fp))
    print(f"Fingerprint bits: {len(fp.ToBitString())}")

    print("\n" + "=" * 60)
    print("Molecular Visualization → aspirin.png")
    print("=" * 60)
    save_molecule_image(aspirin, "aspirin.png")
    print("Saved aspirin.png")

    print("\n" + "=" * 60)
    print("3D conformer → aspirin_3d.sdf")
    print("=" * 60)
    generate_3d_conformer(aspirin, "aspirin_3d.sdf")
    print("Saved aspirin_3d.sdf")

    print("\n" + "=" * 60)
    print("Docking prep → aspirin_ligand.sdf")
    print("=" * 60)
    prepare_ligand_for_docking(aspirin, "aspirin_ligand.sdf")
    print("Saved aspirin_ligand.sdf")