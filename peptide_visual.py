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

import csv
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

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
# CSV export (peptide library assets)
# ============================================================

_FILENAME_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def _sanitize_export_id(raw_id: str, fallback: str = "peptide") -> str:
    cleaned = _FILENAME_SAFE.sub("_", raw_id.strip())
    return cleaned or fallback


def _resolve_smiles_from_row(
    row: Dict[str, str],
    *,
    hydrolyze_sidechains: bool = True,
) -> tuple[Optional[str], List[str]]:
    """
    Return (smiles, errors). Uses ``smiles`` column first, else builds from ``sequence``.
    """
    errors: List[str] = []
    smiles = (row.get("smiles") or "").strip()
    if smiles:
        if not validate_smiles(smiles):
            errors.append("invalid SMILES in CSV")
            return None, errors
        canonical = canonicalize_smiles(smiles)
        return canonical, errors

    sequence = (row.get("sequence") or "").strip()
    if not sequence:
        errors.append("missing smiles and sequence")
        return None, errors

    from peptide_parser import sequence_to_smiles

    built = sequence_to_smiles(sequence, hydrolyze_sidechains=hydrolyze_sidechains)
    if built is None:
        errors.append("could not build SMILES from sequence")
        return None, errors
    return built, errors


def export_peptides_from_csv(
    csv_path: Union[str, Path],
    output_dir: Union[str, Path],
    *,
    hydrolyze_sidechains: bool = True,
    skip_3d: bool = False,
    write_grid: bool = False,
    write_descriptor_summary: bool = True,
) -> Dict[str, Any]:
    """
    Read a peptide CSV and write per-row assets (same pattern as the aspirin demo).

    For each valid row writes under ``output_dir``:

    - ``{id}.png`` — 2D depiction
    - ``{id}_3d.sdf`` — 3D conformer (UFF), unless ``skip_3d=True``
    - ``{id}_ligand.sdf`` — docking-style 3D (MMFF), unless ``skip_3d=True``

    Also writes ``export_summary.csv`` in ``output_dir`` with per-row status.

    CSV columns (header row):

    - **id** or **name** (required): used for filenames
    - **smiles** (recommended): canonical or valid SMILES
    - **sequence** (optional): used when ``smiles`` is empty; built via ``sequence_to_smiles``
    """
    csv_path = Path(csv_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    results: List[Dict[str, Any]] = []
    smiles_for_grid: List[str] = []

    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("CSV has no header row")

        fieldnames_lower = {f.lower(): f for f in reader.fieldnames}
        if "id" not in fieldnames_lower and "name" not in fieldnames_lower:
            raise ValueError("CSV must include an 'id' or 'name' column")
        if "smiles" not in fieldnames_lower and "sequence" not in fieldnames_lower:
            raise ValueError("CSV must include 'smiles' and/or 'sequence' column")

        id_key = fieldnames_lower.get("id") or fieldnames_lower.get("name")
        for index, row in enumerate(reader):
            raw_id = row.get(id_key, "") if id_key else ""
            if not raw_id:
                raw_id = f"row_{index + 1}"
            pep_id = _sanitize_export_id(str(raw_id), fallback=f"row_{index + 1}")

            entry: Dict[str, Any] = {
                "id": pep_id,
                "smiles": "",
                "png": False,
                "sdf_3d": False,
                "ligand_sdf": False,
                "errors": [],
            }

            smiles, resolve_errors = _resolve_smiles_from_row(
                row,
                hydrolyze_sidechains=hydrolyze_sidechains,
            )
            entry["errors"].extend(resolve_errors)

            if smiles is None:
                results.append(entry)
                continue

            entry["smiles"] = smiles
            png_path = output_dir / f"{pep_id}.png"
            sdf_3d_path = output_dir / f"{pep_id}_3d.sdf"
            ligand_path = output_dir / f"{pep_id}_ligand.sdf"

            entry["png"] = save_molecule_image(smiles, str(png_path))
            if not entry["png"]:
                entry["errors"].append("failed to save PNG")

            if not skip_3d:
                entry["sdf_3d"] = generate_3d_conformer(smiles, str(sdf_3d_path))
                if not entry["sdf_3d"]:
                    entry["errors"].append("3D embed/UFF failed")
                entry["ligand_sdf"] = prepare_ligand_for_docking(
                    smiles, str(ligand_path)
                )
                if not entry["ligand_sdf"]:
                    entry["errors"].append("docking prep/MMFF failed")
            else:
                entry["sdf_3d"] = None
                entry["ligand_sdf"] = None

            if entry["png"]:
                smiles_for_grid.append(smiles)

            results.append(entry)

    if write_grid and smiles_for_grid:
        save_grid_image(smiles_for_grid, str(output_dir / "library_grid.png"))

    if write_descriptor_summary:
        _write_descriptor_summary_csv(results, output_dir / "descriptors_summary.csv")

    summary_path = output_dir / "export_summary.csv"
    _write_export_summary_csv(results, summary_path)

    exported = sum(
        1
        for r in results
        if r["png"] and (skip_3d or (r["sdf_3d"] and r["ligand_sdf"]))
    )
    partial = sum(1 for r in results if r["png"] and r["errors"])
    failed = len(results) - exported - partial

    return {
        "csv_path": str(csv_path),
        "output_dir": str(output_dir),
        "results": results,
        "summary": {
            "total": len(results),
            "fully_exported": exported,
            "partial": partial,
            "failed": failed,
            "export_summary_csv": str(summary_path),
        },
    }


def _write_export_summary_csv(results: List[Dict[str, Any]], path: Path) -> None:
    fieldnames = [
        "id",
        "smiles",
        "png",
        "sdf_3d",
        "ligand_sdf",
        "errors",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(
                {
                    "id": row["id"],
                    "smiles": row.get("smiles", ""),
                    "png": row["png"],
                    "sdf_3d": row.get("sdf_3d"),
                    "ligand_sdf": row.get("ligand_sdf"),
                    "errors": "; ".join(row.get("errors") or []),
                }
            )


def _write_descriptor_summary_csv(results: List[Dict[str, Any]], path: Path) -> None:
    rows_with_smiles = [r for r in results if r.get("smiles")]
    if not rows_with_smiles:
        return

    desc_rows = batch_calculate_descriptors([r["smiles"] for r in rows_with_smiles])
    for meta, desc in zip(rows_with_smiles, desc_rows):
        desc["id"] = meta["id"]

    fieldnames = ["id", "smiles", "valid"] + sorted(
        k for k in desc_rows[0].keys() if k not in ("id", "smiles", "valid")
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(desc_rows)


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