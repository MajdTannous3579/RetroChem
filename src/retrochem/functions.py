import os
from rdkit import Chem
from urllib.request import urlopen
from urllib.parse import quote


def name_to_smiles(name: str) -> str:
    """
    Convert a chemical name to its SMILES representation using the NCI Chemical Identifier Resolver API.

    Parameters
    ----------
    name : str
        The common or IUPAC name of the molecule.

    Returns
    -------
    str
        The corresponding SMILES string.

    Raises
    ------
    ValueError
        If the conversion fails or the API is unreachable.
    """
    try:
        url = f"https://cactus.nci.nih.gov/chemical/structure/{quote(name)}/smiles"
        response = urlopen(url)
        smiles = response.read().decode('utf-8').strip()
        if not smiles:
            raise ValueError(f"Empty response for name: {name}")
        return smiles
    except Exception as e:
        raise ValueError(f"Could not convert name '{name}' to SMILES: {e}")



def structure_to_smiles(
    structure_input: Chem.Mol | str,
    input_format: str | None = None,
) -> str:
    """
    Convert a molecular structure (file path, MolBlock string, or RDKit Mol object)
    to its SMILES representation.

    Parameters
    ----------
    structure_input : str or rdkit.Chem.Mol
        If str:
            - Path to a structure file (e.g., .mol, .sdf).
            - A MolBlock string (containing 'M  END').
        If rdkit.Chem.Mol:
            - An RDKit Mol object.
    input_format : str, optional
        If provided and structure_input is a file path, forces the file format (e.g., 'sdf', 'mol').

    Returns
    -------
    str
        The SMILES string for the input structure (non-canonical).

    Raises
    ------
    ValueError
        If the input cannot be parsed or conversion fails.
    """
    # Load molecule
    if isinstance(structure_input, Chem.Mol):
        mol = structure_input
    elif isinstance(structure_input, str):
        if os.path.isfile(structure_input):
            fmt = input_format or os.path.splitext(structure_input)[1].lstrip('.')
            if fmt.lower() in ('sdf', 'mol'):
                mol = Chem.MolFromMolFile(structure_input, sanitize=True, removeHs=False)
            else:
                raise ValueError(f"Unsupported file format: {fmt}")
        elif 'M  END' in structure_input:
            mol = Chem.MolFromMolBlock(structure_input, sanitize=True, removeHs=False)
        else:
            raise ValueError("String input is neither a valid file path nor a MolBlock.")
    else:
        raise TypeError("structure_input must be an RDKit Mol or a string path/MolBlock.")

    if mol is None:
        raise ValueError("Failed to parse molecular structure.")

    # Generate non-canonical SMILES
    smiles = Chem.MolToSmiles(mol, canonical=False)
    return smiles


def canonicalize_smiles(smiles: str) -> str:
    """
    Convert a SMILES string to its canonical form using RDKit.

    Parameters
    ----------
    smiles : str
        The input SMILES string (canonical or non-canonical).

    Returns
    -------
    str
        The canonical SMILES representation.

    Raises
    ------
    ValueError
        If the input SMILES cannot be parsed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)
