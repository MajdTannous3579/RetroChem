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
