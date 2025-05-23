from typing import *
from rdkit import Chem 
from rdkit.Chem import rdChemReactions as Reac
import json
import os

# Set working directory to the script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

SmartsConditionsPair = Tuple[str, dict]
"""
A tuple of a smart reaction rule, and a dictionnary of corresponding conditions
"""
SmilesConditionsPair= Tuple[str,dict]
"""
A tuple of a smile, and a dictionnary of corresponding conditions, needed to get to this smile
"""

def reverse_reaction_generator(reaction_smart: SmartsConditionsPair)->Callable[[str], SmilesConditionsPair | None]:
    """
    Creates a callable that performs a single reverse reaction transformation on a SMILES string.

    Parameters
    ----------
    reaction_smart : tuple[str, dict]
        A tuple where the first element is a forward reaction SMARTS string,
        and the second is a dictionary of associated reaction conditions.

    Returns
    -------
    Callable[[str], tuple[str, dict] or None]
        A function that takes a SMILES string of a product,
        and returns a dot-separated SMILES string of reactants with the same conditions, or None if no match found.

    Raises
    ------
    ValueError
        If the input SMILES string is invalid.
    """
    rxn = Reac.ReactionFromSmarts(reaction_smart[0])
    cond = reaction_smart[1]
    def reverser_to_smiles(smiles: str) -> str | None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        try:
            prods = rxn.RunReactants((mol,))
        except Exception:
            return None
        if not prods: # if no products could be generated
            return None
        first = prods[0]
        first_smiles = [Chem.MolToSmiles(m, canonical=True) for m in first]
        combo = ".".join(first_smiles)
        return (combo, cond)
    return reverser_to_smiles

def load_database(path: str)->List[SmartsConditionsPair] | None:
    """
    Load a JSON-based reaction database from disk and convert each entry to a tuple.

    Parameters
    ----------
    path : str
        Path to the `.db` file containing a list of SMARTS-condition entries.

    Returns
    -------
    list of (str, dict) or None
        A list of tuples containing reaction SMARTS and associated conditions,
        or None if the file is missing or malformed.
    """
    try:
        with open(path, "r") as file:
            ret = json.load(file)
            return [tuple(pair) for pair in ret]
    except:
        return None

REACTION_DATABASES: dict[str, list[SmartsConditionsPair]] = {}
"""
A dictionary mapping database names to their corresponding list of (SMARTS, conditions) reaction rules.
Used during UI selection and query operations.
"""
REACTION_REVERSERS: dict[str, List[Callable[[str], SmartsConditionsPair | None]]] = {}
"""
A dictionary mapping database names to a list of compiled reverse-reaction callables.
These are pre-compiled at registration time for performance.
"""


def register_database(values: List[SmartsConditionsPair], database: str)->None:
    """
    Registers a list of SMARTS reaction rules and compiles reverse reaction functions.

    Parameters
    ----------
    values : list of (str, dict)
        List of (SMARTS, conditions) reaction definitions to register.
    database : str
        The name under which to store and retrieve this database.
    """
    REACTION_DATABASES[database] = values
    REACTION_REVERSERS[database] = [
        reverse_reaction_generator(i) for i in values
    ]

def clear_registered_databases():
    """
    Clears all previously registered databases and reverser functions.
    Useful during UI resets or database refresh operations.
    """
    REACTION_DATABASES.clear()
    REACTION_REVERSERS.clear()

def list_reactants(smiles: str, database: str)->List[SmilesConditionsPair] | None:
    """
    Apply all reverse reaction functions from a specified database to a given product SMILES.

    Parameters
    ----------
    smiles : str
        The target molecule in SMILES format.
    database : str
        The name of the database whose reversers should be used.

    Returns
    -------
    list of (str, dict) or None
        A list of tuples with dot-separated reactants and associated conditions.
        Returns None if the database is not registered.
    """
    get = REACTION_REVERSERS.get(database)
    if get is None:
        return None
    ret = []
    for fn in get:
        val = fn(smiles)
        if val is not None:
            ret.append(val)
    return ret


def add_new_smart(database_name: str, product: str, reactants: list[str], conditions: dict[str, str] = dict())->None:
    """
    Adds a new reaction SMARTS rule to a named database and saves it to disk.

    Parameters
    ----------
    database_name : str
        Name of the database file (without `.db` extension).
    product : str
        SMILES string of the reaction product.
    reactants : list of str
        List of SMILES strings for all reactants in the reaction.
    conditions : dict[str, str], optional
        Dictionary of reaction conditions (solvent, temperature, catalyst, etc.).

    Notes
    -----
    - Existing reactions in the database are preserved and the new one is appended.
    - Stereochemistry markers ('\\', '/') are stripped for simplicity.
    """
    file_path = f'{database_name}.db'
    previous = load_database(file_path) or []
    product_smarts = Chem.MolToSmarts(Chem.MolFromSmiles(product)).replace('\\', '-').replace('/', '-')
    reactants_smarts = '.'.join([Chem.MolToSmarts(Chem.MolFromSmiles(i)).replace('\\', '-').replace('/', '-') for i in reactants])
    previous.append((f'{product_smarts}>>{reactants_smarts}', conditions))
    register_database(previous, database_name)
    with open(file_path, 'w') as file:
        file.write('[\n')
        if len(previous) != 0:
            file.write(f'  [\n    "{previous[0][0]}",\n    {json.dumps(previous[0][1])}\n  ]')
        for i in range(1, len(previous)):
            file.write(f',  [\n    "{previous[i][0]}",\n    {json.dumps(previous[i][1])}\n  ]')
        file.write('\n]\n')