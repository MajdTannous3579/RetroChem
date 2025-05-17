from typing import *
from rdkit import Chem 
from rdkit.Chem import rdChemReactions as Reac
import json
import os

# Set working directory to the script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

SmartsConditionsPair = Tuple[str, dict]
SmilesConditionsPair= Tuple[str,dict]

def reverse_reaction_generator(reaction_smart: SmartsConditionsPair)->Callable[[str], SmilesConditionsPair | None]:
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
        first_smiles = sorted(Chem.MolToSmiles(m, canonical=True) for m in first)
        combo = ".".join(first_smiles)
        return (combo, cond)
    return reverser_to_smiles

def load_database(path: str)->List[SmartsConditionsPair] | None:
    try:
        with open(path, "r") as file:
            ret = json.load(file)
            return ret
    except:
        return None

REACTION_DATABASES: dict[str, list[SmartsConditionsPair]] = {}
REACTION_REVERSERS: dict[str, List[Callable[[str], SmartsConditionsPair | None]]] = {}


def register_database(values: List[SmartsConditionsPair], database: str)->None:
    REACTION_DATABASES[database] = values
    REACTION_REVERSERS[database] = [
        reverse_reaction_generator(i) for i in values
    ]

def list_reactants(smiles: str, database: str)->List[SmilesConditionsPair] | None:
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
    file_path = f'{database_name}.db'
    previous = load_database(file_path) or []
    product_smarts = Chem.MolToSmarts(Chem.MolFromSmiles(product)).replace('\\', '-').replace('/', '-')
    reactants_smarts = '.'.join([Chem.MolToSmarts(Chem.MolFromSmiles(i)).replace('\\', '-').replace('/', '-') for i in reactants])
    previous.append((f'{product_smarts}>>{reactants_smarts}', conditions))
    with open(file_path, 'w') as file:
        file.write('[\n')
        if len(previous) != 0:
            file.write(f'  [\n    "{previous[0][0]}",\n    {json.dumps(previous[0][1])}\n  ]')
        for i in range(1, len(previous)):
            file.write(f'  [\n    "{previous[i][0]}",\n    {json.dumps(previous[i][1])}\n  ]')
        file.write('\n]\n')