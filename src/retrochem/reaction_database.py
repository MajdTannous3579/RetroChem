from typing import *
from rdkit import Chem 
from rdkit.Chem import rdChemReactions as Reac

def reverse_reaction_generator(reaction_smart: str)->Callable[[str], str | None]:
    rxn = Reac.ReactionFromSmarts(reaction_smart)
    def reverser_to_smiles(smiles: str) -> str | None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        prods = rxn.RunReactants((mol,))
        if not prods: # if no products could be generated
            return None
        first = prods[0] # FIXME: recheck if all tuples
        first_smiles = []
        for mol in first:
            first_smiles.append(Chem.MolToSmiles(mol, canonical = True))
        combo = ".".join(first_smiles)
        return combo
    return reverser_to_smiles

REACTION_DATABASE: List[str] = [
    "[C:1](-O[*:2])(-O[*:3])-[!O:4]>>[C:1](=O)-[!O:4].O[*:2].O[*:3]", #reverse acetilisation
]

REACTION_REVERSERS: List[Callable[[str], str | None]] = [
    reverse_reaction_generator(i) for i in REACTION_DATABASE
]

def list_reactants(smiles: str)->list[str]:
    return list(filter(lambda x: x is not None, [fn(smiles) for fn in REACTION_REVERSERS]))


if __name__ == "__main__":
    print(fn("CC(C)(OC)OC") for fn in REACTION_REVERSERS)