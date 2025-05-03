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
        combo = Chem.CombineMols( first[0], Chem.CombineMols(first[1], first[2]))
        return Chem.MolToSmiles(combo, canonical=True)
    return reverser_to_smiles

REACTION_DATABASE: List[str] = [
    "[C:1](-O[*:2])(-O[*:3])-[!O:4]>>[C:1](=O)-[!O:4].O[*:2].O[*:3]",
]

REACTION_REVERSERS: List[Callable[[str], str | None]] = [
    reverse_reaction_generator(i) for i in REACTION_DATABASE
]

if __name__ == '__main__':
    test = "CC(C)(OC)OC"
    print(list(filter(lambda x: x is not None, [fn(test) for fn in REACTION_REVERSERS])))