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
    "[C:1]([OX2H:2])([OX2H:3])([!O:4])([!O:5])>>[C:1](=O)([!O:4])([!O:5]).O[*:2].O[*:3]", #ketone hydration
    "[C:1]([OX2H:2])([!O:3])([!O:4])>>[C:1](=O)([!O:3])([!O:4])", #carbonyl reduction with hydride 
    "[C:1]([OX2H:2])([S:3](=O)(=O)[O-])([!O:4])>>[C:1](=O)([!O:4]).[Na+].[O-]S(=O)(=O)O", #bisulfide addition
    "[CX3:1]([#6:2])([#6:3])=[NX2:4]([#6:5])([#6:6])>>[C:1](=O)[#6:3][#6:4].[N:4]([#6:5])([#6:6])", #imine formation
    "[CX4H2:1][NX3:2]([#6:3])([#6:4])>>[C:1](=O)[NX3:2]([#6:3])([#6:4])" #amide reduction
    "[C:1]([CN:2])([#6:3])([#6:4])([OX2H:5])>>[C:1](=O)([#6:3])([#6:4])" #Cyanohydride formation
    "[CX4H:1]([#6:2])([#6:3])-[NX3;H1;X3:4]([#6:5])>>[C:1](=O)([#6:2])([#6:3]).[N]([#6:4])([#6:5])" #carbonyl to amine conversion
    "[CX4H:1]([#6:2])([#6:3])>>[C:1](=O)([#6:2])([#6:3].[NH2(NH2)])" #wolf kischner carbonyl to alkane 
    "[CX4H:1]([#6:2])([#6:3])>>[C:1](=O)([#6:2])([#6:3])"

]

REACTION_REVERSERS: List[Callable[[str], str | None]] = [
    reverse_reaction_generator(i) for i in REACTION_DATABASE
]

def list_reactants(smiles: str)->list[str]:
    return list(filter(lambda x: x is not None, [fn(smiles) for fn in REACTION_REVERSERS]))

