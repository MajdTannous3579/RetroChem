import pytest 
from rdkit import Chem 
from rdkit.Chem import Draw 
from rdkit.Chem import rdChemReactions as Reac

# SMARTS that matches a tetra-valent carbon attached to 2 hetero-O groups
# plus any single "R" (non-oxygen) substituent.
ACETAL_SMARTS = "[C:1](-O[*:2])(-O[*:3])-[!O:4]"

# 2⃣  Reverse-acetalisation reaction:
#     keep the R-group (atom 4) and carbon atom 1, convert it to a carbonyl.
#     The OR groups become separate alcohol fragments (dropped here).
REACTION_SMARTS = (
    "[C:1](-O[*:2])(-O[*:3])-[!O:4]>>"
    "[C:1](=O)-[!O:4]"                # product 0  (aldehyde/ketone fragment)
    ".O[*:2].O[*:3]"                  # product 1+2 (two alcohols)
)
rxn = Reac.ReactionFromSmarts(REACTION_SMARTS)

def reverse_acetalisation(smiles: str) -> str | None:
    """
    Given a SMILES, return a SMILES with each acetal centre converted
    to a carbonyl (plus two alcohol by-products).  If no acetal is found,
    return None.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    # Run the reaction (returns tuple of product tuples)
    prods = rxn.RunReactants((mol,))
    if not prods:
        return None

    # Take the first product tuple, merge into one combined molecule
    first = prods[0]                  # (prod0_mol, prod1_mol, prod2_mol)
    combo = Chem.CombineMols( first[0], Chem.CombineMols(first[1], first[2]))
    return Chem.MolToSmiles(combo, canonical=True)


if __name__ == "__main__":
    test = "CC(C)(OC)OC"   # tert-butyl dimethyl acetal
    print("input :", test)
    print("output:", reverse_acetalisation(test))
