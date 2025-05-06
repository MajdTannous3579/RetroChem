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
    "[CX3:1]([#6:2])([#6:3])=[NX2:4]([#6:5])>>[C:1](=O)[#6:3][#6:4].[N:4]([#6:5])([#6:6])", #imine formation
    "[CX4H:1][NX3:2]([#6:3])([#6:4])>>[C:1](=O)[NX3:2]([#6:3])([#6:4])", #amide reduction
    "[C:1]([CN:2])([#6:3])([#6:4])([OX2H:5])>>[C:1](=O)([#6:3])([#6:4])", #Cyanohydride formation
    "[CX4H:1]([#6:2])([#6:3])-[NX3;H1;X3:4]([#6:5])>>[C:1](=O)([#6:2])([#6:3]).[N]([#6:4])([#6:5])", #carbonyl to amine conversion 
    "[CX4H2:1]([#6:2])([#6:3])>>[C:1](=O)([#6:2])([#6:3])", #clemensen rxn
    "[CX4H2:1]([#6:2])([#6:3])>>[C:1](=O)([#6:2])([#6:3]).NN", #wolf kischner rxn 
    "[#6:1][C:2]([O:6][H])[#6:3]C(=O)[#6:5]>>[#6:1][CH2:2]C(=O)[#6:3].[#6:4]C(=O)[#6:5]", #aldol reaction
    "[#6:1][C:2]=[C:3]C(=O)[#6:4]>>[#6:1][C:2]([O:5][H])[C:3]C(=O)[#6:4]", #aldol condensation
    "[#6:1][CH2:2][OX2H:3]>>[#6:1][CX2H:2]=[O:3]", #aldehyde reduction (not ketone)
    "[a:1]C(=O)C(O)[a:5]>>[a:1]C=O.[a:5]C=O", #benzoin reaction 
    "[a:1][C:2](=O)[C:3]([O:4][Si:5])[a:6]>>[a:1]C=O.[a:6]C([O:4][Si:5])", #cross benzoin reaction
    "[#6:3][CH2:1][OX2H:2]>>[CH2:1]=[O:2].[#6:3][Mg:4][Br:5]", #grignard primary alcohol
    "[C:1]([#6:2])([#6:4])[OX2H:3]>>[C:1]([#6:2])=[O:3].[#6:4][Mg:5][Br:6]", #grignard secondary 
    "[C:1]([#6:2])([#6:3])([#6:5])[OX2H:4]>>[C:1]([#6:2])([#6:3])=[O:4].[#6:5][Mg:6][Br:7]", #grignard tertiary
    "[C:1]([Cl:2])[C:3]>>[C:1]=[C:3].[Cl:2]", #acid addition to alkene (chloride)
    "[C:1]([Br:2])[C:3]>>[C:1]=[C:3].[Br:2]", #acid addition to alkene (bromide)
    "[C:1]([I:2])[C:3]>>[C:1]=[C:3].[I:2]", #acid addition to alkene (iodine)
    "[C:1]([F:2])[C:3]>>[C:1]=[C:3].[F:2]", #acid addition to alkene (floride)
    "[C:1]([OX2H:2])[C:3]>>[C:1]=[C:3].[O]", #hydration of alkene 
    "[C:1]([Br:2])[C:3]([Br:4])>>[C:1]=[C:3].[Br:2][Br:4]", #addition of hallide to alkene
    "[C:1]([Cl:2])[C:3]([Cl:4])>>[C:1]=[C:3].[Cl:2][Cl:4]", #addition of hallide to alkene
    "[C:1]([I:2])[C:3]([I:4])>>[C:1]=[C:3].[I:2][I:4]", #addition of hallide to alkene
    "[C:1]1O[C:2]1>>[C:1]=[C:2]", #epoxydation
    "[#6:1]=[#6:2][OX2H]>>[#6:1]([CX2]#C)[#6:2].[OH2]", #alkyne hydration 
    "[#6:1]=C=C>>[#6:1]([CX2]#C)[#6:2]", #alkyne+base --> allene pH<30 
    "[#6:1][#6:2][C:3]#[CH:4]>>[#6:1][C:3]#[C:4][#6:2]", #alkyne+base --> allene pH>30-35
    "[#6:1][C:2]=[C:3][#6:4]>>[#6:1][C:2]#[C:3][#6:4].[H][H]" #lindlar hydrogenation alkynes 
    
      
]


from typing import List

REACTION_CONDITIONS: List[dict] = [
    {"medium": "acidic"},  # reverse acetalization
    {"medium": "acidic", "solvent": "water"},  # ketone hydration
    {
        "step1": "hydride (e.g. LAH) attacks central carbon",
        "step2": "acidic medium protonates the oxyanion",
    },  # carbonyl reduction with hydride
    {},  # bisulfite addition
    {"medium": "acidic"},  # imine formation
    {"presence": "hydride (e.g. LAH) attacks and reduces central carbon"},  # amide reduction
    {},  # cyanohydrin formation
    {"solvent": "NaBH3CN", "pH": "4-6", "condition": "one‐pot reductive amination"},  # one‐pot imine reduction
    {"solvent": "HCl", "state": "reflux", "catalyst": "Zn"},  # Clemmensen reduction
    {
        "Solvent": "Water and NaOH", 
        "Temperature" : "453.15K", 
        "Presence": "Hydrazine"
    }, #wolfkischner rxn
    {
        "option1": "acidic medium, water solvent: proton adds to carbonyl → enolization",
        "option2": "basic medium: enolate forms via base attack",
        "products": "racemic"
    },  # aldol addition
    {
        "option1": "acidic medium, water solvent: proton adds to carbonyl → enol → dehydration",
        "option2": "basic medium: enolate → condensation → dehydration",
        "products": "racemic"
    },  # aldol condensation

    {
        "reductor": "hydride : NaBH4 (reduces aldehydes but not ketones)", 
        "Temperature": "195.15 K", 
        "Solvent": "EtOH/DCM 1:1"
    }, #aldehyde reduction but not ketone 
    {
        "presence" : "CN- , H2O", 
        "Solvent" : "EtOH"
    }, #benzoin rxn
    {
        "presence" : "CN- , H2O", 
        "Solvent" : "EtOH"
    }, #cross benzoin rxn
    {""}, #grignard primary 
    {""}, #grignard secondary 
    {""}, #grignard tertiary 
    {"note": "according to markovnikov rules, the hallide will be placed on the most substituted position to form a more stable carbocation"}, #Chloride addition to alkene 
    {"note": "according to markovnikov rules, the hallide will be placed on the most substituted position to form a more stable carbocation"}, #Brooride addition to alkene 
    {"note": "according to markovnikov rules, the hallide will be placed on the most substituted position to form a more stable carbocation"}, #Iodide addition to alkene 
    {"note": "according to markovnikov rules, the hallide will be placed on the most substituted position to form a more stable carbocation"}, #Fluoride addition to alkene
    {"note": "according to markovnikov rules, the hydroxyl will be placed on the most substituted position to form a more stable carbocation"}, #hydration of alkene
    {"note": "not a radical addition but an ionic one"}, #hallide addition to alkene 
    {"note": "not a radical addition but an ionic one"}, #hallide addition to alkene 
    {"note": "not a radical addition but an ionic one"}, #hallide addition to alkene 
    {"catalyst" : "m-CPBA"}, #epoxydation
    {"solvent": "water" , "catalyst":"Hg+"}, #alkyne hydration
    {"medium":"base", "pH":"< 30"}, #alkyne+base --> allene pH<30 
    {"medium":"base", "pH":"> 30-35"}, #alkyne+base --> alkyne pushed pH>30-35
    {"catalyst":"Pd"}, #lindlar hydrogenation alkynes 


]

REACTION_REVERSERS: List[Callable[[str], str | None]] = [
    reverse_reaction_generator(i) for i in REACTION_DATABASE
]

def list_reactants(smiles: str)->list[str]:
    return list(filter(lambda x: x is not None, [fn(smiles) for fn in REACTION_REVERSERS]))