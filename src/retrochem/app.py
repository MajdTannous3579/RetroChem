import streamlit as st
import pandas as pd
from streamlit_ketcher import st_ketcher  # type: ignore
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

from retrochem.functions import name_to_smiles, canonicalize_smiles
import retrochem.reaction_database as rd

# ─── PAGE CONFIG & TITLE ─────────────────────────────────────────────────────
st.set_page_config(
    page_title="RetroChem",
    layout="wide",
    page_icon="🧪",
)
st.title("RetroChem – Your Organic Chemistry Guide")

# ─── RESET CALLBACK ───────────────────────────────────────────────────────────
def reset_all():
    """Clear retrosynthesis state to start fresh"""
    st.session_state.selected_smiles = None
    st.session_state.reactant_list = None
    st.session_state.combos = []
    st.session_state.history = []

# Show a Start Over button at the top
st.button(
    "🧹 Start Over",
    key="reset",
    on_click=reset_all,
)

# ─── OTHER CALLBACKS ──────────────────────────────────────────────────────────
def start_retro(smi):
    if not smi:
        st.warning("⚠️ Provide a molecule first.")
        return
    reset_all()  # clear any previous state
    st.session_state.selected_smiles = smi
    canon = canonicalize_smiles(smi)
    st.session_state.combos = rd.list_reactants(canon)

def choose_combo(idx):
    combo = [i[0] for i in st.session_state.combos][idx] #changed the variable "combos" that is here a tuple into the "true" combo like it is defined in the database file, kept it like this until here to keep track of history for both conditions and smarts and not have to separate them from the beginning
    st.session_state.reactant_list = combo.split(".")

def choose_reactant(part_smi):
    st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = part_smi
    st.session_state.reactant_list = None
    canon = canonicalize_smiles(part_smi)
    st.session_state.combos = rd.list_reactants(canon)

# ─── SESSION STATE SETUP ───────────────────────────────────────────────────────
if 'selected_smiles' not in st.session_state:
    st.session_state.selected_smiles = None
if 'reactant_list' not in st.session_state:
    st.session_state.reactant_list = None
if 'combos' not in st.session_state:
    st.session_state.combos = []
if 'history' not in st.session_state:
    st.session_state.history = []

# ─── DISPLAY CURRENT TARGET & BACK BUTTON ────────────────────────────────────
if st.session_state.selected_smiles:
    st.subheader("🔄 Current molecule for retrosynthesis")
    smi = st.session_state.selected_smiles
    st.write("SMILES:", smi)
    mol = Chem.MolFromSmiles(smi)
    if mol:
        st.image(MolToImage(mol, size=(200,200)))

    if st.session_state.history:
        if st.button("⬅️ Back", key="back"):
            prev = st.session_state.history.pop()
            st.session_state.selected_smiles = prev
            st.session_state.reactant_list = None
            canon_prev = canonicalize_smiles(prev)
            st.session_state.combos = rd.list_reactants(canon_prev)
    st.markdown("---")

# ─── INPUT MODE (BEFORE FIRST RETRO) ──────────────────────────────────────────
st.subheader("🔬 Input Molecule")
mode = st.radio("Choose input mode:", ["Name", "Draw structure"], horizontal=True)
smiles_input = ""

if mode == "Name":
    raw_name = st.text_input("Name of molecule", key="name_input")
    if raw_name:
        try:
            smiles_input = name_to_smiles(raw_name)
            st.write("SMILES:", smiles_input)
            mol0 = Chem.MolFromSmiles(smiles_input)
            if mol0:
                st.image(MolToImage(mol0, size=(200,200)))
        except Exception as e:
            st.error(f"❌ Name → SMILES failed: {e}")
else:
    drawn = st_ketcher("", height=400, key="draw_input")
    if drawn:
        smiles_input = drawn
        st.write("SMILES:", smiles_input)
        mol0 = Chem.MolFromSmiles(smiles_input)
        if mol0:
            st.image(MolToImage(mol0, size=(200,200)))

st.markdown("---")

# ─── INITIAL RETROSYNTHESIS TRIGGER ───────────────────────────────────────────
if st.session_state.selected_smiles is None:
    st.button(
        "🔄 Retrosynthesis",
        key="start",
        on_click=start_retro,
        args=(smiles_input,),
    )

# ─── SHOW RETROSYNTHESIS OPTIONS ─────────────────────────────────────────────
elif st.session_state.reactant_list is None:
    st.subheader("🧩 Retrosynthesis Options")
    combos = st.session_state.combos
    if not combos:
        st.info("Your product is too simple to be retrosynthesized, or not in our database")
    else:
        cond = [i[1] for i in combos]
        combos = [i[0] for i in combos] #changed the variable "combos" that is here a tuple into the "true" combo like it is defined in the database file, kept it like this until here to keep track of history for both conditions and smarts and not have to separate them from the beginning
        cols = st.columns(len(combos))
        for i, combo in enumerate(combos):
            with cols[i]:
                mol = Chem.MolFromSmiles(combo)
                if mol:
                    st.image(MolToImage(mol, size=(200,200)))
                st.button(
                    f"Option {i+1}",
                    key=f"opt_{i}",
                    on_click=choose_combo,
                    args=(i,),
                )
                keys   = list([i.capitalize() for i in cond[i].keys()])
                values = list(cond[i].values())
                df = pd.DataFrame({
                    "Conditions": values,
                }, index=keys)
                st.table(df)

# ─── SPLIT & SELECT NEXT REACTANT ─────────────────────────────────────────────
else:
    st.subheader("🔹 Select a reactant to continue")
    parts = st.session_state.reactant_list
    cols2 = st.columns(len(parts))
    for j, part in enumerate(parts):
        with cols2[j]:
            st.write(part)
            mol_p = Chem.MolFromSmiles(part)
            if mol_p:
                st.image(MolToImage(mol_p, size=(200,200)))
            st.button(
                f"Retrosynthesize Reactant {j+1}",
                key=f"sel_{j}",
                on_click=choose_reactant,
                args=(part,),
            )
