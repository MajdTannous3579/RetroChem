import streamlit as st
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

# ─── SESSION STATE INITIALIZATION ────────────────────────────────────────────
if 'selected_smiles' not in st.session_state:
    st.session_state['selected_smiles'] = None
if 'reactant_list' not in st.session_state:
    st.session_state['reactant_list'] = None
if 'combos' not in st.session_state:
    st.session_state['combos'] = []
if 'history' not in st.session_state:
    st.session_state['history'] = []

# ─── DISPLAY CURRENT TARGET & BACK BUTTON ────────────────────────────────────
if st.session_state.selected_smiles:
    st.subheader("🔄 Current molecule for retrosynthesis")
    smi = st.session_state.selected_smiles
    st.write("SMILES:", smi)
    mol = Chem.MolFromSmiles(smi)
    if mol:
        st.image(MolToImage(mol, size=(200,200)))

    # Back button
    if st.session_state.history:
        if st.button("⬅️ Back"):
            prev = st.session_state.history.pop()
            st.session_state.selected_smiles = prev
            st.session_state.reactant_list = None
            # recompute options for prev
            canon = canonicalize_smiles(prev)
            st.session_state.combos = rd.list_reactants(canon)
    st.markdown("---")

# ─── INPUT MODE (ONLY BEFORE FIRST RETRO) ────────────────────────────────────
st.subheader("🔬 Input Molecule")
mode = st.radio("Choose input mode:", ["Name", "Draw structure"], horizontal=True)
smiles_input = ""

if mode == "Name":
    raw_name = st.text_input("Name of molecule")
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
    drawn = st_ketcher("", height=400)
    if drawn:
        smiles_input = drawn
        st.write("SMILES:", smiles_input)
        mol0 = Chem.MolFromSmiles(smiles_input)
        if mol0:
            st.image(MolToImage(mol0, size=(200,200)))

st.markdown("---")

# ─── INITIAL RETROSYNTHESIS BUTTON ────────────────────────────────────────────
if st.session_state.selected_smiles is None:
    if st.button("🔄 Retrosynthesis"):
        if not smiles_input:
            st.warning("⚠️ Provide a molecule first.")
        else:
            # kick off the first retrosynthesis
            st.session_state.history = []
            st.session_state.selected_smiles = smiles_input
            st.session_state.reactant_list = None
            canon = canonicalize_smiles(smiles_input)
            st.session_state.combos = rd.list_reactants(canon)

# ─── SHOW RETROSYNTHESIS OPTIONS ──────────────────────────────────────────────
elif st.session_state.reactant_list is None:
    st.subheader("🧩 Retrosynthesis Options")
    combos = st.session_state.combos
    if not combos:
        st.info("No disconnection rules matched.")
    else:
        cols = st.columns(len(combos))
        for i, combo in enumerate(combos):
            with cols[i]:
                mol = Chem.MolFromSmiles(combo)
                if mol:
                    st.image(MolToImage(mol, size=(200,200)))
                if st.button(f"Option {i+1}"):
                    parts = combo.split(".")
                    st.session_state.reactant_list = parts

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
            if st.button(f"Retrosynthesize Reactant {j+1}"):
                # push current to history
                st.session_state.history.append(st.session_state.selected_smiles)
                # set new target
                st.session_state.selected_smiles = part
                st.session_state.reactant_list = None
                # compute its options
                canon = canonicalize_smiles(part)
                st.session_state.combos = rd.list_reactants(canon)
