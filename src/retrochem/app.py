import streamlit as st
from streamlit_ketcher import st_ketcher  # type: ignore
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

from retrochem.functions import name_to_smiles, canonicalize_smiles 
import retrochem.reaction_database as rd

# ─── 1) PAGE CONFIG & TITLE ──────────────────────────────────────────────────
st.set_page_config(
    page_title="RetroChem",
    layout="wide",
    page_icon="🧪",
)
st.title("RetroChem - Your Organic Chemistry Guide")
st.subheader("🔬 Input Molecule")
mode = st.radio("Choose input mode:", ["Name", "Draw structure"], horizontal=True)

smiles_input = ""
# ─── 2) INPUT SECTION ─────────────────────────────────────────────────────────
if mode == "Name":
    st.subheader("🔬 Input Molecule by Name")
    raw_name = st.text_input("Name of molecule", key="name_input")
    if raw_name:
        try:
            smiles_input = name_to_smiles(raw_name)
            st.write("SMILES:", smiles_input)
            mol = Chem.MolFromSmiles(smiles_input)
            if mol:
                img = MolToImage(mol, size=(200, 200))
                st.image(img)
        except Exception as e:
            st.error(f"❌ Name → SMILES conversion failed:\n{e}")


else:
    st.subheader("🔬 Draw Structure")
    smiles_input = st_ketcher("", height=450)
    if smiles_input:
        st.write("SMILES:", smiles_input)
        mol2 = Chem.MolFromSmiles(smiles_input)
        if mol2:
            img2 = MolToImage(mol2, size=(200, 200))
            st.image(img2)

st.markdown("---")  # divider before run

# ─── 3) RUN & OUTPUT ───────────────────────────────────────────────────────────
if st.button("🔄 Retrosynthesis"):
    # Choose input preference: drawn > named
    if not smiles_input:
        st.warning("⚠️ Please provide a molecule via name or drawing.")
    else:
        canon = canonicalize_smiles(smiles_input)
        result = rd.list_reactants(canon)
        st.markdown("---")
        st.header("🧩 Retrosynthesis Output")
        if not result:
            st.info("No disconnection rules matched.")
        for products in result:
            st.image(MolToImage(Chem.MolFromSmiles(products)))