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
st.title("RetroChem - Your Organic Chemistry Guide")

# Initialize session state
if 'reactant_list' not in st.session_state:
    st.session_state['reactant_list'] = None
if 'selected_smiles' not in st.session_state:
    st.session_state['selected_smiles'] = None

# Display currently selected reactant for next retro step
if st.session_state['selected_smiles']:
    st.subheader("🔄 Current molecule for retrosynthesis")
    sel = st.session_state['selected_smiles']
    st.write("SMILES:", sel)
    mol_sel = Chem.MolFromSmiles(sel)
    if mol_sel:
        st.image(MolToImage(mol_sel, size=(200,200)))
    st.markdown("---")

# ─── INPUT MODE SELECTION ─────────────────────────────────────────────────────
st.subheader("🔬 Input Molecule")
mode = st.radio("Choose input mode:", ["Name", "Draw structure"], horizontal=True)

smiles_input = ""

# ─── NAME INPUT ───────────────────────────────────────────────────────────────
if mode == "Name":
    st.subheader("By Name")
    raw_name = st.text_input("Name of molecule", key="name_input")
    if raw_name:
        try:
            smiles_input = name_to_smiles(raw_name)
            st.write("SMILES:", smiles_input)
            mol = Chem.MolFromSmiles(smiles_input)
            if mol:
                st.image(MolToImage(mol, size=(200, 200)))
        except Exception as e:
            st.error(f"❌ Name → SMILES conversion failed:\n{e}")

# ─── DRAW STRUCTURE INPUT ─────────────────────────────────────────────────────
else:
    st.subheader("Draw Structure")
    drawn = st_ketcher("", height=450, key="draw_input")
    if drawn:
        smiles_input = drawn
        st.write("SMILES:", smiles_input)
        mol2 = Chem.MolFromSmiles(smiles_input)
        if mol2:
            st.image(MolToImage(mol2, size=(200, 200)))

st.markdown("---")  # divider before run

# ─── RUN & OUTPUT ────────────────────────────────────────────────────────────
if st.button("🔄 Retrosynthesis", key="retro_btn"):
    if not smiles_input:
        st.warning("⚠️ Please provide a molecule via the selected input mode.")
    else:
        # clear previous selections
        st.session_state['reactant_list'] = None
        st.session_state['selected_smiles'] = None
        canon = canonicalize_smiles(smiles_input)
        combos = rd.list_reactants(canon)
        st.subheader("🧩 Retrosynthesis Options")
        if not combos:
            st.info("No disconnection rules matched.")
        else:
            # Display each combination option
            cols = st.columns(len(combos))
            for i, combo_smiles in enumerate(combos):
                with cols[i]:
                    mol_prod = Chem.MolFromSmiles(combo_smiles)
                    if mol_prod:
                        st.image(MolToImage(mol_prod, size=(200,200)))
                    if st.button(f"Option {i+1}", key=f"opt_{i}"):
                        # split into individual reactant SMILES
                        parts = combo_smiles.split('.')
                        st.session_state['reactant_list'] = parts

# ─── REACTANT SELECTION ───────────────────────────────────────────────────────
if st.session_state['reactant_list']:
    st.markdown("---")
    st.subheader("🔹 Select a reactant to continue retrosynthesis")
    cols2 = st.columns(len(parts))
    for j, part_smiles in enumerate(parts):
        with cols2[j]:
            st.write(part_smiles)
            mol_part = Chem.MolFromSmiles(part_smiles)
            if mol_part:
                st.image(MolToImage(mol_part, size=(200,200)))
            if st.button(f"Select {j+1}", key=f"sel_{j}"):
                st.session_state['selected_smiles'] = part_smiles
    


# # ─── 1) PAGE CONFIG & TITLE ──────────────────────────────────────────────────
# st.set_page_config(
#     page_title="RetroChem",
#     layout="wide",
#     page_icon="🧪",
# )
# st.title("RetroChem - Your Organic Chemistry Guide")
# st.subheader("🔬 Input Molecule")
# mode = st.radio("Choose input mode:", ["Name", "Draw structure"], horizontal=True)

# smiles_input = ""
# # ─── 2) INPUT SECTION ─────────────────────────────────────────────────────────
# if mode == "Name":
#     st.subheader("🔬 Input Molecule by Name")
#     raw_name = st.text_input("Name of molecule", key="name_input")
#     if raw_name:
#         try:
#             smiles_input = name_to_smiles(raw_name)
#             st.write("SMILES:", smiles_input)
#             mol = Chem.MolFromSmiles(smiles_input)
#             if mol:
#                 img = MolToImage(mol, size=(200, 200))
#                 st.image(img)
#         except Exception as e:
#             st.error(f"❌ Name → SMILES conversion failed:\n{e}")


# else:
#     st.subheader("🔬 Draw Structure")
#     smiles_input = st_ketcher("", height=450)
#     if smiles_input:
#         st.write("SMILES:", smiles_input)
#         mol2 = Chem.MolFromSmiles(smiles_input)
#         if mol2:
#             img2 = MolToImage(mol2, size=(200, 200))
#             st.image(img2)

# st.markdown("---")  # divider before run

# # ─── 3) RUN & OUTPUT ───────────────────────────────────────────────────────────
# if st.button("🔄 Retrosynthesis"):
#     # Choose input preference: drawn > named
#     if not smiles_input:
#         st.warning("⚠️ Please provide a molecule via name or drawing.")
#     else:
#         canon = canonicalize_smiles(smiles_input)
#         result = rd.list_reactants(canon)
#         st.markdown("---")
#         st.header("🧩 Retrosynthesis Output")
#         if not result:
#             st.info("No disconnection rules matched.")
#         for products in result:
#             st.image(MolToImage(Chem.MolFromSmiles(products)))