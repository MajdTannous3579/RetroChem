import streamlit as st
from streamlit_ketcher import st_ketcher  # type: ignore

# Backend functions
from retrochem.functions import name_to_smiles, canonicalize_smiles
from retrochem.example_module import hello_smiles

# Page config
st.set_page_config(page_title="RetroChem", layout="centered")
st.title("RetroChem - Your Organic Chemistry Guide")

# 1) Input choice
input_mode = st.radio("Choose input method:", ["Name", "Draw structure"])

smiles_input = ""
if input_mode == "Name":
    raw_name = st.text_input("Enter a molecule name")
    if raw_name:
        try:
            smiles_input = name_to_smiles(raw_name)
            st.markdown(f"**SMILES:** `{smiles_input}`")
        except Exception as e:
            st.error(f"❌ Could not convert name to SMILES: {e}")
else:
    st.markdown("### Draw your molecule:")
    # st_ketcher(initial_smiles: str = "", height: int = 600, molecule_format: str = "SMILES")
    ketcher_smiles = st_ketcher("", height=400)
    if ketcher_smiles:
        smiles_input = ketcher_smiles
        st.markdown(f"**Detected SMILES:** `{smiles_input}`")

# 2) Retrosynthesis
if st.button("Retrosynthesis"):
    if not smiles_input:
        st.warning("⚠️ Please provide a molecule (by name or drawing).")
    else:
        try:
            canon = canonicalize_smiles(smiles_input)
            result = hello_smiles(canon)
        except Exception as e:
            st.error(f"❌ Error during retrosynthesis: {e}")
        else:
            st.success(result)
