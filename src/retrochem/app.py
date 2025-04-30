import streamlit as st
from retrochem import hello_smiles 
import rdkit  
import re

allowed = re.compile(r"^[A-Za-z0-9\- ]+$")   # simple pattern

st.title("RetroChem - Your Organic Chemistry guide")
molecule = st.text_input("Enter a molecule")

if not molecule.strip():          # empty after removing spaces
        st.warning("⚠️ Please enter a molecule name (text).")
else:
        st.success(f"Got it: {molecule}")



if st.button("Retrosynthesis"):
    st.success( hello_smiles(molecule) )