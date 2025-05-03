import streamlit as st
from streamlit_ketcher import st_ketcher  # type: ignore

from retrochem.functions import name_to_smiles, canonicalize_smiles 
import reaction_database as rd

st.set_page_config(
    page_title="RetroChem",
    layout="wide",
    page_icon="🧪",
)

# ─── 2) MAIN ──────────────────────────────────────────────────────────────────
st.title("RetroChem - Your Organic Chemistry Guide")

# Input mode selector at top (no sidebar)
mode = st.radio("🔬 Select input mode:", ["Name", "Draw structure"], horizontal=True)

smiles_input = ""
run = False

if mode == "Name":
    st.subheader("Enter a molecule name")
    raw_name = st.text_input("Name of molecule")
    if st.button("🔄 Retrosynthesis"):
        run = True
        if raw_name:
            try:
                smiles_input = name_to_smiles(raw_name)
            except Exception as e:
                st.error(f"❌ Name → SMILES conversion failed:\n{e}")
        else:
            st.warning("⚠️ Please enter a molecule name.")
else:
    st.subheader("Draw your molecule")
    sk_smiles = st_ketcher("", height=450)
    if st.button("🔄 Retrosynthesis"):
        run = True
        if sk_smiles:
            smiles_input = sk_smiles
        else:
            st.warning("⚠️ Please draw a molecule and click Apply.")

# ─── 4) RUN & OUTPUT ───────────────────────────────────────────────────────────
if run and smiles_input:
    try:
        # Canonicalize & run your backend
        canon = canonicalize_smiles(smiles_input)
    except Exception as e:
        st.error(f"❌ Error during retrosynthesis:\n{e}")
    else:
        st.markdown("---")
        st.header("🧩 Retrosynthesis Output")
