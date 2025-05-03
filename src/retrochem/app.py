import streamlit as st
from streamlit_ketcher import st_ketcher  # type: ignore

from retrochem.functions import name_to_smiles, canonicalize_smiles 

# ─── 1) PAGE CONFIG ───────────────────────────────────────────────────────────
st.set_page_config(
    page_title="RetroChem",
    layout="wide",     # give the main area full width
    page_icon="🧪",
)

# ─── 2) SIDEBAR ────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🔬 Input Mode")
    mode = st.radio("", ["Name", "Draw structure"])

# ─── 3) MAIN ──────────────────────────────────────────────────────────────────
st.title("RetroChem - Your Organic Chemistry Guide")

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
    # Full‐width sketcher; height chosen so Apply is visible
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
