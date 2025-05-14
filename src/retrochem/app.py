import streamlit as st
import pandas as pd
import os
from streamlit_ketcher import st_ketcher  # type: ignore
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

from retrochem.functions import name_to_smiles, canonicalize_smiles
import retrochem.reaction_database as rd

# ─── PAGE CONFIG & CUSTOM CSS THEME ─────────────────────────────────────────
st.set_page_config(
    page_title="RetroChem",
    layout="wide",
    page_icon="🧪",
)
# Inline CSS: chem-green theme
st.markdown(
    """
    <style>
    h1, h2, h3 { color: #2E7D32 !important; text-align: center; }
    .stButton>button {
        background-color: #388E3C !important;
        border-color: #2E7D32 !important;
        color: white !important;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ─── SESSION STATE SETUP ──────────────────────────────────────────────────────
for key, default in [
    ("page", "home"),
    ("selected_smiles", None),
    ("reactant_list", None),
    ("database", None),
    ("combos", []),
    ("history", []),
]:
    if key not in st.session_state:
        st.session_state[key] = default

# ─── LANDING PAGE ─────────────────────────────────────────────────────────────
if st.session_state.page == "home":
    st.title("🧪 Welcome to RetroChem")
    st.markdown("Your Organic Chemistry Retrosynthesis Assistant")
    st.write("---")
    st.markdown(
        "Click below to begin drawing or naming your target molecule, then choose your database and start retrosynthesis."
    )
    if st.button("🔬 Start Retrosynthesis"):
        st.session_state.page = "main"
        st.experimental_rerun()
    st.stop()

# ─── MAIN APP TITLE ───────────────────────────────────────────────────────────
st.title("RetroChem - Your Organic Chemistry Guide")

# ─── CALLBACKS ─────────────────────────────────────────────────────────────────

def reset_all():
    """Clear all retrosynthesis state for a fresh start."""
    for k in ["selected_smiles", "reactant_list", "database", "combos", "history"]:
        st.session_state[k] = [] if k in ("combos", "history") else None


def start_retro(smi, db):
    if not smi:
        st.warning("⚠️ Provide a molecule first.")
        return
    if not db:
        st.warning("⚠️ Choose a database first.")
        return
    reset_all()
    st.session_state.selected_smiles = smi
    canon = canonicalize_smiles(smi)
    st.session_state.combos = rd.list_reactants(canon, db)


def choose_combo(idx):
    combo = st.session_state.combos[idx][0]
    st.session_state.reactant_list = combo.split('.')


def refresh_databases():
    for path in os.listdir('.'):
        if path.endswith('.db'):
            db = rd.load_database(path)
            if db is None:
                st.warning(f"⚠️ Could not load {path}")
            rd.register_database(db, path.removesuffix('.db'))


def choose_database(db):
    st.session_state.database = db


def choose_reactant(part_smi, db):
    st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = part_smi
    st.session_state.reactant_list = None
    canon = canonicalize_smiles(part_smi)
    st.session_state.combos = rd.list_reactants(canon, db)

# ─── AUTO-LOAD DATABASES ON STARTUP ───────────────────────────────────────────
if not rd.REACTION_DATABASES:
    refresh_databases()

# ─── SIDEBAR: GLOBAL STATUS & CONTROLS ─────────────────────────────────────────
with st.sidebar:
    st.button("🧹 Start Over", on_click=reset_all)
    st.write("---")
    st.header("🔍 Retrosynthesis Status")
    st.markdown(f"**Database:** {st.session_state.database or '_(none selected)_'}")
    if st.session_state.selected_smiles:
        st.markdown(f"**Molecule:** `{st.session_state.selected_smiles}`")
        mol = Chem.MolFromSmiles(st.session_state.selected_smiles)
        if mol:
            st.image(MolToImage(mol, size=(150, 150)))
    else:
        st.markdown("**Molecule:** _(none entered)_")
    if st.session_state.history:
        st.write("---")
        st.subheader("History")
        for prev in reversed(st.session_state.history):
            st.write(prev)

# ─── DATABASE SELECTION ───────────────────────────────────────────────────────
st.markdown("---")
st.markdown("## 🔍 Database Selection")
st.caption("Select which reaction database to use for your retrosynthesis.")
col1, col2 = st.columns([1, 3])
with col1:
    st.button("Load/Refresh DBs", on_click=refresh_databases)
with col2:
    dbs = list(rd.REACTION_DATABASES.keys())
    if dbs:
        cols = st.columns(len(dbs))
        for i, d in enumerate(dbs):
            with cols[i]:
                st.button(d, on_click=choose_database, args=(d,))
    else:
        st.info("No .db files found.")

# ─── INPUT MOLECULE SECTION ───────────────────────────────────────────────────
st.markdown("---")
st.markdown("## 🧬 Input Molecule")
st.caption("Draw your molecule or enter its IUPAC name here, then click Start below.")
mode = st.radio("Mode:", ["Name", "Draw"], horizontal=True)
smiles_input = None
if mode == "Name":
    name = st.text_input("Molecule name")
    if name:
        try:
            smiles_input = name_to_smiles(name)
            st.write("SMILES:", smiles_input)
            mol0 = Chem.MolFromSmiles(smiles_input)
            if mol0:
                st.image(MolToImage(mol0, size=(200, 200)))
        except Exception as e:
            st.error(f"❌ Name→SMILES failed: {e}")
else:
    drawn = st_ketcher("", key="draw_input")
    if drawn:
        smiles_input = drawn
        st.write("SMILES:", smiles_input)
        mol0 = Chem.MolFromSmiles(smiles_input)
        if mol0:
            st.image(MolToImage(mol0, size=(200, 200)))

# ─── RETROSYNTHESIS TRIGGER ───────────────────────────────────────────────────
st.markdown("---")
st.markdown("## 🔄 Retrosynthesis")
st.caption("When you’re ready, press to compute possible reactants.")
if st.session_state.selected_smiles is None:
    st.button("🔄 Start", on_click=start_retro, args=(smiles_input, st.session_state.database))

# ─── RETROSYNTHESIS OPTIONS ──────────────────────────────────────────────────
elif st.session_state.reactant_list is None:
    st.markdown("## 🧩 Retrosynthesis Options")
    st.caption("Select one of the reactant options below.")
    combos = st.session_state.combos
    if not combos:
        st.info("No routes found.")
    else:
        conds = [c[1] for c in combos]
        smis = [c[0] for c in combos]
        cols = st.columns(len(smis))
        for i, smi in enumerate(smis):
            with cols[i]:
                m = Chem.MolFromSmiles(smi)
                if m:
                    st.image(MolToImage(m, size=(200, 200)))
                st.button(f"Option {i+1}", on_click=choose_combo, args=(i,))
                df = pd.DataFrame({"Conditions": list(conds[i].values())}, index=[k.capitalize() for k in conds[i].keys()])
                st.table(df)

# ─── NEXT REACTANT ────────────────────────────────────────────────────────────
else:
    st.markdown("## 🔹 Next Reactant")
    st.caption("Choose a fragment for further retrosynthesis.")
    parts = st.session_state.reactant_list
    cols = st.columns(len(parts))
    for j, p in enumerate(parts):
        with cols[j]:
            st.write(p)
            m = Chem.MolFromSmiles(p)
            if m:
                st.image(MolToImage(m, size=(200, 200)))
            st.button(f"Reactant {j+1}", on_click=choose_reactant, args=(p, st.session_state.database))
