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

# ─── RERUN COMPATIBILITY ─────────────────────────────────────────────────────
def rerun():
    if hasattr(st, "rerun"):
        st.rerun()

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

# ─── DATABASE REFRESH ─────────────────────────────────────────────────────────
def refresh_databases():
    for path in os.listdir('.'):
        if path.endswith('.db'):
            db = rd.load_database(path)
            if db is None:
                st.warning(f"⚠️ Could not load {path}")
            else:
                rd.register_database(db, path.removesuffix('.db'))

# Auto-load on startup
if not rd.REACTION_DATABASES:
    refresh_databases()

# ─── CALLBACKS ─────────────────────────────────────────────────────────────────
def reset_all():
    """Clear all retrosynthesis state for a fresh start."""
    for k in ['selected_smiles', 'reactant_list', 'database', 'combos', 'history']:
        st.session_state[k] = [] if k in ('combos', 'history') else None

def go_back():
    """Return to the previous molecule in history."""
    if st.session_state.history:
        prev = st.session_state.history.pop()
        st.session_state.selected_smiles = prev
        st.session_state.reactant_list = None
        canon = canonicalize_smiles(prev)
        st.session_state.combos = rd.list_reactants(canon, st.session_state.database)
    else:
        st.warning("⚠️ No previous molecule to go back to.")

def start_retro(smi):
    """Run retrosynthesis on a new target, preserving history."""
    db = st.session_state.database
    if not smi:
        st.warning("⚠️ Provide a molecule first.")
        return
    if not db:
        st.warning("⚠️ Choose a database first.")
        return
    # if we're already on a target, push it to history
    if st.session_state.selected_smiles:
        st.session_state.history.append(st.session_state.selected_smiles)
    # set up new target
    st.session_state.selected_smiles = smi
    st.session_state.reactant_list = None
    canon = canonicalize_smiles(smi)
    st.session_state.combos = rd.list_reactants(canon, db)

def choose_combo(idx):
    combo = st.session_state.combos[idx][0]
    st.session_state.reactant_list = combo.split('.')

def choose_reactant(part_smi):
    db = st.session_state.database
    st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = part_smi
    st.session_state.reactant_list = None
    canon = canonicalize_smiles(part_smi)
    st.session_state.combos = rd.list_reactants(canon, db)

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
        rerun()
    st.stop()

if st.session_state.page == "main":
# ─── MAIN APP TITLE ───────────────────────────────────────────────────────────
    st.title("RetroChem - Your Organic Chemistry Guide")

# ─── SIDEBAR: GLOBAL CONTROLS & STATUS ────────────────────────────────────────
    with st.sidebar:
        st.button("🧹 Start Over", on_click=reset_all)
        if st.session_state.history:
            st.button("🔙 Back", on_click=go_back)
        st.write("---")
        st.header("🔍 Retrosynthesis Status")
        st.markdown(f"**Database:** {st.session_state.database or '_(none selected)_'}")
        dbs = list(rd.REACTION_DATABASES.keys())
        dbs_display = dbs + ["➕ Add or edit your own database"]
        if dbs_display:
            chosen = st.selectbox(
                "Choose Database", dbs_display,
                index=dbs_display.index(st.session_state.database) if st.session_state.database in dbs_display else 0,
                key="db_selector",
            )
            if chosen == "➕ Add or edit your own database":
                st.session_state.page = "builder"
                rerun()
            else:
                st.session_state.database = chosen

        st.write("---")
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

# ─── INPUT MOLECULE SECTION ───────────────────────────────────────────────────
    st.markdown("---")
    st.markdown("## 🧬 Input Molecule")
    st.caption("Draw your molecule or enter its IUPAC name here, then click Retrosynthesize.")
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

# ─── RETROSYNTHESIS NAVIGATION ────────────────────────────────────────────────
    st.markdown("---")
    st.markdown("## 🔄 Retrosynthesis")
    # always allow retrosynthesis on the current input
    if smiles_input:
        st.button("🔄 Retrosynthesize", on_click=start_retro, args=(smiles_input,))

# show options or next fragments
    if st.session_state.selected_smiles and st.session_state.reactant_list is None:
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
                    df = pd.DataFrame(
                        {"Conditions": list(conds[i].values())},
                        index=[k.capitalize() for k in conds[i].keys()],
                    )
                    st.table(df)

elif st.session_state.reactant_list:
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
                st.button(f"Reactant {j+1}", on_click=choose_reactant, args=(p,))

# ─── DATABASE BUILDER ─────────────────────────────────────────────────────────

if st.session_state.page == "builder":
    st.title("🛠️ Custom Database Builder")

# ─── Session State Init ─────────────────────────────────────────────
    if "builder_reactants" not in st.session_state:
        st.session_state.builder_reactants = []
    if "builder_product" not in st.session_state:
        st.session_state.builder_product = None

    db_name = st.text_input("Database Name (will be saved as .db if new, or edit if already existing)")

    st.markdown("## ➕ Molecule Entry")
    input_mode = st.radio("Input Molecule By:", ["Name", "Draw"], horizontal=True)
    role = st.radio("Molecule Role", ["Product", "Reactant"], horizontal=True)

    mol_smiles = None
    if input_mode == "Name":
        name = st.text_input("Molecule name")
        if name:
            try:
                mol_smiles = name_to_smiles(name)
                st.success(f"SMILES: {mol_smiles}")
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    st.image(MolToImage(mol, size=(200, 200)))
            except Exception as e:
                st.error(f"Error: {e}")
    else:
        drawn = st_ketcher("", key="builder_draw")
        if drawn:
            mol_smiles = drawn
            st.success(f"SMILES: {mol_smiles}")
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                st.image(MolToImage(mol, size=(200, 200)))

    if st.button("✅ Add Molecule to Reaction"):
        if not mol_smiles:
            st.warning("No molecule to add.")
        elif role == "Reactant":
            st.session_state.builder_reactants.append(mol_smiles)
        else:
            st.session_state.builder_product = mol_smiles

# ─── Preview Reaction ───────────────────────────────────────────────
    if st.session_state.builder_product or st.session_state.builder_reactants:
        st.markdown("## 👁️ Reaction Preview")

        cols = st.columns(len(st.session_state.builder_reactants) + 2)
        for i, smi in enumerate(st.session_state.builder_reactants):
            with cols[i]:
                st.image(MolToImage(Chem.MolFromSmiles(smi), size=(150, 150)))
                st.caption("Reactant")

        with cols[len(st.session_state.builder_reactants)]:
            st.markdown("### ➡️")

        if st.session_state.builder_product:
            with cols[-1]:
                st.image(MolToImage(Chem.MolFromSmiles(st.session_state.builder_product), size=(150, 150)))
                st.caption("Product")
        else:
            with cols[-1]:
                st.info("Product not yet defined")

# ─── Conditions ─────────────────────────────────────────────────────
    st.markdown("### ⚙️ Conditions")
    conds = st.text_area("Conditions (key: value per line)")
    parsed_conds = {}
    if conds:
        for line in conds.splitlines():
            if ':' in line:
                k, v = map(str.strip, line.split(':', 1))
                parsed_conds[k] = v

# ─── Final Add Button ───────────────────────────────────────────────
    if st.button("🧪 Save Reaction to Database"):
        if not db_name:
            st.error("Please enter a database name.")
        elif not st.session_state.builder_product:
            st.error("No product specified.")
        elif not st.session_state.builder_reactants:
            st.error("No reactants specified.")
        else:
            try:
                rd.add_new_smart(
                    db_name,
                    st.session_state.builder_product,
                    st.session_state.builder_reactants,
                    parsed_conds
                )
                refresh_databases()
                st.success("✅ Reaction saved!")

                # Reset builder
                st.session_state.builder_product = None
                st.session_state.builder_reactants = []
            except Exception as e:
                st.error(f"Failed to save: {e}")

# ─── Navigation ─────────────────────────────────────────────────────
    if st.button("⬅️ Back to Main App"):
        st.session_state.page = "main"
        rerun()