import streamlit as st
import pandas as pd
import os
from streamlit_ketcher import st_ketcher  # type: ignore
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

from retrochem.preliminary_functions import name_to_smiles, canonicalize_smiles
import retrochem.reaction_database as rd

# ─── PAGE CONFIG & THEME ─────────────────────────────────────────────────────
st.set_page_config(page_title="RetroChem", layout="wide", page_icon="🧪")
st.markdown(
    """
    <style>
    h1, h2, h3 { color: #2E7D32 !important; text-align: center; }
    .stButton>button { background-color: #388E3C !important; border-color: #2E7D32 !important; color: white !important; }
    </style>
    """,
    unsafe_allow_html=True,
)

# ─── GLOBAL UTILS ─────────────────────────────────────────────────────────────
def rerun():
    st.experimental_rerun()

def select_database():
    sel = st.session_state.db_selector
    if sel == "➕ Add or edit your own database":
        st.session_state.page = "builder"
    else:
        st.session_state.database = sel
    # auto re-render on session_state change

# ─── SESSION STATE SETUP ──────────────────────────────────────────────────────
defaults = {
    "page": "home",
    "selected_smiles": None,
    "reactant_list": None,
    "database": None,
    "combos": [],
    "history": [],
    "builder_reactants": [],
    "builder_product": None,
    "builder_mol_name": "",
    "builder_conds": "",
    "builder_role": "Product",
    "builder_input_mode": "Name",
    "builder_draw": "",
}
for key, val in defaults.items():
    if key not in st.session_state:
        st.session_state[key] = val

# ─── LOAD DATABASES ──────────────────────────────────────────────────────────
def refresh_databases():
    """Scan for SQLite .db files and register them, skipping invalid or double-suffixed files silently."""
    for path in os.listdir('.'):
        # only load single-suffix .db files
        if not path.endswith('.db') or path.endswith('.db.db'):
            continue
        try:
            db = rd.load_database(path)
            if db:
                rd.register_database(db, path.removesuffix('.db'))
        except Exception:
            # skip any files that fail to load
            continue

# Auto-load on startup (only once)
if not rd.REACTION_DATABASES:
    refresh_databases()
if not rd.REACTION_DATABASES:
    refresh_databases()

# ─── RETROSYNTHESIS CALLBACKS ─────────────────────────────────────────────────
def reset_all():
    for k in ['selected_smiles','reactant_list','database','combos','history']:
        st.session_state[k] = [] if isinstance(st.session_state[k], list) else None

def go_back():
    if st.session_state.history:
        prev = st.session_state.history.pop()
        st.session_state.selected_smiles = prev
        st.session_state.reactant_list = None
        st.session_state.combos = rd.list_reactants(canonicalize_smiles(prev), st.session_state.database)
    else:
        st.warning("⚠️ No previous molecule to go back to.")

def start_retro(smi):
    if not smi:
        st.warning("⚠️ Provide a molecule first.")
        return
    if not st.session_state.database:
        st.warning("⚠️ Choose a database first.")
        return
    if st.session_state.selected_smiles:
        st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = smi
    st.session_state.reactant_list = None
    st.session_state.combos = rd.list_reactants(canonicalize_smiles(smi), st.session_state.database)

def choose_combo(idx):
    combo = st.session_state.combos[idx][0]
    st.session_state.reactant_list = combo.split('.')

def choose_reactant(part):
    st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = part
    st.session_state.reactant_list = None
    st.session_state.combos = rd.list_reactants(canonicalize_smiles(part), st.session_state.database)

# ─── PAGE: HOME ───────────────────────────────────────────────────────────────
if st.session_state.page == "home":
    st.title("🧪 Welcome to RetroChem")
    st.markdown("Your Organic Chemistry Retrosynthesis Assistant")
    st.write("---")
    st.markdown("Click below to begin drawing or naming your target, then choose your database and start retrosynthesis.")
    if st.button("🔬 Start Retrosynthesis"):
        st.session_state.page = "main"
        rerun()
    st.stop()

# ─── PAGE: MAIN ───────────────────────────────────────────────────────────────
elif st.session_state.page == "main":
    st.title("RetroChem - Your Organic Chemistry Guide")
    with st.sidebar:
        st.button("🧹 Start Over", on_click=reset_all)
        if st.session_state.history:
            st.button("🔙 Back", on_click=go_back)
        st.write("---")
        st.header("🔍 Retrosynthesis Status")
        dbs = list(rd.REACTION_DATABASES.keys())
        st.selectbox("Choose Database", dbs+ ["➕ Add or edit your own database"], key="db_selector", on_change=select_database)
        st.write(f"**Current DB:** {st.session_state.database or '_(none selected)_'}")
        st.write("---")
        if st.session_state.selected_smiles:
            st.markdown(f"**Molecule:** `{st.session_state.selected_smiles}`")
            mol = Chem.MolFromSmiles(st.session_state.selected_smiles)
            if mol: st.image(MolToImage(mol, size=(150,150)))
        if st.session_state.history:
            st.write("---")
            st.subheader("History")
            for prev in reversed(st.session_state.history): st.write(prev)

    st.markdown("---")
    st.markdown("## 🧬 Input Molecule")
    mode = st.radio("Mode:", ["Name","Draw"], horizontal=True)
    smiles = None
    if mode == "Name":
        name = st.text_input("Molecule name", key="main_name")
        if name:
            try:
                smiles = name_to_smiles(name)
                st.write("SMILES:", smiles)
                m = Chem.MolFromSmiles(smiles)
                if m: st.image(MolToImage(m, size=(200,200)))
            except Exception as e:
                st.error(f"❌ Name→SMILES failed: {e}")
    else:
        drawn = st_ketcher("", key="main_draw")
        if drawn:
            smiles = drawn
            st.write("SMILES:", smiles)
            m = Chem.MolFromSmiles(smiles)
            if m: st.image(MolToImage(m, size=(200,200)))

    st.markdown("---")
    if smiles and st.button("🔄 Retrosynthesize", on_click=start_retro, args=(smiles,)):
        pass
    if st.session_state.selected_smiles and st.session_state.reactant_list is None:
        st.markdown("## 🧩 Options")
        combos = st.session_state.combos
        if not combos: st.info("No routes found.")
        else:
            cols = st.columns(len(combos))
            for i,(smi, cond) in enumerate(combos):
                with cols[i]:
                    m = Chem.MolFromSmiles(smi)
                    if m: st.image(MolToImage(m, size=(200,200)))
                    st.button(f"Option {i+1}", on_click=choose_combo, args=(i,))
                    df = pd.DataFrame({'Conditions': list(cond.values())}, index=[k.capitalize() for k in cond.keys()])
                    st.table(df)
    elif st.session_state.reactant_list:
        st.markdown("## 🔹 Next Reactant")
        parts = st.session_state.reactant_list
        cols = st.columns(len(parts))
        for j,p in enumerate(parts):
            with cols[j]:
                st.write(p)
                m = Chem.MolFromSmiles(p)
                if m: st.image(MolToImage(m, size=(150,150)))
                st.button(f"Reactant {j+1}", on_click=choose_reactant, args=(p,))

# ─── PAGE: BUILDER ────────────────────────────────────────────────────────────
elif st.session_state.page == "builder":
    st.title("🛠️ Custom Database Builder")

    # Callback to clear only reaction fields (keep db name)
    def clear_builder_fields():
        st.session_state.builder_reactants = []
        st.session_state.builder_product = None
        st.session_state.builder_mol_name = ""
        st.session_state.builder_draw = ""
        st.session_state.builder_conds = ""
        st.session_state.just_saved = False

    # Database name input (persistent)
    raw_name = st.text_input('Database Name (no ".db" suffix)', key='builder_db_name')
    db_base = raw_name.removesuffix('.db') if raw_name else None

    # Show new-reaction prompt if just saved
    if st.session_state.get('just_saved', False):
        st.success("✅ Reaction added! What next?")
        col1, col2 = st.columns([1,1])
        with col1:
            if st.button("➕ Add New Reaction"):
                clear_builder_fields()
                rerun()
        with col2:
            if st.button("⬅️ Back to Main App"):
                st.session_state.page = "main"
                rerun()
        # skip rest of form until user clicks Add New Reaction
        st.stop()

    # Molecule entry
    role = st.radio("Molecule Role", ["Product", "Reactant"], key='builder_role', horizontal=True)
    mode = st.radio("Input by", ["Name", "Draw"], key='builder_input_mode', horizontal=True)

    mol_smiles = None
    if mode == "Name":
        name = st.text_input("Molecule name", key='builder_mol_name')
        if name:
            try:
                mol_smiles = name_to_smiles(name)
                st.success(f"SMILES: {mol_smiles}")
            except Exception as e:
                st.error(f"Error: {e}")
    else:
        drawn = st_ketcher("", key='builder_draw')
        if drawn:
            mol_smiles = drawn
            st.success(f"SMILES: {mol_smiles}")

    # Add molecule to reaction
    if mol_smiles and st.button("✅ Add Molecule to Reaction"):
        if role == "Reactant":
            st.session_state.builder_reactants.append(mol_smiles)
        else:
            st.session_state.builder_product = mol_smiles

    # Reaction preview
    if st.session_state.builder_product or st.session_state.builder_reactants:
        st.markdown("## 👁️ Reaction Preview")
        cols = st.columns(len(st.session_state.builder_reactants) + 2)
        for i, smi in enumerate(st.session_state.builder_reactants):
            with cols[i]:
                st.image(MolToImage(Chem.MolFromSmiles(smi), size=(150, 150)))
                st.caption("Reactant")
        cols[len(st.session_state.builder_reactants)].markdown("➡️")
        if st.session_state.builder_product:
            cols[-1].image(MolToImage(Chem.MolFromSmiles(st.session_state.builder_product), size=(150, 150)))
            cols[-1].caption("Product")
        else:
            cols[-1].info("Product not set")

    # Conditions input
    conds = st.text_area("Conditions (key: value per line)", key='builder_conds')
    parsed = {k:v for k,v in (line.split(':',1) for line in conds.splitlines() if ':' in line)}

    # Save reaction button
    if st.button("🧪 Save Reaction to Database"):
        if not db_base:
            st.error("Please enter a database name (without .db suffix).")
        elif not st.session_state.builder_product:
            st.error("No product specified.")
        elif not st.session_state.builder_reactants:
            st.error("No reactants specified.")
        else:
            try:
                rd.add_new_smart(db_base,
                                 st.session_state.builder_product,
                                 st.session_state.builder_reactants,
                                 parsed)
                # register new DB if not already
                new_file = f"{db_base}.db"
                try:
                    new_db = rd.load_database(new_file)
                    if new_db:
                        rd.register_database(new_db, db_base)
                except:
                    pass
                # flag for showing Add New Reaction
                st.session_state.just_saved = True
                rerun()
            except Exception as e:
                st.error(f"Failed to save: {e}")
    
    # Always allow navigation back to main app
    if st.button("⬅️ Back to Main App"):
        st.session_state.page = "main"
        rerun()
