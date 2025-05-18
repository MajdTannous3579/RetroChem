import streamlit as st
import pandas as pd
import os
from streamlit_ketcher import st_ketcher  # type: ignore
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

from retrochem.functions import name_to_smiles, canonicalize_smiles
import retrochem.reaction_database as rd

# ─── PAGE CONFIG & THEME ─────────────────────────────────────────────────────
st.set_page_config(page_title="RetroChem", layout="wide", page_icon="🧪")
st.markdown(
    """
    <style>
    h1, h2, h3 { color: #2E7D32 !important; text-align: center; }
    .stButton>button { background-color: #388E3C !important; border-color: #2E7D32 !important; color: white !important; }
    .builder-panel { background-color: #f5f5f5; padding: 16px; border-radius: 8px; }
    </style>
    """,
    unsafe_allow_html=True,
)

# ─── GLOBAL UTILITIES ─────────────────────────────────────────────────────────
def rerun():
    st.experimental_rerun()

def select_database():
    sel = st.session_state.db_selector
    if sel == "➕ Add or edit your own database":
        st.session_state.page = "builder"
    else:
        st.session_state.database = sel
# ─── SESSION STATE SETUP ──────────────────────────────────────────────────────
def init_state():
    defaults = {
        "page": "home",
        "selected_smiles": None,
        "reactant_list": None,
        "database": None,
        "combos": [],
        "history": [],
        # builder fields
        "builder_db_name": "",
        "builder_product": None,
        "builder_reactants": [],
        "builder_conditions": "",
        "builder_input_mode": "Name",
        "builder_mol_name": "",
        "builder_draw": "",
        "just_saved": False,
    }
    for k,v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v
init_state()

# ─── LOAD DATABASES ──────────────────────────────────────────────────────────
def refresh_databases():
    for fn in os.listdir('.'):
        if fn.endswith('.db') and not fn.endswith('.db.db'):
            try:
                db = rd.load_database(fn)
                if db:
                    rd.register_database(db, fn.removesuffix('.db'))
            except:
                pass
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
        st.warning("⚠️ No previous molecule.")

def start_retro(smi):
    if not smi:
        st.warning("Provide a molecule first.")
        return
    if not st.session_state.database:
        st.warning("Choose a database first.")
        return
    if st.session_state.selected_smiles:
        st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = smi
    st.session_state.reactant_list = None
    st.session_state.combos = rd.list_reactants(canonicalize_smiles(smi), st.session_state.database)

def choose_combo(idx):
    st.session_state.reactant_list = st.session_state.combos[idx][0].split('.')

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
    st.markdown("Click below to begin drawing or naming your target molecule.")
    if st.button("🔬 Start Retrosynthesis"):
        st.session_state.page = "main"
        rerun()
    st.stop()

# ─── PAGE: MAIN ───────────────────────────────────────────────────────────────
elif st.session_state.page == "main":
    st.title("RetroChem - Your Guide")
    with st.sidebar:
        st.button("🧹 Start Over", on_click=reset_all)
        if st.session_state.history:
            st.button("🔙 Back", on_click=go_back)
        st.write("---")
        st.header("🔍 Retrosynthesis Status")
        dbs = list(rd.REACTION_DATABASES.keys())
        st.selectbox("Choose Database", dbs+ ["➕ Add or edit your own database"], key="db_selector", on_change=select_database)
        st.markdown(f"**Current DB:** {st.session_state.database or '_(none)_'}")
    st.markdown("---")
    st.subheader("🧬 Input Molecule")
    mode = st.radio("Mode", ["Name","Draw"], horizontal=True)
    smiles = None
    if mode == "Name":
        name = st.text_input("Molecule name", key="main_name")
        if name:
            try:
                smiles = name_to_smiles(name)
                st.write(f"SMILES: `{smiles}`")
                m = Chem.MolFromSmiles(smiles)
                if m: st.image(MolToImage(m, size=(200,200)))
            except Exception as e:
                st.error(e)
    else:
        drawn = st_ketcher("", key="main_draw")
        if drawn:
            smiles = drawn
            st.write(f"SMILES: `{smiles}`")
            m = Chem.MolFromSmiles(smiles)
            if m: st.image(MolToImage(m, size=(200,200)))
    st.markdown("---")
    if smiles and st.button("🔄 Retrosynthesize", on_click=start_retro, args=(smiles,)):
        pass
    if st.session_state.selected_smiles and st.session_state.reactant_list is None:
        st.subheader("🧩 Options")
        combos = st.session_state.combos
        if not combos:
            st.info("No routes found.")
        else:
            cols = st.columns(len(combos))
            for i,(smi,cond) in enumerate(combos):
                with cols[i]:
                    m = Chem.MolFromSmiles(smi)
                    if m: st.image(MolToImage(m, size=(200,200)))
                    st.button(f"Option {i+1}", on_click=choose_combo, args=(i,))
                    df = pd.DataFrame({'Conditions': list(cond.values())}, index=[k.capitalize() for k in cond.keys()])
                    st.table(df)
    elif st.session_state.reactant_list:
        st.subheader("🔹 Next Fragment")
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
    # Background panel for builder
    st.markdown("""
    <div class='builder-panel'>
    """, unsafe_allow_html=True)
    st.header("🛠️ Custom Database Builder")

    # Wizard progress based on steps completed
    total_steps = 4
    done = int(bool(st.session_state.builder_product)) + int(bool(st.session_state.builder_reactants)) + int(bool(st.session_state.builder_conditions)) + int(st.session_state.just_saved)
    progress = done / total_steps
    st.progress(progress)
    st.markdown(f"Step {done+1} of {total_steps+1}")

    # Step 1: Database Name (sticky header)
    st.subheader("Step 1: Choose Database")
    st.text_input("Database Name", key="builder_db_name", help="Name of your custom reaction database (no .db suffix)")
    db_base = st.session_state.builder_db_name.removesuffix('.db') if st.session_state.builder_db_name else None
    st.write("---")

    # Step 2: Add Product
    with st.expander("Step 2: Add Product", expanded=True):
        st.subheader("Step 2: Add Product")
        mode_p = st.radio("Input by", ["Name","Draw"], key="builder_input_mode", help="How to specify your product molecule", horizontal=True)
        p_smiles = None
        if mode_p == "Name":
            n = st.text_input("Product name", key="builder_mol_name")
            if n:
                try:
                    p_smiles = name_to_smiles(n)
                    st.success(f"SMILES: {p_smiles}")
                except Exception as e:
                    st.error(e)
        else:
            d = st_ketcher("", key="builder_draw")
            if d:
                p_smiles = d
                st.success(f"SMILES: {p_smiles}")
        if p_smiles and st.button("Add Product"):
            st.session_state.builder_product = p_smiles
            rerun()
        if st.session_state.builder_product:
            st.image(MolToImage(Chem.MolFromSmiles(st.session_state.builder_product), size=(150,150)))
            st.markdown("✅ Product added")
    st.write("---")

    # Step 3: Add Reactants
    with st.expander("Step 3: Add Reactants", expanded=True):
        st.subheader("Step 3: Add Reactants")
        multi = st.text_area("Or paste reactant names, comma-separated", key="builder_reactants_batch", help="Quickly add multiple reactants by name")
        if multi:
            for r in (x.strip() for x in multi.split(',')):
                try:
                    s = name_to_smiles(r)
                    st.session_state.builder_reactants.append(s)
                except:
                    continue
            st.success(f"Added {len(multi.split(','))} reactants")
            st.session_state.builder_reactants_batch = ""
            rerun()
        mode_r = st.radio("Input reactant by", ["Name","Draw"], key="builder_input_mode_r", help="Specify a single reactant", horizontal=True)
        r_smiles = None
        if mode_r == "Name":
            nr = st.text_input("Reactant name")
            if nr:
                try:
                    r_smiles = name_to_smiles(nr)
                except:
                    st.error("Conversion failed")
        else:
            dr = st_ketcher("")
            if dr:
                r_smiles = dr
        if r_smiles and st.button("Add Reactant"):
            st.session_state.builder_reactants.append(r_smiles)
            rerun()
        if st.session_state.builder_reactants:
            cols = st.columns(len(st.session_state.builder_reactants))
            for i,smi in enumerate(st.session_state.builder_reactants):
                with cols[i]: st.image(MolToImage(Chem.MolFromSmiles(smi), size=(100,100))); st.markdown("✅")
    st.write("---")

    # Step 4: Enter Conditions
    with st.expander("Step 4: Enter Conditions", expanded=True):
        st.subheader("Step 4: Enter Conditions")
        st.text_area("Conditions (key: value per line)", key="builder_conditions", help="E.g. catalyst: Pd/C")
    st.write("---")

    # Step 5: Save Reaction
    with st.expander("Step 5: Save Reaction", expanded=True):
        st.subheader("Step 5: Save Reaction")
        if st.button("🧪 Save Reaction to Database"):
            if not db_base:
                st.error("Enter a database name.")
            elif not st.session_state.builder_product:
                st.error("Add a product.")
            elif not st.session_state.builder_reactants:
                st.error("Add reactants.")
            else:
                parsed = {k:v for k,v in (line.split(':',1) for line in st.session_state.builder_conditions.splitlines() if ':' in line)}
                rd.add_new_smart(db_base, st.session_state.builder_product, st.session_state.builder_reactants, parsed)
                st.success("Reaction saved!")
                st.session_state.just_saved = True
                rerun()
        if st.session_state.just_saved:
            if st.button("➕ Add New Reaction"):
                # clear fields but keep DB name
                for f in ['builder_product','builder_reactants','builder_conditions']: st.session_state[f] = None if f=='builder_product' else [] if f=='builder_reactants' else ""
                st.session_state.just_saved = False
                rerun()
            if st.button("⬅️ Back to Main App"):
                st.session_state.page = "main"
                rerun()
    st.markdown("</div>", unsafe_allow_html=True)
