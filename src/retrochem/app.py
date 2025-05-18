import streamlit as st
import os
import sqlite3
import json
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
from streamlit_ketcher import st_ketcher  # type: ignore
from preliminary_functions import name_to_smiles, canonicalize_smiles
import reaction_database as rd

# ─── PAGE CONFIG & THEME ─────────────────────────────────────────────────────
st.set_page_config(page_title="RetroChem", layout="wide", page_icon="🧪")
st.markdown(
    """
    <style>
    h1, h2, h3 { color: #2E7D32 !important; text-align: center; }
    .stButton>button {
        background-color: #388E3C !important;
        border-color: #2E7D32 !important;
        color: white !important;
    }
    .builder-panel {
        background-color: #f5f5f5;
        padding: 16px;
        border-radius: 8px;
        margin-bottom: 16px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# ─── UTILS ───────────────────────────────────────────────────────────────────
def save_reaction_to_sqlite(db_name, product_smiles, reactant_smiles, conditions):
    # strip .db suffix
    if db_name.lower().endswith('.db'):
        db_name = db_name[:-3]

    # canonicalize before saving
    prod = canonicalize_smiles(product_smiles)
    reacs = [canonicalize_smiles(r) for r in reactant_smiles]

    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(script_dir, f"{db_name}.db")

    # if exists but not SQLite, delete
    if os.path.exists(db_path):
        with open(db_path, 'rb') as f:
            header = f.read(16)
        if not header.startswith(b"SQLite format 3"):
            os.remove(db_path)

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS reactions (
            product    TEXT,
            reactants  TEXT,
            conditions TEXT
        )
    """)
    cur.execute(
        "INSERT INTO reactions (product, reactants, conditions) VALUES (?, ?, ?)",
        (prod, ".".join(reacs), json.dumps(conditions))
    )
    conn.commit()
    conn.close()

# ─── STATE HELPERS ──────────────────────────────────────────────────────────
def reset_all():
    for key in ('selected_smiles','reactant_list','database','combos','history'):
        st.session_state[key] = [] if key in ('combos','history') else None
    st.session_state.page = 'home'
    st.session_state.main_name = ''
    st.session_state.main_draw = ''

def go_back():
    if st.session_state.history:
        prev = st.session_state.history.pop()
        st.session_state.selected_smiles = prev
        st.session_state.reactant_list = None
        st.session_state.combos = rd.list_reactants(
            canonicalize_smiles(prev),
            st.session_state.database
        )
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
    st.session_state.combos = rd.list_reactants(
        canonicalize_smiles(smi),
        st.session_state.database
    )

def choose_combo(idx):
    st.session_state.reactant_list = st.session_state.combos[idx][0].split('.')

def choose_reactant(part):
    st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = part
    st.session_state.reactant_list = None
    st.session_state.combos = rd.list_reactants(
        canonicalize_smiles(part),
        st.session_state.database
    )

# ─── SESSION INIT ────────────────────────────────────────────────────────────
def init_state():
    defaults = {
        'page':'home','selected_smiles':None,'reactant_list':None,
        'database':None,'combos':[],'history':[],
        'main_name':'','main_draw':'',
        'builder_db_name':'','builder_product':None,
        'builder_reactants':[],'builder_conditions':'',
        'builder_input_mode':'Name','builder_mol_name':'','builder_draw':'',
        'builder_input_mode_r':'Name','builder_mol_name_r':'','builder_draw_r':'',
        'clear_reactant_name':False,'just_saved':False
    }
    for k,v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

init_state()

# ─── DATABASE LOADING ─────────────────────────────────────────────────────────
def refresh_databases():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for fn in os.listdir(script_dir):
        if fn.endswith('.db') and not fn.endswith('.db.db'):
            try:
                db = rd.load_database(os.path.join(script_dir, fn))
                if db:
                    rd.register_database(db, fn.removesuffix('.db'))
            except:
                pass

if not rd.REACTION_DATABASES:
    refresh_databases()

# ─── PAGE DISPATCH ───────────────────────────────────────────────────────────
if st.session_state.page == 'home':
    st.title('🧪 Welcome to RetroChem')
    st.markdown('Your Organic Chemistry Retrosynthesis Assistant')
    st.write('---')
    if st.button('🔬 Start Retrosynthesis', key='home_start'):
        st.session_state.page = 'main'
        st.experimental_rerun()

elif st.session_state.page == 'main':
    # ensure new DBs are loaded before showing main
    refresh_databases()

    st.title('RetroChem - Your Guide')
    with st.sidebar:
        st.button('🧹 Start Over', on_click=reset_all, key='sidebar_reset')
        if st.session_state.history:
            st.button('🔙 Back', on_click=go_back, key='sidebar_back')
        st.write('---')
        st.header('🔍 Retrosynthesis Status')

        dbs = list(rd.REACTION_DATABASES.keys())
        opts = dbs + ['➕ Add or edit your own database']
        idx = opts.index(st.session_state.database) if st.session_state.database in opts else 0
        sel = st.selectbox('Choose Database', opts, index=idx, key='sidebar_dbselect')
        if sel == '➕ Add or edit your own database':
            st.session_state.page = 'builder'
            st.experimental_rerun()
        else:
            st.session_state.database = sel

        st.markdown(f"**Current DB:** {st.session_state.database or '_(none)_'}")

    st.write('---')
    # ─── Input Molecule ────────────────────────────────────────────────────
    st.subheader('🧬 Input Molecule')
    mode = st.radio('Mode', ['Name','Draw'], horizontal=True)
    smiles = None
    if mode == 'Name':
        name = st.text_input('Molecule name', key='main_name')
        if name:
            try:
                smiles = name_to_smiles(name)
                st.markdown(f'SMILES: `{smiles}`')
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    st.image(MolToImage(mol, (200,200)))
            except Exception as e:
                st.error(e)
    else:
        drawn = st_ketcher('', key='main_draw')
        if drawn:
            smiles = drawn
            st.markdown(f'SMILES: `{smiles}`')
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                st.image(MolToImage(mol, (200,200)))

    if smiles and st.button('🔄 Retrosynthesize', on_click=start_retro, args=(smiles,), key='main_retro'):
        pass

    # ─── Show results ───────────────────────────────────────────────────────
    if st.session_state.selected_smiles and not st.session_state.reactant_list:
        st.subheader('🧩 Options')
        combos = st.session_state.combos
        if not combos:
            st.info('No routes found.')
        else:
            cols = st.columns(len(combos))
            for i, (s, _) in enumerate(combos):
                with cols[i]:
                    st.image(MolToImage(Chem.MolFromSmiles(s), (150,150)))
                    st.button(f'Option {i+1}', on_click=choose_combo, args=(i,), key=f'opt_{i}')

    elif st.session_state.reactant_list:
        st.subheader('🔹 Next Fragment')
        parts = st.session_state.reactant_list
        cols = st.columns(len(parts))
        for j, p in enumerate(parts):
            with cols[j]:
                st.image(MolToImage(Chem.MolFromSmiles(p), (150,150)))
                st.button(f'Reactant {j+1}', on_click=choose_reactant, args=(p,), key=f'react_{j}')

elif st.session_state.page == 'builder':
    st.markdown("<div class='builder-panel'>", unsafe_allow_html=True)
    st.header('🛠️ Custom Database Builder')

    # top nav only when not just saved
    if not st.session_state.just_saved:
        nav1, nav2 = st.columns(2)
        if nav1.button('🏠 Return to Home', key='nav_return_home'):
            st.session_state.page = 'home'
            st.experimental_rerun()
        if nav2.button('⬅️ Back to Main App', key='nav_back_main'):
            st.session_state.page = 'main'
            st.experimental_rerun()

    # after save, show only these two
    if st.session_state.just_saved:
        c1, c2 = st.columns(2)
        if c1.button('➕ Add New Reaction', key='saved_add_reaction'):
            # reset builder but keep DB name
            st.session_state.builder_product = None
            st.session_state.builder_reactants = []
            st.session_state.builder_conditions = ''
            st.session_state.clear_reactant_name = False
            st.session_state.just_saved = False
            st.experimental_rerun()
        if c2.button('⬅️ Back to Main App', key='saved_back_main'):
            st.session_state.page = 'main'
            st.experimental_rerun()
        st.markdown('</div>', unsafe_allow_html=True)
        st.stop()

    # progress bar
    total = 5
    done = sum(bool(st.session_state[k]) for k in [
        'builder_db_name','builder_product','builder_reactants','builder_conditions'
    ])
    st.progress(done/total)

    # Step 1: DB name
    st.subheader('Step 1: Choose Database')
    st.text_input('Database Name (no .db suffix)', key='builder_db_name')
    st.write('---')

    # Step 2: Product
    with st.expander('Step 2: Add Product', expanded=True):
        mode_p = st.radio('Product Input', ['Name','Draw'], key='builder_input_mode', horizontal=True)
        ps = None
        if mode_p == 'Name':
            pn = st.text_input('Product name', key='builder_mol_name')
            if pn:
                try:
                    ps = name_to_smiles(pn)
                except:
                    pass
        else:
            ps = st_ketcher('', key='builder_draw')
        if ps and st.button('Add Product', key='add_prod'):
            st.session_state.builder_product = ps
            st.experimental_rerun()
        if st.session_state.builder_product:
            st.image(MolToImage(Chem.MolFromSmiles(st.session_state.builder_product), (100,100)))
            st.success('✅ Product added')
    st.write('---')

    # Step 3: Reactants (single only)
    with st.expander('Step 3: Add Reactants', expanded=True):
        if st.session_state.clear_reactant_name:
            st.session_state.builder_mol_name_r = ''
            st.session_state.builder_draw_r = ''
            st.session_state.clear_reactant_name = False

        mode_r = st.radio('Reactant Input', ['Name','Draw'], key='builder_input_mode_r', horizontal=True)
        rs = None
        if mode_r == 'Name':
            rn = st.text_input('Reactant name', key='builder_mol_name_r')
            if rn:
                try:
                    rs = name_to_smiles(rn)
                except:
                    pass
        else:
            rs = st_ketcher('', key='builder_draw_r')

        if rs and st.button('Add Reactant', key='add_react'):
            st.session_state.builder_reactants.append(rs)
            st.session_state.clear_reactant_name = True
            st.experimental_rerun()

        if st.session_state.builder_reactants:
            cols = st.columns(len(st.session_state.builder_reactants))
            for i, s in enumerate(st.session_state.builder_reactants):
                with cols[i]:
                    st.image(MolToImage(Chem.MolFromSmiles(s), (80,80)))
                    st.markdown('✅')
    st.write('---')

    # Step 4: Conditions
    with st.expander('Step 4: Enter Conditions', expanded=True):
        st.text_area('Conditions (k: v per line)', key='builder_conditions')
    st.write('---')

    # Step 5: Save Reaction
    with st.expander('Step 5: Save Reaction', expanded=True):
        if st.button('🧪 Save Reaction to Database', key='save_rxn'):
            dbn  = st.session_state.builder_db_name.strip()
            prod = st.session_state.builder_product
            reac = st.session_state.builder_reactants
            cond = st.session_state.builder_conditions

            if not dbn:
                st.error('Enter DB name')
            elif not prod:
                st.error('Add product')
            elif not reac:
                st.error('Add reactants')
            elif not cond.strip():
                st.error('Enter conditions')
            else:
                cd = {
                    k.strip(): v.strip()
                    for k, v in (
                        line.split(':', 1)
                        for line in cond.splitlines()
                        if ':' in line
                    )
                }
                save_reaction_to_sqlite(dbn, prod, reac, cd)
                refresh_databases()
                st.session_state.database = dbn
                st.success('Reaction saved!')
                st.session_state.just_saved = True
                st.experimental_rerun()

    st.markdown('</div>', unsafe_allow_html=True)
