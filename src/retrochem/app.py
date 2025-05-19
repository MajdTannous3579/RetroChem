import streamlit as st
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
from streamlit_ketcher import st_ketcher  # type: ignore
from preliminary_functions import name_to_smiles, canonicalize_smiles
import reaction_database as rd

# ─── PAGE CONFIG & THEME ─────────────────────────────────────────────────────
st.set_page_config(page_title="RetroChem", layout="wide", page_icon="🧪")
st.markdown("""
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
""", unsafe_allow_html=True)

# ─── CALLBACKS ────────────────────────────────────────────────────────────────
def go_home():
    for k in list(st.session_state.keys()):
        if k != 'page':
            del st.session_state[k]
    st.session_state.page = 'home'

def go_main():
    st.session_state.page = 'main'

def clear_history():
    st.session_state.history = []

def do_retrosynthesis(smiles, db):
    if st.session_state.get('selected_smiles'):
        st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = smiles
    st.session_state.reactant_list = None
    st.session_state.combos = rd.list_reactants(canonicalize_smiles(smiles), db)

def select_option(idx):
    s, _ = st.session_state.combos[idx]
    st.session_state.reactant_list = s.split('.')

def select_fragment(idx):
    p = st.session_state.reactant_list[idx]
    st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = p
    st.session_state.reactant_list = None
    st.session_state.combos = rd.list_reactants(canonicalize_smiles(p), st.session_state.database)

def reset_builder():
    st.session_state.builder_product = None
    st.session_state.builder_reactants = []
    st.session_state.builder_conditions = ''
    st.session_state.just_saved = False

def handle_db_select():
    sel = st.session_state.main_db_select
    add_label = '➕ Add or edit your own database'
    placeholder = 'Select Database'
    if sel == add_label:
        st.session_state.page = 'builder'
    elif sel != placeholder:
        st.session_state.database = sel

# ─── SESSION STATE HELPERS ───────────────────────────────────────────────────
def init_state():
    defaults = {
        'page': 'home',
        'main_name': '', 'main_draw': '',
        'selected_smiles': None, 'reactant_list': None,
        'database': None, 'combos': [], 'history': [],
        'builder_db_name': '', 'builder_product': None,
        'builder_reactants': [], 'builder_conditions': '',
        'just_saved': False
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

init_state()

# ─── DATABASE MANAGEMENT ─────────────────────────────────────────────────────
def refresh_databases():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for fn in os.listdir(script_dir):
        if fn.endswith('.db') and not fn.endswith('.db.db'):
            path = os.path.join(script_dir, fn)
            db = rd.load_database(path)
            if db is not None:
                rd.register_database(db, fn.removesuffix('.db'))

if not rd.REACTION_DATABASES:
    refresh_databases()

# ─── PAGE DISPATCH ───────────────────────────────────────────────────────────
if st.session_state.page == 'home':
    st.title('🧪 Welcome to RetroChem')
    st.markdown('Your Organic Chemistry Retrosynthesis Assistant')
    st.write('---')
    st.button('🔬 Start Retrosynthesis',
              key='home_start_retrosynthesis',
              on_click=go_main)

elif st.session_state.page == 'main':
    refresh_databases()
    st.title('RetroChem - Your Guide')

    # ─── Sidebar ──────────────────────────────────────────────────────────────
    with st.sidebar:
        st.button('🏠 Return to Home',
                  key='main_return_home',
                  on_click=go_home)
        st.button('🧹 Start Over',
                  key='main_start_over',
                  on_click=clear_history)
        if st.session_state.history:
            st.button('🔙 Back',
                      key='main_back',
                      on_click=lambda: st.session_state.history.pop())

        if st.session_state.history:
            st.markdown('---')
            st.subheader('🔄 History')
            for i, smi in enumerate(st.session_state.history, start=1):
                st.markdown(f"{i}. `{smi}`")

        st.markdown('---')
        st.header('🔍 Retrosynthesis Status')

        # ─── Selectbox with on_change ───────────────────────────────────────
        dbs     = list(rd.REACTION_DATABASES.keys())
        placeholder   = 'Select Database'
        add_label     = '➕ Add or edit your own database'
        opts    = [placeholder] + dbs + [add_label]
        st.selectbox('Choose Database',
                     opts,
                     index=0,
                     key='main_db_select',
                     on_change=handle_db_select)
        st.markdown(f"**Current DB:** {st.session_state.database or '_(none)_'}")

    # ─── Main Panel ────────────────────────────────────────────────────────────
    st.write('---')
    st.subheader('🧬 Input Molecule')
    mode = st.radio('Mode', ['Name','Draw'], horizontal=True)
    smiles = None

    if mode == 'Name':
        name = st.text_input('Molecule name', key='main_name')
        if name:
            try:
                smiles = name_to_smiles(name)
                st.markdown(f'SMILES: `{smiles}`')
                m = Chem.MolFromSmiles(smiles)
                if m:
                    st.image(MolToImage(m, (200,200)))
            except Exception as e:
                st.error(e)
    else:
        drawn = st_ketcher('', key='main_draw')
        if drawn:
            smiles = drawn
            st.markdown(f'SMILES: `{smiles}`')
            m = Chem.MolFromSmiles(smiles)
            if m:
                st.image(MolToImage(m, (200,200)))

    if smiles:
        st.button('🔄 Retrosynthesize',
                  key='main_retrosynthesize',
                  on_click=do_retrosynthesis,
                  args=(smiles, st.session_state.database))

    if st.session_state.selected_smiles and not st.session_state.reactant_list:
        st.subheader('🧩 Options')
        combos = st.session_state.combos or []
        if not combos:
            st.info('No routes found.')
        else:
            cols = st.columns(len(combos))
            for i, (s, cond) in enumerate(combos):
                with cols[i]:
                    st.image(MolToImage(Chem.MolFromSmiles(s), (150,150)))
                    if cond:
                        df = pd.DataFrame.from_dict(cond, orient='index', columns=['Value'])
                        df.index.name = 'Parameter'
                        st.table(df)
                    st.button(f'Option {i+1}',
                              key=f'opt{i}',
                              on_click=select_option,
                              args=(i,))

    elif st.session_state.reactant_list:
        st.subheader('🔹 Next Fragment')
        parts = st.session_state.reactant_list
        cols = st.columns(len(parts))
        for i, p in enumerate(parts):
            with cols[i]:
                st.image(MolToImage(Chem.MolFromSmiles(p), (150,150)))
                st.button(f'Reactant {i+1}',
                          key=f'react{i}',
                          on_click=select_fragment,
                          args=(i,))

elif st.session_state.page == 'builder':
    st.button('🔙 Back to Retrosynthesis',
              key='builder_return_to_main',
              on_click=go_main)

    st.markdown("<div class='builder-panel'>", unsafe_allow_html=True)
    st.header('🛠️ Custom Database Builder')

    total_steps = 5
    completed = (
        int(bool(st.session_state.builder_db_name.strip())) +
        int(st.session_state.builder_product is not None) +
        int(len(st.session_state.builder_reactants) > 0) +
        int(bool(st.session_state.builder_conditions.strip())) +
        int(st.session_state.just_saved)
    )
    # show fraction text if you like
    st.write(f"Progress: **{completed} / {total_steps}**")
    st.progress(completed / total_steps)

    if st.session_state.just_saved:
        st.success('✅ Reaction saved!')
        c1, _ = st.columns(2)
        c1.button('➕ Add New Reaction',
                  key='builder_add_new_rxn',
                  on_click=reset_builder)
        st.markdown('</div>', unsafe_allow_html=True)
        st.stop()

    # Step 1: Database Name
    st.subheader('Step 1: Database Name')
    st.text_input('Name (no .db)', key='builder_db_name')
    st.write('---')

    # Step 2: Add Product
    st.subheader('Step 2: Add Product')
    mode_p = st.radio('Input', ['Name','Draw'], key='builder_mode_p', horizontal=True)
    ps = None
    if mode_p == 'Name':
        pn = st.text_input('Product name', key='builder_name_p')
        if pn:
            ps = name_to_smiles(pn)
    else:
        ps = st_ketcher('', key='builder_draw_p')

    if ps:
        st.button('Add Product',
                  key='builder_add_product',
                  on_click=lambda: st.session_state.update({'builder_product': ps}))

    if st.session_state.builder_product:
        mol = Chem.MolFromSmiles(st.session_state.builder_product)
        st.image(MolToImage(mol, (100,100)))
        st.success('✅ Product added')
    st.write('---')

    # Step 3: Add Reactants
    st.subheader('Step 3: Add Reactants')
    mode_r = st.radio('Input', ['Name','Draw'], key='builder_mode_r', horizontal=True)
    rs = None
    if mode_r == 'Name':
        rn = st.text_input('Reactant name', key='builder_name_r')
        if rn:
            rs = name_to_smiles(rn)
    else:
        rs = st_ketcher('', key='builder_draw_r')

    if rs:
        st.button('Add Reactant',
                  key='builder_add_reactant',
                  on_click=lambda: st.session_state.builder_reactants.append(rs))

    if st.session_state.builder_reactants:
        for r in st.session_state.builder_reactants:
            mol = Chem.MolFromSmiles(r)
            st.image(MolToImage(mol, (80,80)))
    st.write('---')

    # Step 4: Conditions
    st.subheader('Step 4: Conditions')
    st.markdown(
        "Enter each condition on its own line in `key: value` format.  \n"
        "For example:  \n"
        "`solvent: THF`  \n"
        "`temperature: 25°C`  \n"
        "`catalyst: Pd/C`"
    )
    st.text_area('Conditions', key='builder_conditions', height=100)
    st.write('---')

    # Step 5: Save Reaction
    st.subheader('Step 5: Save Reaction')
    st.button('🧪 Save Reaction',
              key='builder_save_reaction',
              on_click=lambda: (
                  rd.add_new_smart(
                      st.session_state.builder_db_name.strip(),
                      st.session_state.builder_product,
                      st.session_state.builder_reactants,
                      {
                          k.strip(): v.strip()
                          for line in st.session_state.builder_conditions.splitlines()
                          if ':' in line
                          for k, v in [line.split(':', 1)]
                      }
                  ),
                  refresh_databases(),
                  st.session_state.update({
                      'database': st.session_state.builder_db_name.strip(),
                      'just_saved': True
                  })
              ))

    st.markdown('</div>', unsafe_allow_html=True)
