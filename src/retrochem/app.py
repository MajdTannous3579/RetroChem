import streamlit as st
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage
from streamlit_ketcher import st_ketcher
from preliminary_functions import name_to_smiles, canonicalize_smiles
import reaction_database as rd

# -------------- PAGE CONFIG & THEME --------------
st.set_page_config(page_title="RetroChem", layout="wide", page_icon="./logo.png")
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

# -------------- CALLBACKS --------------
def go_home():
    """
    Reset the session state and navigate to the Home page.

    All keys except 'page' are cleared from st.session_state.
    The default page is set to 'home', and the database registry is refreshed.
    """
    # Clear all session state except page, then reset to home and reload databases
    for k in list(st.session_state.keys()):
        if k != 'page':
            del st.session_state[k]
    st.session_state.page = 'home'
    refresh_databases()


def go_main():
    """
    Set the current session state page to 'main'.

    This is used to transition from the Home page to the Retrosynthesis interface.
    """
    st.session_state.page = 'main'


def clear_history():
    """
    Clears all retrosynthesis history and related state while keeping the user on the main page.

    Also refreshes the reaction database registry from disk.
    """
    # Full reset within main: clear all state except page, keep page as main, reload databases
    for k in list(st.session_state.keys()):
        if k != 'page':
            del st.session_state[k]
    st.session_state.page = 'main'
    refresh_databases()


def do_retrosynthesis(smiles, db):
    """
    Executes a retrosynthesis step for the given SMILES string using the selected database. Uses backend function list_reactants.

    Parameters
    ----------
    smiles : str
        The SMILES representation of the molecule to analyze.
    db : str
        The database name (must match a key in REACTION_DATABASES).

    - Pushes current molecule to history.
    """
    if st.session_state.get('selected_smiles'):
        st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = smiles
    st.session_state.reactant_list = None
    st.session_state.combos = rd.list_reactants(smiles, db)


def select_option(idx):
    """
    Selects a retrosynthesis route by index from `st.session_state.combos`.

    Parameters
    ----------
    idx : int
        Index of the route option selected.

    - Parses the dot-separated reactants into a list.
    - Stores them in `st.session_state.reactant_list`.
    """
    s, _ = st.session_state.combos[idx]
    st.session_state.reactant_list = s.split('.')


def select_fragment(idx):
    """
    Recursively apply retrosynthesis to a selected fragment from the current reactant list.

    Parameters
    ----------
    idx : int
        Index of the reactant fragment to analyze.
    ----------
    - Pushes current molecule to history.
    - Updates `selected_smiles` to the selected fragment.
    - Triggers a new retrosynthesis search for that fragment.
    """
    p = st.session_state.reactant_list[idx]
    st.session_state.history.append(st.session_state.selected_smiles)
    st.session_state.selected_smiles = p
    st.session_state.reactant_list = None
    st.session_state.combos = rd.list_reactants(canonicalize_smiles(p), st.session_state.database)


def reset_builder():
    """
    Clears all fields in the custom database builder.

    This includes product, reactants, and conditions.
    Resets the 'just_saved' flag to allow a fresh entry.
    """
    st.session_state.builder_product = None
    st.session_state.builder_reactants = []
    st.session_state.builder_conditions = ''
    st.session_state.just_saved = False

def remove_reactant(idx):
    """
    Removes a single reactant by index from the builder list and re-renders the app.

    Parameters
    ----------
    idx : int
        Index of the reactant to be removed.
    """
    st.session_state.builder_reactants.pop(idx)
    st.rerun()

def handle_db_select():
    """
    Handles user selection from the database selectbox in the sidebar.

    - If 'Add or edit your own database' is selected, redirects to the builder page.
    - Otherwise sets the selected database in `st.session_state.database`.
    """
    sel = st.session_state.main_db_select
    add_label = '➕ Add or edit your own database'
    placeholder = 'Select Database'
    if sel == add_label:
        st.session_state.page = 'builder'
    elif sel != placeholder:
        st.session_state.database = sel


# -------------- SESSION STATE HELPERS --------------
def init_state():
    """
    Initializes default values in Streamlit's session state to ensure a stable initial application state.

    This function should be called once at startup.
    """
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


# -------------- DATABASE MANAGEMENT --------------
def refresh_databases():
    """
    Scans the current directory for `.db` files, loads them, and registers each one.

    - Uses the backend `load_database` and `register_database` functions.
    - Populates the global reaction database registry (REACTION_DATABASES and REACTION_REVERSERS).

    This ensures that all valid databases are available for retrosynthesis after any changes.
    """
    rd.clear_registered_databases()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for fn in os.listdir(script_dir):
        if fn.endswith('.db') and not fn.endswith('.db.db'):
            path = os.path.join(script_dir, fn)
            db = rd.load_database(path)
            if db is not None:
                rd.register_database(db, fn.removesuffix('.db'))

# Register databases on startup
if not rd.REACTION_DATABASES:
    refresh_databases()


# -------------- PAGE DISPATCH --------------
if st.session_state.page == 'home':
    st.title('Welcome to Retrochem')
    st.markdown("<p style='text-align: center;'>Your Organic Chemistry Retrosynthesis Assistant</p>", unsafe_allow_html=True)
    col1, col2, col3 = st.columns([3.8, 5, 2])
    with col2:
        st.image("./logo.png", width = 400)
    col1, col2, col3 = st.columns([2.3, 1, 2])
    with col2:
        st.button('🔬 Start Retrosynthesis',
                key='home_start_retrosynthesis',
                on_click=go_main)

elif st.session_state.page == 'main':
    refresh_databases()
    st.title('RetroChem - Your Retro-Organic Chemistry Guide')

    # -------------- Sidebar --------------
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

        if st.session_state.selected_smiles:
            st.subheader('🧪 Current Molecule')
            st.markdown(f'`{st.session_state.selected_smiles}`')
            mol = Chem.MolFromSmiles(st.session_state.selected_smiles)
            if mol:
                st.image(MolToImage(mol, (100,100)))

        # -------------- Selectbox with on_change --------------
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

    # -------------- Main Panel --------------
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
        st.subheader('📊 Options')
    # Require a database first
        if not st.session_state.database:
            st.warning('🚨 Please select a database before retrosynthesizing.')
        else:
            combos = st.session_state.combos or []
            if not combos:
                st.info('Your product is too simple, or not in the chosen database.')
            else:
                cols = st.columns(len(combos))
                for i, (s, cond) in enumerate(combos):
                    with cols[i]:
                        st.image(
                            MolToImage(Chem.MolFromSmiles(s), (150,150))
                        )
                        if cond:
                            df = pd.DataFrame.from_dict(cond, orient='index', columns=['Value'])
                            df.index.name = 'Parameter'
                            st.table(df)
                        st.button(
                            f'Option {i+1}',
                            key=f'opt{i}',
                            on_click=select_option,
                            args=(i,),
                        )

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

# DATABASE BUILDER
elif st.session_state.page == 'builder':
    c1, c2 = st.columns([1,1])
    with c1:
        st.button('🔙 Back to Retrosynthesis',
                  key='builder_return_to_main',
                  on_click=go_main)
    with c2:
        st.button('🧹 Start Over',
                  key='builder_builder_startover',
                  on_click=reset_builder)

    st.header('🛠️ Custom Database Builder')
    st.write('build or edit your databases here, input your reaction paterns and conditions!')
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
    st.text_input('Database name, will be saved as db if new, or edited if already existing', key='builder_db_name')
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
        if st.button('🗑️', key='remove_product'):
            st.session_state.builder_product = None
            st.rerun()

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
    # for each reactant, show it + a remove button
        cols = st.columns(len(st.session_state.builder_reactants))
        for idx, smi in enumerate(st.session_state.builder_reactants):
            with cols[idx]:
                mol = Chem.MolFromSmiles(smi)
                st.image(MolToImage(mol, (80,80)))
                # remove this reactant button
                if st.button('🗑️', key=f'remove_r{idx}'):
                    remove_reactant(idx)

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

    # Reaction preview
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
