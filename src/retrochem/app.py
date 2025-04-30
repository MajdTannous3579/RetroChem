import streamlit as st
import streamlit.components.v1 as components

# import your backend utilities
from retrochem.functions import (
    name_to_smiles,
    structure_to_smiles,
    canonicalize_smiles,
)
from retrochem.example_module import hello_smiles

st.set_page_config(page_title="RetroChem", layout="centered")
st.title("RetroChem - Your Organic Chemistry Guide")

# 1) Choose input mode
input_mode = st.radio("Choose input method:", ["Name", "Draw structure"])

# 2) Capture raw input (name or MolBlock)
raw_input = ""
if input_mode == "Name":
    raw_input = st.text_input("Enter a molecule name")
else:
    # Embed JSME sketcher for drawing
    jsme_html = """
    <script src="https://unpkg.com/jsme-editor@latest/jsme.nocache.js"></script>
    <div id="sketcher" style="width:600px; height:350px;"></div>
    <script>
      const jsme = new JSApplet.JSME("sketcher","600","350");
      function sendMol() {
        const mol = jsme.molfile();
        // stash MolBlock in the URL query under ?mol=...
        window.history.replaceState(null, null, '?mol=' + encodeURIComponent(mol));
      }
      jsme.listen('change', sendMol);
    </script>
    """
    components.html(jsme_html, height=380)

    # Read back MolBlock from the URL
    params = st.query_params
    raw_input = params.get("mol", [""])[0]
    if raw_input:
        st.markdown("**Structure captured.**")

# 3) Retrosynthesis button
if st.button("Retrosynthesis"):
    if not raw_input:
        st.warning("⚠️ Please provide a molecule (name or drawing).")
    else:
        try:
            # Convert raw input → SMILES
            if input_mode == "Name":
                smiles = name_to_smiles(raw_input)
            else:
                smiles = structure_to_smiles(raw_input)

            # Canonicalize
            smiles = canonicalize_smiles(smiles)

        except Exception as e:
            st.error(f"❌ Conversion to SMILES failed:\n{e}")
        else:
            # Run your existing retrosynthesis logic
            result = hello_smiles(smiles)
            st.success(result)
