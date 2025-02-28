import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

# Title of the app
st.title("🧬 Molecular Visualization App")
st.write("Enter a **SMILES string** to visualize the molecule!")

# Input for SMILES
smiles = st.text_input("SMILES Input", "CCO")  # Default: ethanol

if smiles:
    try:
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:
            # Display the molecule
            st.subheader("Molecular Structure")
            mol_image = Draw.MolToImage(mol, size=(400, 400))
            st.image(mol_image)
        else:
            st.error("❌ Invalid SMILES string. Please try again.")
            
    except Exception as e:
        st.error(f"⚠️ Error: {e}")
