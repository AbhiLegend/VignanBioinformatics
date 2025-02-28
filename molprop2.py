import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

# App title
st.title("🧬 Molecular Property Calculator with Visualization")

# Input field for SMILES
smiles = st.text_input("Enter a SMILES string:", "CCO")

# Process the input
if smiles:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # Display molecular properties
        st.write("### 🧪 Molecular Properties:")
        st.write(f"- **Molecular Weight:** {Descriptors.MolWt(mol):.2f} g/mol")
        st.write(f"- **LogP (Partition Coefficient):** {Descriptors.MolLogP(mol):.2f}")
        st.write(f"- **Number of Hydrogen Donors:** {Descriptors.NumHDonors(mol)}")
        st.write(f"- **Number of Hydrogen Acceptors:** {Descriptors.NumHAcceptors(mol)}")

        # Visualize the molecule
        st.write("### 🖼️ Molecular Structure:")
        mol_image = Draw.MolToImage(mol, size=(400, 400))
        st.image(mol_image, caption="2D Structure")

    else:
        st.error("❌ Invalid SMILES string! Please enter a valid one.")
