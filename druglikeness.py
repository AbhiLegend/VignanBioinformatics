import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import numpy as np

st.title("💊 Drug-Likeness Predictor with Visualization")

def featurize(mol):
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
    ]

# Input field for SMILES string
smiles = st.text_input("Enter a SMILES string:", value="CCO")

# Process input when entered
if smiles:
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        # Visualize the molecule
        st.subheader("🖼️ Molecular Structure:")
        mol_image = Draw.MolToImage(mol, size=(300, 300))
        st.image(mol_image, caption="2D Structure")

        # Calculate features
        features = np.array(featurize(mol)).reshape(1, -1)

        # Placeholder logic for prediction
        prediction = "✅ Likely Drug-like" if Descriptors.MolWt(mol) < 500 else "🚫 Not Drug-like"

        # Display prediction
        st.subheader("🔎 Prediction:")
        st.write(f"### **{prediction}**")

        # Display calculated features
        st.subheader("📊 Molecular Properties:")
        st.write(f"- **Molecular Weight:** {Descriptors.MolWt(mol):.2f} g/mol")
        st.write(f"- **LogP:** {Descriptors.MolLogP(mol):.2f}")
        st.write(f"- **H-Bond Donors:** {Descriptors.NumHDonors(mol)}")
        st.write(f"- **H-Bond Acceptors:** {Descriptors.NumHAcceptors(mol)}")

    else:
        st.error("❌ Invalid SMILES string! Please try again.")
