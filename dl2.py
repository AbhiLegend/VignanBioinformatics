import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import numpy as np

st.title("💊 Advanced Drug-Likeness Predictor")

def featurize(mol):
    """Extract molecular descriptors"""
    return {
        "Molecular Weight": Descriptors.MolWt(mol),
        "LogP (Hydrophobicity)": Descriptors.MolLogP(mol),
        "Hydrogen Bond Donors": Descriptors.NumHDonors(mol),
        "Hydrogen Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Molar Refractivity": Descriptors.MolMR(mol),
    }

def evaluate_drug_likeness(properties):
    """Apply multiple drug-likeness filters"""
    lipinski_criteria = sum([
        properties["Molecular Weight"] < 500,
        properties["LogP (Hydrophobicity)"] < 5,
        properties["Hydrogen Bond Donors"] <= 5,
        properties["Hydrogen Bond Acceptors"] <= 10
    ])

    ghose_criteria = (
        160 <= properties["Molecular Weight"] <= 480 and
        40 <= properties["Molar Refractivity"] <= 130 and
        -0.4 <= properties["LogP (Hydrophobicity)"] <= 5.6
    )

    veber_criteria = (
        properties["Rotatable Bonds"] <= 10 and
        (properties["Hydrogen Bond Donors"] + properties["Hydrogen Bond Acceptors"]) <= 12
    )

    molecular_weight_limit = properties["Molecular Weight"] <= 900  # Hard cap

    # Final decision: Must pass at least Lipinski + (either Ghose or Veber) and not exceed 900 MW
    if (
        lipinski_criteria >= 3 and  # Allow one Lipinski violation
        (ghose_criteria or veber_criteria) and
        molecular_weight_limit
    ):
        return "✅ Likely Drug-like"
    else:
        return "🚫 Not Drug-like (Fails Multiple Criteria)"

# Input field for SMILES
smiles = st.text_input("Enter a SMILES string:", value="CCO")

# Process input
if smiles:
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        # Visualize the molecule
        st.subheader("🖼️ Molecular Structure:")
        mol_image = Draw.MolToImage(mol, size=(300, 300))
        st.image(mol_image, caption="2D Structure")

        # Calculate features
        properties = featurize(mol)

        # Apply improved drug-likeness prediction
        prediction = evaluate_drug_likeness(properties)

        # Display prediction
        st.subheader("🔎 Drug-Likeness Prediction:")
        if "Not Drug-like" in prediction:
            st.warning(prediction)
        else:
            st.success(prediction)

        # Display detailed molecular properties
        st.subheader("📊 Molecular Properties:")
        for key, value in properties.items():
            st.write(f"- **{key}:** {value:.2f}")

    else:
        st.error("❌ Invalid SMILES string! Please try again.")
