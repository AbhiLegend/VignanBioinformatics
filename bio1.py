import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

st.title("🧬 Molecular Similarity for Drug Screening")

smiles1 = st.text_input("Enter first molecule SMILES:", "CCO")
smiles2 = st.text_input("Enter second molecule SMILES:", "CCN")

if smiles1 and smiles2:
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 and mol2:
        # Generate fingerprints
        fp1 = FingerprintMols.FingerprintMol(mol1)
        fp2 = FingerprintMols.FingerprintMol(mol2)

        # Calculate similarity
        similarity = DataStructs.FingerprintSimilarity(fp1, fp2)

        # Display molecules and similarity
        st.write("### Molecular Similarity Score:")
        st.write(f"🔎 Similarity: **{similarity:.2f}**")
        st.image(Draw.MolsToGridImage([mol1, mol2], molsPerRow=2, subImgSize=(200, 200)))

    else:
        st.error("Invalid SMILES entered.")
