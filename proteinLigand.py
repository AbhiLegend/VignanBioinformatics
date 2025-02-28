import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

st.title("🔬 Protein-Ligand Preparation Tool")

smiles = st.text_input("Enter Ligand SMILES:", "CCO")

if smiles:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        
        st.write("### Prepared Ligand:")
        st.image(Draw.MolToImage(mol, size=(400, 400)))

        # Export SDF (for docking tools like AutoDock Vina)
        sdf_data = Chem.MolToMolBlock(mol)
        st.download_button("Download SDF File", data=sdf_data, file_name="ligand.sdf")

    else:
        st.error("Invalid SMILES entered.")
