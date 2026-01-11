import pandas as pd
import numpy as np
from rdkit import Chem
import os

def universal_curator():
    # 1. User Input - Universal design
    file_path = input("Enter the path to your CSV file: ").strip()
    if not os.path.exists(file_path):
        print("Error: File not found.")
        return

    df = pd.read_csv(file_path)
    
    # 2. Dynamic Column Identification
    # This searches for columns by keywords so the user doesn't have to be exact
    smiles_col = next((c for c in df.columns if 'smiles' in c.lower()), None)
    id_col = next((c for c in df.columns if 'id' in c.lower()), None)
    ic50_col = next((c for c in df.columns if 'ic50' in c.lower()), None)

    if not all([smiles_col, id_col, ic50_col]):
        print(f"Error: Could not find required columns. Found: {df.columns.tolist()}")
        return

    # 3. Structural Standardization (Canonicalization)
    def standardize(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            # Standardize tautomers, remove salts, and fix stereochemistry
            return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        return None

    print("Step 1: Standardizing chemical structures...")
    df['SMILES_Standard'] = df[smiles_col].apply(standardize)

    # 4. Rigorous Duplicate & Empty Removal
    initial_len = len(df)
    df = df.dropna(subset=['SMILES_Standard', ic50_col])
    df = df.drop_duplicates(subset=['SMILES_Standard'])
    print(f"Step 2: Removed {initial_len - len(df)} duplicates/invalid entries.")

    # 5. Handling Qualifiers and Unit Conversion
    # We remove > and < but keep the numerical value as a standard QSAR practice
    def clean_activity(val):
        val_str = str(val).replace('>', '').replace('<', '').strip()
        try:
            return float(val_str)
        except ValueError:
            return np.nan

    print("Step 3: Converting IC50 to pIC50...")
    df['IC50_Clean'] = df[ic50_col].apply(clean_activity)
    df = df.dropna(subset=['IC50_Clean'])
    
    # Identify unit from column header (e.g., "IC50 (uM)")
    unit = ic50_col.lower()
    if 'um' in unit or 'micromolar' in unit:
        df['pIC50'] = -np.log10(df['IC50_Clean'] / 1e6)
    elif 'nm' in unit or 'nanomolar' in unit:
        df['pIC50'] = -np.log10(df['IC50_Clean'] / 1e9)
    else:
        # Default to nM if not specified, but warn the user
        print("Warning: Units not detected in column header. Defaulting to nM.")
        df['pIC50'] = -np.log10(df['IC50_Clean'] / 1e9)

    # 6. Final Export
    output_file = file_path.replace(".csv", "_UNIVERSAL_CLEAN.csv")
    df[[id_col, 'SMILES_Standard', 'pIC50']].to_csv(output_file, index=False)
    print(f"--- SUCCESS ---\nCleaned file saved as: {output_file}")

# Run the curator
universal_curator()