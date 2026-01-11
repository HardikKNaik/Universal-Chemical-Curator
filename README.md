# Universal Chemical Curator 

A robust Python tool for the standardized curation of chemical datasets, specifically designed for QSAR and Stereochem-aware modeling.

## Features
- **Intelligent Column Detection**: Automatically finds ID, SMILES, and IC50 columns regardless of naming.
- **RDKit Standardization**: Performs canonicalization to remove "hidden" duplicates.
- **Stereo-Aware**: Preserves (R)/(S) chirality while standardizing the molecular backbone.
- **Activity Processing**: Automatically filters non-quantitative qualifiers (>, <) and converts IC50 to pIC50.

## How to Use
1. Clone the repo: `git clone https://github.com/your-username/Universal-Chemical-Curator.git`
2. Install requirements: `pip install -r requirements.txt`
3. Run the script: `python curate.py`
4. Enter your CSV path when prompted!

## Thank You!
