import pandas as pd
import qcelemental as qcel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import tempfile
import os
import argparse

def get_smiles_and_aromaticity(mol_obj):
    """Converts a qcelemental Molecule object to a SMILES string and checks for aromaticity."""
    try:
        tmp_path = 'temp.xyz'
        mol_obj.to_file(tmp_path)
        charge = int(mol_obj.molecular_charge)
        mol = Chem.MolFromXYZFile(tmp_path)
        rdDetermineBonds.DetermineConnectivity(mol, charge=charge)
        # os.remove(tmp_path)

        if mol is None:
            return "N/A", False

        smiles = Chem.MolToSmiles(mol)
        is_aromatic = False
        if mol.GetNumBonds() > 0:
            for bond in mol.GetBonds():
                if bond.GetIsAromatic():
                    is_aromatic = True
                    break
        if not is_aromatic:
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic():
                    is_aromatic = True
                    break

    except Exception as e:
        print(f"Error processing molecule: {e}")
        return "Error", False
            
    return smiles, is_aromatic

def main(input_file, output_file):
    df = pd.read_pickle(input_file)

    # Calculate error
    df['error'] = (df['SAPT0 TOTAL ENERGY atqz'] - df['SAPT2+3(CCD) TOTAL ENERGY atqz']) * 627.5

    results = []
    for index, row in df.iterrows():
        # Create qcelemental molecule with fragments
        dimer = qcel.models.Molecule.from_data({
            'symbols': row['atomic_numbers'],
            'geometry': row['coordinates'] / qcel.constants.bohr2angstroms,
            'molecular_charge': row['dimer_charge'],
            'molecular_multiplicity': row['dimer_multiplicity'],
            'fragments': [row['monAs'], row['monBs']],
            'fragment_charges': [row['monA_charge'], row['monB_charge']],
            'fragment_multiplicities': [row['monA_multiplicity'], row['monB_multiplicity']],
        })

        # Get monomers
        mon_a = dimer.get_fragment(0)
        mon_b = dimer.get_fragment(1)

        # Get SMILES and aromaticity
        smiles_a, is_aromatic_a = get_smiles_and_aromaticity(mon_a)
        smiles_b, is_aromatic_b = get_smiles_and_aromaticity(mon_b)

        results.append({
            'error': row['error'],
            'smiles_a': smiles_a,
            'is_aromatic_a': is_aromatic_a,
            'smiles_b': smiles_b,
            'is_aromatic_b': is_aromatic_b,
        })

    results_df = pd.DataFrame(results)
    results_df['contains_aromatic'] = results_df['is_aromatic_a'] | results_df['is_aromatic_b']
    
    print(results_df.head())
    print(results_df['contains_aromatic'].value_counts())

    # Visualize the results
    plt.figure(figsize=(10, 7))
    ax = sns.violinplot(x='contains_aromatic', y='error', data=results_df)
    ax.set_xticklabels(['Non-Aromatic', 'Aromatic'])
    plt.xlabel('System Type')
    plt.ylabel('Error (kcal/mol)')
    plt.title('Error Distribution for Aromatic vs. Non-Aromatic Systems')
    
    # Annotate with counts
    n_aromatic = results_df['contains_aromatic'].sum()
    n_non_aromatic = len(results_df) - n_aromatic
    ax.text(0, ax.get_ylim()[1]*0.9, f'n = {n_non_aromatic}', ha='center')
    if n_aromatic > 0:
        ax.text(1, ax.get_ylim()[1]*0.9, f'n = {n_aromatic}', ha='center')
    
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze SAPT errors for aromatic vs. non-aromatic systems.')
    parser.add_argument('--input_file', type=str, default='combined_df_4569.pkl', help='Path to the input pickle file.')
    parser.add_argument('--output_file', type=str, default='aromaticity_plot.png', help='Path to save the output plot.')
    args = parser.parse_args()
    main(args.input_file, args.output_file)
