import pandas as pd
import qcelemental as qcel
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import tempfile
import os
import argparse

def get_cheminformatics_properties(mol_obj):
    """
    Converts a qcelemental Molecule object to an RDKit molecule,
    and calculates a set of cheminformatics properties.
    """
    properties = {
        'smiles': 'N/A',
        'is_aromatic': False,
        'mol_weight': 0,
        'num_rings': 0,
        'num_rotatable_bonds': 0,
        'logp': 0,
        'tpsa': 0,
    }
    tmp_path = 'temp.xyz'
    mol_obj.to_file(tmp_path)
    charge = int(mol_obj.molecular_charge)
    mol = Chem.MolFromXYZFile(tmp_path)
    rdDetermineBonds.DetermineConnectivity(mol, charge=charge)
    os.remove(tmp_path)

    if mol is None:
        return properties

    properties['smiles'] = Chem.MolToSmiles(mol)
    
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
    properties['is_aromatic'] = is_aromatic
    
    properties['mol_weight'] = Descriptors.MolWt(mol)
    properties['num_rings'] = Descriptors.RingCount(mol)
    properties['num_rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
    properties['logp'] = Descriptors.MolLogP(mol)
    properties['tpsa'] = Descriptors.TPSA(mol)

    return properties

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

        # Get cheminformatics properties
        props_a = get_cheminformatics_properties(mon_a)
        props_b = get_cheminformatics_properties(mon_b)

        results.append({
            'error': row['error'],
            'smiles_a': props_a['smiles'],
            'is_aromatic_a': props_a['is_aromatic'],
            'mol_weight_a': props_a['mol_weight'],
            'num_rings_a': props_a['num_rings'],
            'num_rotatable_bonds_a': props_a['num_rotatable_bonds'],
            'logp_a': props_a['logp'],
            'tpsa_a': props_a['tpsa'],
            'smiles_b': props_b['smiles'],
            'is_aromatic_b': props_b['is_aromatic'],
            'mol_weight_b': props_b['mol_weight'],
            'num_rings_b': props_b['num_rings'],
            'num_rotatable_bonds_b': props_b['num_rotatable_bonds'],
            'logp_b': props_b['logp'],
            'tpsa_b': props_b['tpsa'],
        })

    results_df = pd.DataFrame(results)
    results_df['contains_aromatic'] = results_df['is_aromatic_a'] | results_df['is_aromatic_b']
    
    print(results_df.head())

    # Visualize the results
    pair_plot_vars = ['error', 'mol_weight_a', 'mol_weight_b', 'logp_a', 'logp_b', 'tpsa_a', 'tpsa_b']
    sns.pairplot(results_df, vars=pair_plot_vars, hue='contains_aromatic', diag_kind='kde')
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze SAPT errors and cheminformatics properties.')
    parser.add_argument('--input_file', type=str, default='combined_df_4569.pkl', help='Path to the input pickle file.')
    parser.add_argument('--output_file', type=str, default='cheminformatics_analysis.png', help='Path to save the output plot.')
    args = parser.parse_args()
    main(args.input_file, args.output_file)
