import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse

def main(input_file, output_file):
    df = pd.read_pickle(input_file)

    sapt0_cols = ['SAPT0 TOTAL ENERGY adtz', 'SAPT0 TOTAL ENERGY adz', 'SAPT0 TOTAL ENERGY aqz', 'SAPT0 TOTAL ENERGY atqz', 'SAPT0 TOTAL ENERGY atz', 'SAPT0 TOTAL ENERGY jdz']
    sapt23_cols = ['SAPT2+3(CCD) TOTAL ENERGY adtz', 'SAPT2+3(CCD) TOTAL ENERGY adz', 'SAPT2+3(CCD) TOTAL ENERGY aqz', 'SAPT2+3(CCD) TOTAL ENERGY atqz', 'SAPT2+3(CCD) TOTAL ENERGY atz', 'SAPT2+3(CCD) TOTAL ENERGY jdz']
    benchmark_col = 'benchmark ref energy'
    conversion_factor = 627.5

    plot_data = []
    for col in sapt0_cols:
        basis = col.replace('SAPT0 TOTAL ENERGY ', '')
        error = (df[col] - df[benchmark_col]) * conversion_factor
        for val in error.dropna():
            plot_data.append({'Method': 'SAPT0', 'Basis': basis, 'Error': val})

    for col in sapt23_cols:
        basis = col.replace('SAPT2+3(CCD) TOTAL ENERGY ', '')
        error = (df[col] - df[benchmark_col]) * conversion_factor
        for val in error.dropna():
            plot_data.append({'Method': 'SAPT2+3(CCD)', 'Basis': basis, 'Error': val})

    plot_df = pd.DataFrame(plot_data)

    plt.figure(figsize=(20, 10))
    sns.set(font_scale=1.5)
    ax = sns.violinplot(x='Basis', y='Error', hue='Method', data=plot_df, split=True)

    # Calculate and annotate MAE
    basis_sets = sorted(plot_df['Basis'].unique())
    methods = sorted(plot_df['Method'].unique())
    for i, basis in enumerate(basis_sets):
        for j, method in enumerate(methods):
            subset = plot_df[(plot_df['Basis'] == basis) & (plot_df['Method'] == method)]
            mae = subset['Error'].abs().mean()
            x_pos = i + (j * 0.2) - 0.1
            y_pos = 1.8 
            ax.text(x_pos, y_pos, f'MAE: {mae:.2f}', ha='center', va='center', fontsize=12, rotation=90)


    plt.ylabel('Energy Error (kcal/mol)', fontsize=18)
    plt.xlabel('Basis', fontsize=18)
    plt.title('Comparison of SAPT Energy Errors with Different Basis Sets', fontsize=22)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(-2, 2)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a violin plot of SAPT errors.')
    parser.add_argument('input_file', type=str, help='Path to the input pickle file.')
    parser.add_argument('output_file', type=str, help='Path to save the output plot.')
    args = parser.parse_args()
    main(args.input_file, args.output_file)
