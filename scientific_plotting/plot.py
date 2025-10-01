
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import mean_absolute_error, mean_squared_error
import qcelemental as qcel

def plot_violin(df):
    """Generates a violin plot comparing SAPT0 and SAPT2+3(CCD) errors"""
    plt.figure(figsize=(10, 6))
    sapt0_col = "SAPT0 TOTAL ENERGY adz"
    sapt23_col = "SAPT2+3(CCD) TOTAL ENERGY adz"
    benchmark_col = "Benchmark"

    # Unit conversion
    hartree_to_kcalmol = qcel.constants.conversion_factor("hartree", "kcal/mol")
    sapt0_kcalmol = "SAPT0 (kcal/mol)"
    sapt23_kcalmol = "SAPT2+3(CCD) (kcal/mol)"
    df[sapt0_kcalmol] = df[sapt0_col] * hartree_to_kcalmol
    df[sapt23_kcalmol] = df[sapt23_col] * hartree_to_kcalmol

    # Calculate errors
    sapt0_error_col = "SAPT0 Error (kcal/mol)"
    sapt23_error_col = "SAPT2+3(CCD) Error (kcal/mol)"
    df[sapt0_error_col] = df[sapt0_kcalmol] - df[benchmark_col]
    df[sapt23_error_col] = df[sapt23_kcalmol] - df[benchmark_col]

    sns.violinplot(data=df[[sapt0_error_col, sapt23_error_col]], inner="quartile")
    plt.title("SAPT0 vs SAPT2+3(CCD) Error Distribution (adz basis set)")
    plt.ylabel("Error (kcal/mol)")

    # Calculate statistics
    mae_sapt0 = mean_absolute_error(df[benchmark_col], df[sapt0_kcalmol])
    rmse_sapt0 = np.sqrt(mean_squared_error(df[benchmark_col], df[sapt0_kcalmol]))
    corr_sapt0 = df[benchmark_col].corr(df[sapt0_kcalmol])

    mae_sapt23 = mean_absolute_error(df[benchmark_col], df[sapt23_kcalmol])
    rmse_sapt23 = np.sqrt(mean_squared_error(df[benchmark_col], df[sapt23_kcalmol]))
    corr_sapt23 = df[benchmark_col].corr(df[sapt23_kcalmol])

    # Annotate plot with statistics
    stats_text = (
        f"SAPT0:\n" 
        f"MAE: {mae_sapt0:.2f}\n"
        f"RMSE: {rmse_sapt0:.2f}\n"
        f"Correlation: {corr_sapt0:.2f}\n\n"
        f"SAPT2+3(CCD):\n"
        f"MAE: {mae_sapt23:.2f}\n"
        f"RMSE: {rmse_sapt23:.2f}\n"
        f"Correlation: {corr_sapt23:.2f}"
    )
    plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.savefig("plot.png")

def main():
    """Main function to generate plot"""
    df = pd.read_pickle("combined_df_4569.pkl")
    plot_violin(df)

if __name__ == "__main__":
    main()
