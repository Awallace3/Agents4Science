import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pprint import pprint as pp
import qcelemental as qcel
from qm_tools_aw import tools

h2kcalmol = qcel.constants.conversion_factor("hartree", "kcal/mol")


def qcel_mols(row):
    """
    Convert the row to a qcel molecule
    """
    atomic_numbers = [
        row["atomic_numbers"][row["monAs"]],
        row["atomic_numbers"][row["monBs"]],
    ]
    coords = [row["coordinates"][row["monAs"]], row["coordinates"][row["monBs"]]]
    cm = [
        [row["monA_charge"], row["monA_multiplicity"]],
        [row["monB_charge"], row["monB_multiplicity"]],
    ]
    return tools.convert_pos_carts_to_mol(atomic_numbers, coords, cm)


def isolate_curve(df, system_id_contains: str = "01_Water-Water"):
    """
    Isolate a subset of the DataFrame for a specific system_id.
    """
    df = pd.read_pickle("./combined_df_4569.pkl")
    df_subset = df[df["system_id"].str.contains(system_id_contains)].copy()
    print(df_subset)
    df_subset['qcel_molecule'] = df_subset.apply(qcel_mols, axis=1)
    print(df_subset[['system_id', 'SAPT0 ELST ENERGY adz']])
    return


def main():
    isolate_curve()
    return


if __name__ == "__main__":
    main()


def isolate_curve(df, system_id_contains: str = "01_Water-Water"):
    """
    Isolate a subset of the DataFrame for a specific system_id.
    """
    df_subset = df[df["system_id"].str.contains(system_id_contains)].copy()
    print(df_subset[["R", "ELST"]])
    # Plot Potential Energy Curve
    df_subset["R"] = df_subset["R"] / 0.529177249
    df_subset.sort_values(by="R", inplace=True)
    fig = plt.figure(dpi=300)
    plt.plot(
        df_subset["R"],
        df_subset["total"],
        marker="o",
        c="black",
        label="Total",
    )
    plt.plot(
        df_subset["R"],
        df_subset["ELST"],
        marker="o",
        c="red",
        label="ELST",
    )
    plt.plot(
        df_subset["R"],
        df_subset["IND"],
        marker="o",
        c="blue",
        label="IND",
    )
    plt.plot(
        df_subset["R"],
        df_subset["EXCH"],
        marker="o",
        c="green",
        label="EXCH",
    )
    plt.plot(
        df_subset["R"],
        df_subset["DISP"],
        marker="o",
        c="orange",
        label="DISP",
    )
    plt.legend()
    plt.minorticks_on()
    plt.grid(axis="both")
    plt.xlabel(r"Distance ($\AA$)")
    plt.ylabel(r"Interaction Energy (kcal/mol)")
    plt.savefig(f"plots/PES_{system_id_contains}.png")
    return df_subset


def isolate():
    df_LoS = pd.read_pickle("./combined_df_subset_358.pkl")
    df_LoS["qcel_molecule"] = df_LoS.apply(qcel_mols, axis=1)
    pprint(df_LoS.columns.tolist())
    pd.set_option("display.max_rows", None)
    print(df_LoS[["system_id", "SAPT0 ELST ENERGY adz"]])
    df_LoS["ELST"] = df_LoS["SAPT0 ELST ENERGY adz"] * h2kcalmol
    df_LoS["ELST_HL"] = df_LoS["SAPT2+3(CCD)DMP2 ELST ENERGY aqz"] * h2kcalmol
    df_LoS["IND"] = df_LoS["SAPT0 IND ENERGY adz"] * h2kcalmol
    df_LoS["EXCH"] = df_LoS["SAPT0 EXCH ENERGY adz"] * h2kcalmol
    df_LoS["DISP"] = df_LoS["SAPT0 DISP ENERGY adz"] * h2kcalmol
    df_LoS["total"] = df_LoS["ELST"] + df_LoS["IND"] + df_LoS["EXCH"] + df_LoS["DISP"]
    # Collect df subset for a Potential Energy Curve
    df_water_water = isolate_curve(df_LoS, system_id_contains="01_Water-Water")
    df_water_water.to_pickle("df_water_water.pkl")
    df_bz_meoh = isolate_curve(df_LoS, system_id_contains="55_Benzene-MeOH_OH-pi")
    df_bz_meoh.to_pickle("df_bz_meoh.pkl")

    # get last df_water_water row and print
    last_row = df_water_water.iloc[-1]
    print("Last row of df_water_water:")
    print(last_row)
    print(last_row["qcel_molecule"].to_string("psi4"))
    return
