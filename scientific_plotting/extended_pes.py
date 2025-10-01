import pandas as pd
import qcelemental as qcel
from qm_tools_aw import tools
import numpy as np
import psi4
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import imageio
import os

def min_contact_distance(molecule: qcel.models.Molecule) -> float:
    """
    Calculate the minimum contact distance between two fragments.
    """
    frag1_coords = molecule.geometry[molecule.fragments[0]]
    frag2_coords = molecule.geometry[molecule.fragments[1]]
    return np.min(cdist(frag1_coords, frag2_coords))


def isolate_curve(df: pd.DataFrame, system_id_contains: str = "01_Water-Water"):
    """
    Isolate a subset of the DataFrame for a specific system_id.
    """
    df_subset = df[df["system_id"].str.contains(system_id_contains)].copy()

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

    df_subset["qcel_molecule"] = df_subset.apply(qcel_mols, axis=1)
    return df_subset


def displacement_geometry(
    df_subset: pd.DataFrame,
    starting_distance: float,
    ending_distance: float,
    increment: float,
):
    """
    Generate new geometries by displacing a fragment along a vector.
    """
    geom1 = df_subset.iloc[0]["qcel_molecule"]
    geom2 = df_subset.iloc[1]["qcel_molecule"]
    
    frag1_geom1 = geom1.geometry[geom1.fragments[0]]
    frag2_geom1 = geom1.geometry[geom1.fragments[1]]
    
    frag1_geom2 = geom2.geometry[geom2.fragments[0]]
    frag2_geom2 = geom2.geometry[geom2.fragments[1]]

    dist1 = np.linalg.norm(frag1_geom1 - frag1_geom2)
    dist2 = np.linalg.norm(frag2_geom1 - frag2_geom2)

    if dist1 > dist2:
        moving_frag_indices = geom1.fragments[0]
        static_frag_indices = geom1.fragments[1]
        displacement_vector = np.mean(frag1_geom2 - frag1_geom1, axis=0)
    else:
        moving_frag_indices = geom1.fragments[1]
        static_frag_indices = geom1.fragments[0]
        displacement_vector = np.mean(frag2_geom2 - frag2_geom1, axis=0)

    displacement_vector /= np.linalg.norm(displacement_vector)
    
    new_rows = []
    for dist in np.arange(starting_distance, ending_distance + increment, increment):
        dist_bohr = dist * qcel.constants.conversion_factor("angstrom", "bohr")
        
        new_coords = geom1.geometry.copy()
        new_coords[moving_frag_indices] += displacement_vector * dist_bohr
        
        new_geom = qcel.models.Molecule(
            symbols=geom1.symbols,
            geometry=new_coords,
            fragments=geom1.fragments,
            molecular_charge=geom1.molecular_charge,
            molecular_multiplicity=geom1.molecular_multiplicity,
            real=geom1.real,
            comment=geom1.comment,
            provenance=geom1.provenance.dict(),
            connectivity=geom1.connectivity,
        )
        
        new_row = {col: np.nan for col in df_subset.columns}
        new_row["qcel_molecule"] = new_geom
        new_row["displacement_distance"] = dist
        new_row["min_contact_distance"] = min_contact_distance(new_geom)
        new_rows.append(new_row)

    return pd.concat([df_subset, pd.DataFrame(new_rows)], ignore_index=True)


def run_psi4_sapt(molecule: qcel.models.Molecule) -> dict:
    """
    Run SAPT0 calculation using Psi4.
    """
    psi4.core.set_output_file("psi4_output.dat", False)
    mol_str = molecule.to_string("psi4")
    psi4.geometry(mol_str)
    psi4.set_options({"basis": "aug-cc-pVDZ", "freeze_core": "true"})
    try:
        energy = psi4.energy("sapt0")
        sapt_components = {
            "SAPT0 Total": energy,
            "SAPT ELST": psi4.variable("SAPT0 ELST ENERGY"),
            "SAPT EXCH": psi4.variable("SAPT0 EXCH ENERGY"),
            "SAPT IND": psi4.variable("SAPT0 IND ENERGY"),
            "SAPT DISP": psi4.variable("SAPT0 DISP ENERGY"),
        }
    except psi4.PsiException as e:
        print(f"Psi4 calculation failed: {e}")
        sapt_components = {
            "SAPT0 Total": np.nan,
            "SAPT ELST": np.nan,
            "SAPT EXCH": np.nan,
            "SAPT IND": np.nan,
            "SAPT DISP": np.nan,
        }
    return sapt_components


def plot_pes(df: pd.DataFrame):
    """
    Plot the SAPT0 potential energy surface.
    """
    kcal_per_mol = qcel.constants.conversion_factor("hartree", "kcal / mol")
    
    original_df = df.dropna(subset=["SAPT0 TOTAL ENERGY adz"])
    new_df = df[df["SAPT0 TOTAL ENERGY adz"].isnull()]

    plt.figure(figsize=(10, 6))

    for df_plot, marker in [(original_df, "X"), (new_df, "o")]:
        if not df_plot.empty:
            distances = df_plot["min_contact_distance"]
            plt.plot(
                distances,
                df_plot["SAPT0 Total"] * kcal_per_mol,
                marker=marker,
                linestyle="-",
                label="Total",
                color="black",
            )
            plt.plot(
                distances,
                df_plot["SAPT ELST"] * kcal_per_mol,
                marker=marker,
                linestyle="--",
                label="Elst",
                color="red",
            )
            plt.plot(
                distances,
                df_plot["SAPT EXCH"] * kcal_per_mol,
                marker=marker,
                linestyle="--",
                label="Exch",
                color="green",
            )
            plt.plot(
                distances,
                df_plot["SAPT IND"] * kcal_per_mol,
                marker=marker,
                linestyle="--",
                label="Ind",
                color="blue",
            )
            plt.plot(
                distances,
                df_plot["SAPT DISP"] * kcal_per_mol,
                marker=marker,
                linestyle="--",
                label="Disp",
                color="orange",
            )

    plt.xlabel("Minimum Contact Distance (Å)")
    plt.ylabel("Energy (kcal/mol)")
    plt.title("SAPT0/aug-cc-pVDZ Potential Energy Surface")
    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    
    plt.grid(True)
    plt.savefig("water_pes.png")
    plt.show()

def create_pes_gif(df: pd.DataFrame):
    """
    Create a GIF of the water PES to visualize the waters being pulled apart.
    """
    kcal_per_mol = qcel.constants.conversion_factor("hartree", "kcal / mol")
    
    original_df = df.dropna(subset=["SAPT0 TOTAL ENERGY adz"])
    new_df = df[df["SAPT0 TOTAL ENERGY adz"].isnull()]

    temp_dir = "gif_frames"
    os.makedirs(temp_dir, exist_ok=True)
    
    frames = []
    for i, (index, row) in enumerate(new_df.iterrows()):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))

        # Plot PES
        for df_plot, marker in [(original_df, "X"), (new_df, "o")]:
            if not df_plot.empty:
                distances = df_plot["min_contact_distance"]
                ax1.plot(
                    distances,
                    df_plot["SAPT0 Total"] * kcal_per_mol,
                    marker=marker,
                    linestyle="-",
                    label="Total",
                    color="black",
                )
                ax1.plot(
                    distances,
                    df_plot["SAPT ELST"] * kcal_per_mol,
                    marker=marker,
                    linestyle="--",
                    label="Elst",
                    color="red",
                )
                ax1.plot(
                    distances,
                    df_plot["SAPT EXCH"] * kcal_per_mol,
                    marker=marker,
                    linestyle="--",
                    label="Exch",
                    color="green",
                )
                ax1.plot(
                    distances,
                    df_plot["SAPT IND"] * kcal_per_mol,
                    marker=marker,
                    linestyle="--",
                    label="Ind",
                    color="blue",
                )
                ax1.plot(
                    distances,
                    df_plot["SAPT DISP"] * kcal_per_mol,
                    marker=marker,
                    linestyle="--",
                    label="Disp",
                    color="orange",
                )
        
        # Highlight current point
        ax1.plot(
            row["min_contact_distance"],
            row["SAPT0 Total"] * kcal_per_mol,
            "o",
            color="purple",
            markersize=10,
            label="Current Point",
        )

        ax1.set_xlabel("Minimum Contact Distance (Å)")
        ax1.set_ylabel("Energy (kcal/mol)")
        ax1.set_title("SAPT0/aug-cc-pVDZ Potential Energy Surface")
        handles, labels = ax1.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax1.legend(by_label.values(), by_label.keys())
        ax1.grid(True)

        # Plot molecule
        molecule = row["qcel_molecule"]
        coords = molecule.geometry
        symbols = molecule.symbols
        
        # Get atomic radii for plotting
        atomic_radii = {
            "H": 0.31,  # Covalent radius in Angstroms
            "O": 0.66,  # Covalent radius in Angstroms
        } # This is a simplification, real radii are more complex

        for j, (symbol, coord) in enumerate(zip(symbols, coords)):
            color = "blue" if j in molecule.fragments[0] else "red"
            radius = atomic_radii.get(symbol, 0.5) # Default radius if not found
            circle = plt.Circle(coord[:2], radius=radius, color=color, alpha=0.6)
            ax2.add_patch(circle)
            ax2.text(coord[0], coord[1], symbol, ha="center", va="center", fontsize=10)

        ax2.set_aspect("equal", adjustable="box")
        ax2.set_title(f"Min Contact Dist: {row['min_contact_distance']:.2f} Å")
        ax2.axis("off")
        
        frame_path = os.path.join(temp_dir, f"frame_{i:03d}.png")
        plt.savefig(frame_path)
        plt.close(fig)
        frames.append(imageio.v2.imread(frame_path))

    imageio.mimsave("water_pes_animation.gif", frames, fps=5)
    
    # Clean up temporary frames
    for f in os.listdir(temp_dir):
        os.remove(os.path.join(temp_dir, f))
    os.rmdir(temp_dir)