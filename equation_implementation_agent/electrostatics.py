import numpy as np
import qcelemental as qcel

def calculate_electrostatics(qcel_mol, qA, muA, qB, muB):
    """
    Calculates the electrostatic interaction energy between two molecules (A and B)
    based on Equation (4) from Schriber et al., J. Chem. Phys. 154, 184110 (2021).

    This implementation is for the undamped case (f_1=1 and f_2=1).
    """
    hartree_to_kcalmol = qcel.constants.hartree2kcalmol

    # Get geometry and nuclear charges
    geom = qcel_mol.geometry
    Z = qcel_mol.atomic_numbers
    
    # Assuming first 3 atoms are monomer A, next 3 are monomer B
    coordsA = geom[:3]
    coordsB = geom[3:]
    ZA = np.array(Z[:3])
    ZB = np.array(Z[3:])

    qA_electronic = qA - ZA
    qB_electronic = qB - ZB

    E1 = 0.0  # E_ZA_ZB
    E2 = 0.0  # E_ZA_MB
    E3 = 0.0  # E_ZB_MA
    E4 = 0.0  # MTP_MTP

    for i in range(len(coordsA)):
        for j in range(len(coordsB)):
            r_ij_vec = coordsB[j] - coordsA[i]
            r_ij_sq = np.dot(r_ij_vec, r_ij_vec)
            r_ij = np.sqrt(r_ij_sq)
            u_ij = r_ij_vec / r_ij

            # Term 1: \sum_{i \in A, j \in B} Z_i Z_j / r_{ij}
            E1 += ZA[i] * ZB[j] / r_ij

            # Term 2: \sum_{i \in A, j \in B} Z_i T_{ij}^1 M_j
            E2 += ZA[i] * (qB_electronic[j] / r_ij - np.dot(muB[j], u_ij) / r_ij_sq)

            # Term 3: \sum_{i \in A, j \in B} M_i^T T_{ij}^1 Z_j
            E3 += ZB[j] * (qA_electronic[i] / r_ij + np.dot(muA[i], u_ij) / r_ij_sq)

            # Term 4: \sum_{i \in A, j \in B} M_i^T T_{ij}^2 M_j
            # charge-charge
            E4 += qA_electronic[i] * qB_electronic[j] / r_ij
            # charge-dipole
            E4 += qA_electronic[i] * (-np.dot(muB[j], u_ij) / r_ij_sq)
            # dipole-charge
            E4 += qB_electronic[j] * (np.dot(muA[i], u_ij) / r_ij_sq)
            # dipole-dipole
            r_ij_cub = r_ij_sq * r_ij
            dd_term = (np.dot(muA[i], muB[j]) - 3 * np.dot(muA[i], u_ij) * np.dot(muB[j], u_ij)) / r_ij_cub
            E4 += dd_term

    # Convert from Hartree to kcal/mol
    E1_kcal = E1 * hartree_to_kcalmol
    E2_kcal = E2 * hartree_to_kcalmol
    E3_kcal = E3 * hartree_to_kcalmol
    E4_kcal = E4 * hartree_to_kcalmol
    
    total_elst_kcal = E1_kcal + E2_kcal + E3_kcal + E4_kcal

    return total_elst_kcal, E1_kcal, E2_kcal, E3_kcal, E4_kcal