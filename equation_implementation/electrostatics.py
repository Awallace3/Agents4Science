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

def elst_damping_z_mtp_f1(alpha_j, r):
    """
    # Z-MTP interaction from CLIFF
    lam_1 is for Z-q, lam_3 is for Z-mu
    """
    # Z-MTP interaction from CLIFF
    lam_1 = 1.0 - np.exp(-1.0 * np.multiply(alpha_j, r))
    lam_3 = 1.0 - (1.0 + np.multiply(alpha_j, r)) * np.exp(
        -1.0 * np.multiply(alpha_j, r)
    )
    return lam_1, lam_3

def elst_damping_mtp_mtp_f2(alpha_i, alpha_j, r):
    """
    # MTP-MTP interaction from CLIFF
    lam_1 is for q-q, lam_3 is for q-mu, and lam_5 is for mu-mu
    """
    r2 = r**2
    r3 = r2 * r
    a1_2 = alpha_i * alpha_i
    a2_2 = alpha_j * alpha_j
    a1_3 = a1_2 * alpha_i
    e1r = np.exp(-1.0 * alpha_i * r)
    e2r = np.exp(-1.0 * alpha_j * r)
    lam1, lam3, lam5 = (1.0, 1.0, 1.0)
    if abs(alpha_i - alpha_j) > 1e-6:
        A = a2_2 / (a2_2 - a1_2)
        B = a1_2 / (a1_2 - a2_2)
        lam1 -= A * e1r
        lam1 -= B * e2r
        lam3 -= (1.0 + alpha_i * r) * A * e1r
        lam3 -= (1.0 + alpha_j * r) * B * e2r

        lam5 -= (1.0 + alpha_i * r + (1.0 / 3.0) * a1_2 * r2) * A * e1r
        lam5 -= (1.0 + alpha_j * r + (1.0 / 3.0) * a2_2 * r2) * B * e2r
    else:
        lam1 -= (1.0 + 0.5 * alpha_i * r) * e1r
        lam3 -= (1.0 + alpha_i * r + 0.5 * a1_2 * r2) * e1r
        lam5 -= (1.0 + alpha_i * r + 0.5 * a1_2 * r2 + (1.0 / 6.0) * a1_3 * r3) * e1r

    return lam1, lam3, lam5

def calculate_damped_electrostatics(qcel_mol, qA, muA, qB, muB, alphaA, alphaB):
    """
    Calculates the damped electrostatic interaction energy between two molecules (A and B)
    based on Equation (4) from Schriber et al., J. Chem. Phys. 154, 184110 (2021).
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

    E1 = 0.0  # E_ZA_ZB (undamped)
    E2 = 0.0  # E_ZA_MB (damped)
    E3 = 0.0  # E_ZB_MA (damped)
    E4 = 0.0  # MTP_MTP (damped)

    for i in range(len(coordsA)):
        for j in range(len(coordsB)):
            r_ij_vec = coordsB[j] - coordsA[i]
            r_ij_sq = np.dot(r_ij_vec, r_ij_vec)
            r_ij = np.sqrt(r_ij_sq)
            u_ij = r_ij_vec / r_ij

            # Term 1: Nuclear-Nuclear interaction (undamped)
            E1 += ZA[i] * ZB[j] / r_ij

            # Term 2: Nuclear-Multipole interaction (damped)
            lam1_f1, lam3_f1 = elst_damping_z_mtp_f1(alphaB[j], r_ij)
            E2_q = ZA[i] * (qB_electronic[j] / r_ij) * lam1_f1
            E2_mu = ZA[i] * (-np.dot(muB[j], u_ij) / r_ij_sq) * lam3_f1
            E2 += E2_q + E2_mu

            # Term 3: Multipole-Nuclear interaction (damped)
            lam1_f1_i, lam3_f1_i = elst_damping_z_mtp_f1(alphaA[i], r_ij)
            E3_q = ZB[j] * (qA_electronic[i] / r_ij) * lam1_f1_i
            E3_mu = ZB[j] * (np.dot(muA[i], u_ij) / r_ij_sq) * lam3_f1_i
            E3 += E3_q + E3_mu

            # Term 4: Multipole-Multipole interaction (damped)
            lam1_f2, lam3_f2, lam5_f2 = elst_damping_mtp_mtp_f2(alphaA[i], alphaB[j], r_ij)
            E4_qq = (qA_electronic[i] * qB_electronic[j] / r_ij) * lam1_f2
            E4_qmu = (qA_electronic[i] * (-np.dot(muB[j], u_ij) / r_ij_sq)) * lam3_f2
            E4_muq = (qB_electronic[j] * (np.dot(muA[i], u_ij) / r_ij_sq)) * lam3_f2
            r_ij_cub = r_ij_sq * r_ij
            dd_term = (np.dot(muA[i], muB[j]) - 3 * np.dot(muA[i], u_ij) * np.dot(muB[j], u_ij)) / r_ij_cub
            E4_mumu = dd_term * lam5_f2
            E4 += E4_qq + E4_qmu + E4_muq + E4_mumu

    # Convert from Hartree to kcal/mol
    E1_kcal = E1 * hartree_to_kcalmol
    E2_kcal = E2 * hartree_to_kcalmol
    E3_kcal = E3 * hartree_to_kcalmol
    E4_kcal = E4 * hartree_to_kcalmol
    
    total_elst_kcal = E1_kcal + E2_kcal + E3_kcal + E4_kcal

    return total_elst_kcal, E1_kcal, E2_kcal, E3_kcal, E4_kcal

def elst_damping_mtp_mtp_f2_vectorized(alpha_i, alpha_j, r):
    """
    Vectorized version of elst_damping_mtp_mtp_f2.
    """
    r2 = r**2
    r3 = r2 * r
    a1_2 = alpha_i**2
    a2_2 = alpha_j**2
    a1_3 = a1_2 * alpha_i
    e1r = np.exp(-1.0 * alpha_i * r)
    e2r = np.exp(-1.0 * alpha_j * r)

    cond = np.abs(alpha_i - alpha_j) > 1e-6
    
    # Initialize lam arrays
    lam1 = np.ones_like(r)
    lam3 = np.ones_like(r)
    lam5 = np.ones_like(r)

    # Case 1: alpha_i != alpha_j
    A = a2_2 / (a2_2 - a1_2)
    B = a1_2 / (a1_2 - a2_2)
    
    lam1_case1 = 1.0 - A * e1r - B * e2r
    lam3_case1 = 1.0 - (1.0 + alpha_i * r) * A * e1r - (1.0 + alpha_j * r) * B * e2r
    lam5_case1 = 1.0 - (1.0 + alpha_i * r + (1.0 / 3.0) * a1_2 * r2) * A * e1r - (1.0 + alpha_j * r + (1.0 / 3.0) * a2_2 * r2) * B * e2r

    # Case 2: alpha_i == alpha_j
    lam1_case2 = 1.0 - (1.0 + 0.5 * alpha_i * r) * e1r
    lam3_case2 = 1.0 - (1.0 + alpha_i * r + 0.5 * a1_2 * r2) * e1r
    lam5_case2 = 1.0 - (1.0 + alpha_i * r + 0.5 * a1_2 * r2 + (1.0 / 6.0) * a1_3 * r3) * e1r

    lam1 = np.where(cond, lam1_case1, lam1_case2)
    lam3 = np.where(cond, lam3_case1, lam3_case2)
    lam5 = np.where(cond, lam5_case1, lam5_case2)

    return lam1, lam3, lam5

def calculate_damped_electrostatics_vectorized(qcel_mol, qA, muA, qB, muB, alphaA, alphaB):
    """
    Vectorized implementation of the damped electrostatic interaction energy calculation.
    """
    hartree_to_kcalmol = qcel.constants.hartree2kcalmol

    # Get geometry and nuclear charges
    geom = qcel_mol.geometry
    Z = qcel_mol.atomic_numbers
    
    coordsA = geom[:3]
    coordsB = geom[3:]
    ZA = np.array(Z[:3])
    ZB = np.array(Z[3:])

    qA_electronic = qA - ZA
    qB_electronic = qB - ZB

    # Create pairwise arrays
    r_vec = coordsB[None, :, :] - coordsA[:, None, :]
    r_sq = np.sum(r_vec**2, axis=2)
    r = np.sqrt(r_sq)
    u = r_vec / r[..., None]

    # Reshape parameters for broadcasting
    ZA_p = ZA[:, None]
    ZB_p = ZB[None, :]
    qA_p = qA_electronic[:, None]
    qB_p = qB_electronic[None, :]
    alphaA_p = alphaA[:, None]
    alphaB_p = alphaB[None, :]
    muA_p = muA[:, None, :]
    muB_p = muB[None, :, :]

    # Term 1: Nuclear-Nuclear interaction (undamped)
    E1 = np.sum(ZA_p * ZB_p / r)

    # Term 2: Nuclear-Multipole interaction (damped)
    lam1_f1, lam3_f1 = elst_damping_z_mtp_f1(alphaB_p, r)
    E2_q = ZA_p * (qB_p / r) * lam1_f1
    E2_mu = ZA_p * (-np.sum(muB_p * u, axis=2) / r_sq) * lam3_f1
    E2 = np.sum(E2_q + E2_mu)

    # Term 3: Multipole-Nuclear interaction (damped)
    lam1_f1_i, lam3_f1_i = elst_damping_z_mtp_f1(alphaA_p, r)
    E3_q = ZB_p * (qA_p / r) * lam1_f1_i
    E3_mu = ZB_p * (np.sum(muA_p * u, axis=2) / r_sq) * lam3_f1_i
    E3 = np.sum(E3_q + E3_mu)

    # Term 4: Multipole-Multipole interaction (damped)
    lam1_f2, lam3_f2, lam5_f2 = elst_damping_mtp_mtp_f2_vectorized(alphaA_p, alphaB_p, r)
    E4_qq = (qA_p * qB_p / r) * lam1_f2
    E4_qmu = (qA_p * (-np.sum(muB_p * u, axis=2) / r_sq)) * lam3_f2
    E4_muq = (qB_p * (np.sum(muA_p * u, axis=2) / r_sq)) * lam3_f2
    r_cub = r_sq * r
    dd_term = (np.sum(muA_p * muB_p, axis=2) - 3 * np.sum(muA_p * u, axis=2) * np.sum(muB_p * u, axis=2)) / r_cub
    E4_mumu = dd_term * lam5_f2
    E4 = np.sum(E4_qq + E4_qmu + E4_muq + E4_mumu)

    # Convert from Hartree to kcal/mol
    E1_kcal = E1 * hartree_to_kcalmol
    E2_kcal = E2 * hartree_to_kcalmol
    E3_kcal = E3 * hartree_to_kcalmol
    E4_kcal = E4 * hartree_to_kcalmol
    
    total_elst_kcal = E1_kcal + E2_kcal + E3_kcal + E4_kcal

    return total_elst_kcal, E1_kcal, E2_kcal, E3_kcal, E4_kcal
