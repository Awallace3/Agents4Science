import numpy as np
import qcelemental as qcel

HARTREE_TO_KCAL_MOL = 627.509

def elst_damping_mtp_mtp(alpha_i, alpha_j, r):
    """
    # MTP-MTP interaction from CLIFF
    lam_1 is for q-q, lam_3 is for q-mu, and lam_5 is for mu-mu
    """
    r2 = r**2
    r3 = r2 * r
    r4 = r2**2
    r5 = r4 * r
    a1_2 = alpha_i * alpha_i
    a2_2 = alpha_j * alpha_j
    a1_3 = a1_2 * alpha_i
    a2_3 = a2_2 * alpha_j
    a1_4 = a1_3 * alpha_i
    a2_4 = a2_3 * alpha_j
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


def elst_damping_z_mtp(alpha_j, r):
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


def calculate_electrostatics_energy(qcel_mol, qA, muA, qB, muB, alphaA, alphaB):
    """
    Calculate electrostatics energy using equation 4 from Schriber et al. 2021:
    Eelst = ∑∑ [ Zi*Zj/rij + Zi*Tij^f1*Mj + Mi^T*Tij^f1*Zj + Mi^T*Tij^f2*Mj ]
    
    where:
    - Zi, Zj are nuclear charges
    - Mi, Mj are multipole moments
    - After getting Z from atomic numbers, modify qA and qB: qA[i] -= Z_A[i]
    - This makes qA/qB represent total electronic charges (negative)
    - Tij^f1, Tij^f2 are damped interaction tensors
    """
    # Extract atom coordinates and nuclear charges
    coords = qcel_mol.geometry
    atomic_numbers = qcel_mol.atomic_numbers

    # Assuming the first fragment is A and the second is B
    frag_A_indices = qcel_mol.fragments[0]
    frag_B_indices = qcel_mol.fragments[1]

    coords_A = coords[frag_A_indices]
    coords_B = coords[frag_B_indices]
    
    # Nuclear charges (Z in equation 4)
    ZA = atomic_numbers[frag_A_indices].astype(float)
    ZB = atomic_numbers[frag_B_indices].astype(float)
    
    # Modify qA and qB as per README instructions: qA[i] -= Z_A[i]
    # This converts partial atomic charges to total electronic charges
    qA_modified = qA - ZA
    qB_modified = qB - ZB

    total_energy = 0.0

    for i in range(len(coords_A)):
        for j in range(len(coords_B)):
            r_vec = coords_B[j] - coords_A[i]
            r_ij = np.linalg.norm(r_vec)

            if r_ij < 1e-12:  # Avoid division by zero
                continue
            
            # Unit vector from i to j
            r_hat = r_vec / r_ij

            # Term 1: Zi*Zj/rij (nuclear-nuclear interaction, undamped)
            term_ZZ = ZA[i] * ZB[j] / r_ij
            total_energy += term_ZZ

            # Get damping functions
            lam_z_mtp_A = elst_damping_z_mtp(alphaA[i], r_ij)  # (lam_1, lam_3) for Z_A to MTP_B
            lam_z_mtp_B = elst_damping_z_mtp(alphaB[j], r_ij)  # (lam_1, lam_3) for Z_B to MTP_A
            lam_mtp_mtp = elst_damping_mtp_mtp(alphaA[i], alphaB[j], r_ij)  # (lam_1, lam_3, lam_5)

            # Term 2: Zi*Tij^f1*Mj (nuclear i to multipole j interaction)
            # This includes: Zi*qj and Zi*mu_j interactions
            
            # Zi*qj interaction with lam_1 damping
            term_Zq = ZA[i] * qB_modified[j] / r_ij * lam_z_mtp_B[0]  # lam_1 for Z-q
            total_energy += term_Zq
            
            # Zi*mu_j interaction with lam_3 damping (dipole interaction)
            # T_{ij}^{(1)} * mu_j for dipole gives: mu_j · r_hat / r^2
            term_Zmu = ZA[i] * np.dot(muB[j], r_hat) / (r_ij**2) * lam_z_mtp_B[1]  # lam_3 for Z-mu
            total_energy += term_Zmu

            # Term 3: Mi^T*Tij^f1*Zj (multipole i to nuclear j interaction)
            # This includes: qi*Zj and mu_i*Zj interactions
            
            # qi*Zj interaction with lam_1 damping
            term_qZ = qA_modified[i] * ZB[j] / r_ij * lam_z_mtp_A[0]  # lam_1 for q-Z
            total_energy += term_qZ
            
            # mu_i*Zj interaction with lam_3 damping (dipole interaction)
            # mu_i^T * T_{ij}^{(1)} for dipole gives: -mu_i · r_hat / r^2
            term_muZ = -ZB[j] * np.dot(muA[i], r_hat) / (r_ij**2) * lam_z_mtp_A[1]  # lam_3 for mu-Z
            total_energy += term_muZ

            # Term 4: Mi^T*Tij^f2*Mj (multipole i to multipole j interaction)
            # This includes: qi*qj, qi*mu_j, mu_i*qj, and mu_i*mu_j interactions
            
            # qi*qj interaction with lam_1 damping
            term_qq = qA_modified[i] * qB_modified[j] / r_ij * lam_mtp_mtp[0]  # lam_1 for q-q
            total_energy += term_qq
            
            # qi*mu_j interaction with lam_3 damping
            term_qmu = qA_modified[i] * np.dot(muB[j], r_hat) / (r_ij**2) * lam_mtp_mtp[1]  # lam_3 for q-mu
            total_energy += term_qmu
            
            # mu_i*qj interaction with lam_3 damping
            term_muq = -qB_modified[j] * np.dot(muA[i], r_hat) / (r_ij**2) * lam_mtp_mtp[1]  # lam_3 for mu-q
            total_energy += term_muq
            
            # mu_i*mu_j interaction with lam_5 damping
            # Dipole-dipole interaction: (mu_i · mu_j)/r^3 - 3*(mu_i · r_hat)*(mu_j · r_hat)/r^3
            term_mumu = (np.dot(muA[i], muB[j]) / (r_ij**3) - 
                        3.0 * np.dot(muA[i], r_hat) * np.dot(muB[j], r_hat) / (r_ij**3)) * lam_mtp_mtp[2]  # lam_5
            total_energy += term_mumu

    return total_energy * HARTREE_TO_KCAL_MOL
