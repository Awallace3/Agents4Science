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
    # Extract atom coordinates and nuclear charges
    coords = qcel_mol.geometry
    atomic_numbers = qcel_mol.atomic_numbers

    # Assuming the first fragment is A and the second is B
    frag_A_indices = qcel_mol.fragments[0]
    frag_B_indices = qcel_mol.fragments[1]

    coords_A = coords[frag_A_indices]
    coords_B = coords[frag_B_indices]
    
    ZA_nuclear = atomic_numbers[frag_A_indices]
    ZB_nuclear = atomic_numbers[frag_B_indices]

    # Adjust qA and qB by subtracting nuclear charges to get electronic partial charges
    # This assumes qA and qB provided are total charges (nuclear + electronic)
    qA_electronic = qA - ZA_nuclear
    qB_electronic = qB - ZB_nuclear

    total_energy = 0.0

    for i in range(len(coords_A)):
        for j in range(len(coords_B)):
            r_vec = coords_B[j] - coords_A[i]
            r_ij = np.linalg.norm(r_vec)

            if r_ij < 1e-12: # Avoid division by zero for self-interaction or coincident atoms
                continue

            # Damping for q-q interaction (lam1 from elst_damping_mtp_mtp)
            # Here, q-q refers to electronic partial charge - electronic partial charge interaction
            lam_mtp_mtp_qq = elst_damping_mtp_mtp(alphaA[i], alphaB[j], r_ij)[0] # lam1
            term_qq = (qA_electronic[i] * qB_electronic[j] / r_ij) * lam_mtp_mtp_qq
            total_energy += term_qq

            # Damping for Z-mu interaction (lam_3 from elst_damping_z_mtp)
            # Z is electronic partial charge, mu is dipole
            lam_z_mtp_B_zmu = elst_damping_z_mtp(alphaB[j], r_ij)[1] # lam_3
            term_qd_undamped = qA_electronic[i] * np.dot(muB[j], r_vec) / (r_ij**3)
            term_qd_damped = term_qd_undamped * lam_z_mtp_B_zmu
            total_energy += term_qd_damped

            # Damping for Z-mu interaction (lam_3 from elst_damping_z_mtp)
            # Z is electronic partial charge, mu is dipole
            lam_z_mtp_A_zmu = elst_damping_z_mtp(alphaA[i], r_ij)[1] # lam_3
            term_dq_undamped = -qB_electronic[j] * np.dot(muA[i], r_vec) / (r_ij**3)
            term_dq_damped = term_dq_undamped * lam_z_mtp_A_zmu
            total_energy += term_dq_damped

            # Damping for mu-mu interaction (lam5 from elst_damping_mtp_mtp)
            # mu is dipole
            lam_mtp_mtp_mumu = elst_damping_mtp_mtp(alphaA[i], alphaB[j], r_ij)[2] # lam5
            term_dd_undamped = (np.dot(muA[i], muB[j]) / (r_ij**3)) - \
                               (3.0 * np.dot(muA[i], r_vec) * np.dot(muB[j], r_vec) / (r_ij**5))
            term_dd_damped = term_dd_undamped * lam_mtp_mtp_mumu
            total_energy += term_dd_damped

    return total_energy * HARTREE_TO_KCAL_MOL