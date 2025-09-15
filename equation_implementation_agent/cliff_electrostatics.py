
import numpy as np

def _get_g_derivs(k, r, f_type):
    """
    Calculates the derivatives of the damped interaction function g(r) = f(r)/r.
    """
    if r == 0:
        return 0, 0, 0

    exp_kr = np.exp(-k * r)

    if f_type == 1:  # f1(r) = 1 - exp(-k*r)
        f = 1 - exp_kr
        fp = k * exp_kr
        fpp = -k**2 * exp_kr
    elif f_type == 2:  # f2(r) = 1 - (1 + k*r + 0.5*(k*r)**2)*exp(-k*r)
        p = 1 + k * r + 0.5 * (k * r)**2
        f = 1 - p * exp_kr
        fp = 0.5 * k**3 * r**2 * exp_kr
        fpp = 0.5 * k**3 * (2 * r - k * r**2) * exp_kr
    else:
        raise ValueError("f_type must be 1 or 2")

    g = f / r
    gp = (fp * r - f) / r**2
    gpp = (fpp * r**2 - 2 * fp * r + 2 * f) / r**3

    return g, gp, gpp

def cliff_electrostatics(coords_A, Z_A, M_A, Kelst_A, coords_B, Z_B, M_B, Kelst_B):
    """
    Calculates the electrostatic interaction energy between two molecules A and B
    using the CLIFF force field electrostatics term.

    This implementation includes interactions up to dipole-dipole. The full
    CLIFF model also includes interactions with quadrupoles, which are not
    implemented here due to their complexity.

    Args:
        coords_A (np.ndarray): Coordinates of atoms in molecule A (N_A, 3).
        Z_A (np.ndarray): Nuclear charges of atoms in molecule A (N_A,).
        M_A (np.ndarray): Atomic multipoles (q, mu_x, mu_y, mu_z) for molecule A (N_A, 4).
        Kelst_A (np.ndarray): Atomic fitting parameters for molecule A (N_A,).
        coords_B (np.ndarray): Coordinates of atoms in molecule B (N_B, 3).
        Z_B (np.ndarray): Nuclear charges of atoms in molecule B (N_B,).
        M_B (np.ndarray): Atomic multipoles (q, mu_x, mu_y, mu_z) for molecule B (N_B, 4).
        Kelst_B (np.ndarray): Atomic fitting parameters for molecule B (N_B,).

    Returns:
        float: The electrostatic interaction energy.
    """
    energy = 0.0
    for i in range(len(coords_A)):
        for j in range(len(coords_B)):
            r_vec = coords_B[j] - coords_A[i]
            r = np.linalg.norm(r_vec)
            if r == 0:
                continue

            # Combine Kelst parameters (using arithmetic mean as it's not specified)
            kelst = 0.5 * (Kelst_A[i] + Kelst_B[j])

            # 1. Nuclear-Nuclear interaction (undamped)
            energy += Z_A[i] * Z_B[j] / r

            # Get multipoles for atoms i and j
            q_i, mu_i = M_A[i, 0], M_A[i, 1:]
            q_j, mu_j = M_B[j, 0], M_B[j, 1:]

            # --- Damped Interactions ---
            # T1 for nuclear-multipole, T2 for multipole-multipole
            g1, g1p, g1pp = _get_g_derivs(kelst, r, 1)
            g2, g2p, g2pp = _get_g_derivs(kelst, r, 2)

            # Tensors for T1
            T1_00 = g1
            T1_0a = g1p * r_vec / r
            T1_ab = (g1pp - g1p / r) * np.outer(r_vec, r_vec) / r**2 + (g1p / r) * np.identity(3)

            # Tensors for T2
            T2_00 = g2
            T2_0a = g2p * r_vec / r
            T2_ab = (g2pp - g2p / r) * np.outer(r_vec, r_vec) / r**2 + (g2p / r) * np.identity(3)

            # 2. Nuclear-Electron interaction (Z_i with M_j)
            energy += Z_A[i] * (q_j * T1_00 - np.dot(mu_j, T1_0a))

            # 3. Electron-Nuclear interaction (M_i with Z_j)
            energy += (q_i * T1_00 + np.dot(mu_i, T1_0a)) * Z_B[j]

            # 4. Electron-Electron interaction (M_i with M_j)
            energy += q_i * q_j * T2_00
            energy += q_i * np.dot(mu_j, -T2_0a)
            energy += np.dot(mu_i, T2_0a) * q_j
            energy += np.dot(mu_i, T2_ab @ mu_j)

    return energy

if __name__ == '__main__':
    # Example: a water dimer
    # Molecule A (atom 0: O, atom 1: H, atom 2: H)
    coords_A = np.array([
        [0.0000, 0.0000, 0.1173],
        [0.0000, 0.7572, -0.4692],
        [0.0000, -0.7572, -0.4692]
    ])
    Z_A = np.array([8, 1, 1])
    # Charges and dipoles are zero for this example for simplicity
    M_A = np.array([
        [-0.834, 0.0, 0.0, 0.0],
        [0.417, 0.0, 0.0, 0.0],
        [0.417, 0.0, 0.0, 0.0]
    ])
    # Using Kelst values for O2 and HC from Table I in the paper
    Kelst_A = np.array([3.8700, 3.5982, 3.5982])


    # Molecule B (3 Angstroms away on x-axis)
    coords_B = coords_A + np.array([3.0, 0.0, 0.0])
    Z_B = Z_A
    M_B = M_A
    Kelst_B = Kelst_A

    electrostatic_energy = cliff_electrostatics(coords_A, Z_A, M_A, Kelst_A, coords_B, Z_B, M_B, Kelst_B)
    print(f"Electrostatic interaction energy: {electrostatic_energy:.6f} Hartree")
    # Note: The energy unit depends on the units of input.
    # If coordinates are in Angstrom and charges in elementary charge,
    # the energy will be in Hartree * a0 / Angstrom.
    # For a proper calculation, units must be handled carefully.
    # This example is for demonstration purposes.
