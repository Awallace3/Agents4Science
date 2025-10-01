# Prompts
```txt
# Task
implement the electrostatics equation from @schriber_2021_184110.pdf. Implement the equation using numpy and qcelemental. Use a test driven development approach to implement the correct code to match an electrostatics energy for the given system below of -10.532698 kcal/mol.

# Information for pytest

"""
elst_charges_dipoles = -10.532698 kcal/mol
qcel_mol = qcel.models.Molecule.from_data(
'''
0 1
--
0 1
O                    -1.326958220000    -0.105938540000     0.018788150000
H                    -1.931665230000     1.600174310000    -0.021710520000
H                     0.486644270000     0.079598100000     0.009862480000
--
0 1
O                     3.907523240000     0.052757410000     0.001850160000
H                     4.619234940000    -0.775660840000     1.449615410000
H                     4.611000850000    -0.847154680000    -1.406756420000
units bohr
no_com
no_reorient
'''
)


qA=np.array([-0.902827,  0.452194,  0.450633])
muA=np.array([[-0.120995, -0.19035 ,  0.004972],
       [-0.023432,  0.01783 , -0.000382],
       [ 0.026533, -0.013788,  0.000275]])
qB=np.array([-0.905161,  0.452582,  0.45258 ])
muB=np.array([[-0.142614,  0.17417 , -0.003947],
       [ 0.001574, -0.001082,  0.029477],
       [ 0.001404, -0.002555, -0.029395]])
alphaA = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])
alphaB = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])
"""
```

Once stuck on the damping function implementation, provide the following code snippet:
```txt
Here are damping functions with their partial derivatives for you to use. Also,
 get Z from the fragement atomic_numbers. qA with muA are considered M_A with respect to equation 4. Same for B. Once you have Z_A, qA[i] -= Z_A[i]. Don't forget to evaluate the Z_A*Z_B/r term
"""
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
    lam1, lam3, lam5, lam7, lam9 = (1.0, 1.0, 1.0, 1.0, 1.0)
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
"""
```