import numpy as np
import qcelemental as qcel
from elst import calculate_electrostatics_energy

def test_electrostatics_energy():
    elst_charges_dipoles = -10.532698  # kcal/mol
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

    qA = np.array([-0.902827, 0.452194, 0.450633])
    muA = np.array([[-0.120995, -0.19035, 0.004972],
                    [-0.023432, 0.01783, -0.000382],
                    [0.026533, -0.013788, 0.000275]])
    qB = np.array([-0.905161, 0.452582, 0.45258])
    muB = np.array([[-0.142614, 0.17417, -0.003947],
                    [0.001574, -0.001082, 0.029477],
                    [0.001404, -0.002555, -0.029395]])

    alphaA = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])
    alphaB = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])

    calculated_energy = calculate_electrostatics_energy(qcel_mol, qA, muA, qB, muB, alphaA, alphaB)

    np.testing.assert_allclose(calculated_energy, elst_charges_dipoles, rtol=1e-6)
