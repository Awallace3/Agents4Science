import numpy as np
from numpy import array
import qcelemental as qcel
import pytest
from electrostatics import calculate_electrostatics, calculate_damped_electrostatics, calculate_damped_electrostatics_vectorized

@pytest.fixture
def system_data():
    """Provides the system data for the electrostatics calculation."""
    qcel_mol = qcel.models.Molecule.from_data(
        '''
        0 1
        --
        0 1
        O -1.326958220000 -0.105938540000 0.018788150000
        H -1.931665230000 1.600174310000 -0.021710520000
        H 0.486644270000 0.079598100000 0.009862480000
        --
        0 1
        O 3.907523240000 0.052757410000 0.001850160000
        H 4.619234940000 -0.775660840000 1.449615410000
        H 4.611000850000 -0.847154680000 -1.406756420000
        units bohr
        no_com
        no_reorient
        '''
    )

    qA = array([-0.90282657, 0.45219405, 0.45063254])
    muA = array([[-0.12099517, -0.19035041, 0.00497184],
                 [-0.02343171, 0.01783004, -0.00038207],
                 [0.026533, -0.01378759, 0.0002752]])
    qB = array([-0.90516141, 0.45258175, 0.45257968])
    muB = array([[-0.14261352, 0.17416981, -0.00394689],
                 [0.00157417, -0.00108209, 0.02947683],
                 [0.00140431, -0.00255545, -0.029395]])

    ref_E_ZA_ZB = 12056.9380  # kcal/mol
    ref_E_ZA_MB = -12206.0827  # kcal/mol
    ref_E_ZB_MA = -11880.5922  # kcal/mol
    ref_MTP_MTP = 12022.9127  # kcal/mol
    ref_elst = -6.824148  # kcal/mol

    return {
        "qcel_mol": qcel_mol,
        "qA": qA, "muA": muA, "qB": qB, "muB": muB,
        "ref_E_ZA_ZB": ref_E_ZA_ZB,
        "ref_E_ZA_MB": ref_E_ZA_MB,
        "ref_E_ZB_MA": ref_E_ZB_MA,
        "ref_MTP_MTP": ref_MTP_MTP,
        "ref_elst": ref_elst
    }

def test_electrostatics_terms(system_data):
    """Tests the individual terms of the electrostatic calculation."""
    total_elst, E1, E2, E3, E4 = calculate_electrostatics(
        system_data["qcel_mol"], system_data["qA"], system_data["muA"],
        system_data["qB"], system_data["muB"]
    )

    assert np.isclose(E1, system_data["ref_E_ZA_ZB"], atol=1e-4)
    assert np.isclose(E2, system_data["ref_E_ZA_MB"], atol=1e-4)
    assert np.isclose(E3, system_data["ref_E_ZB_MA"], atol=1e-4)
    assert np.isclose(E4, system_data["ref_MTP_MTP"], atol=1e-4)
    assert np.isclose(total_elst, system_data["ref_elst"], atol=1e-4)

def test_damped_electrostatics(system_data):
    """Tests the damped electrostatic calculation."""
    alphaA = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])
    alphaB = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])
    ref_damped_elst = -10.532698

    total_damped_elst, _, _, _, _ = calculate_damped_electrostatics(
        system_data["qcel_mol"], system_data["qA"], system_data["muA"],
        system_data["qB"], system_data["muB"], alphaA, alphaB
    )

    assert np.isclose(total_damped_elst, ref_damped_elst, atol=1e-2)

def test_vectorized_vs_original(system_data):
    """Tests that the vectorized and original damped functions give the same result."""
    alphaA = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])
    alphaB = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])

    original_results = calculate_damped_electrostatics(
        system_data["qcel_mol"], system_data["qA"], system_data["muA"],
        system_data["qB"], system_data["muB"], alphaA, alphaB
    )

    vectorized_results = calculate_damped_electrostatics_vectorized(
        system_data["qcel_mol"], system_data["qA"], system_data["muA"],
        system_data["qB"], system_data["muB"], alphaA, alphaB
    )

    for i in range(len(original_results)):
        assert np.isclose(original_results[i], vectorized_results[i], atol=1e-6)