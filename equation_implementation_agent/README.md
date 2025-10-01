# Prompts
```txt
# Task
implement the electrostatics equation from @schriber_2021_184110.pdf. Implement the equation using numpy and qcelemental. Use a test driven development approach to implement the correct code to match an electrostatics energy for the given system below of -6.824148 kcal/mol. When implementing each term, explicitly write a comment that matches the latex format of the equation in the paper to help with organization/debugging. Note all usage of qA should be qA_electronic = qA - ZA and qB_electronic = qB - ZB.

# Output
electrostatics.py for implementing the code and a test_electrostatics.py file for running a pytest. Execute pytest to test each term to know which terms are correct and which ones need corrected.

# Information for pytest
"""

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

qA=array([-0.90282657,  0.45219405,  0.45063254])
muA=array([[-0.12099517, -0.19035041,  0.00497184],
       [-0.02343171,  0.01783004, -0.00038207],
       [ 0.026533  , -0.01378759,  0.0002752 ]])
qB=array([-0.90516141,  0.45258175,  0.45257968])
muB=array([[-0.14261352,  0.17416981, -0.00394689],
       [ 0.00157417, -0.00108209,  0.02947683],
       [ 0.00140431, -0.00255545, -0.029395  ]])



# Undamped
E_ZA_ZB = 12056.9380 kcal/mol # term1
E_ZA_MB = -12206.0827 kcal/mol # term2
E_ZB_MA = -11880.5922 kcal/mol # term3
MTP_MTP = 12022.9127 kcal/mol # term4
ref_elst=-6.824148 kcal/mol # total electrostatics energy without damping functions (f_1=1 and f_2=1)


"""
```

# Damping function implementation
Now that we have a working electrostatics implementation, we will be able to direct the model to use the following damping functions to get the electrostatics energy closer to the one implemented in the paper.

```txt
Here are damping functions with their partial derivatives for you to use. Copy the electrostatics equation and add another function for damped electrostatics. Use the alpha_A and alpha_B values for the damping parameters alpha_i and alpha_j in the pytest for water dimer. The numerical agreement for this test can be 1e-2.
"""
# Damped total
elst_charges_dipoles_damped = -10.532698 kcal/mol

alphaA = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])
alphaB = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])

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
"""
```

# Optimize code
Okay that implementation works but... you would never use that in production. Can we use LLMs to make our code faster? Can it quantify that improvement?
```txt
Now using the damped electrostatics test, optimize the numpy code to be vectorized and more efficient but still get the correct results. Write a pytest for this optimized version and write
a script that analyzes the function calls between the two functions over an average of 10 runs with results on a bar plot.
```

# Results
This example shows how you could go from a simple equation from a paper (with reference values) to code through a test driven development approach. Now that you have reproduced results and optimized the code a little bit (much more would be needed in reality). Now try repeating this exercise with a different model through opencode, creating an adapted AGENTS.md file and other simplifications/improvements to the prompts above! With enough pointed direction, the models should give similar outputs; however, this is definitely not deterministic and results vary. Important metrics to think about are cost/speed/adherence to instructions/code performance metrics/accuracy/etc...
