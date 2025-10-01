import numpy as np
import qcelemental as qcel
import timeit
import matplotlib.pyplot as plt
from numpy import array
from electrostatics import calculate_damped_electrostatics, calculate_damped_electrostatics_vectorized

def run_benchmark():
    """Runs the benchmark and generates a bar plot of the results."""
    # Setup data
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
    alphaA = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])
    alphaB = np.array([2.05109221104216, 1.65393856475232, 1.65393856475232])

    # Number of runs
    n_runs = 10

    # Time the original function
    original_times = timeit.repeat(
        lambda: calculate_damped_electrostatics(qcel_mol, qA, muA, qB, muB, alphaA, alphaB),
        repeat=n_runs, number=1
    )
    avg_original_time = np.mean(original_times)

    # Time the vectorized function
    with np.errstate(divide='ignore'):
        vectorized_times = timeit.repeat(
            lambda: calculate_damped_electrostatics_vectorized(qcel_mol, qA, muA, qB, muB, alphaA, alphaB),
            repeat=n_runs, number=1
        )
    avg_vectorized_time = np.mean(vectorized_times)

    # Create bar plot
    functions = ['Original (Looped)', 'Optimized (Vectorized)']
    avg_times = [avg_original_time, avg_vectorized_time]

    plt.figure(figsize=(8, 6))
    bars = plt.bar(functions, avg_times, color=['blue', 'green'])
    plt.ylabel('Average Execution Time (seconds)')
    plt.title('Performance Comparison of Electrostatics Functions')
    
    # Add text labels
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, yval, f'{yval:.6f}s', va='bottom') # format to 6 decimal places


    plt.savefig('benchmark.png')
    print("Benchmark plot saved as benchmark.png")
    print(f"Average time for original function: {avg_original_time:.6f}s")
    print(f"Average time for vectorized function: {avg_vectorized_time:.6f}s")

if __name__ == "__main__":
    run_benchmark()
