"""
Functions for computing execution time complexity.
"""

import time
import numpy as np
from networks import build_erdos_renyi_network
from triad_motif_profile import count_undirected_triad_motifs


def compute_undirected_count_execution_time(n=100, p=0.5, num_iterations=100):
    """
    Computes the average execution time (in seconds) of the undirected triad motif counting 
    function on a random undirected network of size n and connection probability p.

    Result is returned a printable text report.
    """

    # Define a list of execution times.
    execution_times = []

    # Run through simulations.
    for _ in range(num_iterations):

        # Build a standard random Erdos-Renyi network with the inputted size and connection
        # probability (NOTE: this should not be included in the execution time).
        network = build_erdos_renyi_network(n, p, directed=False)

        # Start the timer.
        start_time = time.time()

        # Execute the function.
        count_undirected_triad_motifs(network)

        # Append the final time.
        execution_times.append(time.time() - start_time)

    # Build a printable report to show the results.
    report = "Network size: {0}\nConnection probability: {1}\nNumber of iterations: {2}\nAverage execution time: {3:.3f} sec.".format(n, p, num_iterations, np.mean(execution_times))

    return report

