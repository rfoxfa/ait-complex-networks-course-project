"""
Functions for computing execution time complexity.
"""

import time
import numpy as np
import networkx as nx
from triad_motif_profile import count_triad_motifs


def compute_count_execution_time(n=300, p=0.5, directed=False, num_iterations=10):
    """
    Computes the average execution time (in seconds) of the triad motif counting 
    function on a random network of size n and connection probability p.

    Arguments:
        n => The number of nodes in the network.
        p => The probability that any two nodes are connected.
        directed => Whether or not the network is directed.
        num_iterations => The number of function executions used to determine the
        average execution time.

    Returns:
        A printable report showing average execution time and other network specs.
    """

    # Define a list of execution times.
    execution_times = []

    # Run through simulations.
    for _ in range(num_iterations):

        # Build a standard random Erdos-Renyi network with the inputted size and connection
        # probability (NOTE: this should not be included in the execution time).
        network = build_erdos_renyi_network(n, p, directed=directed)

        # Start the timer.
        start_time = time.time()

        # Execute the function.
        count_triad_motifs(network, directed=directed)

        # Append the final time.
        execution_times.append(time.time() - start_time)

    # Build a printable report to show the results.
    report = "Network size: {0}\nConnection probability: {1}\nNumber of iterations: {2}\nDirected: {3}\nAverage execution time: {4:.3f} sec.".format(n, p, num_iterations, directed, np.mean(execution_times))

    return report

