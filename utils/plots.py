"""
Plotting functions.
"""

import matplotlib.pyplot as plt
from triad_motif_profile import extract_triad_motif_significance_profile


def plot_triad_motif_significance_profile(network, title=None, num_runs=10):
    """
    Plots the triad motif significance profile for the input network.

    Arguments:
        network => The input network (directed or undirected).
        num_runs => The number of extraction iterations performed on the network.

    Returns:
        Nothing but produces a plot showing the triad motif significance profile
        of several extraction runs on the network such that the results are 
        overlayed on each other.
    """

    # Label the plot.
    plt.xlabel("Triad Motifs")
    plt.ylabel("Normalized Z-Score")
    plt.title(title or "Triad Motif Significance Profile")

    # Perform the initial run.
    Y = extract_triad_motif_significance_profile(network)

    # Set up the plot based on the first run.
    X = range(len(Y))
    plt.xticks(X)
    plt.plot(X, Y)

    # Perform the other n - 1 runs and plot those too.
    for _ in range(num_runs - 1):
        Y = extract_triad_motif_significance_profile(network)
        plt.plot(X, Y)

    # Show the resulting plot.
    plt.show()

