"""
Plotting functions.
"""

import matplotlib.pyplot as plt
from triad_motif_profile import extract_triad_motif_significance_profile


def plot_triad_motif_significance_profile(motif_profile, directed=False):
    """
    Given a triad motif significance profile, plot the data in a readable
    fashion.
    """

    plt.xlabel("Triad Motifs")
    plt.ylabel("Normalized Z-score")
    plt.title("Triad Motif Significance Profile")
    plt.show()

    raise NotImplementedError
