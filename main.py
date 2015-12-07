#!/usr/bin/env python

"""
Executes and demonstrates the capabilities of the triad motif extraction algorithm.
"""

__author__  = "Ross Flieger-Allison"
__date__    = "06-12-2015"
__version__ = "1.0.7"

import time
import networkx as nx
from utils import networks, plots, time_complexity
from utils import triad_motif_profile as tmp


def main():
    """
    This is the user-interface for the project.

    Anything the user should ever have to type should go here.
    """

    # Protein network.
    # network = networks.load_protein_network()
    # plots.plot_triad_motif_significance_profile(network, title="Traid Motif Significance Profile: Protein Network")
    # plots.plot_motif_changes_over_randomization_steps(network, title="Motif Expression Over Randomization Steps: Protein Network")

    # Social networks.
    network1 = networks.load_s420_electronic_circuit_network()
    network2 = networks.load_s208_electronic_circuit_network()
    # plots.plot_triad_motif_significance_profile(network1, title="Traid Motif Significance Profile: Leader Social Network")
    # plots.plot_motif_changes_over_randomization_steps(network1, title="Motif Expression Over Randomization Steps: Prison Social Network")
    plots.plot_triad_motif_significance_profiles([network1, network2], ["s420", "s208"])

if __name__ == "__main__":
    main()
