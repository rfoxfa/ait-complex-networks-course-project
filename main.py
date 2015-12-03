#!/usr/bin/env python

"""
Executes and demonstrates the capabilities of the triad motif extraction algorithm.
"""

__author__  = "Ross Flieger-Allison"
__date__    = "02-12-2015"
__version__ = "1.0.1"

from utils import networks, plots, time_complexity, triad_motif_profile


def main():
    """
    This is the user-interface for the project.

    Anything the user should ever have to type should go here.
    """

    # Do stuff.
    print(time_complexity.compute_undirected_count_execution_time())
    # print(triad_motif_profile.extract_triad_motif_significance_profile(networks.build_random_network()))


if __name__ == "__main__":
    main()
