"""
This program extracts the triad motif significance profile
of an input network.
"""

import networkx as nx
import numpy as np
from itertools import combinations
from networks import build_random_network


def count_undirected_triad_motifs(network):
    """
    Given an undirected network, returns a dictionary mapping triad motifs
    to their number of occurences within the network.
    """
    
    # Initialize motif counts to 0.
    motif_counts = {1: 0, # Motif 1 (full triangle)
                    2: 0} # Motif 2 (3-node chain)
    
    # Compute the number of unique node triplets in the network.
    triplets = combinations(nx.nodes_iter(network), 3)
    
    # Iterate through the triplets.
    for (a, b, c) in triplets:
        
        # a <-> b <-> c <-> a
        if b in network[a] and c in network[b] and a in network[c]:
            motif_counts[1] += 1
            
        # a <-> b <-> c
        elif b in network[a] and c in network[b]:
            motif_counts[2] += 1
            
        # b <-> c <-> a
        elif c in network[b] and a in network[c]:
            motif_counts[2] += 1
            
        # c <-> a <-> b
        elif a in network[c] and b in network[a]:
            motif_counts[2] += 1
            
    return motif_counts


def count_directed_triad_motifs(network):
    """
    Given a directed network, returns a dictionary mapping triad motifs
    to their number of occurences within the network.
    """
            
    raise NotImplementedError


def compute_triad_motif_z_scores(network):
    """
    Given any input network, returns a dictionary mapping triad motifs to
    their network expression z-scores.
    """

    # Determine the size of the input network.
    n, m = nx.number_of_nodes(network), nx.number_of_edges(network)

    if nx.is_directed(network):

        # Count the number of occurences of each triad motif.
        motif_counts = count_directed_triad_motifs(network)

        # Count the number of occurences of each triad motif on several RANDOM graphs of identical size.
        random_runs = [count_directed_triad_motifs(build_random_network(n, m, directed=True)) for _ in range(10)]

    else:

        # Count the number of occurences of each triad motif.
        motif_counts = count_undirected_triad_motifs(network)

        # Count the triad motifs in 10 random networks of equal size.
        random_runs = [count_undirected_triad_motifs(build_random_network(n, m, directed=False)) for _ in range(10)]

    # Initialize an empty dictionary to populate with our z-score values.
    z_scores = {}

    # Iterate through each motif.
    for motif in motif_counts:

        # Isolate the values for that motif, across all random runs.
        motif_values = map(lambda x: x[motif], random_runs)

        # Compute the z-score for that motif by subtracting the expected value and dividing by 
        # the standard deviation.
        z_scores[motif] = (motif_counts[motif] - np.mean(motif_values)) / np.std(motif_values)

    return z_scores
        

def normalize_z_scores(motif_z_scores):
    """
    Normalizes network expression z-scores to make significance profile scores.
    """

    # Define the normalization factor as the sum of the squares of all the z-scores.
    normalization_factor = float(sum(filter(lambda x: x ** 2, motif_z_scores.values())))

    # Iterate through the motifs and normalize their scores.
    for motif in motif_z_scores.keys():
        motif_z_scores[motif] = motif_z_scores[motif] / normalization_factor

    return motif_z_scores


def extract_triad_motif_significance_profile(network):
    """
    Given any input network, returns a triad motif significance profile.
    """

    # Make sure the network labels are encoded as integer IDs (makes things easier).
    network = nx.convert_node_labels_to_integers(network)

    # Build a dictionary of motifs and their corresponding normalized motif expression z-scores.
    significance_profile = normalize_z_scores(compute_triad_motif_z_scores(network))

    return significance_profile

