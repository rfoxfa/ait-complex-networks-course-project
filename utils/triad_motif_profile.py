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
    Given an undirected network, returns a fixed-size array with indices representing
    unique triad motifs and the values representing their their number of occurences 
    within the network.
    """
    
    # Initialize motif counts to 0.
    # NOTE: for undirected networks, there are only 2 connected non-isomorphic triads.
    motif_counts = np.zeros(shape=(2,), dtype=np.int)
    
    # Compute the number of unique node triplets in the network.
    triplets = combinations(nx.nodes_iter(network), 3)
    
    # Iterate through the triplets.
    for (a, b, c) in triplets:
        
        # a <-> b <-> c <-> a
        if b in network[a] and c in network[b] and a in network[c]:
            motif_counts[0] += 1
            
        # a <-> b <-> c
        elif b in network[a] and c in network[b]:
            motif_counts[1] += 1
            
        # b <-> c <-> a
        elif c in network[b] and a in network[c]:
            motif_counts[1] += 1
            
        # c <-> a <-> b
        elif a in network[c] and b in network[a]:
            motif_counts[1] += 1
            
    return motif_counts


def count_directed_triad_motifs(network):
    """
    Given a directed network, returns a fixed-size array with indices representing
    unique triad motifs and the values representing their their number of occurences 
    within the network.
    """
            
    raise NotImplementedError


def compute_triad_motif_z_scores(network):
    """
    Given any input network, returns a fixed-size array with indices representing
    unique triad motifs and the values representing their network expression z-scores.
    """

    # Determine the size of the input network.
    n, m = nx.number_of_nodes(network), nx.number_of_edges(network)

    if nx.is_directed(network):

        # Count the number of occurences of each triad motif.
        motif_counts = count_directed_triad_motifs(network)

        # Count the number of occurences of each triad motif on several RANDOM graphs of identical size.
        random_runs = np.array([count_directed_triad_motifs(build_random_network(n, m, directed=True)) for _ in range(10)])

    else:

        # Count the number of occurences of each triad motif.
        motif_counts = count_undirected_triad_motifs(network)

        # Count the triad motifs in 10 random networks of equal size.
        random_runs = np.array([count_undirected_triad_motifs(build_random_network(n, m, directed=False)) for _ in range(10)])

    # Initialize an empty array to populate with our z-score values.
    motif_z_scores = np.ndarray(shape=(len(motif_counts),), dtype=np.float)

    # Iterate through each motif.
    for i in range(len(motif_counts)):

        # Isolate the values for that motif, across all random runs.
        motif_values = random_runs[i, :]

        # Compute the z-score for that motif by subtracting the expected value and dividing by 
        # the standard deviation.
        motif_z_scores[i] = (motif_counts[i] - np.mean(motif_values)) / np.std(motif_values)

    return motif_z_scores
        

def normalize_z_scores(motif_z_scores):
    """
    Normalizes network expression z-scores to make significance profile scores.
    """

    # Define the normalization factor as the square root of the sum of the squares of 
    # all the z-scores.
    normalization_factor = np.sqrt(np.ndarray.sum(motif_z_scores ** 2, dtype=np.float))

    # Normalize the motif scores.
    motif_z_scores /= normalization_factor

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

