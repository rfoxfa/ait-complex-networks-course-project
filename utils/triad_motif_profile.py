"""
This program extracts the triad motif significance profile
of an input network.
"""

import networkx as nx
import numpy as np
from itertools import combinations
from networks import randomize


def count_triad_motifs(network, directed=None):
    """
    Counts the occurences of triad motifs in a network.

    Arguments:
        network => The input network.
        directed => Whether or not the network is directed.

    Returns:
        A fixed-size array with indices representing unique triad motifs and the values 
        representing their number of occurences within the network.
    """

    if (directed or nx.is_directed(network)):

        # Initialize an array for storing our motif counts.
        motif_counts = np.zeros(shape=(13,), dtype=np.int)

        raise NotImplementedError

    else:

        # Initialize an array for storing our motif counts. Each index represents the following motif:
        # 0 => triangle
        # 1 => chain
        motif_counts = np.zeros(shape=(2,), dtype=np.int)
    
        # Iterate through the edges in the network.
        for source, target in nx.edges_iter(network):
            
            # Find all source neighbors such that their node ID is greater than the target (which, by the sorted 
            # nature of nx.edges_iter, will always be greater than source).
            source_neighbors = [neighbor for neighbor in nx.all_neighbors(network, source) if neighbor > target]

            # Find all target neighbors such that their node ID is greater than the target. The prevents repeated
            # consideration of node triplets.
            target_neighbors = [neighbor for neighbor in nx.all_neighbors(network, target) if neighbor > target]

            # Determine the number of unshared neighbors between source and target.
            # Given a set and a list, symmetric_difference will return an XOR of two lists.
            num_unshared_neighbors = len(set(source_neighbors).symmetric_difference(target_neighbors))
            
            # The number of triangle motifs is the number of common neighbors of source and target.
            motif_counts[0] += len(source_neighbors) - num_unshared_neighbors
            
            # The number of chain motifs is the number of unshared neighbors or source and target.
            motif_counts[1] += num_unshared_neighbors
    
    return motif_counts


def compute_normalized_triad_motif_z_scores(network, num_rand_instances=10, num_rewirings=None):
    """
    Computes the normalized triad motif z-score for each connected non-isomorphic triadic subgraph
    in the input network.

    Arguments:
        network => The input network (can be directed or undirected).
        num_rand_instances => The number of randomly-rewired network instances used when computing
        z-score values.
        num_rewirings => The number of edge rewirings performed when randomizing the network.

    Returns:
        A fixed-size numpy array where each index corresponds to a predefined unique triad motif
        and where the value at each index represents the normalized z-score for the average
        over- or underexpression of a triad motif in the network.
    """

    # Determine if the network is directed or not (store to avoid recalculation).
    directed = nx.is_directed(network)

    # Store the number of connected non-isomorphic motifs to avoid hard-coding these values.
    num_motifs = 13 if directed else 2

    # Count the number of occurences of each triad motif.
    original_motif_counts = count_triad_motifs(network, directed=directed)

    # Initialize an array for storing motif counts in randomized instances.
    rand_motif_counts = np.zeros(shape=(num_motifs,), dtype=np.int)

    # Initialize an array for storing squared motif counts in randomized instances (for standard deviation).
    rand_motif_squared_counts = np.zeros(shape=(num_motifs,), dtype=np.int)

    # Iterate through random instances.
    for _ in range(num_rand_instances):

        # Randomize the network with either the inputted number of rewirings or m rewirings.
        rand_network = randomize(network, num_rewirings=num_rewirings, directed=directed)

        # Store the number of occurences of each motif in the randomized instance.
        counts = count_triad_motifs(rand_network, directed=directed)

        # Add the counts to the total random motif counts.
        rand_motif_counts += counts

        # Add the squared counts to the total squared counts.
        rand_motif_squared_counts += counts * counts

    # Divide the random motif counts by the number of instances to make them into average counts.
    rand_motif_counts = rand_motif_counts / float(num_rand_instances)

    # Divide the random motif squared counts by the number of instances to make them into averages.
    rand_motif_squared_counts = rand_motif_squared_counts / float(num_rand_instances)

    # Compute the random motif standard deviation.
    rand_motif_std_dev = np.sqrt(rand_motif_squared_counts - (rand_motif_counts * rand_motif_counts))

    # Compute the z-scores.
    motif_z_scores = (original_motif_counts - rand_motif_counts) / rand_motif_std_dev

    return motif_z_scores


def extract_triad_motif_significance_profile(network, num_rand_instances=10, num_rewirings=None):
    """
    Computes the triad motif significance profile of the input network.

    Arguments:
        network => The input network (can be directed or undirected).
        num_rand_instances => The number of randomly-rewired network instances used when computing
        z-score values.
        num_rewirings => The number of edge rewirings performed when randomizing the network.

    Returns:
        A fixed-size numpy array where each index corresponds to a predefined unique triad motif
        and where the value at each index represents the normalized z-score for the average
        over- or underexpression of a triad motif in the network.
    """

    # Make sure the network labels are encoded as integer IDs (makes everything easier).
    network = nx.convert_node_labels_to_integers(network)

    # Build an array of normalized motif expression z-scores (indices indicate unique motifs).
    significance_profile = compute_normalized_triad_motif_z_scores(network, 
                                                                   num_rand_instances=num_rand_instances, 
                                                                   num_rewirings=num_rewirings)

    return significance_profile

