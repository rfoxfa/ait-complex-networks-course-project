"""
This program extracts the triad motif significance profile
of an input network.
"""

import networkx as nx
import itertools
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

    if directed or nx.is_directed(network):

        # Initialize an array for storing our motif counts. Each index represents the following motif:
        # 0 => a <- b -> c
        # 1 => a -> b <- c
        # 2 => b -> a -> c
        # 3 => a -> b <-> c
        # 4 => c <-> a -> b
        # 5 => a <-> b <-> c
        # 6 => a <- b -> c <- a
        # 7 => a -> b -> c -> a
        # 8 => c <-> a -> b <- c
        # 9 => c <- a -> b <-> c
        # 10 => a <-> c -> b -> a
        # 11 => a <-> b <-> c -> a
        # 12 => a <-> b <-> c <-> a
        motif_counts = np.zeros(shape=(13,), dtype=np.int)
    
        # Store a set of visited node combinations so repeats don't occur.
        visited_triplets = []

        # Iterate through the valid edges.
        for a, b in sorted(nx.edges_iter(network)):
                                            
            # Take all unique c nodes that form valid a, b, c triplets.
            # By ensuring that all neighbors have node ID greater than b we prevent any unique
            # node triplets from being repeated.
            c_neighbors = set([neighbor for neighbor in
                               itertools.chain(nx.all_neighbors(network, a), nx.all_neighbors(network, b))
                               if neighbor != a and neighbor != b])

            # Iterate through the valid a, b, c triplets.
            for c in c_neighbors:

                # Make sure unique node triplets aren't repeated.
                sorted_abc = sorted((a, b, c))
                if sorted_abc in visited_triplets:
                    continue
                else:
                    visited_triplets.append(sorted_abc)
            
                # a <-------> b
                # a     c     b
                if network.has_edge(b, a):

                    # a <-------> b
                    # a --> c     b
                    if network.has_edge(a, c):

                        # a <-------> b
                        # a <-> c     b
                        if network.has_edge(c, a):

                            # a <-------> b
                            # a <-> c --> b
                            if network.has_edge(c, b):

                                # a <-------> b
                                # a <-> c <-> b
                                if network.has_edge(b, c):
                                    motif_counts[12] += 1

                                # a <-------> b
                                # a <-> c --> b
                                else:
                                    motif_counts[11] += 1

                            # a <-------> b
                            # a <-> c <-- b
                            elif network.has_edge(b, c):
                                motif_counts[11] += 1

                            # a <-------> b
                            # a <-> c     b
                            else:
                                motif_counts[5] += 1

                        # a <-------> b
                        # a --> c     b 
                        else:

                            # a <-------> b
                            # a --> c --> b
                            if network.has_edge(c, b):

                                # a <-------> b
                                # a --> c <-> b
                                if network.has_edge(b, c):
                                    motif_counts[11] += 1

                                # a <-------> b
                                # a --> c --> b
                                else:
                                    motif_counts[10] += 1

                            # a <-------> b
                            # a --> c <-- b
                            elif network.has_edge(b, c):
                                motif_counts[8] += 1

                            # a <-------> b
                            # a --> c     b
                            else:
                                motif_counts[4] += 1

                    # a <-------> b
                    # a <-- c     b
                    elif network.has_edge(c, a):

                        # a <-------> b
                        # a <-- c --> b
                        if network.has_edge(c, b):

                            # a <-------> b
                            # a <-- c <-> b
                            if network.has_edge(b, c):
                                motif_counts[11] += 1

                            # a <-------> b
                            # a <-- c --> b
                            else:
                                motif_counts[9] += 1

                        # a <-------> b
                        # a <-- c <-- b
                        elif network.has_edge(b, c):
                            motif_counts[10] += 1

                        # a <-------> b
                        # a <-- c     b
                        else:
                            motif_counts[3] += 1

                    # a <-------> b
                    # a     c <-- b
                    elif network.has_edge(b, c):

                        # a <-------> b
                        # a     c <-> b
                        if network.has_edge(c, b):
                            motif_counts[5] += 1

                        # a <-------> b
                        # a     c <-- b
                        else:
                            motif_counts[4] += 1

                    # a <-------> b
                    # a     c --> b
                    else:
                        motif_counts[3] += 1

                # a --------> b
                # a     c     b
                else:

                    # a --------> b
                    # a --> c     b
                    if network.has_edge(a, c):

                        # a --------> b
                        # a <-> c     b
                        if network.has_edge(c, a):

                            # a --------> b
                            # a <-> c --> b
                            if network.has_edge(c, b):

                                # a --------> b
                                # a <-> c <-> b
                                if network.has_edge(b, c):
                                    motif_counts[11] += 1

                                # a --------> b
                                # a <-> c --> b
                                else:
                                    motif_counts[8] += 1

                            # a --------> b
                            # a <-> c <-- b
                            elif network.has_edge(b, c):
                                motif_counts[10] += 1

                            # a --------> b
                            # a <-> c     b
                            else:
                                motif_counts[4] += 1

                        # a --------> b
                        # a --> c     b 
                        else:

                            # a --------> b
                            # a --> c --> b 
                            if network.has_edge(c, b):

                                # a --------> b
                                # a --> c <-> b 
                                if network.has_edge(b, c):
                                    motif_counts[9] += 1

                                # a --------> b
                                # a --> c --> b 
                                else:
                                    motif_counts[6] += 1

                            # a --------> b
                            # a --> c <-- b 
                            elif network.has_edge(b, c):
                                motif_counts[6] += 1

                            # a --------> b
                            # a --> c     b 
                            else:
                                motif_counts[0] += 1

                    # a --------> b
                    # a <-- c     b
                    elif network.has_edge(c, a):

                        # a --------> b
                        # a <-- c --> b
                        if network.has_edge(c, b):

                            # a --------> b
                            # a <-- c <-> b
                            if network.has_edge(b, c):
                                motif_counts[10] += 1

                            # a --------> b
                            # a <-- c --> b
                            else:
                                motif_counts[6] += 1

                        # a --------> b
                        # a <-- c <-- b
                        elif network.has_edge(b, c):
                            motif_counts[7] += 1

                        # a --------> b
                        # a <-- c     b
                        else:
                            motif_counts[2] += 1

                    # a --------> b
                    # a     c <-- b
                    elif network.has_edge(b, c):

                        # a --------> b
                        # a     c <-> b
                        if network.has_edge(c, b):
                            motif_counts[3] += 1

                        # a --------> b
                        # a     c <-- b
                        else:
                            motif_counts[2] += 1

                    # a --------> b
                    # a     c --> b
                    else:
                        motif_counts[1] += 1

    else:

        # Initialize an array for storing our motif counts. Each index represents the following motif:
        # 0 => a - b - c - a
        # 1 => a - b - c
        motif_counts = np.zeros(shape=(2,), dtype=np.int)
    
        # Iterate through the edges in the network.
        for a, b in sorted(nx.edges_iter(network)):
            
            # Find all node a neighbors such that their node ID is greater than node b (which, by the sorted 
            # nature of nx.edges_iter, will always be greater than node a).
            a_neighbors = set([neighbor for neighbor in network[a] if neighbor > b])

            # Find all node b neighbors such that their node ID is greater than node b. The prevents repeated
            # consideration of node triplets.
            b_neighbors = set([neighbor for neighbor in network[b] if neighbor > b])
            
            # The number of triangle motifs is the number of common neighbors of a and b.
            motif_counts[0] += len(a_neighbors.intersection(b_neighbors))
            
            # The number of chain motifs is the number of unshared neighbors or a and b.
            motif_counts[1] += len(a_neighbors.symmetric_difference(b_neighbors))
    
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

    # Count the number of occurences of each triad motif.
    original_motif_counts = count_triad_motifs(network, directed=directed)

    # Initialize an array for storing motif counts in randomized instances.
    rand_motif_counts = []

    # Iterate through random instances.
    for i in range(num_rand_instances):

        # Update the user.
        print("Randomizing network...")

        # Randomize the network.
        rand_network = randomize(network, num_rewirings=num_rewirings)

        # Update the user.
        print("Counting motifs...")

        # Store the number of occurences of each motif in the randomized instance.
        rand_motif_counts.append(count_triad_motifs(rand_network, directed=directed))

        # Tell the user that we've finished one random instance.
        print("Counted triad motifs for randomized network instance #{}.".format(i + 1))

    # Update the user.
    print("Computing z-scores...")

    # Stack the counts as an array.
    rand_motif_counts = np.vstack(rand_motif_counts)

    # Divide the random motif counts by the number of instances to make them into average counts.
    avg_rand_motif_counts = np.mean(rand_motif_counts, axis=0)

    # Compute the random motif standard deviation.
    rand_motif_std_dev = np.std(rand_motif_counts, axis=0)

    # Compute the z-scores (ignoring division by 0 warnings and place the resulting NaNs/infs with 0s).
    with np.errstate(divide='ignore', invalid='ignore'):
        motif_z_scores = (original_motif_counts - avg_rand_motif_counts) / rand_motif_std_dev
        motif_z_scores[motif_z_scores == np.inf] = 0
        motif_z_scores = np.nan_to_num(motif_z_scores)

    # Update the user.
    print("Done.\n")

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

    # Notify the user that we've started.
    print("Extracting significance profile...")

    # Build an array of normalized motif expression z-scores (indices indicate unique motifs).
    significance_profile = compute_normalized_triad_motif_z_scores(network, 
                                                                   num_rand_instances=num_rand_instances, 
                                                                   num_rewirings=num_rewirings)

    return significance_profile

