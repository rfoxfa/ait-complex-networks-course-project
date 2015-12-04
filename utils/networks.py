"""
Builds and returns various networkx networks.
"""

import networkx as nx
import random


def randomize(network, num_rewirings=None, directed=False):
    """
    Randomizes the network such that the degree sequence is preserved. 

    Arguments:
        network => The input network.
        num_rewirings => The number of rewirings done to make the 
    """

    # Set a counter for the number of rewirings successfully completed.
    rewirings_completed = 0

    # Set the required number of rewirings (either inputted value or the total number of edges).
    num_rewirings = num_rewirings or nx.number_of_edges(network)

    # Don't terminate until all the necessary rewirings have been performed.
    while rewirings_completed < num_rewirings:

        # Randomly selected a link from the network.
        link1 = (source1, target1) = random.choice(nx.edges(network))

        # Randomly selected a DIFFERENT link from the network.
        link2 = (source2, target2) = random.choice([link for link in nx.edges(network) if link != link1])

        # If the two links do not share any nodes...
        if not any(node in link1 for node in link2):

            # If the graph is directed, there is only one option.
            # If the graph is undirected, there are two options, each with a 50-50 chance.
            if not directed and random.random() < 0.5:

                # Rewire links A-B and C-D to A-C and B-D.
                new_link1 = (source1, source2)
                new_link2 = (target1, target2)

            else:

                # Rewire links A-B and C-D to A-D and C-B.
                new_link1 = (source1, target2)
                new_link2 = (source2, target1)

            # If the new links aren't in the network already, replace the old links with the new links.
            if not network.has_edge(*new_link1) and not network.has_edge(*new_link2):

                # Remove the old links.
                network.remove_edges_from([link1, link2])

                # Add the new links.
                network.add_edges_from([new_link1, new_link2])

                # Incrememnt our rewiring counter for each individual rewired edge.
                rewirings_completed += 2

    return network



def build_random_network(n=100, m=300, directed=False):
    """
    Builds a random network with n nodes (default n=100) and m edges
    (default m=300).
    """

    return nx.gnm_random_graph(n, m, directed=directed)


def build_erdos_renyi_network(n=100, p=0.5, directed=False):
    """
    Builds a random Erdos-Renyi network with n nodes and connection
    probability p (default p=0.5).
    """

    return nx.erdos_renyi_graph(n, p, directed=directed)


def build_watts_strogatz_network(n=100, k=5, p=0.5, directed=False):
    """
    Builds a Watts-Strogatz network with n nodes (default n=100) where each
    node is joined with its k nearest neighbors in a ring topology (default k=5)
    and the probability of rewiring each edge is p (default p=0.5).
    """

    return nx.watts_strogatz_graph(n, k, p, directed=directed)


def build_barabasi_albert_network(n=100, m=300, directed=False):
    """
    Builds a Barabasi-Albert network with n nodes (default n=100) and 
    m edges (default m=300).
    """

    return nx.barabase_albert_graph(n, m, directed=directed)

