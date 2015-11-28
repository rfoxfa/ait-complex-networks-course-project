"""
Builds and returns various networkx networks.
"""

import networkx as nx


def build_random_network(n=100, m=300):
    """
    Builds a random network with n nodes (default n=100) and m edges
    (default m=300).
    """

    return nx.gnm_random_graph(n, m)


def build_erdos_renyi_network(n=100, p=0.5):
    """
    Builds a random Erdos-Renyi network with n nodes and connection
    probability p (default p=0.5).
    """

    return nx.erdos_renyi_graph(n, p)


def build_watts_strogatz_network(n=100, k=5, p=0.5):
    """
    Builds a Watts-Strogatz network with n nodes (default n=100) where each
    node is joined with its k nearest neighbors in a ring topology (default k=5)
    and the probability of rewiring each edge is p (default p=0.5).
    """

    return nx.watts_strogatz_graph(n, k, p)


def build_barabasi_albert_network(n=100, m=300):
    """
    Builds a Barabasi-Albert network with n nodes (default n=100) and 
    m edges (default m=300).
    """

    return nx.barabase_albert_graph(n, m)

