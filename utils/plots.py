"""
Plotting functions.
"""

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import matplotlib.ticker as ticker
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib._png import read_png
from triad_motif_profile import extract_triad_motif_significance_profile, count_triad_motifs
from networks import random_rewiring


def plot_triad_motif_significance_profile(network, num_runs=20, title=None):
    """
    Plots the triad motif significance profile for the input network.

    Arguments:
        network => The input network (directed or undirected).
        num_runs => The number of extraction iterations performed on the network.
        title => A custom title for the plot (let's you specify which plot you're looking at).

    Returns:
        Nothing but produces a plot showing the triad motif significance profile
        of several extraction runs on the network such that the results are 
        overlayed on each other.
    """

    # Determine if the network is directed or not.
    directed = nx.is_directed(network)
    
    # Create a figure with hardcoded dimensions and one subplot.
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(1, 1, 1)

    # Set some parameters specific to the directionality of the network.
    if directed:
        num_motifs = 13
        ax.margins(0.05, 0.05)
        motif_images = [OffsetImage(read_png("../project/static/directed_motif_1.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_2.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_3.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_4.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_5.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_6.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_7.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_8.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_9.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_10.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_11.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_12.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_13.png"), zoom=0.6)]
    else:
        num_motifs = 2
        ax.margins(0.5, 0.5)
        motif_images = [OffsetImage(read_png("../project/static/undirected_motif_1.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/undirected_motif_2.png"), zoom=0.6)]

    # Define x values (integer for each motif).
    X = np.arange(1, num_motifs + 1)
    plt.xticks(X)

    # Iterate through the extraction instances (synonymous with each color instance).
    for _ in range(num_runs):

        # Run the extraction code.
        Y = extract_triad_motif_significance_profile(network)

        # Plot the y-values with straight lines connecting them.
        ax.plot(X, Y)
        
    # Define the x-y coordinate offsets for the images on the plot.
    y_offset = -50
    x_offset = 0
        
    # Add each image to the plot.
    for image, x in zip(motif_images, X):
        ax.add_artist(AnnotationBbox(image, (x, ax.get_ylim()[0]), xybox=(x_offset, y_offset), xycoords='data',
                                     boxcoords='offset points', frameon=False))  


    # Title and label the plot.
    plt.title(title or "Triad Motif Significance Profile")
    plt.xlabel("Triad Motifs", labelpad=75)
    plt.ylabel("Normalized Z-Score")

    # Add the y = 0 line for readability.
    ax.axhline(y=0, color="black")

    # Set the y-axis tick marks.
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    
    # Add extra spacing at the bottom of the plot for the images.
    plt.subplots_adjust(bottom=0.2)  

    # Show the resulting plot.
    plt.show()


def plot_motif_changes_over_randomization_steps(network, colormap="rainbow", rewiring_limit=None, title=None):
    """
    Plots the convergence of motif counts over several randomization steps.

    Arguments:
        network => The input network.
        colormap => A matplotlib colormap for giving different motif lines different colors.
        rewiring_limit => An upper bound on the number of allowable edge rewirings (default set to 3 times
        the number of edges in the network).
        title => A custom plot title.

    Returns:
        Nothing, but shows a plot with motif counts converging over time.
    """

    # Determine if the network is directed or not.
    directed = nx.is_directed(network)
    
    # Create a figure with hardcoded dimensions and one subplot.
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(1, 1, 1)

    # Set some parameters specific to the directionality of the network.
    if directed:
        num_motifs = 13
    else:
        num_motifs = 2

    # Count the motif before randomization.
    starting_counts = count_triad_motifs(network)

    # Keep a history of motif counts over time.
    motif_count_history = starting_counts

    # Keep a history of average motif counts over time (these are our y-values).
    avg_motif_count_history = starting_counts

    # Keep track of how many iterations have been performed.
    num_rewirings_performed = 0

    # Set a hard limit on the number than can be performed.
    rewiring_limit = rewiring_limit or 3 * nx.number_of_edges(network)

    while num_rewirings_performed < rewiring_limit:

        # Randomly rewire a pair of edges in the network.
        network = random_rewiring(network)

        # Increment the counter.
        num_rewirings_performed += 1

        # Count the motifs in the slightly altered network.
        new_motif_counts = count_triad_motifs(network)

        # Add those counts to the motif history.
        motif_count_history = np.vstack((motif_count_history, new_motif_counts))

        # Determine the average motif counts over the current history.
        avg_motif_counts = np.mean(motif_count_history, axis=0)   
   
        # Add that to the avg_motif_count_history stack of averages.
        avg_motif_count_history = np.vstack((avg_motif_count_history, avg_motif_counts))

    # Define x values.
    X = np.arange(num_rewirings_performed + 1)

    # Iterate through the extraction instances (synonymous with each color instance).
    for motif in range(num_motifs):

        # Parse all motif counts for a specific motif over time.
        Y = avg_motif_count_history[:, motif]

        # Plot the y-values with straight lines connecting them.
        ax.plot(X, Y, label="Motif {}".format(motif + 1), linewidth=2.0)

    # Add a legend for the lines in the plot.
    ax.legend()

    # Set the axis tick marks.
    ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(5))
        
    # Title and label the plot.
    plt.title(title or "Motif Expression During Random Rewiring")
    plt.xlabel("Number of Random Edge Rewirings")
    plt.ylabel("Average Motif Counts")

    # Show the resulting plot.
    plt.show()


def plot_triad_motif_significance_profiles(networks, labels, title=None):
    """
    Plots the triad motif significance profile for the input network.

    Arguments:
        networks => A list of input networks (directed or undirected).
        labels => A list of network labels for the legend.
        title => A custom title for the plot (let's you specify which plot you're looking at).

    Returns:
        Nothing but produces a plot showing the triad motif significance profile
        of a list of networks such that the results are overlayed on each other.
    """

    # Determine if the network is directed or not.
    directed = nx.is_directed(networks[0])
    
    # Create a figure with hardcoded dimensions and one subplot.
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(1, 1, 1)

    # Set some parameters specific to the directionality of the network.
    if directed:
        num_motifs = 13
        ax.margins(0.05, 0.05)
        motif_images = [OffsetImage(read_png("../project/static/directed_motif_1.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_2.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_3.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_4.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_5.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_6.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_7.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_8.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_9.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_10.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_11.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_12.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/directed_motif_13.png"), zoom=0.6)]
    else:
        num_motifs = 2
        ax.margins(0.5, 0.5)
        motif_images = [OffsetImage(read_png("../project/static/undirected_motif_1.png"), zoom=0.6),
                        OffsetImage(read_png("../project/static/undirected_motif_2.png"), zoom=0.6)]

    # Define x values (integer for each motif).
    X = np.arange(1, num_motifs + 1)
    plt.xticks(X)

    # Iterate through the extraction instances (synonymous with each color instance).
    for network, label in zip(networks, labels):

        # Run the extraction code.
        Y = extract_triad_motif_significance_profile(network)

        # Plot the y-values with straight lines connecting them.
        ax.plot(X, Y, label=label, linewidth=2.0)

    # Add a legend for the lines in the plot.
    ax.legend()
        
    # Define the x-y coordinate offsets for the images on the plot.
    y_offset = -50
    x_offset = 0
        
    # Add each image to the plot.
    for image, x in zip(motif_images, X):
        ax.add_artist(AnnotationBbox(image, (x, ax.get_ylim()[0]), xybox=(x_offset, y_offset), xycoords='data',
                                     boxcoords='offset points', frameon=False))  


    # Title and label the plot.
    plt.title(title or "Triad Motif Significance Profiles")
    plt.xlabel("Triad Motifs", labelpad=75)
    plt.ylabel("Normalized Z-Score")

    # Add the y = 0 line for readability.
    ax.axhline(y=0, color="black")

    # Set the y-axis tick marks.
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    
    # Add extra spacing at the bottom of the plot for the images.
    plt.subplots_adjust(bottom=0.2)  

    # Show the resulting plot.
    plt.show()

