# Extracting Triad Motif Significance Profiles from Complex Networks

__Institution__: AIT Budapest

__Instructor__: Dániel Ábel

__Project Description__: Write a program that can extract the triad motif significance profile of an input network. Compare the motif profile of the networks used in the practical session. Give a 5-8 minute presentation (beside the motif significance profiles) showing how the frequency of a given motif is changing when randomizing the network, and include snapshots of the networks before and after the randomization.

### Background:

A __triad motif__ is a 3-node, connected subgraph in which the configuration of the links is predefined.

For undirected networks, we have only 2 types of triadic motifs, while for directed networks we have 13.

The __motif discovery problem__ comprises two main steps: 

1. Calculating the number of occurrences of a sub-graph.
2. Evaluating the sub-graph significance. The recurrence is "significant" if it is detectable far more than expected. Roughly speaking, the expected number of appearances of a sub-graph can be determined by a Null-model, which is defined by an ensemble of random networks with some of the same properties as the original network.

