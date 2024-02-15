"""
This package computes free splittings of graph braid groups.

The package implements computational methods developed in the mathematical paper 'Graph of groups decompositions of graph braid groups' by Daniel Berlyne (the author of this package).
Given a graph braid group, the splitter will search for a free splitting of the group, returning the splitting if one is found. 
The factors of this splitting are expressed as either: 
- free groups; 
- simpler graph braid groups; 
- graphs of groups where the vertex groups and edge groups are simpler graph braid groups. 
Next, the splitter expresses each of these free factors as a direct product if it is directly decomposable. 
The splitter will then attempt to further split these direct factors, alternating between free products and direct products until a splitting with freely indecomposable factors is acquired. 
The splitter will then evaluate any known graph braid groups that appear in the splitting.

Package contents: 
- `splitter` contains the main code for detecting free splittings of graph braid groups.
- `graph` contains methods for performing graph computations.
- `graph_braid_group` contains methods for performing graph braid group computations.
- `graph_of_groups` contains methods for performing graph of groups computations.
"""