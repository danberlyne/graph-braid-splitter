# About this tool
This tool implements computational methods developed in the mathematical paper ['Graph of groups decompositions of graph braid groups'](https://www.worldscientific.com/doi/10.1142/S0218196723500583), written by Daniel Berlyne (the author of this package). 

Given a graph braid group, this tool will search for a free splitting of the group, returning the splitting if one is found. The factors of this splitting are expressed as either free groups, simpler graph braid groups, or graphs of groups where the vertex groups and edge groups are simpler graph braid groups. Next, the tool expresses each of these free factors as a direct product if it is directly decomposable. The tool will then attempt to further split these direct factors, alternating between free products and direct products until a splitting with freely indecomposable factors is acquired. The tool will then evaluate any known graph braid groups that appear in the splitting. A list of known graph braid groups is given below.

Please note: it has not yet been proven whether all free splittings of graph braid groups can be constructed via the methods in ['Graph of groups decompositions of graph braid groups'](https://www.worldscientific.com/doi/10.1142/S0218196723500583). Thus, if this tool cannot find a free splitting, do not assume that one does not exist. Likewise, the factors in a splitting returned by this tool are freely indecomposable via the methods of ['Graph of groups decompositions of graph braid groups'](https://www.worldscientific.com/doi/10.1142/S0218196723500583), but may turn out to be freely decomposable via some other method.

# Instructions
1. Enter the data for your graph braid group in the `gbg_data.txt` file. Detailed instructions are given in the file, but briefly, this should include the following.
- The number of particles/strands.
- The adjacency matrix for your graph. Note that [this handy online graph calculator](https://www.mas.ncl.ac.uk/graph-curvature/)[^1] can easily produce an adjacency matrix from a hand-drawn graph, which can then be copy-pasted into `gbg_data.txt`.
- An initial configuration of the particles, if required.

[^1]: *The Graph Curvature Calculator and the curvatures of cubic graphs, Experimental Mathematics, 2019
(arXiv:1712.03033 [math.CO])*

2. Run `splitter.py` and follow the on-screen instructions.

# How to read output data in `splitting.txt`
If a free splitting is found, the factors will be displayed in the terminal. Factors will either be: 
- free groups, displayed in the form Z or F_k;
- other known groups, displayed in the form given in `known_gbgs.txt`;
- graph braid groups, displayed in the form B_n(Gamma_i) or RB_n(Gamma_i);
- graphs of groups, displayed in the form G_i.

Further information about the graphs Gamma_i and the graphs of groups G_i can be found in `splitting.txt`. 
1. For Gamma_i, its adjacency matrix is given.
2. For G_i, first the adjacency matrix of its graph is given. Then:
- for each vertex (numbered according to its row in the adjacency matrix), the data of its associated graph braid group is given, including the type of braid group and number of particles (displayed in the form B_n or RB_n), the adjacency matrix, and the initial configuration;
- for each edge (given as a 2-tuple of the row numbers of the vertices it connects), the data of its associated graph braid group is given.

Tip: Click 'load' in the bottom-right corner of [this handy online graph calculator](https://www.mas.ncl.ac.uk/graph-curvature/)[^1] and paste the adjacency matrix to quickly obtain a picture of the graph.

# Known graph braid groups
Below is a list of known graph braid groups that are non-cyclic and directly and freely indecomposable, with citations in the literature. These are used to aid in computations. The author would appreciate any further contributions to this list. The data for these is included in `known_gbgs.txt` (note: in order to make these graphs canonical, any vertices of degree 2 must be removed before adding the data to this list).
1. B_2(K_5) = pi_1(N_6), where K_5 is the complete graph on 5 vertices and N_6 is the closed non-orientable surface of genus 6.
- Example 3.18 of ['Graph of groups decompositions of graph braid groups'](https://www.worldscientific.com/doi/10.1142/S0218196723500583).
- Example 5.1 of ['Configuration spaces and braid groups of graphs'](https://www.proquest.com/docview/304583880).
2. B_2(K_{3,3}) = pi_1(N_4), where K_{3,3} is the complete bipartite graph on two sets of 3 vertices and N_4 is the closed non-orientable surface of genus 4.
- Example 3.19 of ['Graph of groups decompositions of graph braid groups'](https://www.worldscientific.com/doi/10.1142/S0218196723500583).
- Example 5.2 of ['Configuration spaces and braid groups of graphs'](https://www.proquest.com/docview/304583880).
3. B_3(\theta_4) = pi_1(S_3), where \theta_n is a generalised theta graph and S_3 is the closed orientable surface of genus 3.
- Example 4.3 of ['Characteristics of graph braid groups'](https://arxiv.org/pdf/1101.2648.pdf).