# About this tool
This tool implements computational methods developed in the mathematical paper ['Graph of groups decompositions of graph braid groups'](https://arxiv.org/pdf/2209.03860.pdf), written by Daniel Berlyne (the author of this package). 

Given a graph braid group, this tool will search for a free splitting of the group, returning the splitting if one is found. The factors of this splitting are expressed as either free groups, simpler graph braid groups, or graphs of groups where the vertex groups and edge groups are simpler graph braid groups. Next, the tool expresses each of these free factors as a direct product if it is directly decomposable. The tool will then attempt to further split these direct factors, alternating between free products and direct products until a splitting with freely indecomposable factors is acquired. The tool will then evaluate any known graph braid groups that appear in the splitting. A list of known graph braid groups is given below.

Please note: it has not yet been proven whether all free splittings of graph braid groups can be constructed via the methods in ['Graph of groups decompositions of graph braid groups'](https://arxiv.org/pdf/2209.03860.pdf). Thus, if this tool cannot find a free splitting, do not assume that one does not exist. Likewise, the factors in a splitting returned by this tool are freely indecomposable via the methods of ['Graph of groups decompositions of graph braid groups'](https://arxiv.org/pdf/2209.03860.pdf), but may turn out to be freely decomposable via some other method.

# Instructions
1. Enter the data for your graph braid group in the `gbg_data.txt` file. This includes:
- the number of particles/strands;
- the adjacency matrix for your graph;
- an initial configuration of the particles, if required.
Detailed instructions are given in the file.
2. Run `gbg_splitter.py` and follow the on-screen instructions.

# Known graph braid groups
Below is a list of known graph braid groups that are non-cyclic and directly and freely indecomposable, with citations in the literature. These are used to aid in computations. The author would appreciate any further contributions to this list. The data for these is included in `known_gbgs.txt` (note: in order to make these graphs canonical, any vertices of degree 2 must be removed before adding the data to this list).
1. B_2(K_5) = pi_1(N_6), where K_5 is the complete graph on 5 vertices and N_6 is the closed non-orientable surface of genus 6.
- Example 3.18 of ['Graph of groups decompositions of graph braid groups'](https://arxiv.org/pdf/2209.03860.pdf).
- Example 5.1 of ['Configuration spaces and braid groups of graphs'](https://www.proquest.com/docview/304583880).
2. B_2(K_{3,3}) = pi_1(N_4), where K_{3,3} is the complete bipartite graph on two sets of 3 vertices and N_4 is the closed non-orientable surface of genus 4.
- Example 3.19 of ['Graph of groups decompositions of graph braid groups'](https://arxiv.org/pdf/2209.03860.pdf).
- Example 5.2 of ['Configuration spaces and braid groups of graphs'](https://www.proquest.com/docview/304583880).
3. B_3(\theta_4) = pi_1(S_3), where \theta_n is a generalised theta graph and S_3 is the closed orientable surface of genus 3.
- Example 4.3 of ['Characteristics of graph braid groups'](https://arxiv.org/pdf/1101.2648.pdf).