#!/usr/bin/env python3
# gbg_splitter.py - Detects free splittings of graph braid groups

import math, copy, sys, ast, itertools
from collections import defaultdict

MatrixDimensionException = Exception('Adjacency matrix must be square.')
MatrixFormatException = Exception('Adjacency matrix must only contain positive integers and spaces.')
MatrixSymmetryException = Exception('Adjacency matrix must be symmetric.')
ParticleFormatException = Exception('Number of particles must be a positive integer.')
ConfigDimensionException = Exception('Number of integers in initial configuration must equal number of particles.')
ConfigFormatException = Exception('Initial configuration must only contain spaces and integers between 1 and the number of rows in the adjacency matrix.')
VertexException = Exception('Each connected component of the graph must have at least as many vertices as particles starting in that component.')


#############################
# General purpose functions #
#############################

# Returns all the different partitions of a non-negative integer as a sum of a fixed number of non-negative integers.
# The partitions are returned in lexicographic order.
def integer_partitions(total_sum, num_parts):
    if num_parts == 0:
        yield tuple()
    elif num_parts == 1:
        yield (total_sum,)
    else:
        for value in range(total_sum + 1):
            for partition in integer_partitions(total_sum - value, num_parts - 1):
                yield (value,) + partition

def start_exit_sequence():
    while True:
        exit_sequence = input('Press [Enter] to exit.')
        if not exit_sequence:
            sys.exit()


###########
# Classes #
###########

# Class for performing graph computations.
class Graph:
    # A graph is defined by its adjacency matrix.
    def __init__(self, adj_matrix):
        # The adjacency matrix should be a tuple of tuples or a list of lists.
        # Rows are numbered using `i` and columns are numbered using `j`.
        self.adj_matrix = adj_matrix
        # The vertices of the graph are numbered and encoded in a list.
        self.vertices = [i for i in range(len(adj_matrix))]
        self.num_vertices = len(self.vertices)
        # The edges of the graph are encoded as a list of 2-tuples (i, j) where i <= j.
        # That is, we only use the upper triangle of the matrix.
        self.edges = [(i, j) for j in range(self.num_vertices) for i in range(j+1) for n in range(adj_matrix[i][j])]
        self.num_edges = len(self.edges)
        self.essential_vertices = [v for v in self.vertices if self.get_degree(v) != 2]

    # Returns dictionary of connected components of the graph, where:
    # the keys are 2-tuples of vertex sets and edge sets of the connected components of the graph, considered as subgraphs;
    # the values are the components as standalone 'Graph' objects.
    def get_connected_components(self, subgraph = None):
        if subgraph == None:
            subgraph = (self.vertices, self.edges)
        # The empty subgraph has no connected components.
        if subgraph[0] == []:
            return {}
        # Otherwise, start with the component containing first vertex of `subgraph`
        components = {self.get_component(subgraph[0][0], subgraph)[0]: self.get_component(subgraph[0][0], subgraph)[1]}
        for v in subgraph[0]:
            if v not in [w for component in components for w in component[0]]:
                components.update({self.get_component(v, subgraph)[0]: self.get_component(v, subgraph)[1]})
        return components

    # Returns a 2-tuple where the first entry is the 2-tuple (vertex set, edge set) of the connected component of `subgraph` containing `vertex` 
    # (`subgraph` defaults to the original graph) and the second entry is the connected component as a standalone Graph object.
    # Note, `subgraph` should be a (vertex set, edge set) 2-tuple of sublists of `self.vertices` and `self.edges`
    def get_component(self, vertex, subgraph = None):
        if subgraph == None:
            subgraph = (self.vertices, self.edges)
        component_vertices = [vertex]
        component_edges = []
        self.iterate_component(component_vertices, component_edges, vertex, subgraph)
        # Sort vertices and edges to avoid KeyErrors when trying to access components in dictionary.
        # Vertex values in the new Graph object will be their indices in this list.
        component_vertices.sort()
        component_edges.sort()
        component_adj_matrix = [[0 for j in range(len(component_vertices))] for i in range(len(component_vertices))]
        for edge in component_edges:
            if edge[0] != edge[1]:
                # Component vertices must be referred to by their index in the list when constructing the adjacency matrix
                component_adj_matrix[component_vertices.index(edge[0])][component_vertices.index(edge[1])] += 1
                component_adj_matrix[component_vertices.index(edge[1])][component_vertices.index(edge[0])] += 1
            else:
                component_adj_matrix[component_vertices.index(edge[0])][component_vertices.index(edge[1])] += 1
        return ((tuple(component_vertices), tuple(component_edges)), Graph(component_adj_matrix))
    
    def iterate_component(self, component_vertices, component_edges, current_vertex, subgraph = None):
        if subgraph == None:
            subgraph = (self.vertices, self.edges)
        for edge in subgraph[1]:
            if current_vertex in edge and edge not in component_edges:
                component_edges.append(edge)
                adjacent_vertex = edge[edge.index(current_vertex) - 1]
                if adjacent_vertex not in component_vertices:    
                    component_vertices.append(adjacent_vertex)
                    self.iterate_component(component_vertices, component_edges, adjacent_vertex, subgraph)

    def get_num_connected_components(self):
        return len(self.get_connected_components())

    # Returns a 2-tuple where: 
    # the first entry is the 2-tuple of vertices and edges of the graph minus `removed_edges`, considered as a subgraph of the original graph;
    # the second entry is the graph minus `removed_edges` as an object of class Graph.
    def get_graph_minus_open_edges(self, removed_edges):
        adj_matrix_minus_open_edges = copy.deepcopy(self.adj_matrix)
        subgraph_vertices = [v for v in self.vertices]
        subgraph_edges = [e for e in self.edges]
        for (i, j) in removed_edges:
            if i == j:
                adj_matrix_minus_open_edges[i][j] -= 1
            else:
                adj_matrix_minus_open_edges[i][j] -= 1
                adj_matrix_minus_open_edges[j][i] -= 1
            subgraph_edges.remove((i, j))
        return ((subgraph_vertices, subgraph_edges), Graph(adj_matrix_minus_open_edges))

    # Same as above, but we remove a single closed edge
    # Note, this is equivalent to removing this edge's two vertices together with all edges containing one of these vertices
    def get_graph_minus_closed_edge(self, removed_edge):
        (i, j) = removed_edge
        subgraph_vertices = [v for v in self.vertices if v not in (i, j)]
        subgraph_edges = [e for e in self.edges if i not in e and j not in e]
        # Removing edge (i, j) corresponds to removing the ith and jth rows and columns from the adjacency matrix
        adj_matrix_minus_closed_edge = [[self.adj_matrix[k][l] for l in range(len(self.adj_matrix[k])) if l not in (i, j)] for k in range(len(self.adj_matrix)) if k not in (i, j)]
        return ((subgraph_vertices, subgraph_edges), Graph(adj_matrix_minus_closed_edge))
    
    # We remove a single vertex from our graph, together with all edges containing the vertex.
    def get_graph_minus_vertex(self, removed_vertex):
        subgraph_vertices = [v for v in self.vertices if v != removed_vertex]
        subgraph_edges = [e for e in self.edges if removed_vertex not in e]
        # Removing the vertex corresponds to removing its row and column from the adjacency matrix
        adj_matrix_minus_vertex = [[self.adj_matrix[k][l] for l in range(len(self.adj_matrix[k])) if l != removed_vertex] for k in range(len(self.adj_matrix)) if k != removed_vertex]
        return ((subgraph_vertices, subgraph_edges), Graph(adj_matrix_minus_vertex))

    # Goes through the list `edges` one by one and checks if the edge is non-separating. If so, removes that (open) edge and repeats on the new graph.
    # Returns the pruned graph and a list of the edges that were removed.
    def prune(self, edges):
        modified_graph = self
        removed_edges = []
        for edge in edges:
            if not modified_graph.is_separating(edge):
                removed_edges.append(edge)
                modified_graph = modified_graph.get_graph_minus_open_edges([edge])[1]
        return (modified_graph, removed_edges)

    # Returns True if `edge` separates the graph into two connected components.
    def is_separating(self, edge):
        if self.get_graph_minus_open_edges([edge])[1].get_num_connected_components() > 1:
            return True
        else:
            return False
    
    def get_degree(self, vertex):
        degree = 0
        for edge in self.edges:
            if edge == (vertex, vertex):
                degree += 2
            elif vertex in edge:
                degree += 1
        return degree
    
    # Returns a list of all vertices that can be reached from `vertex` via a path of length at least 1 and at most `distance`.
    # The ball is centreless, meaning if `vertex` appears in the list, this means it is contained in a cycle.
    def get_centreless_ball(self, vertex, distance):
        ball = []
        visited_edges = []
        if distance > 0:
            self.iterate_ball(ball, visited_edges, vertex, distance)
        return ball
    
    def iterate_ball(self, ball, visited_edges, current_vertex, remaining_distance):
        for edge in self.edges:
            if current_vertex in edge and edge not in visited_edges:
                adjacent_vertex = edge[edge.index(current_vertex) - 1]
                if adjacent_vertex not in ball:
                    ball.append(adjacent_vertex)
                    visited_edges.append(edge)
                    if remaining_distance >= 2:
                        self.iterate_ball(ball, visited_edges, adjacent_vertex, remaining_distance - 1)

    # Removes all vertices of degree 2 from the graph so that only essential vertices remain.
    def make_essential(self):
        non_essential_vertices = [v for v in self.vertices if v not in self.essential_vertices]
        # If all vertices have degree 2, then we have a cycle graph and must therefore designate one degree 2 vertex as essential.
        if self.essential_vertices == []:
            non_essential_vertices.pop()
        new_edges = []
        removed_edges = []
        for v in non_essential_vertices:
            adjacent_edges = [edge for edge in self.edges if edge not in removed_edges and v in edge] + [edge for edge in new_edges if v in edge]
            edge_1 = adjacent_edges[0]
            edge_2 = adjacent_edges[1]
            removed_edges.append(edge_1)
            removed_edges.append(edge_2)
            new_edges = [edge for edge in new_edges if edge not in removed_edges]
            new_edges.append((edge_1[edge_1.index(v) - 1], edge_2[edge_2.index(v) - 1]))
        modified_adj_matrix = [[self.adj_matrix[i][j] for j in self.vertices] for i in self.vertices]
        for (i, j) in new_edges:
            modified_adj_matrix[i][j] += 1
            modified_adj_matrix[j][i] += 1
        essential_adj_matrix = [[modified_adj_matrix[i][j] for j in self.essential_vertices] for i in self.essential_vertices]
        return Graph(essential_adj_matrix)

# Class for processing graphs of groups
# Note, monomorphisms are omitted since we do not need this data to detect free splittings.
# For this reason, we do not need our graphs to be directed, either.
# See [Graph of groups decompositions of graph braid groups](https://arxiv.org/pdf/2209.03860.pdf) for a description of the monomorphisms.
class GraphOfGroups:
    def __init__(self, graph, vertex_groups, edge_groups):
        # `graph` should be an object of class `Graph`.
        self.graph = graph
        # `vertex_groups` should be a dictionary, with keys given by `graph.vertices`.
        self.vertex_groups = vertex_groups
        # `edge_groups` should be a dictionary, with keys given by `graph.edges`. Note, `graph.edges` may contain repeat entries for multi-edges. 
        # In this case, extra entries should be added to each edge 2-tuple to make it a unique key.
        self.edge_groups = edge_groups

    # Returns a splitting of the fundamental group of the graph of groups as a free product of non-trivial graphs of groups, given as a tuple of the factors.
    # Otherwise, returns the original graph of groups.
    def get_free_splitting(self):
        # Note: `is_trivial` method must be implemented for the class that the edge group is in.
        trivial_edges = [edge for edge in self.edge_groups if self.edge_groups[edge].is_trivial()]
        # Turns edges in the graph of groups (which are labelled to distinguish multi-edges) into actual edges in `graph.edges`.
        trivial_edges_unlabelled = [(edge[0], edge[1]) for edge in trivial_edges]
        if len(trivial_edges) == 0:
            return (self,)
        # First, remove non-separating edges with trivial edge group, one by one.
        (modified_graph, removed_edges_unlabelled) = self.graph.prune(trivial_edges_unlabelled)
        removed_edges = [trivial_edges[trivial_edges_unlabelled.index(edge)] for edge in removed_edges_unlabelled]
        # We record the edges we removed, which each correspond to a free Z factor, then continue our analysis on the graph minus these edges.
        num_free_Z = len(removed_edges)
        modified_edge_groups = {edge: self.edge_groups[edge] for edge in self.edge_groups if edge not in removed_edges}
        modified_gog = GraphOfGroups(modified_graph, self.vertex_groups, modified_edge_groups)
        # If a separating edge has trivial edge group and all the vertex groups of a connected component are trivial, then we throw away that commponent.
        trivial_sep_edges = [edge for edge in trivial_edges if edge not in removed_edges]
        # This is a 3-tuple of a subgraph, a list of the associated vertex groups, and a list of the associated edge groups
        reduced_gog_data = modified_gog.reduce(trivial_sep_edges)
        reduced_trivial_sep_edges = [edge for edge in trivial_sep_edges if edge in reduced_gog_data[2]]
        # If we threw away everything when reducing, that means all of our remaining factors are trivial.
        if reduced_gog_data[0][0] == []:
            factor_gogs = []
        # Otherwise, remove all remaining separating edges with trivial edge group. The fundamental groups of the connected components will be our factors.
        else:
            split_graph = ([v for v in reduced_gog_data[0][0]], [(e[0], e[1]) for e in reduced_gog_data[2] if e not in reduced_trivial_sep_edges])
            factor_components = {self.trim(component): self.graph.get_connected_components(self.trim(component))[self.trim(component)] for component in self.graph.get_connected_components(split_graph)}
            # List of (subgraph, vertex groups, edge groups) tuples for each connected component of `split_graph`.
            factor_data = [(component, {v: self.vertex_groups[v] for v in component[0]}, {e: self.edge_groups[e] for e in reduced_gog_data[2] if (e[0], e[1]) in component[1]}) for component in factor_components]
            # Now turn this data into standalone graphs of groups.
            factor_gogs = []
            for data in factor_data:
                factor_graph = factor_components[data[0]]
                factor_vertex_groups = {v: data[1][data[0][0][v]] for v in factor_graph.vertices}
                # Need to carefully reindex edges due to identical multi-edges in the Graph object that may have different edge groups.
                multi_edges = [e for e in data[2] if len(e) > 2]
                non_multi_edges = [e for e in data[2] if len(e) == 2]
                old_edges = multi_edges + non_multi_edges
                reindexed_multi_edges = [(data[0][0].index(e[0]), data[0][0].index(e[1])) + e[2:] for e in multi_edges]
                reindexed_non_multi_edges = [(data[0][0].index(e[0]), data[0][0].index(e[1])) for e in non_multi_edges]
                reindexed_edges = reindexed_multi_edges + reindexed_non_multi_edges
                factor_edge_groups = {e: data[2][old_edges[reindexed_edges.index(e)]] for e in reindexed_edges}
                factor_gogs.append(GraphOfGroups(factor_graph, factor_vertex_groups, factor_edge_groups))
        return (f'F_{num_free_Z}',) + tuple(factor_gogs)
    
    # Given a collection of separating edges with trivial edge groups, checks for triviality of the fundamental groups of the connected components of the 
    # graph minus the edge, for each edge. If any are trivial, removes these components. Then removes degree 1 vertex groups whose incident edge group is
    # equal to the vertex group. Returns the remainder, given as a 3-tuple subgraph of groups.
    def reduce(self, trivial_sep_edges):
        if trivial_sep_edges == []:
            return ((self.graph.vertices, self.graph.edges), self.vertex_groups, self.edge_groups)
        trivial_components = []
        reduced_graph = self.iterate_reduction((self.graph.vertices, self.graph.edges), trivial_sep_edges, trivial_components)
        reduced_vertex_groups = {v: self.vertex_groups[v] for v in self.vertex_groups if v in reduced_graph[0]}
        # May have multiple copies of the same edge in `reduced_graph`; need to reapply labelling to distinguish them before assigning edge groups.
        all_edges = [edge for edge in self.edge_groups]
        all_edges_unlabelled = [(edge[0], edge[1]) for edge in all_edges]
        reduced_edges = [all_edges[all_edges_unlabelled.index(edge)] for edge in reduced_graph[1]]
        reduced_edge_groups = {e: self.edge_groups[e] for e in reduced_edges}
        # Need to convert above from subgraph form to Graph form.
        return (reduced_graph, reduced_vertex_groups, reduced_edge_groups)

    def iterate_reduction(self, subgraph, trivial_sep_edges, trivial_components):
        edge = trivial_sep_edges[0]
        graph_minus_edge = (subgraph[0], [e for e in subgraph[1] if e != edge])
        components = self.graph.get_connected_components(graph_minus_edge)
        for component in components:
            trivial = True
            for v in component[0]:
                if self.vertex_groups[v].is_trivial():
                    continue
                trivial = False
                break
            if trivial:
                trivial_components.append(component)
        trivial_component_vertices = [v for component in trivial_components for v in component[0]]
        trivial_component_edges = [e for component in trivial_components for e in component[1]]
        reduced_graph = ([v for v in subgraph[0] if v not in trivial_component_vertices], [e for e in subgraph[1] if e[0] not in trivial_component_vertices and e[1] not in trivial_component_vertices])
        reduced_sep_edges = [e for e in trivial_sep_edges if e not in trivial_component_edges and e != edge]
        if reduced_sep_edges == []:
            return reduced_graph
        else:
            return self.iterate_reduction(reduced_graph, reduced_sep_edges, trivial_components)
        
    # Iteratively removes degree 1 vertices from `subgraph` whose vertex groups are the same as the incident edge group.
    def trim(self, subgraph):
        # We trim after all trivial edges have been removed, so any remaining edge should be labelled with a non-trivial edge group.
        non_trivial_edges = [edge for edge in self.edge_groups if not self.edge_groups[edge].is_trivial()]
        non_trivial_edges_unlabelled = [(edge[0], edge[1]) for edge in non_trivial_edges]
        reduced_graph = ([v for v in subgraph[0]], [e for e in subgraph[1]])
        for v in subgraph[0]:
            incident_edges = []
            for e in subgraph[1]:
                if e == (v, v):
                    incident_edges.append(e)
                    incident_edges.append(e)
                elif v in e:
                    incident_edges.append(e)
            if len(incident_edges) == 1:
                incident_edge_group = self.edge_groups[non_trivial_edges[non_trivial_edges_unlabelled.index(incident_edges[0])]]
                vertex_group = self.vertex_groups[v]
                if incident_edge_group.is_same(vertex_group):
                    reduced_graph = ([vertex for vertex in reduced_graph[0] if vertex != v], [edge for edge in reduced_graph[1] if edge != e])
        if reduced_graph == subgraph:
            return (tuple(reduced_graph[0]), tuple(reduced_graph[1]))
        else:
            return self.trim(reduced_graph)

# Class for processing graph braid groups.
# Note: This handles reduced graph braid groups, which are only isomorphic to the graph braid group if `is_reduced` returns True.
class GraphBraidGroup:
    def __init__(self, graph, num_particles, initial_config = None):
        # `graph` should be a simple graph (i.e. no loops or multi-edges), expressed as an object of class Graph.
        self.graph = graph
        self.num_particles = num_particles
        self.graph_components = self.graph.get_connected_components()
        # `initial_config` should be a list of vertices of `graph` of length `num_particles`.
        # The initial configuration serves as a base point in case the configuration space is not connected.
        # If the config space is connected, `initial_config` can be left blank and defaults to the first `num_particles` vertices of the graph.
        if initial_config == None:
            self.initial_config = [self.graph.vertices[i] for i in range(num_particles)]
        else:
            self.initial_config = initial_config
        # Dictionary whose keys are the connected components of the graph containing vertices in `initial_config`, as (vertex_set, edge_set)-tuples,
        # and whose values are the number of particles that start in that component.
        self.num_initial_particles_per_component = self.get_num_particles_per_component(self.initial_config)
        # Uses Corollary 2.16 of [Graph of groups decompositions of graph braid groups](https://arxiv.org/pdf/2209.03860.pdf).
        self.num_connected_components = math.comb(self.num_particles + self.graph.get_num_connected_components() - 1, self.graph.get_num_connected_components() - 1)

    # Returns the graph of groups decomposition given by splitting along the edges of `graph` in list `edges`.
    # The edges in `edges` must all share a common vertex.
    # Uses Theorem 3.5 of [Graph of groups decompositions of graph braid groups](https://arxiv.org/pdf/2209.03860.pdf).
    def get_graph_of_groups(self, edges):
        # `graph_minus_open_edges` is a 2-tuple consisting of:
        # - a (vertex set, edge set) 2-tuple as a subgraph of the original graph;
        # - the graph as a standalone Graph object.
        graph_minus_open_edges = self.graph.get_graph_minus_open_edges(edges)
        # List of dictionaries, where each dictionary assigns a number of particles to each component of `graph_minus_open_edges`.
        # Components are given as (vertex set, edge set) 2-tuples.
        # This list consists of all such assignments that are compatible with `initial_config`.
        compatible_particles_per_component = self.get_compatible_particles_per_component(edges)
        # The vertex groups are the different braid groups of the graph minus the open edges that arise from different initial configurations.
        # These correspond to all the different assignments of particles to the connected components of `graph_minus_open_edges` that are 
        # compatible with `initial_config`.
        # We therefore enumerate the particle assignments by putting them in a list and construct a dictionary where:
        # - the keys are the indices of the particle assignments in the list; 
        # - the values are the braid groups corresponding to the particle assignments.
        vertex_particle_assignments = [particle_assignment for particle_assignment in compatible_particles_per_component]
        vertex_groups = {i: GraphBraidGroup(graph_minus_open_edges[1], self.num_particles, self.generate_initial_config(vertex_particle_assignments[i])) 
                         for i in range(len(vertex_particle_assignments))}
        graphs_minus_closed_edges = {edge: self.graph.get_graph_minus_closed_edge(edge) for edge in edges}
        # The edge groups are braid groups of the graph minus a single closed edge in `edges`, with one fewer particle.
        # An edge group connects two vertex groups if the corresponding particle assignments differ by moving a single particle (the 'active particle') 
        # along the closed edge.
        # The initial configuration of the edge group is an initial configuration of either of the vertex groups with the active particle removed.
        # Note, removing the closed edge may further disconnect the graph, thus there may be multiple possible configurations, 
        # corresponding to the different compatible assignments of particles to connected components of the graph minus the closed edge.
        # There is one edge group for each compatible assignment.
        # We therefore construct a dictionary where:
        # - the keys are 4-tuples consisting of
        #  -- the indices of the particle assignments of the two vertex groups
        #  -- the edge by which the assignments differ
        #  -- the index of a compatible assignment for the edge group in the list of all compatible assignments for that edge group;
        # - the values are the corresponding graph braid groups.
        edge_groups = {(sorted([vertex_particle_assignments.index(assignment_1), vertex_particle_assignments.index(assignment_2)])[0], 
                        sorted([vertex_particle_assignments.index(assignment_1), vertex_particle_assignments.index(assignment_2)])[1], 
                        edge, 
                        self.get_compatible_assignments(edges, assignment_1, assignment_2, edge).index(compatible_assignment)): 
                        GraphBraidGroup(graphs_minus_closed_edges[edge][1], self.num_particles - 1, self.reindex(self.generate_initial_config(compatible_assignment), [edge[0], edge[1]])) 
                        for (assignment_1, assignment_2, edge) in self.get_adjacent_assignments(edges, compatible_particles_per_component) 
                        for compatible_assignment in self.get_compatible_assignments(edges, assignment_1, assignment_2, edge)}
        # Adjacency matrix of the graph of groups.
        gog_adj_matrix = [[0 for v in vertex_groups] for w in vertex_groups]
        for e in edge_groups:
            if e[0] != e[1]:
                gog_adj_matrix[e[0]][e[1]] += 1
                gog_adj_matrix[e[1]][e[0]] += 1
            else:
                gog_adj_matrix[e[0]][e[1]] += 1
        gog_graph = Graph(gog_adj_matrix)
        return GraphOfGroups(gog_graph, vertex_groups, edge_groups)
  
    # Returns a list of dictionaries, where each dictionary assigns a number of particles to each component of the graph minus the open edges in `edges`.
    # This list consists of all such assignments that are compatible with `initial_config`.
    # Note, the edges in `edges` must share a common vertex.
    def get_compatible_particles_per_component(self, edges):
        if len(edges) > 1:
            edges_as_sets = [set(edge) for edge in edges]
            (common_vertex,) = set.intersection(*edges_as_sets)
        else:
            # If we are splitting over a single edge, just choose an arbitrary vertex of that edge.
            common_vertex = edges[0][0]
        graph_minus_open_edges = self.graph.get_graph_minus_open_edges(edges)
        # Tuple of the connected components of `graph_minus_open_edges`, where each component is expressed as a (vertex set, edge set) 2-tuple.
        gmoe_components = tuple([component for component in self.graph.get_connected_components(graph_minus_open_edges[0])])
        # Dictionary where:
        # - keys are components of `graph` that do not contain the common vertex of `edges`, i.e. remain unchanged in `graph_minus_open_edges`, as (vertex set, edge set) 2-tuples;
        # - values are the number of particles in that component in `initial_config`.
        old_particles_per_component = {component: self.num_initial_particles_per_component[component] for component in self.graph_components if common_vertex not in component[0]}
        old_partition = tuple([old_particles_per_component[component] for component in old_particles_per_component])
        num_old_components = len(old_partition)
        if num_old_components == 0:
            total_old_particles = 0
        else:
            total_old_particles = sum(i for i in old_partition)
        old_component_vertex_sets = [set(component[0]) for component in old_particles_per_component]
        new_components = [component for component in gmoe_components if set(component[0]) not in old_component_vertex_sets]
        new_partitions = integer_partitions(self.num_particles - total_old_particles, graph_minus_open_edges[1].get_num_connected_components() - num_old_components)
        new_particles_per_component = [{new_component: new_partition[new_components.index(new_component)] for new_component in new_components} for new_partition in new_partitions if self.has_sufficient_capacity(new_components, new_partition)]
        # List of dictionaries of particle assignments consisting of `old_particles_per_component` augmented with each of the assignments in `new_particles_per_component`.
        compatible_particles_per_component = [{component: old_particles_per_component[component] for component in old_particles_per_component} for i in range(len(new_particles_per_component))]
        for i in range(len(new_particles_per_component)):    
            compatible_particles_per_component[i].update(new_particles_per_component[i])
        return compatible_particles_per_component
    
    # Takes as input a list of connected components of a graph as (vertex set, edge set) tuples, and an integer partition of the same length.
    # Returns True if each component has at least as many vertices as the corresponding integer in the partition.
    def has_sufficient_capacity(self, components, partition):
        for component in components:
            if len(component[0]) < partition[components.index(component)]:
                return False
        return True
    
    def get_num_particles_per_component(self, config):
        num_particles_per_component = defaultdict(int)
        for v in config:
                num_particles_per_component[self.graph.get_component(v)[0]] += 1
        return num_particles_per_component

    # Given the connected components of a subgraph and the number of particles in each component, 
    # returns a configuration of particles on the subgraph with that number of particles in each component.
    # This configuration is given as a list of vertices.
    # `particles_per_component` should be a dictionary where:
    # - the keys are connected components as (vertex set, edge set) 2-tuples;
    # - the values are the number of particles in that component.
    def generate_initial_config(self, particles_per_component):
        initial_config_by_component = [[component[0][i] for i in range(particles_per_component[component])] for component in particles_per_component]
        return [vertex for component_config in initial_config_by_component for vertex in component_config]
    
    # Takes as input:
    # - a list of edges;
    # - a list of assignments of particles to the connected components of the graph minus the open edges in `edges`.
    # Each assignment in the list is given as a dictionary, where:
    # - keys are the connected components of the graph minus the open edges, given as (vertex set, edge set) 2-tuples;
    # - values are the number of particles assigned to that component.
    # Returns a list of 3-tuples, where each 3-tuple consists of: 
    # - two assignments that differ by moving a single particle along one of the open edges;
    # - the edge that they differ by.
    # Note, if an edge connects a connected component to itself, then we get (assignment, assignment, edge) for each assignment.
    # Uses Proposition 3.7 of [Graph of groups decompositions of graph braid groups](https://arxiv.org/pdf/2209.03860.pdf).
    def get_adjacent_assignments(self, edges, particle_assignments):
        gmoe = self.graph.get_graph_minus_open_edges(edges)
        adjacent_assignments = []
        for edge in edges:
            for i_1, assignment_1 in enumerate(particle_assignments):
                for assignment_2 in particle_assignments[i_1:]:
                    differences = {component: assignment_2[component] - assignment_1[component] for component in assignment_1 if assignment_2[component] - assignment_1[component] != 0}
                    # Checks if the assignments are equal and `edge` connects the same component to itself
                    if len(differences) == 0 and gmoe[1].get_component(edge[0])[0] == gmoe[1].get_component(edge[1])[0]:
                        adjacent_assignments.append((assignment_1, assignment_2, edge))
                    # Checks if `assignment_2` differs from `assignment_1` by a single particle along an edge.
                    elif sum(differences[component] for component in differences if differences[component] >= 1) == 1:
                        # First vertex set is for the component in which `assignment_1` has more particles.
                        differing_component_vertex_sets = [component[0] for component in differences if differences[component] < 0] + [component[0] for component in differences if differences[component] > 0]
                        # Order the assignments in the 3-tuple according to the order of the vertices in `edge`.
                        if edge[0] in differing_component_vertex_sets[0] and edge[1] in differing_component_vertex_sets[1]:
                            # Particle moves along `edge` from `assignment_1` to `assignment_2`.
                            adjacent_assignments.append((assignment_1, assignment_2, edge))
                        elif edge[0] in differing_component_vertex_sets[1] and edge[1] in differing_component_vertex_sets[0]:
                            # Particle moves along `edge` from `assignment_2` to `assignment_1`.
                            adjacent_assignments.append((assignment_2, assignment_1, edge))
        return adjacent_assignments

    # Takes as input:
    # - a list `edges` of the open edges being removed in the vertex spaces;
    # - two particle assignments that differ by moving a single particle (the 'active particle') along one of the edges in `edges`;
    # - the edge `active_edge` that the particle assignments differ by.
    # Returns a list of dictionaries, where each dictionary assigns a number of particles to each component of the graph minus the closed edge `active_edge`.
    # This list consists of all such assignments that are compatible with both `assignment_1` and `assignment_2`.
    # Note, the total number of particles is one less than `self.num_particles`; 
    # if the missing particle is placed at the vertices of `active_edge`, then we get `assignment_1` and `assignment_2`, respectively.
    # Uses Proposition 3.8 of [Graph of groups decompositions of graph braid groups](https://arxiv.org/pdf/2209.03860.pdf).
    def get_compatible_assignments(self, edges, assignment_1, assignment_2, active_edge):
        gmoe = self.graph.get_graph_minus_open_edges(edges)
        # Dictionary of connected components of `gmoe`, where keys are (vertex set, edge set) tuples and values are Graph objects.
        gmoe_components = self.graph.get_connected_components(gmoe[0])
        # The connected components of the graph minus the open edges that the active particle starts and ends in, as (vertex set, edge set) tuples.
        # Note, these components may be the same. ('Start' and 'end' are with respect to the order of the vertices in `active_edge`.)
        (initial_component,) = [component for component in gmoe_components if active_edge[0] in component[0]]
        (terminal_component,) = [component for component in gmoe_components if active_edge[1] in component[0]]
        # Particle assignments for all components of the graph minus the open edges except the ones that the active particle is moving between.
        old_assignment = {component: assignment_1[component] for component in assignment_1 if component not in (initial_component, terminal_component)}
        # Recall, the value of a vertex in a standalone Graph object is equal to its index in the list of vertices as a subgraph.
        initial_vertex_reindexed = initial_component[0].index(active_edge[0])
        terminal_vertex_reindexed = terminal_component[0].index(active_edge[1])
        
        if initial_component != terminal_component:
            # If `active_edge` connects two different components of gmoe, then we must consider all assignments of particles in `initial_component` and
            # `terminal_component` to components of the initial/terminal component minus the initial/terminal vertex, respectively.
            initial_component_minus_vertex = gmoe_components[initial_component].get_graph_minus_vertex(initial_vertex_reindexed)
            # This is `initial_component_minus_vertex` as a subgraph of the original graph, i.e. a (vertex set, edge set) tuple.
            icmv_subgraph = ([initial_component[0][v] for v in initial_component_minus_vertex[0][0]], [(initial_component[0][e[0]], initial_component[0][e[1]]) for e in initial_component_minus_vertex[0][1]])
            initial_subcomponents = [component for component in self.graph.get_connected_components(icmv_subgraph)]
            terminal_component_minus_vertex = gmoe_components[terminal_component].get_graph_minus_vertex(terminal_vertex_reindexed)
            tcmv_subgraph = ([terminal_component[0][v] for v in terminal_component_minus_vertex[0][0]], [(terminal_component[0][e[0]], terminal_component[0][e[1]]) for e in terminal_component_minus_vertex[0][1]])
            terminal_subcomponents = [component for component in self.graph.get_connected_components(tcmv_subgraph)]
            subcomponents = initial_subcomponents + terminal_subcomponents
            # We use `assignment_2` for the partitions of the initial component because we want to exclude the active particle,
            # which is in the terminal component in `assignment_2`. Similarly, we choose `assignment_1` for partitions of the terminal component.
            new_partitions_initial = integer_partitions(assignment_2[initial_component], len(initial_subcomponents))
            new_partitions_terminal = integer_partitions(assignment_1[terminal_component], len(terminal_subcomponents))
            new_partitions = [initial_partition + terminal_partition for initial_partition in new_partitions_initial for terminal_partition in new_partitions_terminal]
            new_assignments = [{component: partition[subcomponents.index(component)] for component in subcomponents} for partition in new_partitions if self.has_sufficient_capacity(subcomponents, partition)]
        else:
            # If `active_edge` connects a component of gmoe to itself, then we must consider all assignments of particles in `initial_component` to components
            # of `initial_component` minus the initial and terminal vertices. 
            initial_component_minus_vertices = gmoe_components[initial_component].get_graph_minus_closed_edge((initial_vertex_reindexed, terminal_vertex_reindexed))
            icmv_subgraph = ([initial_component[0][v] for v in initial_component_minus_vertices[0][0]], [(initial_component[0][e[0]], initial_component[0][e[1]]) for e in initial_component_minus_vertices[0][1]])
            subcomponents = [component for component in self.graph.get_connected_components(icmv_subgraph)]
            # Need to subtract 1 because we want to exclude the active particle (which is present in the initial component in both assignments this time).
            new_partitions = integer_partitions(assignment_2[initial_component] - 1, len(subcomponents))
            new_assignments = [{component: partition[subcomponents.index(component)] for component in subcomponents} for partition in new_partitions if self.has_sufficient_capacity(subcomponents, partition)]
        
        compatible_assignments = [{component: old_assignment[component] for component in old_assignment} for i in range(len(new_assignments))]
        for i in range(len(new_assignments)):    
            compatible_assignments[i].update(new_assignments[i])
        return compatible_assignments

    # `old_vertices` is a list of vertices of the original graph that doesn't contain any vertices in `removed_vertices`.
    # Returns the values of `old_vertices` as vertices of the Graph object consisting of the original graph with `removed_vertices` removed.
    def reindex(self, old_vertices, removed_vertices):
        new_vertices = []
        for v in old_vertices:
            shift = 0
            for w in removed_vertices:
                if v > w:
                    shift += 1
            v -= shift
            new_vertices.append(v)
        return new_vertices

    # Returns True if the graph braid group is the trivial group.
    def is_trivial(self):
        components = self.graph_components
        # Check if any connected components of the graph have non-trivial braid group.
        # We loop over `num_initial_particles_per_component` because that only includes components that contain at least one particle.
        for comp in self.num_initial_particles_per_component:
            # If the graph contains a cycle and has more vertices than particles, then its braid group is non-trivial.
            for edge in components[comp].edges:
                if not components[comp].is_separating(edge) and components[comp].num_vertices > self.num_initial_particles_per_component[comp]:
                    return False
            # If the graph contains at least two particles and a vertex of degree at least three, and
            # has at least two more vertices than particles, then its braid group is non-trivial.
            if self.num_initial_particles_per_component[comp] >= 2 and components[comp].num_vertices >= self.num_initial_particles_per_component[comp] + 2:
                for vertex in components[comp].vertices:
                    if components[comp].get_degree(vertex) > 2:
                        return False
        return True
    
    # Returns True if the braid group of the graph is isomorphic to the reduced braid group of the graph.
    def is_reduced(self):
        components = self.graph_components
        for comp in self.num_initial_particles_per_component:
            for v in components[comp].vertices:
                # All cycles in the component must have length at least n+1, where n is the number of particles in the component.
                ball_cycle = components[comp].get_centreless_ball(v, self.num_initial_particles_per_component[comp])
                if v in ball_cycle:
                    return False
                # All paths between essential vertices of the component must have length at least n-1.
                if v in components[comp].essential_vertices:
                    ball_essential = components[comp].get_centreless_ball(v, self.num_initial_particles_per_component[comp] - 2)
                    for w in ball_essential:
                        if w in components[comp].essential_vertices:
                            return False
        return True
    
    # Returns `final_splitting`, a list of direct factors of the graph braid group.
    # Each direct factor in the list is expressed as a list of its free factors.
    # Each free factor in this list is expressed as a list of its direct factors.
    # Each of these direct factors is expressed as a list of its free factors, and so on, terminating in freely indecomposable factors.
    def get_splitting(self):
        final_splitting = []
        self.iterate_splitting(final_splitting)
        return final_splitting

    def iterate_splitting(self, current_splitting):
        # Filter out the trivial factors of the graph braid group.
        non_trivial_direct_factors = self.factorise()
        # Split each non-trivial factor.
        direct_factor_splittings = {direct_factor: (direct_factor,) for direct_factor in non_trivial_direct_factors}
        for direct_factor in non_trivial_direct_factors:
            splittings = [direct_factor.get_graph_of_groups([edge]).get_free_splitting() for edge in direct_factor.graph.edges]
            # If there is a splitting that gives a free group, choose this.
            for splitting in splittings:
                if len(splitting) == 1 and type(splitting[0]) == str:
                    direct_factor_splittings[direct_factor] = splitting
                    break
            # Otherwise, choose any other non-trivial splitting.
            if type(direct_factor_splittings[direct_factor][0]) != str:
                for splitting in splittings:
                    if len(splitting) > 1:
                        direct_factor_splittings[direct_factor] = splitting
                        break
        for direct_factor in direct_factor_splittings:
            # If a direct factor is freely indecomposable, append it as-is.
            if direct_factor_splittings[direct_factor] == (direct_factor,):
                current_splitting.append(direct_factor)
            # Otherwise, append a list of the free factors, where each free factor is either of the form F_n or a graph of groups.
            # If a free factor is a non-trivial single-vertex graph of groups, split the graph braid group a direct product (if possible).
            # Keep alternating between direct product splittings and free product splittings like this.
            else:
                next_free_splitting = []
                for free_factor in direct_factor_splittings[direct_factor]:
                    if type(free_factor) == str:
                        next_free_splitting.append(free_factor)
                    elif len(free_factor.graph.edges) != 0:
                        next_free_splitting.append(free_factor)
                    elif not free_factor.vertex_groups[free_factor.graph.vertices[0]].is_trivial(): 
                        next_direct_splitting = []
                        free_factor.vertex_groups[free_factor.graph.vertices[0]].iterate_splitting(next_direct_splitting)
                        next_free_splitting.append(next_direct_splitting)
                current_splitting.append(next_free_splitting)

    # Returns a list of the non-trivial direct factors of the graph braid group, as graph braid groups.
    def factorise(self):
        non_trivial_factors = []
        for comp in self.graph_components:
            # Checks if each component has particles in it.
            if comp in self.num_initial_particles_per_component:
                # Checks if each component has non-trivial braid group.
                gbg_factor = GraphBraidGroup(self.graph_components[comp], self.num_initial_particles_per_component[comp])
                if not gbg_factor.is_trivial():
                    non_trivial_factors.append(gbg_factor)
        return non_trivial_factors
    
    # Returns True if the non-trivial factors of `self` have the same graphs and same numbers of particles as those of the GraphBraidGroup object `gbg`.
    # Note, if both graph braid groups are reduced then it is sufficient for their graphs to be homeomorphic.
    def is_same(self, gbg):
        non_trivial_factors_1 = self.factorise()
        non_trivial_factors_2 = gbg.factorise()
        has_same_factor = False
        for factor_1 in non_trivial_factors_1:
            for factor_2 in non_trivial_factors_2:
                if factor_1.num_particles == factor_2.num_particles:
                    if factor_1.is_reduced() and factor_2.is_reduced():
                        if is_same(factor_1.graph.make_essential().adj_matrix, factor_2.graph.make_essential().adj_matrix):
                            has_same_factor = True
                            break
                    elif is_same(factor_1.graph.adj_matrix, factor_2.graph.adj_matrix):
                        has_same_factor = True
                        break
            if has_same_factor == False:
                return False
        has_same_factor = False
        for factor_2 in non_trivial_factors_2:
            for factor_1 in non_trivial_factors_1:
                if factor_2.num_particles == factor_1.num_particles:
                    if factor_2.is_reduced() and factor_1.is_reduced():
                        if is_same(factor_2.graph.make_essential().adj_matrix, factor_1.graph.make_essential().adj_matrix):
                            has_same_factor = True
                            break
                    elif is_same(factor_2.graph.adj_matrix, factor_1.graph.adj_matrix):
                        has_same_factor = True
                        break
            if has_same_factor == False:
                return False
        return True


#########################
# Specialised functions #
#########################

# Prompt user to enter data manually if it has not been entered in 'gbg_data.txt'.
def enter_data_manually():
    while True:
        print('No data detected in gbg_data.txt. Enter data manually? (y/n)')
        response = input().lower()
        if response == 'n':
            while True:
                exit_sequence = input('Please enter your data in gbg_data.txt then run the script again. Press [Enter] to exit.')
                if not exit_sequence:
                    sys.exit()
        elif response == 'y':
            # Loop for manual entry of number of particles.
            while True:
                particle_response = input('Please enter the number of particles/strands in the graph braid group:\n')
                if particle_response.isdecimal():
                    if int(particle_response) > 0:
                        num_particles = int(particle_response)
                        break
                    else:
                        print('The number of particles must be a positive integer.')    
                else:
                    print('The number of particles must be a positive integer.')
            # Loop for manual entry of adjacency matrix.
            while True:
                correct_formatting = True
                matrix_response_1 = input('Please enter the first row of the adjacency matrix of the graph, with entries separated by single spaces:\n').split()
                if matrix_response_1 == []:
                    correct_formatting = False
                    print('Incorrect formatting.')
                    continue
                for char in matrix_response_1:
                    if not char.isdecimal():
                        correct_formatting = False
                        print('Incorrect formatting.')
                        break
                    elif int(char) < 0:
                        correct_formatting = False
                        print('Integers cannot be negative.')
                        break
                if correct_formatting:
                    matrix_response = [matrix_response_1]
                    for i in range(2, len(matrix_response_1) + 1):
                        while True:
                            correct_formatting = True
                            current_matrix_response = input(f'Please enter row {i} of the adjacency matrix of the graph:\n').split()
                            for char in current_matrix_response:
                                if not char.isdecimal():
                                    correct_formatting = False
                                    print('Incorrect formatting.')
                                    break
                                elif int(char) < 0:
                                    correct_formatting = False
                                    print('Integers must be positive.')
                                    break
                            if len(current_matrix_response) != len(matrix_response_1):
                                correct_formatting = False
                                print('Incorrect formatting.')
                            else:
                                for j in range(i-1):
                                    if current_matrix_response[j] != matrix_response[j][i-1]:
                                        correct_formatting = False
                                        print('Matrix must be symmetric.')
                            if correct_formatting:
                                matrix_response.append(current_matrix_response)
                                break
                    adj_matrix = [[int(entry) for entry in matrix_response[i]] for i in range(len(matrix_response_1))]
                    break
            # Loop for manual entry of initial configuration.
            if Graph(adj_matrix).get_num_connected_components() == 1:
                initial_config = None
            else:
                while True:
                    correct_formatting = True
                    config_response = input(f'Please enter the initial configuration of the {num_particles} particles on the vertices of the graph. Enter this data in the form of {num_particles} integers separated by single spaces, corresponding to the row numbers of the vertices in the adjacency matrix:\n').split()
                    for char in config_response:
                        if not char.isdecimal():
                            correct_formatting = False
                            print('Incorrect formatting.')
                            break
                        elif int(char) < 1 or int(char) > len(adj_matrix):
                            correct_formatting = False
                            print(f'Integers must be between 1 and {len(adj_matrix)}.')
                            break
                    if len(config_response) != num_particles:
                        correct_formatting = False
                        print('Number of integers must match number of particles.')
                    if correct_formatting:
                        initial_config = [int(entry) - 1 for entry in config_response]
                        break
            break
    return(num_particles, adj_matrix, initial_config)

# Takes as input a list of strings corresponding from the lines of `gbg_data.txt` below the dotted line in the file. 
# Returns the number of particles, adjacency matrix of the graph, and the initial configuration in the correct format.
def get_data_from_file(gbg_data):
    particle_data = gbg_data[0].lstrip('Number of particles:').strip()
    if particle_data.isdecimal():
        if int(particle_data) > 0:
            num_particles = int(particle_data)
        else:
            raise ParticleFormatException
    else: 
        raise ParticleFormatException

    matrix_size = len(gbg_data[2].split())
    for i in range(matrix_size):
        if len(gbg_data[2 + i]) >= 21:
            if gbg_data[2 + i][:21] == 'Initial configuration':
                raise MatrixDimensionException
        for entry in gbg_data[2 + i].split():
            if not entry.isdecimal():
                raise MatrixFormatException
            elif int(entry) < 0:
                raise MatrixFormatException
        if len(gbg_data[2 + i].split()) != matrix_size:
            raise MatrixDimensionException
    adj_matrix = [[int(entry) for entry in gbg_data[2 + i].split()] for i in range(matrix_size)]
    for i in range(matrix_size):
        for j in range(i+1):
            if adj_matrix[i][j] != adj_matrix[j][i]:
                raise MatrixSymmetryException

    config_data = gbg_data[2 + matrix_size].lstrip('Initial configuration:').split()
    if config_data == []:
        if Graph(adj_matrix).get_num_connected_components() == 1:
            initial_config = None
        else:
            raise ConfigDimensionException
    elif len(config_data) != num_particles:
        raise ConfigDimensionException
    else:
        for entry in config_data:
            if not entry.isdecimal():
                raise ConfigFormatException
            elif int(entry) < 1 or int(entry) > len(adj_matrix):
                raise ConfigFormatException

    return (num_particles, adj_matrix, initial_config)

# Replaces all factors in a splitting with string versions and adds detailed data to `splitting.txt` where appropriate.
def stringify_factors(splitting, braid_factor_counter, gog_factor_counter):
    known_gbgs = get_known_gbgs_from_file()
    for i, factor in enumerate(splitting):
        if type(factor) == GraphBraidGroup:
            if factor.is_reduced():
                essential_graph = factor.graph.make_essential()
                essential_adj_matrix_hashable = tuple(tuple(row) for row in essential_graph.adj_matrix)
                if (essential_adj_matrix_hashable, factor.num_particles) in known_gbgs: 
                    splitting[i] = known_gbgs[(essential_adj_matrix_hashable, factor.num_particles)]
                else:
                    file = open('splitting.txt', 'r')
                    splitting_data = file.readlines()
                    file.close()
                    previous_matrices = [ast.literal_eval(line.partition('adjacency matrix: ')[2].strip()) for line in splitting_data if line.startswith('Ga')]
                    is_new_graph = True
                    for previous_matrix in previous_matrices:
                        if is_same(essential_graph.adj_matrix, previous_matrix):
                            splitting[i] = f'B_{factor.num_particles}(Gamma_{previous_matrices.index(previous_matrix) + 1})'
                            is_new_graph = False
                            break
                    if is_new_graph:
                        splitting[i] = f'B_{factor.num_particles}(Gamma_{braid_factor_counter})'
                        file = open('splitting.txt', 'a')
                        file.write(f'Gamma_{braid_factor_counter} adjacency matrix: {essential_graph.adj_matrix} \n')
                        file.close()
                        braid_factor_counter += 1
            else:
                file = open('splitting.txt', 'r')
                splitting_data = file.readlines()
                file.close()
                previous_matrices = [ast.literal_eval(line.partition('adjacency matrix: ')[2].strip()) for line in splitting_data if line.startswith('Ga')]
                is_new_graph = True
                for previous_matrix in previous_matrices:
                    if is_same(factor.adj_matrix, previous_matrix):
                        splitting[i] = f'RB_{factor.num_particles}(Gamma_{previous_matrices.index(previous_matrix) + 1})'
                        is_new_graph = False
                        break
                if is_new_graph:
                    splitting[i] = f'RB_{factor.num_particles}(Gamma_{braid_factor_counter})'
                    file = open('splitting.txt', 'a')
                    file.write(f'Gamma_{braid_factor_counter} adjacency matrix: {factor.adj_matrix} \n')
                    file.close()
                    braid_factor_counter += 1
        elif type(factor) == GraphOfGroups:
            splitting[i] = f'G_{gog_factor_counter}'
            file = open('splitting.txt', 'a')
            file.write(f'G_{gog_factor_counter} adjacency matrix: {factor.graph.adj_matrix} \n') 
            for v in factor.vertex_groups:
                if factor.vertex_groups[v].is_reduced():
                    file.write(f'Vertex group {v+1}: B_{factor.vertex_groups[v].num_particles}, adjacency matrix: {factor.vertex_groups[v].graph.adj_matrix}, initial configuration: {[i+1 for i in factor.vertex_groups[v].initial_config]} \n')
                else:
                    file.write(f'Vertex group {v+1}: RB_{factor.vertex_groups[v].num_particles}, adjacency matrix: {factor.vertex_groups[v].graph.adj_matrix}, initial configuration: {[i+1 for i in factor.vertex_groups[v].initial_config]} \n')
            for e in factor.edge_groups:
                if factor.edge_groups[e].is_reduced():
                    file.write(f'Edge group {(e[0]+1, e[1]+1)}: B_{factor.edge_groups[e].num_particles}, adjacency matrix: {factor.edge_groups[e].graph.adj_matrix}, initial configuration: {[i+1 for i in factor.edge_groups[e].initial_config]} \n')
                else:
                    file.write(f'Edge group {(e[0]+1, e[1]+1)}: RB_{factor.edge_groups[e].num_particles}, adjacency matrix: {factor.edge_groups[e].graph.adj_matrix}, initial configuration: {[i+1 for i in factor.edge_groups[e].initial_config]} \n')
            file.close()  
            gog_factor_counter += 1      
        elif type(factor) == list:
            stringify_factors(factor, braid_factor_counter, gog_factor_counter)

# Gets list of known graph braid groups from `known_gbgs.txt` and puts them in a dictionary.
def get_known_gbgs_from_file():
    file = open('known_gbgs.txt')
    known_gbg_data = [line.strip() for line in file.readlines()]
    file.close()
    # Dictionary where keys are (adj_matrix, num_particles) tuples and values are the braid groups as strings.
    known_gbgs = {(ast.literal_eval(known_gbg_data[3*i + 1]), int(known_gbg_data[3*i + 2])): known_gbg_data[3*i] for i in range(int(len(known_gbg_data) / 3))}
    # Adds all possible permutations of the vertex labels of the graph to the dictionary.
    known_gbgs_unordered = {braid_data: known_gbgs[braid_data] for braid_data in known_gbgs}
    for (matrix, particles) in known_gbgs:
        for perm in itertools.permutations(range(len(matrix)), len(matrix)):
            permuted_matrix = tuple(tuple(matrix[perm[i]][perm[j]] for i in range(len(matrix))) for j in range(len(matrix)))
            known_gbgs_unordered.setdefault((permuted_matrix, particles), known_gbgs[(matrix, particles)])
    return known_gbgs_unordered

# Takes as input two adjacency matrices and returns True if they define the same graph.
# Tip: To check if two Graph objects are homeomorphic, first use `make_essential` on both and then use `is_same` on their adjacency matrices.
def is_same(adj_matrix_1, adj_matrix_2):
    matrix_2_permutations = [[[adj_matrix_2[perm[i]][perm[j]] for i in range(len(adj_matrix_2))] for j in range(len(adj_matrix_2))] for perm in itertools.permutations(range(len(adj_matrix_2)), len(adj_matrix_2))]
    if adj_matrix_1 in matrix_2_permutations:
        return True
    else:
        return False

# Converts `stringified_splitting` into a list of strings, where each string represents the splitting of a direct factor of the original group. 
# These strings have been converted from list form into single strings with 'x' and '*' characters as appropriate.
def combine_strings(stringified_splitting, is_direct = True):
    for i, factor in enumerate(stringified_splitting):
        # If a factor is a list of non-lists, then replace the factor in `stringified_splitting` with its string representation.
        if type(factor) == list:
            for subfactor in factor:
                # If we find a list in factor's subfactors, then feed factor back into `combine_strings`.
                if type(subfactor) == list:
                    combine_strings(factor, not is_direct)
                    break
            # If `stringified splitting` is a direct splitting, then `factor` is a free splitting.
            if is_direct:
                # Collect together F_m terms and Z terms.
                free_rank = sum(int(subfactor[2:]) for subfactor in factor if subfactor.startswith('F_')) + len([subfactor for subfactor in factor if subfactor == 'Z'])
                factor = [f'F_{k}' for k in [free_rank] if free_rank > 1] + ['Z' for k in [free_rank] if free_rank == 1] + [subfactor for subfactor in factor if not subfactor.startswith('F_') and not subfactor == 'Z']
                # Collect together like terms and express in simplified notation.
                factor_count = defaultdict(int)
                for subfactor in factor:
                    factor_count[subfactor] += 1
                factor = [f'*^{factor_count[subfactor]}(' + subfactor + ')' for subfactor in factor_count if factor_count[subfactor] > 1] + [subfactor for subfactor in factor_count if factor_count[subfactor] == 1]
                if len(factor) > 1:
                    stringified_splitting[i] = '(' + ' * '.join(factor) + ')'
                else:
                    stringified_splitting[i] = factor[0]
            else:
                # Collect together like terms and express as powers.
                factor_count = defaultdict(int)
                for subfactor in factor:
                    factor_count[subfactor] += 1
                factor = [subfactor + f'^{factor_count[subfactor]}' for subfactor in factor_count if factor_count[subfactor] > 1 and len(subfactor) == 1] + ['(' + subfactor + ')' + f'^{factor_count[subfactor]}' for subfactor in factor_count if factor_count[subfactor] > 1 and len(subfactor) > 1] + [subfactor for subfactor in factor_count if factor_count[subfactor] == 1]
                if len(factor) > 1:
                    stringified_splitting[i] = '(' + ' x '.join(factor) + ')'
                else:
                    stringified_splitting[i] = factor[0]


#############
# Main code #
#############

def main():
    print('Checking gbg_data.txt...')
    file = open('gbg_data.txt')
    file_as_list = file.readlines()
    file.close()
    for i, line in enumerate(file_as_list):
        file_as_list[i] = line.strip()
    non_empty_lines = [line for line in file_as_list if line != '']
    dotted_line = non_empty_lines.index('-' * 50)
    gbg_data = non_empty_lines[dotted_line + 1:]

    # If data has not been entered in `gbg_data.txt`, prompt user to enter it manually.
    if len(gbg_data) == 3:
        (num_particles, adj_matrix, initial_config) = enter_data_manually()
    # Otherwise, get data from `gbg_data.txt` and verify it is formatted correctly.
    else: 
        (num_particles, adj_matrix, initial_config) = get_data_from_file(gbg_data)

    gbg = GraphBraidGroup(Graph(adj_matrix), num_particles, initial_config)
    if not gbg.is_reduced():
        print('WARNING: In order to perform computations for B_n(\Gamma), the graph \Gamma must satisfy the following conditions:')
        print('1. All cycles must have length at least n+1.')
        print('2. All paths between vertices of degree not equal to 2 must have length at least n-1.')
        print('At least one of these conditions is not satisfied by your graph.')
        print('If you choose to continue, any results obtained will only be true for the reduced braid group RB_n(\Gamma).')
        print('Do you wish to continue? (y/n)')
        while True:
            response = input().lower()
            if response == 'n':
                while True:
                    exit_sequence = input('Please amend your data then run the script again. Press [Enter] to exit.')
                    if not exit_sequence:
                        sys.exit()
            elif response == 'y':
                for comp in gbg.num_initial_particles_per_component:
                    if len(comp[0]) < gbg.num_initial_particles_per_component[comp]:
                        raise VertexException
                break

    print('Verified successfully.')
    print('Searching for free splittings...')

    if gbg.is_trivial():
        if gbg.is_reduced():
            print(f'B_{num_particles} = 1')
            start_exit_sequence()
        else:
            print(f'RB_{num_particles} = 1')
            start_exit_sequence()

    # Returns a nested list, where the odd levels of the list correspond to direct factors and the even levels correspond to free factors.
    splitting = gbg.get_splitting()

    print('Search complete.')

    if len(splitting) == 1:
        if type(splitting[0]) != str:
            if type(splitting[0]) != list:
                print('No splittings found.')
                start_exit_sequence()
            elif len(splitting[0]) == 1:
                print('No splittings found.')
                start_exit_sequence()

    print('Splitting found. Converting data to readable format...')

    # Makes a fresh copy of `splitting.txt`, the file containing detailed splitting information (for the graph of group factors and the graph braid group factors).
    file = open('splitting.txt', 'w')
    file.write('')
    file.close()

    # Turns all factors in the splitting into strings and adds detailed data to `splitting.txt` where appropriate.
    stringify_factors(splitting, 1, 1)

    combine_strings(splitting)
    final_string = ' x '.join(splitting)

    # Prints the splitting and adds it to the beginning of `splitting.txt`
    with open('splitting.txt','r') as contents:
        save = contents.read()
    if gbg.is_reduced():
        with open('splitting.txt','w') as contents:
            contents.write(f'B_{num_particles} = ' + final_string + '\n\n')
        with open('splitting.txt','a') as contents:
            contents.write(save)
        print(f'B_{num_particles} = ' + final_string)
    else:
        with open('splitting.txt','w') as contents:
            contents.write(f'RB_{num_particles} = ' + final_string + '\n\n')
        with open('splitting.txt','a') as contents:
            contents.write(save)
        print(f'RB_{num_particles} = ' + final_string)

    if 'G_' in final_string or 'Gamma_' in final_string:
        print('In the above splitting, Gamma_i are graphs and G_i are fundamental groups of graphs of groups. More detailed data can be found in splitting.txt.')
    start_exit_sequence()

if __name__ == '__main__':
    main()