#!/usr/bin/env python3
# graph_braid_group.py - GraphBraidGroup class file.

import math
import itertools
from collections import defaultdict
from graph import Graph
from graph_of_groups import GraphOfGroups

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

# Takes as input two adjacency matrices and returns True if they define the same graph.
# Tip: To check if two Graph objects are homeomorphic, first use `make_essential` on both and then use `is_same` on their adjacency matrices.
def is_same(adj_matrix_1, adj_matrix_2):
    for perm in itertools.permutations(range(len(adj_matrix_2)), len(adj_matrix_2)):
        if adj_matrix_1 == [[adj_matrix_2[perm[i]][perm[j]] for i in range(len(adj_matrix_2))] for j in range(len(adj_matrix_2))]:
            return True
    return False

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
        if initial_config is None:
            self.initial_config = [self.graph.vertices[i] for i in range(num_particles)]
        else:
            self.initial_config = initial_config
        # Dictionary whose keys are the connected components of the graph containing vertices in `initial_config`, as (vertex_set, edge_set)-tuples,
        # and whose values are the number of particles that start in that component.
        self.num_initial_particles_per_component = self.get_num_particles_per_component(self.initial_config)
        # Uses Corollary 2.16 of [Graph of groups decompositions of graph braid groups](https://arxiv.org/pdf/2209.03860.pdf).
        self.num_connected_components = math.comb(self.num_particles + self.graph.get_num_connected_components() - 1, self.graph.get_num_connected_components() - 1)

    def __eq__(self, other):
        return self.is_same(other)

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
    
    # Note: `num_particle_per_component` does not list components that contain zero particles. This is an important feature; see e.g. `is_trivial()`.
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
            if len(initial_subcomponents) == 0:
                new_partitions = [terminal_partition for terminal_partition in new_partitions_terminal]
            elif len(terminal_subcomponents) == 0:
                new_partitions = [initial_partition for initial_partition in new_partitions_initial]
            else:
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
                if len(splitting) == 1 and isinstance(splitting[0], str):
                    direct_factor_splittings[direct_factor] = splitting
                    break
            # Otherwise, choose any other non-trivial splitting.
            if not isinstance(direct_factor_splittings[direct_factor][0], str):
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
                    if isinstance(free_factor, str):
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
            if not has_same_factor:
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
            if not has_same_factor:
                return False
        return True
