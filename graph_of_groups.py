#!/usr/bin/env python3
# graph_of_groups.py - GraphOfGroups class file.

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