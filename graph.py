#!/usr/bin/env python3
# graph.py - Graph class file.

import copy

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

    def __eq__(self, other):
        return self.adj_matrix == other.adj_matrix

    # Returns dictionary of connected components of the graph, where:
    # the keys are 2-tuples of vertex sets and edge sets of the connected components of the graph, considered as subgraphs;
    # the values are the components as standalone 'Graph' objects.
    def get_connected_components(self, subgraph = None):
        if subgraph is None:
            subgraph = (self.vertices, self.edges)
        # The empty subgraph has no connected components.
        if subgraph[0] == []:
            return {}
        # Otherwise, start with the component containing first vertex of `subgraph`
        components = {self.get_component(subgraph[0][0], subgraph)[0]: self.get_component(subgraph[0][0], subgraph)[1]}
        for v in subgraph[0]:
            if v not in {w for component in components for w in component[0]}:
                components.update({self.get_component(v, subgraph)[0]: self.get_component(v, subgraph)[1]})
        return components

    # Returns a 2-tuple where the first entry is the 2-tuple (vertex set, edge set) of the connected component of `subgraph` containing `vertex` 
    # (`subgraph` defaults to the original graph) and the second entry is the connected component as a standalone Graph object.
    # Note, `subgraph` should be a (vertex set, edge set) 2-tuple of sublists of `self.vertices` and `self.edges`
    def get_component(self, vertex, subgraph = None):
        if subgraph is None:
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
        if subgraph is None:
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
        adj_matrix_minus_closed_edge = [[self.adj_matrix[k][m] for m in range(len(self.adj_matrix[k])) if m not in (i, j)] for k in range(len(self.adj_matrix)) if k not in (i, j)]
        return ((subgraph_vertices, subgraph_edges), Graph(adj_matrix_minus_closed_edge))
    
    # We remove a single vertex from our graph, together with all edges containing the vertex.
    def get_graph_minus_vertex(self, removed_vertex):
        subgraph_vertices = [v for v in self.vertices if v != removed_vertex]
        subgraph_edges = [e for e in self.edges if removed_vertex not in e]
        # Removing the vertex corresponds to removing its row and column from the adjacency matrix
        adj_matrix_minus_vertex = [[self.adj_matrix[k][m] for m in range(len(self.adj_matrix[k])) if m != removed_vertex] for k in range(len(self.adj_matrix)) if k != removed_vertex]
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
        # If all vertices have degree 2, then our graph consists of cycles and we must therefore designate one degree 2 vertex as essential in each cycle.
        if self.essential_vertices == []:
            for component in self.get_connected_components():
                non_essential_vertices.remove(component[0][0])
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