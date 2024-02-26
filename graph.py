#!/usr/bin/env python3
# graph.py - Graph class file.
"""
A module for graph functionality.

Module contents: 
- The `Graph` class, which implements methods for performing graph computations.
"""

import copy
from collections import defaultdict
from types import NotImplementedType

class Graph:
    """Class for performing graph computations."""

    def __init__(self, adj_matrix: list[list[int]] | tuple[tuple[int, ...], ...]) -> None:
        """
        Initialises graph attributes.

        `adj_matrix` is the adjacency matrix; it should be a tuple of tuples or a list of lists. 
        
        The vertices of the graph are numbered according to their row/column in the matrix and encoded in a list.
        The edges of the graph are encoded as a list of 2-tuples (i, j) where i <= j.
        """
        self.adj_matrix = adj_matrix
        if self.adj_matrix == [[]] or []:
            self.vertices = []
        else:
            self.vertices = [i for i in range(len(adj_matrix))]
        self.num_vertices = len(self.vertices)
        if self.num_vertices == 0:
            self.edges = []
        else:
            self.edges = [(i, j) for i in range(self.num_vertices) for j in range(i, self.num_vertices) for n in range(adj_matrix[i][j])]
        self.num_edges = len(self.edges)
        self.essential_vertices = [v for v in self.vertices if self.get_degree(v) != 2]

    def __eq__(self, other: object) -> bool | NotImplementedType:
        """Two graphs are equal if they have the same adjacency matrix or they both have trivial adjacency matrix."""
        trivial = ([], [[]])
        if not isinstance(other, Graph):
            return NotImplemented
        return self.adj_matrix == other.adj_matrix or (self.adj_matrix in trivial and other.adj_matrix in trivial)
    
    def __hash__(self) -> int:
        """Graphs are hashed via their adjacency matrix."""
        trivial = ([], [[]])
        if self.adj_matrix in trivial:
            return hash(tuple())
        else:
            return hash(tuple(tuple(row) for row in self.adj_matrix))

    def get_connected_components(self, subgraph: tuple[list[int], list[tuple[int, int]]] | tuple[tuple[int, ...], tuple[tuple[int, int], ...]] = None) -> dict[tuple[tuple[int, ...], tuple[tuple[int, int], ...]], 'Graph']:
        """
        Returns dictionary of connected components of the (sub)graph.

        If a subgraph is specified, returns the connected components of the subgraph.
        `subgraph` should be a (vertex set, edge set) 2-tuple of sublists of `self.vertices` and `self.edges`.
        Vertices are represented by integers corresponding to the row number in the adjacency matrix, and edges are 2-tuples of integers.
        If left blank, `subgraph` defaults to the entire graph.
        
        In the returned dictionary:
        - the keys are 2-tuples of vertex sets and edge sets of the connected components of the graph, considered as subgraphs;
        - the values are the components as standalone Graph objects.
        """
        if subgraph is None:
            subgraph = (self.vertices, self.edges)
        if subgraph[0] == []:
            return {} # The empty subgraph has no connected components.
        initial_component = self.get_component(subgraph[0][0], subgraph) # Otherwise, start with the component containing first vertex of `subgraph`.
        components = {initial_component[0]: initial_component[1]}
        for v in subgraph[0]:
            if v not in {w for component in components for w in component[0]}:
                new_component = self.get_component(v, subgraph)
                components.update({new_component[0]: new_component[1]})
        return components

    def get_component(self, vertex: int, subgraph: tuple[list[int], list[tuple[int, int]]] | tuple[tuple[int, ...], tuple[tuple[int, int], ...]] = None) -> tuple[tuple[tuple[int, ...], tuple[tuple[int, int], ...]], 'Graph']:
        """
        Returns a 2-tuple representing the connected component of the (sub)graph containing the given vertex.
        
        `subgraph` should be a (vertex set, edge set) 2-tuple of sublists of `self.vertices` and `self.edges`.
        Vertices are represented by integers corresponding to the row number in the adjacency matrix, and edges are 2-tuples of integers.
        If left blank, `subgraph` defaults to the entire graph.
        
        The first entry in the returned 2-tuple is the (vertex set, edge set) 2-tuple of the connected component of `subgraph` containing `vertex`.
        The second entry in the returned 2-tuple is the connected component as a standalone `Graph` object.
        """
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
            # Component vertices must be referred to by their index in the list when constructing the adjacency matrix
            v0, v1 = component_vertices.index(edge[0]), component_vertices.index(edge[1])
            if edge[0] != edge[1]:
                component_adj_matrix[v0][v1] += 1
                component_adj_matrix[v1][v0] += 1
            else:
                component_adj_matrix[v0][v1] += 1
        return ((tuple(component_vertices), tuple(component_edges)), Graph(component_adj_matrix))
    
    def iterate_component(self, component_vertices: list[int], component_edges: list[tuple[int, int]], current_vertex: int, subgraph: tuple[list[int], list[tuple[int, int]]] | tuple[tuple[int, ...], tuple[tuple[int, int], ...]] = None) -> None:
        if subgraph is None:
            subgraph = (self.vertices, self.edges)
        for edge in subgraph[1]:
            if current_vertex in edge and edge not in component_edges:
                for e in [e for e in subgraph[1] if e == edge]:
                    component_edges.append(e)
                adjacent_vertex = edge[edge.index(current_vertex) - 1]
                if adjacent_vertex not in component_vertices:    
                    component_vertices.append(adjacent_vertex)
                    self.iterate_component(component_vertices, component_edges, adjacent_vertex, subgraph)

    def get_num_connected_components(self) -> int:
        """Returns the number of connected components of the graph."""
        return len(self.get_connected_components())

    def get_graph_minus_open_edges(self, removed_edges: list[tuple[int, int]]) -> tuple[tuple[list[int], list[tuple[int, int]]], 'Graph']:
        """
        Returns a 2-tuple representing the graph minus the given edges, considered as open edges.

        `removed_edges` should be a list of edges in `self.edges` to be removed from the graph.
        Edges are 2-tuples of integers corresponding to the row number and column number in the adjacency matrix.

        In the returned 2-tuple:
        - the first entry is the 2-tuple of vertices and edges of the graph minus `removed_edges`, considered as a subgraph of the original graph;
        - the second entry is the graph minus `removed_edges` as an object of class `Graph`.
        """
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

    def get_graph_minus_closed_edge(self, removed_edge: tuple[int, int]) -> tuple[tuple[list[int], list[tuple[int, int]]], 'Graph']:
        """
        Returns a 2-tuple representing the graph minus the given edge, considered as a closed edge.

        `removed_edge` should be an edge in `self.edges` to be removed from the graph.
        Edges are 2-tuples of integers corresponding to the row number and column number in the adjacency matrix.

        Note, this is equivalent to removing this edge's two vertices together with all edges containing one of these vertices.
        
        In the returned 2-tuple:
        - the first entry is the 2-tuple of vertices and edges of the graph minus `removed_edge`, considered as a subgraph of the original graph;
        - the second entry is the graph minus `removed_edge` as an object of class `Graph`.
        """
        (i, j) = removed_edge
        subgraph_vertices = [v for v in self.vertices if v not in (i, j)]
        subgraph_edges = [e for e in self.edges if i not in e and j not in e]
        # Removing edge (i, j) corresponds to removing the ith and jth rows and columns from the adjacency matrix
        adj_matrix_minus_closed_edge = [[self.adj_matrix[k][m] for m in range(len(self.adj_matrix[k])) if m not in (i, j)] for k in range(len(self.adj_matrix)) if k not in (i, j)]
        return ((subgraph_vertices, subgraph_edges), Graph(adj_matrix_minus_closed_edge))


    def get_graph_minus_vertex(self, removed_vertex: int) -> tuple[tuple[list[int], list[tuple[int, int]]], 'Graph']:
        """
        Returns a 2-tuple representing the graph minus the given vertex.

        `removed_vertex` should be a vertex in `self.vertices` to be removed from the graph.
        Vertices are integers corresponding to the row number in the adjacency matrix.

        Note, this is equivalent to removing the vertex together with all edges containing it.

        In the returned 2-tuple:
        - the first entry is the 2-tuple of vertices and edges of the graph minus `removed_vertex`, considered as a subgraph of the original graph;
        - the second entry is the graph minus `removed_vertex` as an object of class `Graph`.
        """
        subgraph_vertices = [v for v in self.vertices if v != removed_vertex]
        subgraph_edges = [e for e in self.edges if removed_vertex not in e]
        # Removing the vertex corresponds to removing its row and column from the adjacency matrix
        adj_matrix_minus_vertex = [[self.adj_matrix[k][m] for m in range(len(self.adj_matrix[k])) if m != removed_vertex] for k in range(len(self.adj_matrix)) if k != removed_vertex]
        return ((subgraph_vertices, subgraph_edges), Graph(adj_matrix_minus_vertex))

    def prune(self, edges: list[tuple[int, int]]) -> tuple['Graph', list[tuple[int, int]]]:
        """
        Iteratively checks whether the given edges are non-separating and removes them if so, returning the pruned graph and a list of removed edges.
        
        Goes through the list `edges` one by one and checks if the edge is non-separating. 
        If so, removes that (open) edge and repeats on the new graph.

        `edges` should be a list of candidate edges in `self.edges` for pruning.
        Edges are 2-tuples of integers corresponding to the row number and column number in the adjacency matrix.

        Returns a 2-tuple consisting of: 
        - a `Graph` object representing the pruned graph;
        - a list of the removed edges.
        """
        modified_graph = self
        removed_edges = []
        num_components = self.get_num_connected_components()
        for edge in edges:
            modified_graph_minus_edge = modified_graph.get_graph_minus_open_edges([edge])[1]
            if modified_graph_minus_edge.get_num_connected_components() == num_components:
                removed_edges.append(edge)
                modified_graph = modified_graph_minus_edge
        return (modified_graph, removed_edges)

    def is_separating(self, edge: tuple[int, int]) -> bool:
        """
        Returns True if the given edge separates the graph into more connected components.
        
        `edge` should be an edge in `self.edges`.
        Edges are 2-tuples of integers corresponding to the row number and column number in the adjacency matrix.
        """
        if self.get_graph_minus_open_edges([edge])[1].get_num_connected_components() > self.get_num_connected_components():
            return True
        else:
            return False
    
    def get_degree(self, vertex: int) -> int:
        """
        Returns the degree of the given vertex in the graph.
        
        `vertex` should be a vertex in `self.vertices`.
        Vertices are integers corresponding to the row number in the adjacency matrix.
        """
        degree = 0
        for edge in self.edges:
            if edge == (vertex, vertex):
                degree += 2
            elif vertex in edge:
                degree += 1
        return degree
    
    def get_centreless_ball(self, vertex: int, distance: int) -> list[int]:
        """
        Returns a list of all vertices that can be reached from the given vertex via a path of length at least 1 and at most the given distance.

        `vertex` should be a vertex in `self.vertices`.
        Vertices are integers corresponding to the row number in the adjacency matrix.
        `distance` should be an integer representing the radius of the ball.

        The ball is centreless, meaning if `vertex` appears in the list, this means it is contained in a cycle.
        """
        ball = []
        remaining_edges = defaultdict(int)
        for edge in self.edges:
            remaining_edges[edge] += 1
        if distance > 0:
            self.iterate_ball(ball, remaining_edges, vertex, distance)
        return ball
    
    def iterate_ball(self, ball: list[int], remaining_edges: dict[tuple[int, int], int], current_vertex: int, remaining_distance: int) -> None:
        for edge in self.edges:
            if current_vertex in edge and remaining_edges[edge] > 0:
                adjacent_vertex = edge[edge.index(current_vertex) - 1]
                if adjacent_vertex not in ball:
                    ball.append(adjacent_vertex)
                    remaining_edges[edge] -= 1
                    if remaining_distance >= 2:
                        self.iterate_ball(ball, remaining_edges, adjacent_vertex, remaining_distance - 1)

    def make_essential(self) -> 'Graph':
        """
        Removes vertices of degree 2 from the graph until only essential vertices remain, returning the essential graph.
        
        Returns the essential graph as a `Graph` object.
        """
        non_essential_vertices = [v for v in self.vertices if v not in self.essential_vertices]
        updated_essential_vertices = self.essential_vertices
        # If all vertices have degree 2, then our graph consists of cycles and we must therefore designate one degree 2 vertex as essential in each cycle.
        if self.essential_vertices == []:
            for component in self.get_connected_components():
                non_essential_vertices.remove(component[0][0])
            updated_essential_vertices = [v for v in self.vertices if v not in non_essential_vertices]
        new_edges = []
        removed_edges = []
        # For each non-essential vertex, remove its adjacent edges and connect the two vertices on either side.
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
            if i == j:
                modified_adj_matrix[i][j] += 1
            else:
                modified_adj_matrix[i][j] += 1
                modified_adj_matrix[j][i] += 1
        essential_adj_matrix = [[modified_adj_matrix[i][j] for j in updated_essential_vertices] for i in updated_essential_vertices]
        return Graph(essential_adj_matrix)