from graph import Graph

# TO DO: Create a class for each type of graph to be used in testing:
# - 100 cycles each of length 100

class TestGraphVertex:
    adj_matrix = [[0]]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {(tuple(self.test_graph.vertices), tuple(self.test_graph.edges)) : self.test_graph}
    def test_get_component(self):
        assert self.test_graph.get_component(0) == ((tuple(self.test_graph.vertices), tuple(self.test_graph.edges)), self.test_graph)
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 1
    def test_get_graph_minus_vertex(self):
        assert self.test_graph.get_graph_minus_vertex(0) == (([], []), Graph([[]]))
    def test_get_degree(self):
        assert self.test_graph.get_degree(0) == 0
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(0, 1) == []
    def test_make_essential(self):
        assert self.test_graph.make_essential() == self.test_graph

class TestGraphEdge:
    adj_matrix = [[0,1],[1,0]]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {(tuple(self.test_graph.vertices), tuple(self.test_graph.edges)) : self.test_graph}
    def test_get_component(self):
        assert self.test_graph.get_component(0) == ((tuple(self.test_graph.vertices), tuple(self.test_graph.edges)), self.test_graph)
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 1
    def test_get_graph_minus_open_edges(self):
        assert self.test_graph.get_graph_minus_open_edges([(0,1)]) == (([0,1], []), Graph([[0,0],[0,0]]))
    def test_get_graph_minus_closed_edge(self):
        assert self.test_graph.get_graph_minus_closed_edge((0,1)) == (([], []), Graph([[]]))
    def test_get_graph_minus_vertex(self):
        assert self.test_graph.get_graph_minus_vertex(0) == (([1], []), Graph([[0]]))
    def test_prune(self):
        assert self.test_graph.prune([(0,1)]) == (Graph([[0,1],[1,0]]), [])
    def test_is_separating(self):
        assert self.test_graph.is_separating((0,1))
    def test_get_degree(self):
        assert self.test_graph.get_degree(0) == 1
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(0, 1) == [1]
    def test_make_essential(self):
        assert self.test_graph.make_essential() == self.test_graph

class TestGraphLoop:
    adj_matrix = [[1]]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {(tuple(self.test_graph.vertices), tuple(self.test_graph.edges)) : self.test_graph}
    def test_get_component(self):
        assert self.test_graph.get_component(0) == ((tuple(self.test_graph.vertices), tuple(self.test_graph.edges)), self.test_graph)
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 1
    def test_get_graph_minus_open_edges(self):
        assert self.test_graph.get_graph_minus_open_edges([(0,0)]) == (([0], []), Graph([[0]]))
    def test_get_graph_minus_closed_edge(self):
        assert self.test_graph.get_graph_minus_closed_edge((0,0)) == (([], []), Graph([[]]))
    def test_get_graph_minus_vertex(self):
        assert self.test_graph.get_graph_minus_vertex(0) == (([], []), Graph([[]]))
    def test_prune(self):
        assert self.test_graph.prune([(0,0)]) == (Graph([[0]]), [(0,0)])
    def test_is_separating(self):
        assert not self.test_graph.is_separating((0,0))
    def test_get_degree(self):
        assert self.test_graph.get_degree(0) == 2
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(0, 2) == [0]
    def test_make_essential(self):
        assert self.test_graph.make_essential() == self.test_graph

class TestGraphSegment:
    # Segment of length 3.
    adj_matrix = [[0,1,0,0],[1,0,1,0],[0,1,0,1],[0,0,1,0]]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {(tuple(self.test_graph.vertices), tuple(self.test_graph.edges)) : self.test_graph}
    def test_get_component(self):
        assert self.test_graph.get_component(0) == ((tuple(self.test_graph.vertices), tuple(self.test_graph.edges)), self.test_graph)
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 1
    def test_get_graph_minus_open_edges(self):
        assert self.test_graph.get_graph_minus_open_edges([(0,1), (2,3)]) == (([0,1,2,3], [(1,2)]), Graph([[0,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,0]]))
    def test_get_graph_minus_closed_edge(self):
        assert self.test_graph.get_graph_minus_closed_edge((1,2)) == (([0,3], []), Graph([[0,0],[0,0]]))
    def test_get_graph_minus_vertex(self):
        assert self.test_graph.get_graph_minus_vertex(0) == (([1,2,3], [(1,2), (2,3)]), Graph([[0,1,0],[1,0,1],[0,1,0]]))
    def test_prune(self):
        assert self.test_graph.prune([(0,1), (1,2)]) == (self.test_graph, [])
    def test_is_separating(self):
        assert self.test_graph.is_separating((1,2))
    def test_get_degree(self):
        assert self.test_graph.get_degree(1) == 2
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(1, 2) == [0,2,3]
    def test_make_essential(self):
        assert self.test_graph.make_essential() == Graph([[0,1],[1,0]])

class TestGraph3Vertex:
    # Three disjoint vertices.
    adj_matrix = [[0,0,0],[0,0,0],[0,0,0]]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {((0,), tuple()) : Graph([[0]]), ((1,), tuple()) : Graph([[0]]), ((2,), tuple()) : Graph([[0]])}
    def test_get_component(self):
        assert self.test_graph.get_component(0) == (((0,), tuple()), Graph([[0]]))
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 3
    def test_get_graph_minus_vertex(self):
        assert self.test_graph.get_graph_minus_vertex(0) == (([1,2], []), Graph([[0,0],[0,0]]))
    def test_get_degree(self):
        assert self.test_graph.get_degree(2) == 0
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(1, 2) == []
    def test_make_essential(self):
        assert self.test_graph.make_essential() == self.test_graph

class TestGraph3Segment:
    # Three disjoint segments of length 1, 2, 3, respectively.
    adj_matrix = [[0,1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,1,0,1,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0,0],[0,0,0,0,0,1,0,1,0],[0,0,0,0,0,0,1,0,1],[0,0,0,0,0,0,0,1,0]]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {((0,1), ((0,1),)) : Graph([[0,1],[1,0]]), ((2,3,4), ((2,3),(3,4))) : Graph([[0,1,0],[1,0,1],[0,1,0]]), ((5,6,7,8), ((5,6),(6,7),(7,8))) : Graph([[0,1,0,0],[1,0,1,0],[0,1,0,1],[0,0,1,0]])}
    def test_get_component(self):
        assert self.test_graph.get_component(2) == (((2,3,4), ((2,3),(3,4))), Graph([[0,1,0],[1,0,1],[0,1,0]]))
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 3
    def test_get_graph_minus_open_edges(self):
        assert self.test_graph.get_graph_minus_open_edges([(6,7)]) == (([0,1,2,3,4,5,6,7,8], [(0,1),(2,3),(3,4),(5,6),(7,8)]), Graph([[0,1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,1,0,1,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0,0],[0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,1,0]]))
    def test_get_graph_minus_closed_edge(self):
        assert self.test_graph.get_graph_minus_closed_edge((2,3)) == (([0,1,4,5,6,7,8], [(0,1),(5,6),(6,7),(7,8)]), Graph([[0,1,0,0,0,0,0],[1,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,1,0,0],[0,0,0,1,0,1,0],[0,0,0,0,1,0,1],[0,0,0,0,0,1,0]]))
    def test_get_graph_minus_vertex(self):
        assert self.test_graph.get_graph_minus_vertex(6) == (([0,1,2,3,4,5,7,8], [(0,1),(2,3),(3,4),(7,8)]), Graph([[0,1,0,0,0,0,0,0],[1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,1,0,1,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0]]))
    def test_prune(self):
        assert self.test_graph.prune([(0,1), (2,3)]) == (self.test_graph, [])
    def test_is_separating(self):
        assert self.test_graph.is_separating((2,3))
    def test_get_degree(self):
        assert self.test_graph.get_degree(8) == 1
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(3, 1) == [2,4]
    def test_make_essential(self):
        assert self.test_graph.make_essential() == Graph([[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,1,0,0],[0,0,1,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]])

class TestGraphStar:
    # Three-pronged star graph
    adj_matrix = [[0,1,1,1],[1,0,0,0],[1,0,0,0],[1,0,0,0]]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {(tuple(self.test_graph.vertices), tuple(self.test_graph.edges)) : self.test_graph}
    def test_get_component(self):
        assert self.test_graph.get_component(2) == ((tuple(self.test_graph.vertices), tuple(self.test_graph.edges)), self.test_graph)
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 1
    def test_get_graph_minus_open_edges(self):
        assert self.test_graph.get_graph_minus_open_edges([(0,1), (0,2), (0,3)]) == (([0,1,2,3], []), Graph([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]))
    def test_get_graph_minus_closed_edge(self):
        assert self.test_graph.get_graph_minus_closed_edge((0,3)) == (([1,2], []), Graph([[0,0],[0,0]]))
    def test_get_graph_minus_vertex(self):
        assert self.test_graph.get_graph_minus_vertex(0) == (([1,2,3], []), Graph([[0,0,0],[0,0,0],[0,0,0]]))
    def test_prune(self):
        assert self.test_graph.prune([(0,1)]) == (self.test_graph, [])
    def test_is_separating(self):
        assert self.test_graph.is_separating((0,3))
    def test_get_degree(self):
        assert self.test_graph.get_degree(0) == 3
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(0, 1) == [1,2,3]
    def test_make_essential(self):
        assert self.test_graph.make_essential() == self.test_graph

class TestGraphTheta4:
    adj_matrix = [[0,4],[4,0]]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {(tuple(self.test_graph.vertices), tuple(self.test_graph.edges)) : self.test_graph}
    def test_get_component(self):
        assert self.test_graph.get_component(0) == ((tuple(self.test_graph.vertices), tuple(self.test_graph.edges)), self.test_graph)
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 1
    def test_get_graph_minus_open_edges(self):
        assert self.test_graph.get_graph_minus_open_edges([(0,1), (0,1)]) == (([0,1], [(0,1), (0,1)]), Graph([[0,2],[2,0]]))
    def test_get_graph_minus_closed_edge(self):
        assert self.test_graph.get_graph_minus_closed_edge((0,1)) == (([], []), Graph([[]]))
    def test_get_graph_minus_vertex(self):
        assert self.test_graph.get_graph_minus_vertex(0) == (([1], []), Graph([[0]]))
    def test_prune(self):
        assert self.test_graph.prune([(0,1), (0,1), (0,1), (0,1)]) == (Graph([[0,1],[1,0]]), [(0,1), (0,1), (0,1)])
    def test_is_separating(self):
        assert not self.test_graph.is_separating((0,1))
    def test_get_degree(self):
        assert self.test_graph.get_degree(0) == 4
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(0, 2) == [1,0]
    def test_make_essential(self):
        assert self.test_graph.make_essential() == self.test_graph

class TestGraphK100:
    adj_matrix = [[1 if i != j else 0 for i in range(100)] for j in range(100)]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {(tuple(self.test_graph.vertices), tuple(self.test_graph.edges)) : self.test_graph}
    def test_get_component(self):
        assert self.test_graph.get_component(99) == ((tuple(self.test_graph.vertices), tuple(self.test_graph.edges)), self.test_graph)
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 1
    def test_get_graph_minus_open_edges(self):
        assert self.test_graph.get_graph_minus_open_edges([(i,99) for i in range(99)]) == (([i for i in range(100)], [(i,j) for i in range(100) for j in range(99) if i < j]), Graph([[1 if i != j and i != 99 and j != 99 else 0 for i in range(100)] for j in range(100)]))
    def test_get_graph_minus_closed_edge(self):
        assert self.test_graph.get_graph_minus_closed_edge((98,99)) == (([i for i in range(98)], [(i,j) for i in range(98) for j in range(98) if i < j]), Graph([[1 if i != j else 0 for i in range(98)] for j in range(98)]))
    def test_get_graph_minus_vertex(self):
        assert self.test_graph.get_graph_minus_vertex(99) == (([i for i in range(99)], [(i,j) for i in range(99) for j in range(99) if i < j]), Graph([[1 if i != j else 0 for i in range(99)] for j in range(99)]))
    def test_prune(self):
        assert self.test_graph.prune([(i,99) for i in range(99)]) == (Graph([[1 if (i == 98 and j == 99) or (i != j and i != 99 and j != 99) else 0 for i in range(100)] for j in range(100)]), [(i,99) for i in range(98)])
    def test_is_separating(self):
        assert not self.test_graph.is_separating((41,92))
    def test_get_degree(self):
        assert self.test_graph.get_degree(33) == 99
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(99, 1) == [i for i in range(99)]
    def test_make_essential(self):
        assert self.test_graph.make_essential() == self.test_graph

class TestGraph100Cycles:
    adj_matrix = [[1 if (i % 100 == (j+1) % 100 or j % 100 == (i+1) % 100) and i // 100 == j // 100 else 0 for i in range(100*100)] for j in range(100*100)]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {(tuple([100*k + i for i in range(100)]), tuple([(100*k + i, 100*k + i+1) for i in range(99)] + [(100*k, 100*k + 99)])) : 
                                                              Graph([[1 if (i == (j+1) % 100 or j == (i+1) % 100) else 0 for i in range(100)] for j in range(100)]) 
                                                              for k in range(100)}
    def test_get_component(self):
        assert self.test_graph.get_component(0) == ((tuple([i for i in range(100)]), tuple([(0, 1)] + [(0, 99)] + [(i, i+1) for i in range(1, 99)])), 
                                                    Graph([[1 if (i == (j+1) % 100 or j == (i+1) % 100) else 0 for i in range(100)] for j in range(100)]))
    def test_get_num_connected_components(self):
        assert self.test_graph.get_num_connected_components() == 100
    def test_get_graph_minus_open_edges(self):
        assert self.test_graph.get_graph_minus_open_edges([]) == ((self.test_graph.vertices, self.test_graph.edges), self.test_graph)
    def test_get_graph_minus_closed_edge(self):
        edges = [(i, i+1) for i in range(97)]
        for k in range(1, 100):
            edges += [(100*k + i, 100*k + i+1) for i in range(99)] + [(100*k, 100*k + 99)]
        assert self.test_graph.get_graph_minus_closed_edge((98,99)) == (([i for i in range(100*100) if i not in [98,99]], edges), 
                                                                        Graph([[1 if (i == j+1 or j == i+1) and (j < 97 or i < 97) else 0 for i in range(98)] + [1 if (i % 100 == (j+1) % 100 or j % 100 == (i+1) % 100) and (i-98) // 100 == (j-98) // 100 else 0 for i in range(98, 100*100 - 2)] for j in range(100*100 - 2)]))
    def test_get_graph_minus_vertex(self):
        edges = [(i, i+1) for i in range(98)]
        for k in range(1, 100):
            edges += [(100*k + i, 100*k + i+1) for i in range(99)] + [(100*k, 100*k + 99)]
        assert self.test_graph.get_graph_minus_vertex(99) == (([i for i in range(100*100) if i != 99], edges), 
                                                              Graph([[1 if (i == j+1 or j == i+1) and (j < 98 or i < 98) else 0 for i in range(99)] + [1 if (i % 100 == (j+1) % 100 or j % 100 == (i+1) % 100) and (i-99) // 100 == (j-99) // 100 else 0 for i in range(99, 100*100 - 1)] for j in range(100*100 - 1)]))
    def test_prune(self):
        assert self.test_graph.prune([(100*k, 100*k + 99) for k in range(100)]) == (Graph([[1 if (i == j+1 or j == i+1) and (j % 100 not in [99,0] or i % 100 not in [99,0]) else 0 for i in range(100*100)] for j in range(100*100)]), 
                                                                                    [(100*k, 100*k + 99) for k in range(100)])
    def test_is_separating(self):
        assert not self.test_graph.is_separating((9905,9906))
    def test_get_degree(self):
        assert self.test_graph.get_degree(9990) == 2
    def test_get_centreless_ball(self):
        assert self.test_graph.get_centreless_ball(0, 100) == [i for i in range(1, 100)] + [0]
    def test_make_essential(self):
        assert self.test_graph.make_essential() == Graph([[1 if i == j else 0 for i in range(100)] for j in range(100)])