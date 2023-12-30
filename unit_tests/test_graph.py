from graph import Graph

# TO DO: Create a class for each type of graph to be used in testing:
# - Single vertex with single edge
# - Segment graph
# - Three single vertices
# - Three segments
# - Star graph
# - Two star graphs, one subdivided
# - Star and single vertex
# - Star and segment
# - K_5 not subdivided
# - K_3,3 subdivided
# - Theta_4
# - Lambda graph
# - 100-vertex complete graph
# - 100-prong star graph
# - 100 cycles each of length 100

class TestGraphVertex:
    adj_matrix = [[0]]
    test_graph = Graph(adj_matrix)

    def test_get_connected_components(self):
        assert self.test_graph.get_connected_components() == {(self.test_graph.vertices, self.test_graph.edges) : self.test_graph}
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
        assert self.test_graph.get_connected_components() == {(self.test_graph.vertices, self.test_graph.edges) : self.test_graph}
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