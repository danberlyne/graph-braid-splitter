from graph_of_groups import GraphOfGroups
from graph import Graph
from graph_braid_group import GraphBraidGroup

class TestGraphOfGroups:
    adj_matrix = [[0]]
    test_graph = Graph(adj_matrix)
    test_vg = {0 : GraphBraidGroup(Graph([[0]]), 1)}
    test_eg = {}
    test_gog = GraphOfGroups(test_graph, test_vg, test_eg)

    def test_get_free_splitting(self):
        assert True
    def test_reduce(self):
        assert True
    def test_trim():
        assert True