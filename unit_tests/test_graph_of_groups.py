from graph_of_groups import GraphOfGroups
from graph import Graph
from graph_braid_group import GraphBraidGroup

# TO DO: Create a class for each of the outputs of the tests for get_graph_of_groups in the gbg section, plus these:
# - Single vertex
# - Single edge
# - Single loop
# - K5
# - Segment of length 100
# - Large binary tree

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