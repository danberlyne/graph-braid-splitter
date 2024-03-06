from graph_of_groups import GraphOfGroups
from graph import Graph
from graph_braid_group import GraphBraidGroup

# TO DO: Create a class for each of:
# - Single vertex
# - Single edge
# - Single loop
# - K5
# - Segment of length 100
# - Large binary tree

class TestGOGVertex:
    adj_matrix = [[0]]
    test_graph = Graph(adj_matrix)
    test_vg = {0: GraphBraidGroup(Graph([[0,1,1],[1,0,1],[1,1,0]]), 1)} # Vertex group is Z
    test_eg = {}
    test_gog = GraphOfGroups(test_graph, test_vg, test_eg)

    def test_get_free_splitting(self):
        assert self.test_gog.get_free_splitting() == (self.test_gog,)
    def test_reduce(self):
        assert True
    def test_trim(self):
        assert True
    def test_is_same(self):
        new_gog = GraphOfGroups(self.test_graph, {0: GraphBraidGroup(Graph([[0,2],[2,0]]), 1)}, self.test_eg)
        assert self.test_gog.is_same(new_gog)

class TestGOGEdge:
    adj_matrix = [[0,1],[1,0]]
    test_graph = Graph(adj_matrix)
    test_vg = {0: GraphBraidGroup(Graph([[0,1,1],[1,0,1],[1,1,0]]), 1), 
               1: GraphBraidGroup(Graph([[0,1,1],[1,0,1],[1,1,0]]), 1)} # Vertex groups are both Z
    test_eg = {(0,1): GraphBraidGroup(Graph([[0,1,0],[1,0,1],[0,1,0]]), 1)} # Edge group is trivial
    test_gog = GraphOfGroups(test_graph, test_vg, test_eg)

    def test_get_free_splitting(self):
        gph = Graph([[0]])
        vg = {0: GraphBraidGroup(Graph([[0,1,1],[1,0,1],[1,1,0]]), 1)}
        eg = {}
        splitting = ('F_0', GraphOfGroups(gph, vg, eg), GraphOfGroups(gph, vg, eg))
        assert self.test_gog.get_free_splitting()[0] == splitting[0]
        for i in range(1, len(splitting)):
            assert self.test_gog.get_free_splitting()[i].is_same(splitting[i])
    def test_reduce(self):
        assert self.test_gog.reduce([(0,1)]) == (([0,1], [(0,1)]), self.test_vg, self.test_eg)
    def test_trim(self):
        assert self.test_gog.trim(((0,1), tuple())) == ((0,1), tuple())
    def test_is_same(self):
        new_gog = GraphOfGroups(Graph([[0]]), {0: GraphBraidGroup(Graph([[0,2],[2,0]]), 1)}, {})
        assert not self.test_gog.is_same(new_gog)