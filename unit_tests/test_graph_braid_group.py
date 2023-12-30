from graph_braid_group import integer_partitions, is_same, GraphBraidGroup
from graph import Graph

def test_integer_partitions():
    assert True

def test_is_same():
    assert True

class TestGraphBraidGroup:
    adj_matrix = [[0]]
    test_graph = Graph(adj_matrix)
    test_particles = 1
    test_config = [0]
    test_gbg = GraphBraidGroup(test_graph, test_particles, test_config)

    def test_get_graph_of_groups(self):
        assert True
    def test_get_compatible_particles_per_component(self):
        assert True
    def test_has_sufficient_capacity(self):
        assert True
    def test_get_num_particles_per_component(self):
        assert True
    def test_generate_initial_config(self):
        assert True
    def test_get_adjacent_assignments(self):
        assert True
    def test_get_compatible_assignments(self):
        assert True
    def test_reindex(self):
        assert True
    def test_is_trivial(self):
        assert True
    def test_is_reduced(self):
        assert True
    def test_get_splitting(self):
        assert True
    def test_factorise(self):
        assert True
    def test_is_same(self):
        assert True