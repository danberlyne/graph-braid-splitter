from graph_braid_group import integer_partitions, is_same, GraphBraidGroup
from graph import Graph

# TO DO: Create a class for each type of graph to be used in testing:
# - Single vertex, 1 particle
# - Segment of length 4, 3 particles
# - Three segments of length 2, two particles
# - 3-prong star graph and single vertex, (2,1) particles
# - Two 3-prong star graphs, one subdivided, (3,3) particles
# - K_5 subdivided, 2 particles
# - Delta graph, 2 particles
# - 100-prong star graph, 10 particles
# - 100 cycles each of length 100, 100 particles

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