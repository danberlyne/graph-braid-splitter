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

def test_integer_partitions_0_0():
    assert list(integer_partitions(0, 0)) == [tuple()]

def test_integer_partitions_5_0():
    assert list(integer_partitions(5, 0)) == [tuple()]

def test_integer_partitions_2_1():
    assert list(integer_partitions(2, 1)) == [(2,)]

def test_integer_partitions_4_2():
    assert list(integer_partitions(4, 2)) == [(0,4), (1,3), (2,2), (3,1), (4,0)]

def test_integer_partitions_10_3():
    assert list(integer_partitions(10, 3)) == [(0,0,10), (0,1,9), (0,2,8), (0,3,7), (0,4,6), (0,5,5), (0,6,4), (0,7,3), (0,8,2), (0,9,1), (0,10,0),
                                               (1,0,9), (1,1,8), (1,2,7), (1,3,6), (1,4,5), (1,5,4), (1,6,3), (1,7,2), (1,8,1), (1,9,0),
                                               (2,0,8), (2,1,7), (2,2,6), (2,3,5), (2,4,4), (2,5,3), (2,6,2), (2,7,1), (2,8,0),
                                               (3,0,7), (3,1,6), (3,2,5), (3,3,4), (3,4,3), (3,5,2), (3,6,1), (3,7,0),
                                               (4,0,6), (4,1,5), (4,2,4), (4,3,3), (4,4,2), (4,5,1), (4,6,0),
                                               (5,0,5), (5,1,4), (5,2,3), (5,3,2), (5,4,1), (5,5,0),
                                               (6,0,4), (6,1,3), (6,2,2), (6,3,1), (6,4,0),
                                               (7,0,3), (7,1,2), (7,2,1), (7,3,0),
                                               (8,0,2), (8,1,1), (8,2,0),
                                               (9,0,1), (9,1,0),
                                               (10,0,0)]

def test_is_same_vertex():
    assert is_same([[0]], [[0]])

def test_is_same_dumbell():
    adj_matrix_1 = [[0,1,0,1,1,0],[1,0,1,0,0,1],[0,1,0,0,0,1],[1,0,0,0,1,0],[1,0,0,1,0,0],[0,1,1,0,0,0]]
    adj_matrix_2 = [[0,1,1,1,0,0],[1,0,1,0,0,0],[1,1,0,0,0,0],[1,0,0,0,1,1],[0,0,0,1,0,1],[0,0,0,1,1,0]]
    assert is_same(adj_matrix_1, adj_matrix_2)

def test_is_same_two_cycles():
    adj_matrix_11 = [[0,1,0,0,0,1],[1,0,1,0,0,0],[0,1,0,1,0,0],[0,0,1,0,1,0],[0,0,0,1,0,1],[1,0,0,0,1,0]]
    adj_matrix_12 = [[0,1,1,0,0,0],[1,0,0,1,0,0],[1,0,0,0,1,0],[0,1,0,0,0,1],[0,0,1,0,0,1],[0,0,0,1,1,0]]
    adj_matrix_1 = [row + [0 for i in range(6)] for row in adj_matrix_11] + [[0 for i in range(6)] + row for row in adj_matrix_12]
    adj_matrix_21 = [[0,1,0,1,0,0],[1,0,1,0,0,0],[0,1,0,0,0,1],[1,0,0,0,1,0],[0,0,0,1,0,1],[0,0,1,0,1,0]]
    adj_matrix_22 = [[0,1,1,0,0,0],[1,0,0,0,1,0],[1,0,0,1,0,0],[0,0,1,0,0,1],[0,1,0,0,0,1],[0,0,0,1,1,0]]
    adj_matrix_2 = [row + [0 for i in range(6)] for row in adj_matrix_21] + [[0 for i in range(6)] + row for row in adj_matrix_22]
    assert is_same(adj_matrix_1, adj_matrix_2)

def test_is_same_cycle_dumbell():
    adj_matrix_1 = [[0,1,0,1,1,0],[1,0,1,0,0,1],[0,1,0,0,0,1],[1,0,0,0,1,0],[1,0,0,1,0,0],[0,1,1,0,0,0]]
    adj_matrix_2 = [[0,1,0,0,0,1],[1,0,1,0,0,0],[0,1,0,1,0,0],[0,0,1,0,1,0],[0,0,0,1,0,1],[1,0,0,0,1,0]]
    assert not is_same(adj_matrix_1, adj_matrix_2)

class TestGBGVertex_1:
    adj_matrix = [[0]]
    test_graph = Graph(adj_matrix)
    test_particles = 1
    test_config = [0]
    test_gbg = GraphBraidGroup(test_graph, test_particles, test_config)

    def test_has_sufficient_capacity(self):
        assert self.test_gbg.has_sufficient_capacity([([0], [])], (1,))
    def test_get_num_particles_per_component(self):
        assert self.test_gbg.get_num_particles_per_component([0]) == {((0,), tuple()): 1}
    def test_generate_initial_config(self):
        assert self.test_gbg.generate_initial_config({((0,), tuple()): 1}) == [0]
    def test_reindex(self):
        assert self.test_gbg.reindex([0], []) == [0]
    def test_is_trivial(self):
        assert self.test_gbg.is_trivial()
    def test_is_reduced(self):
        assert self.test_gbg.is_reduced()
    def test_get_splitting(self):
        assert self.test_gbg.get_splitting() == []
    def test_factorise(self):
        assert self.test_gbg.factorise() == []
    def test_is_same(self):
        assert self.test_gbg.is_same(GraphBraidGroup(Graph([[0,1,0,0,0],[1,0,1,0,0],[0,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]]), 3))

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