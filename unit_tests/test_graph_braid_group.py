from graph_braid_group import integer_partitions, is_same, GraphBraidGroup
from graph import Graph
from graph_of_groups import GraphOfGroups
from math import comb
from itertools import combinations

def split_star(prongs, particles):
    if prongs > particles + 1 and particles > 2:
        return [[f'F_{comb(prongs-1, particles-1) - 1}', split_star(prongs-1, particles-1), split_star(prongs-1, particles)]]
    elif prongs == particles + 1 and particles > 2:
        return [[f'F_{comb(prongs-1, particles-1) - 1}', split_star(prongs-1, particles-1)]]
    elif prongs > particles + 1 and particles == 2:
        return [[f'F_{comb(prongs-1, particles-1) - 1}', split_star(prongs-1, particles)]]
    else:
        return [['F_1']]

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

class TestGBGSegment_5:
    # Segment of length 4
    adj_matrix = [[0,1,0,0,0],[1,0,1,0,0],[0,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]]
    test_graph = Graph(adj_matrix)
    test_particles = 5
    test_config = None
    test_gbg = GraphBraidGroup(test_graph, test_particles, test_config)

    def test_get_graph_of_groups(self):
        assert self.test_gbg.get_graph_of_groups([(0,1), (1,2)]).is_same(GraphOfGroups(Graph([[0]]), {0: GraphBraidGroup(Graph([[0,0,0,0,0],[0,0,0,0,0],[0,0,0,1,0],[0,0,1,0,1],[0,0,0,1,0]]), 5)}, {}))
    def test_get_compatible_particles_per_component(self):
        assert self.test_gbg.get_compatible_particles_per_component([(1,2)]) == [{((0,1), ((0,1),)): 2, ((2,3,4), ((2,3), (3,4))): 3}]
    def test_has_sufficient_capacity(self):
        assert self.test_gbg.has_sufficient_capacity([(self.test_graph.vertices, self.test_graph.edges)], (5,))
    def test_get_num_particles_per_component(self):
        assert self.test_gbg.get_num_particles_per_component(self.test_gbg.initial_config) == {(tuple(self.test_graph.vertices), tuple(self.test_graph.edges)): 5}
    def test_generate_initial_config(self):
        assert self.test_gbg.generate_initial_config({(tuple(self.test_graph.vertices), tuple(self.test_graph.edges)): 5}) == [0,1,2,3,4]
    def test_get_adjacent_assignments(self):
        assert self.test_gbg.get_adjacent_assignments([(1,2)], [{((0,1), ((0,1),)): 2, ((2,3,4), ((2,3), (3,4))): 3}]) == []
    def test_get_compatible_assignments(self):
        assert True
    def test_reindex(self):
        assert self.test_gbg.reindex([0,2,4], [1,3]) == [0,1,2]
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

class TestGBG3Segment_2_2_0:
    # Three segments of length 2
    adj_matrix = [[0,1,0,0,0,0,0,0,0],[1,0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0],
                  [0,0,0,0,1,0,0,0,0],[0,0,0,1,0,1,0,0,0],[0,0,0,0,1,0,0,0,0],
                  [0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,1,0,1],[0,0,0,0,0,0,0,1,0]]
    test_graph = Graph(adj_matrix)
    test_particles = 4
    test_config = [0,1,3,4]
    test_gbg = GraphBraidGroup(test_graph, test_particles, test_config)

    def test_get_graph_of_groups(self):
        assert self.test_gbg.get_graph_of_groups([(0,1)]).is_same(GraphOfGroups(Graph([[0,1],[1,0]]), 
                                                                           {0: GraphBraidGroup(Graph([[0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,1,0,1,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,1,0,1],[0,0,0,0,0,0,0,1,0]]), 4, [1,2,3,4]), 
                                                                            1: GraphBraidGroup(Graph([[0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,1,0,1,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,1,0,1],[0,0,0,0,0,0,0,1,0]]), 4, [0,2,3,4])}, 
                                                                           {(0,1,(0,1),0): GraphBraidGroup(Graph([[0,0,0,0,0,0,0],[0,0,1,0,0,0,0],[0,1,0,1,0,0,0],[0,0,1,0,0,0,0],[0,0,0,0,0,1,0],[0,0,0,0,1,0,1],[0,0,0,0,0,1,0]]), 3, [0,1,2])})
                                                                 )
    def test_get_compatible_particles_per_component(self):
        assert self.test_gbg.get_compatible_particles_per_component([(0,1)]) == [{((0,), tuple()): 0, ((1,2), ((1,2),)): 2, ((3,4,5), ((3,4),(4,5))): 2, ((6,7,8), ((6,7),(7,8))): 0}, 
                                                                                 {((0,), tuple()): 1, ((1,2), ((1,2),)): 1, ((3,4,5), ((3,4),(4,5))): 2, ((6,7,8), ((6,7),(7,8))): 0}]
    def test_has_sufficient_capacity(self):
        assert not self.test_gbg.has_sufficient_capacity([((0,1,2), ((0,1),(1,2))), ((3,4,5), ((3,4),(4,5))), ((6,7,8), ((6,7),(7,8)))], (4,0,0))
    def test_get_num_particles_per_component(self):
        assert self.test_gbg.get_num_particles_per_component(self.test_gbg.initial_config) == {((0,1,2), ((0,1),(1,2))): 2, ((3,4,5), ((3,4),(4,5))): 2}
    def test_generate_initial_config(self):
        assert self.test_gbg.generate_initial_config({((0,1,2), ((0,1),(1,2))): 1, ((3,4,5), ((3,4),(4,5))): 0, ((6,7,8), ((6,7),(7,8))): 3}) == [0,6,7,8]
    def test_get_adjacent_assignments(self):
        assert self.test_gbg.get_adjacent_assignments([(0,1)], 
                                                      [{((0,), tuple()): 1, ((1,2), ((1,2),)): 1, ((3,4,5), ((3,4),(4,5))): 2, ((6,7,8), ((6,7),(7,8))): 0}, 
                                                       {((0,), tuple()): 0, ((1,2), ((1,2),)): 2, ((3,4,5), ((3,4),(4,5))): 2, ((6,7,8), ((6,7),(7,8))): 0}], 
                                                     ) == [({((0,), tuple()): 1, ((1,2), ((1,2),)): 1, ((3,4,5), ((3,4),(4,5))): 2, ((6,7,8), ((6,7),(7,8))): 0}, 
                                                            {((0,), tuple()): 0, ((1,2), ((1,2),)): 2, ((3,4,5), ((3,4),(4,5))): 2, ((6,7,8), ((6,7),(7,8))): 0},
                                                            (0,1))]
    def test_get_compatible_assignments(self):
        assert self.test_gbg.get_compatible_assignments([(0,1)], 
                                                    {((0,), tuple()): 1, ((1,2), ((1,2),)): 1, ((3,4,5), ((3,4),(4,5))): 2, ((6,7,8), ((6,7),(7,8))): 0}, 
                                                    {((0,), tuple()): 0, ((1,2), ((1,2),)): 2, ((3,4,5), ((3,4),(4,5))): 2, ((6,7,8), ((6,7),(7,8))): 0},
                                                    (0,1)
                                                    ) == [{((2,), tuple()): 1, ((3,4,5), ((3,4),(4,5))): 2, ((6,7,8), ((6,7),(7,8))): 0}]
    def test_reindex(self):
        assert self.test_gbg.reindex([1,4,7], [0,3,6]) == [0,2,4]
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

class TestGBGTwo3ProngStar_3_3:
    # Two 3-pronged star graphs, one subdivided
    adj_matrix = [[0,1,1,1,0,0,0,0,0,0,0],
                  [1,0,0,0,0,0,0,0,0,0,0],
                  [1,0,0,0,0,0,0,0,0,0,0],
                  [1,0,0,0,0,0,0,0,0,0,0],
                  [0,0,0,0,0,1,1,1,0,0,0],
                  [0,0,0,0,1,0,0,0,1,0,0],
                  [0,0,0,0,1,0,0,0,0,1,0],
                  [0,0,0,0,1,0,0,0,0,0,1],
                  [0,0,0,0,0,1,0,0,0,0,0],
                  [0,0,0,0,0,0,1,0,0,0,0],
                  [0,0,0,0,0,0,0,1,0,0,0]]
    test_graph = Graph(adj_matrix)
    test_particles = 6
    test_config = [0,1,2,4,5,6]
    test_gbg = GraphBraidGroup(test_graph, test_particles, test_config)

    def test_get_graph_of_groups(self):
        assert self.test_gbg.get_graph_of_groups([(4,5)]).is_same(GraphOfGroups(Graph([[0,2,0],[2,0,3],[0,3,0]]), 
                                                                           {0: GraphBraidGroup(Graph([[0,1,1,1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,1,0,0,0],[0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,1,0,0,0,0,1,0],[0,0,0,0,1,0,0,0,0,0,1],[0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,1,0,0,0]]), 6, [0,1,2,5,7,8]), 
                                                                            1: GraphBraidGroup(Graph([[0,1,1,1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,1,0,0,0],[0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,1,0,0,0,0,1,0],[0,0,0,0,1,0,0,0,0,0,1],[0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,1,0,0,0]]), 6, [0,1,2,5,6,7]),
                                                                            2: GraphBraidGroup(Graph([[0,1,1,1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,1,0,0,0],[0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,1,0,0,0,0,1,0],[0,0,0,0,1,0,0,0,0,0,1],[0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,1,0,0,0]]), 6, [0,1,2,4,6,7])}, 
                                                                           {(0,1,(4,5),0): GraphBraidGroup(Graph([[0,1,1,1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0,0]]), 5, [0,1,2,5,6]),
                                                                            (0,1,(4,5),1): GraphBraidGroup(Graph([[0,1,1,1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0,0]]), 5, [0,1,2,4,6]),
                                                                            (1,2,(4,5),0): GraphBraidGroup(Graph([[0,1,1,1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0,0]]), 5, [0,1,2,5,8]),
                                                                            (1,2,(4,5),1): GraphBraidGroup(Graph([[0,1,1,1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0,0]]), 5, [0,1,2,4,5]),
                                                                            (1,2,(4,5),2): GraphBraidGroup(Graph([[0,1,1,1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0,0]]), 5, [0,1,2,4,7])})
                                                                 )
    def test_get_compatible_particles_per_component(self):
        assert self.test_gbg.get_compatible_particles_per_component([(4,5)]) == [{((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 1, ((5,8), ((5,8),)): 2}, 
                                                                                 {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 2, ((5,8), ((5,8),)): 1},
                                                                                 {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 3, ((5,8), ((5,8),)): 0}]
    def test_has_sufficient_capacity(self):
        assert not self.test_gbg.has_sufficient_capacity([((0,1,2,3), ((0,1),(0,2),(0,3))), ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))), ((5,8), ((5,8),))], (3,0,3))
    def test_get_num_particles_per_component(self):
        assert self.test_gbg.get_num_particles_per_component(self.test_config) == {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,5,6,7,8,9,10), ((4,5),(4,6),(4,7),(5,8),(6,9),(7,10))): 3}
    def test_generate_initial_config(self):
        assert self.test_gbg.generate_initial_config({((0,1,2,3), ((0,1),(0,2),(0,3))): 0, ((4,5,6,7,8,9,10), ((4,5),(4,6),(4,7),(5,8),(6,9),(7,10))): 6}) == [4,5,6,7,8,9]
    def test_get_adjacent_assignments(self):
        assert self.test_gbg.get_adjacent_assignments([(4,5)], 
                                                      [{((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 3, ((5,8), ((5,8),)): 0}, 
                                                       {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 2, ((5,8), ((5,8),)): 1},
                                                       {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 1, ((5,8), ((5,8),)): 2}], 
                                                     ) == [({((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 3, ((5,8), ((5,8),)): 0}, 
                                                            {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 2, ((5,8), ((5,8),)): 1},
                                                            (4,5)),
                                                           ({((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 2, ((5,8), ((5,8),)): 1}, 
                                                            {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 1, ((5,8), ((5,8),)): 2},
                                                            (4,5))]
    def test_get_compatible_assignments(self):
        assert self.test_gbg.get_compatible_assignments([(4,5)], 
                                                        {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 3, ((5,8), ((5,8),)): 0}, 
                                                        {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((4,6,7,9,10), ((4,6),(4,7),(6,9),(7,10))): 2, ((5,8), ((5,8),)): 1},
                                                        (4,5)
                                                       ) == [{((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((6,9), ((6,9),)): 0, ((7,10), ((7,10),)): 2, ((8,), tuple()): 0},
                                                             {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((6,9), ((6,9),)): 1, ((7,10), ((7,10),)): 1, ((8,), tuple()): 0},
                                                             {((0,1,2,3), ((0,1),(0,2),(0,3))): 3, ((6,9), ((6,9),)): 2, ((7,10), ((7,10),)): 0, ((8,), tuple()): 0}]
    def test_reindex(self):
        assert True
    def test_is_trivial(self):
        assert not self.test_gbg.is_trivial()
    def test_is_reduced(self):
        assert not self.test_gbg.is_reduced()
    def test_get_splitting(self):
        assert self.test_gbg.get_splitting() == [['F_3']]
    def test_factorise(self):
        factor_matrix = [[0,1,1,1,0,0,0],
                         [1,0,0,0,1,0,0],
                         [1,0,0,0,0,1,0],
                         [1,0,0,0,0,0,1],
                         [0,1,0,0,0,0,0],
                         [0,0,1,0,0,0,0],
                         [0,0,0,1,0,0,0]]
        assert len(self.test_gbg.factorise()) == 1 and self.test_gbg.factorise()[0].is_same(GraphBraidGroup(Graph(factor_matrix), 3))
    def test_is_same(self):
        assert not self.test_gbg.is_same(GraphBraidGroup(Graph([[0,1,0,0,0],[1,0,1,0,0],[0,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]]), 3))

class TestGBGK5_2:
    # Complete graph on 5 vertices
    adj_matrix = [[0,1,1,1,1],
                  [1,0,1,1,1],
                  [1,1,0,1,1],
                  [1,1,1,0,1],
                  [1,1,1,1,0]]
    test_graph = Graph(adj_matrix)
    test_particles = 2
    test_config = [0,1]
    test_gbg = GraphBraidGroup(test_graph, test_particles, test_config)

    def test_get_graph_of_groups(self):
        assert self.test_gbg.get_graph_of_groups([(3,4)]).is_same(GraphOfGroups(Graph([[1]]), 
                                                                  {0: GraphBraidGroup(Graph([[0,1,1,1,1],[1,0,1,1,1],[1,1,0,1,1],[1,1,1,0,0],[1,1,1,0,0]]), 2, [0,1])}, 
                                                                  {(0,0,(3,4),0): GraphBraidGroup(Graph([[0,1,1],[1,0,1],[1,1,0]]), 1, [0])})
                                                                 )
    def test_get_compatible_particles_per_component(self):
        assert self.test_gbg.get_compatible_particles_per_component([(3,4)]) == [{((0,1,2,3,4), ((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4))): 2}]
    def test_has_sufficient_capacity(self):
        assert self.test_gbg.has_sufficient_capacity([((0,1,2,3,4), ((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4)))], (2,))
    def test_get_num_particles_per_component(self):
        assert self.test_gbg.get_num_particles_per_component(self.test_config) == {((0,1,2,3,4), ((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4))): 2}
    def test_generate_initial_config(self):
        assert self.test_gbg.generate_initial_config({((0,1,2,3,4), ((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4))): 2}) == [0,1]
    def test_get_adjacent_assignments(self):
        assert self.test_gbg.get_adjacent_assignments([(3,4)], 
                                                      [{((0,1,2,3,4), ((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4))): 2}] 
                                                     ) == [({((0,1,2,3,4), ((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4))): 2}, 
                                                            {((0,1,2,3,4), ((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4))): 2},
                                                            (3,4))]
    def test_get_compatible_assignments(self):
        assert self.test_gbg.get_compatible_assignments([(3,4)], 
                                                        {((0,1,2,3,4), ((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4))): 2}, 
                                                        {((0,1,2,3,4), ((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4))): 2},
                                                        (3,4)
                                                       ) == [{((0,1,2), ((0,1),(0,2),(1,2))): 1}]
    def test_reindex(self):
        assert True
    def test_is_trivial(self):
        assert not self.test_gbg.is_trivial()
    def test_is_reduced(self):
        assert self.test_gbg.is_reduced()
    def test_get_splitting(self):
        assert len(self.test_gbg.get_splitting()) == 1 and self.test_gbg.get_splitting()[0].is_same(self.test_gbg)
    def test_factorise(self):
        assert len(self.test_gbg.factorise()) == 1 and self.test_gbg.factorise()[0].is_same(self.test_gbg)
    def test_is_same(self):
        assert not self.test_gbg.is_same(GraphBraidGroup(Graph([[0,1,0,0,0],[1,0,1,0,0],[0,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]]), 3))

class TestGBG10Prong_5:
    # 10-prong radial tree, not subdivided
    adj_matrix = [[0] + [1]*10] + [[1] + [0]*10 for i in range(10)]
    test_graph = Graph(adj_matrix)
    test_particles = 5
    test_config = [0,1,2,3,4]
    test_gbg = GraphBraidGroup(test_graph, test_particles)

    def test_get_graph_of_groups(self):
        gog_graph = Graph([[0, comb(9,4)], [comb(9,4), 0]])
        graph_minus_open_edge = Graph([[0,0] + [1]*9] + [[0]*11] + [[1] + [0]*10 for i in range(9)])
        graph_minus_closed_edge = Graph([[0] * 9] * 9)
        assert self.test_gbg.get_graph_of_groups([(0,1)]).is_same(GraphOfGroups(gog_graph, 
                                                                  {0: GraphBraidGroup(graph_minus_open_edge, 5, [0,1,2,3,4]), 1: GraphBraidGroup(graph_minus_open_edge, 5, [0,2,3,4,5])}, 
                                                                  {(0,1,(0,1),i): GraphBraidGroup(graph_minus_closed_edge, 4, [0,1,2,3]) for i in range(comb(9,4))})
                                                                 )
    def test_get_compatible_particles_per_component(self):
        assert self.test_gbg.get_compatible_particles_per_component([(0,1)]) == [{((0,) + tuple([i+2 for i in range(9)]), tuple([(0,i+2) for i in range(9)])): 4, ((1,), tuple()): 1},
                                                                                 {((0,) + tuple([i+2 for i in range(9)]), tuple([(0,i+2) for i in range(9)])): 5, ((1,), tuple()): 0}
                                                                                ]
    def test_has_sufficient_capacity(self):
        assert self.test_gbg.has_sufficient_capacity([((0,)+tuple([i+2 for i in range(9)]), ((0,i+2) for i in range(9))), ((1,), tuple())], (4,1))
    def test_get_num_particles_per_component(self):
        assert self.test_gbg.get_num_particles_per_component(self.test_config) == {(tuple([i for i in range(11)]), tuple([(0,i+1) for i in range(10)])): 5}
    def test_generate_initial_config(self):
        assert self.test_gbg.generate_initial_config({((0,) + tuple([i+2 for i in range(9)]), tuple([(0,i+2) for i in range(9)])): 4, ((1,), tuple()): 1}) == [0,2,3,4,1]
    def test_get_adjacent_assignments(self):
        assert self.test_gbg.get_adjacent_assignments([(0,1)], 
                                                      [{((0,) + tuple([i+2 for i in range(9)]), tuple([(0,i+2) for i in range(9)])): 5, ((1,), tuple()): 0},
                                                       {((0,) + tuple([i+2 for i in range(9)]), tuple([(0,i+2) for i in range(9)])): 4, ((1,), tuple()): 1}] 
                                                     ) == [({((0,) + tuple([i+2 for i in range(9)]), tuple([(0,i+2) for i in range(9)])): 5, ((1,), tuple()): 0},
                                                            {((0,) + tuple([i+2 for i in range(9)]), tuple([(0,i+2) for i in range(9)])): 4, ((1,), tuple()): 1},
                                                            (0,1))]
    def test_get_compatible_assignments(self):
        assert sorted(self.test_gbg.get_compatible_assignments([(0,1)], 
                                                               {((0,) + tuple([i+2 for i in range(9)]), tuple([(0,i+2) for i in range(9)])): 5, ((1,), tuple()): 0},
                                                               {((0,) + tuple([i+2 for i in range(9)]), tuple([(0,i+2) for i in range(9)])): 4, ((1,), tuple()): 1},
                                                               (0,1)), 
                      key=lambda d: tuple(d.values())
                     ) == sorted([{((i,), tuple()): 1 if i in c else 0 for i in range(2,11)} for c in combinations(range(2,11),4)], 
                                 key=lambda d: tuple(d.values()))
    def test_reindex(self):
        assert True
    def test_is_trivial(self):
        assert not self.test_gbg.is_trivial()
    def test_is_reduced(self):
        assert not self.test_gbg.is_reduced()
    def test_get_splitting(self):
        print(f'True:\n{self.test_gbg.get_splitting()}')
        print(f'Test:\n{split_star(10,5)}')
        assert self.test_gbg.get_splitting() == split_star(10, 5)
    def test_factorise(self):
        assert len(self.test_gbg.factorise()) == 1 and self.test_gbg.factorise()[0].is_same(self.test_gbg)
    def test_is_same(self):
        assert not self.test_gbg.is_same(GraphBraidGroup(Graph([[0,1,0,0,0],[1,0,1,0,0],[0,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]]), 3))