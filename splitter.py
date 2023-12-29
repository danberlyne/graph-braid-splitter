#!/usr/bin/env python3
# splitter.py - Detects free splittings of graph braid groups

import sys
import ast
import itertools
from collections import defaultdict
from graph import Graph
from graph_of_groups import GraphOfGroups
from graph_braid_group import GraphBraidGroup

MatrixDimensionException = Exception('Adjacency matrix must be square.')
MatrixFormatException = Exception('Adjacency matrix must either only contain positive integers and spaces or use the list of lists format.')
MatrixSymmetryException = Exception('Adjacency matrix must be symmetric.')
GraphFormatException = Exception('Graph cannot contain loops or multi-edges.')
ParticleFormatException = Exception('Number of particles must be a positive integer.')
ConfigDimensionException = Exception('Number of integers in initial configuration must equal number of particles.')
ConfigFormatException = Exception('Initial configuration must only contain spaces and integers between 1 and the number of rows in the adjacency matrix.')
VertexException = Exception('Each connected component of the graph must have at least as many vertices as particles starting in that component.')

def start_exit_sequence():
    while True:
        exit_sequence = input('Press [Enter] to exit.')
        if not exit_sequence:
            sys.exit()

# Prompt user to enter data manually if it has not been entered in 'gbg_data.txt'.
def enter_data_manually():
    while True:
        print('No data detected in gbg_data.txt. Enter data manually? (y/n)')
        response = input().lower()
        if response == 'n':
            prompt_data_entry()
        elif response == 'y':
            num_particles = enter_particles()
            adj_matrix = enter_matrix()
            initial_config = enter_config(num_particles, adj_matrix)
            return(num_particles, adj_matrix, initial_config)

def prompt_data_entry():
    while True:
        exit_sequence = input('Please enter your data in gbg_data.txt then run the script again. Press [Enter] to exit.')
        if not exit_sequence:
            sys.exit()

# Loop for manual entry of number of particles.
def enter_particles():
    while True:
        particle_response = input('Please enter the number of particles/strands in the graph braid group:\n')
        if particle_response.isdecimal():
            if int(particle_response) > 0:
                return int(particle_response)
            else:
                print('The number of particles must be a positive integer.')    
        else:
            print('The number of particles must be a positive integer.')

# Loop for manual entry of adjacency matrix.
def enter_matrix():
    while True:
        matrix_response_1 = input('Please enter the first row of the adjacency matrix of the graph, with entries separated by single spaces:\n').split()
        correct_formatting = verify_matrix_formatting(matrix_response_1, matrix_response_1, 1)
        if correct_formatting:
            matrix_response = [matrix_response_1]
            for row_num in range(2, len(matrix_response_1) + 1):
                while True:
                    current_matrix_response = input(f'Please enter row {row_num} of the adjacency matrix of the graph:\n').split()
                    correct_formatting = verify_matrix_formatting(matrix_response_1, current_matrix_response, row_num)
                    if correct_formatting:
                        is_symmetric = verify_symmetric(matrix_response, current_matrix_response, row_num)
                        if is_symmetric:
                            matrix_response.append(current_matrix_response)
                            break
            return [[int(entry) for entry in matrix_response[i]] for i in range(len(matrix_response_1))]
        
def verify_matrix_formatting(initial_matrix_response, current_matrix_response, row_num):
    if len(initial_matrix_response) == 0 or len(current_matrix_response) != len(initial_matrix_response):
        print('Incorrect formatting.')
        return False
    for char in current_matrix_response:
        if not char.isdecimal():
            print('Incorrect formatting.')
            return False
        elif int(char) < 0:
            print('Integers cannot be negative.')
            return False
        elif int(char) > 1:
            print('Graph cannot contain multi-edges.')
            return False
    if int(current_matrix_response[row_num - 1]) != 0:
        print('Graph cannot contain loops.')
        return False
    return True

def verify_symmetric(matrix_response, current_matrix_response, row_num):
    for j in range(row_num - 1):
        if current_matrix_response[j] != matrix_response[j][row_num - 1]:
            print('Matrix must be symmetric.')
            return False
    return True
        
# Loop for manual entry of initial configuration.
def enter_config(num_particles, adj_matrix):
    if Graph(adj_matrix).get_num_connected_components() == 1:
        return None
    else:
        while True:
            config_response = input(f'Please enter the initial configuration of the {num_particles} particles on the vertices of the graph. Enter this data in the form of {num_particles} integers separated by single spaces, corresponding to the row numbers of the vertices in the adjacency matrix:\n').split()
            correct_formatting = verify_config_formatting(config_response, num_particles, adj_matrix)
            if correct_formatting:
                return [int(entry) - 1 for entry in config_response]
            
def verify_config_formatting(config_response, num_particles, adj_matrix):
    if len(config_response) != num_particles:
        print('Number of integers must match number of particles.')
        return False
    for char in config_response:
        if not char.isdecimal():
            print('Incorrect formatting.')
            return False
        elif int(char) < 1 or int(char) > len(adj_matrix):
            print(f'Integers must be between 1 and {len(adj_matrix)}.')
            return False
    return True

# Takes as input a list of strings corresponding from the lines of `gbg_data.txt` below the dotted line in the file. 
# Returns the number of particles, adjacency matrix of the graph, and the initial configuration in the correct format.
def get_data_from_file(gbg_data):
    num_particles = convert_particle_data(gbg_data)
    if is_list_format(gbg_data):
        adj_matrix = convert_matrix_data_list(gbg_data)
    else:
        adj_matrix = convert_matrix_data_lines(gbg_data)
    initial_config = convert_config_data(gbg_data, num_particles, adj_matrix)
    return (num_particles, adj_matrix, initial_config)

def convert_particle_data(gbg_data):
    particle_data = gbg_data[0].lstrip('Number of particles:').strip()
    if particle_data.isdecimal():
        if int(particle_data) > 0:
            return int(particle_data)
        else:
            raise ParticleFormatException
    else: 
        raise ParticleFormatException

def is_list_format(gbg_data):
    if gbg_data[2].replace(' ', '')[:2] == '[[':
        return True
    else:
        return False

def convert_matrix_data_list(gbg_data):
    # Convert data into a list of rows of the matrix, with each row in the form of a string of integers separated by commas.
    rows = gbg_data[2].replace(' ', '').lstrip('[[').rstrip(']]').split('],[')
    string_matrix = [rows[i].split(',') for i in range(len(rows))]
    for i in range(len(string_matrix)):
        for entry in string_matrix[i]:
            if not entry.isdecimal():
                raise MatrixFormatException
            elif int(entry) < 0:
                raise MatrixFormatException
        if len(string_matrix[i]) != len(string_matrix):
            raise MatrixDimensionException
    for i in range(len(string_matrix)):
        for j in range(i+1):
            if string_matrix[i][j] != string_matrix[j][i]:
                raise MatrixSymmetryException
    return [[int(string_matrix[i][j]) for j in range(len(string_matrix))] for i in range(len(string_matrix))]

def convert_matrix_data_lines(gbg_data):
    matrix_size = len(gbg_data[2].split())
    for i in range(matrix_size):
        if len(gbg_data[2 + i]) >= 21:
            if gbg_data[2 + i][:21] == 'Initial configuration':
                raise MatrixDimensionException
        for entry in gbg_data[2 + i].split():
            if not entry.isdecimal():
                raise MatrixFormatException
            elif int(entry) < 0:
                raise MatrixFormatException
        if len(gbg_data[2 + i].split()) != matrix_size:
            raise MatrixDimensionException
    adj_matrix = [[int(entry) for entry in gbg_data[2 + i].split()] for i in range(matrix_size)]
    for i in range(matrix_size):
        for j in range(i+1):
            if adj_matrix[i][j] != adj_matrix[j][i]:
                raise MatrixSymmetryException
    return adj_matrix

def convert_config_data(gbg_data, num_particles, adj_matrix):
    if is_list_format(gbg_data):
        config_data = gbg_data[3].lstrip('Initial configuration:').split()
    else:
        config_data = gbg_data[2 + len(adj_matrix)].lstrip('Initial configuration:').split()
    if config_data == []:
        if Graph(adj_matrix).get_num_connected_components() == 1:
            return None
        else:
            raise ConfigDimensionException
    elif len(config_data) != num_particles:
        raise ConfigDimensionException
    else:
        for entry in config_data:
            if not entry.isdecimal():
                raise ConfigFormatException
            elif int(entry) < 1 or int(entry) > len(adj_matrix):
                raise ConfigFormatException
    return [int(entry) - 1 for entry in config_data]

# Replaces all factors in a splitting with string versions and adds detailed data to `splitting.txt` where appropriate.
def stringify_factors(splitting, braid_factor_counter, gog_factor_counter):
    known_gbgs = get_known_gbgs_from_file()
    for i, factor in enumerate(splitting):
        if isinstance(factor, GraphBraidGroup):
            if factor.is_reduced():
                essential_graph = factor.graph.make_essential()
                essential_adj_matrix_hashable = tuple(tuple(row) for row in essential_graph.adj_matrix)
                if (essential_adj_matrix_hashable, factor.num_particles) in known_gbgs: 
                    splitting[i] = known_gbgs[(essential_adj_matrix_hashable, factor.num_particles)]
                else:
                    file = open('splitting.txt', 'r')
                    splitting_data = file.readlines()
                    file.close()
                    previous_matrices = [ast.literal_eval(line.partition('adjacency matrix: ')[2].strip()) for line in splitting_data if line.startswith('Ga')]
                    is_new_graph = True
                    for previous_matrix in previous_matrices:
                        if is_same(essential_graph.adj_matrix, previous_matrix):
                            splitting[i] = f'B_{factor.num_particles}(Gamma_{previous_matrices.index(previous_matrix) + 1})'
                            is_new_graph = False
                            break
                    if is_new_graph:
                        splitting[i] = f'B_{factor.num_particles}(Gamma_{braid_factor_counter})'
                        file = open('splitting.txt', 'a')
                        file.write(f'Gamma_{braid_factor_counter} adjacency matrix: {essential_graph.adj_matrix} \n')
                        file.close()
                        braid_factor_counter += 1
            else:
                file = open('splitting.txt', 'r')
                splitting_data = file.readlines()
                file.close()
                previous_matrices = [ast.literal_eval(line.partition('adjacency matrix: ')[2].strip()) for line in splitting_data if line.startswith('Ga')]
                is_new_graph = True
                for previous_matrix in previous_matrices:
                    if is_same(factor.adj_matrix, previous_matrix):
                        splitting[i] = f'RB_{factor.num_particles}(Gamma_{previous_matrices.index(previous_matrix) + 1})'
                        is_new_graph = False
                        break
                if is_new_graph:
                    splitting[i] = f'RB_{factor.num_particles}(Gamma_{braid_factor_counter})'
                    file = open('splitting.txt', 'a')
                    file.write(f'Gamma_{braid_factor_counter} adjacency matrix: {factor.adj_matrix} \n')
                    file.close()
                    braid_factor_counter += 1
        elif isinstance(factor, GraphOfGroups):
            splitting[i] = f'G_{gog_factor_counter}'
            file = open('splitting.txt', 'a')
            file.write(f'G_{gog_factor_counter} adjacency matrix: {factor.graph.adj_matrix} \n') 
            for v in factor.vertex_groups:
                if factor.vertex_groups[v].is_reduced():
                    file.write(f'Vertex group {v+1}: B_{factor.vertex_groups[v].num_particles}, adjacency matrix: {factor.vertex_groups[v].graph.adj_matrix}, initial configuration: {[i+1 for i in factor.vertex_groups[v].initial_config]} \n')
                else:
                    file.write(f'Vertex group {v+1}: RB_{factor.vertex_groups[v].num_particles}, adjacency matrix: {factor.vertex_groups[v].graph.adj_matrix}, initial configuration: {[i+1 for i in factor.vertex_groups[v].initial_config]} \n')
            for e in factor.edge_groups:
                if factor.edge_groups[e].is_reduced():
                    file.write(f'Edge group {(e[0]+1, e[1]+1)}: B_{factor.edge_groups[e].num_particles}, adjacency matrix: {factor.edge_groups[e].graph.adj_matrix}, initial configuration: {[i+1 for i in factor.edge_groups[e].initial_config]} \n')
                else:
                    file.write(f'Edge group {(e[0]+1, e[1]+1)}: RB_{factor.edge_groups[e].num_particles}, adjacency matrix: {factor.edge_groups[e].graph.adj_matrix}, initial configuration: {[i+1 for i in factor.edge_groups[e].initial_config]} \n')
            file.close()  
            gog_factor_counter += 1      
        elif isinstance(factor, list):
            stringify_factors(factor, braid_factor_counter, gog_factor_counter)

# Gets list of known graph braid groups from `known_gbgs.txt` and puts them in a dictionary.
def get_known_gbgs_from_file():
    file = open('known_gbgs.txt')
    known_gbg_data = [line.strip() for line in file.readlines()]
    file.close()
    # Dictionary where keys are (adj_matrix, num_particles) tuples and values are the braid groups as strings.
    known_gbgs = {(ast.literal_eval(known_gbg_data[3*i + 1]), int(known_gbg_data[3*i + 2])): known_gbg_data[3*i] for i in range(int(len(known_gbg_data) / 3))}
    # Adds all possible permutations of the vertex labels of the graph to the dictionary.
    known_gbgs_unordered = {braid_data: known_gbgs[braid_data] for braid_data in known_gbgs}
    for (matrix, particles) in known_gbgs:
        for perm in itertools.permutations(range(len(matrix)), len(matrix)):
            permuted_matrix = tuple(tuple(matrix[perm[i]][perm[j]] for i in range(len(matrix))) for j in range(len(matrix)))
            known_gbgs_unordered.setdefault((permuted_matrix, particles), known_gbgs[(matrix, particles)])
    return known_gbgs_unordered

# Takes as input two adjacency matrices and returns True if they define the same graph.
# Tip: To check if two Graph objects are homeomorphic, first use `make_essential` on both and then use `is_same` on their adjacency matrices.
def is_same(adj_matrix_1, adj_matrix_2):
    matrix_2_permutations = [[[adj_matrix_2[perm[i]][perm[j]] for i in range(len(adj_matrix_2))] for j in range(len(adj_matrix_2))] for perm in itertools.permutations(range(len(adj_matrix_2)), len(adj_matrix_2))]
    if adj_matrix_1 in matrix_2_permutations:
        return True
    else:
        return False

# Converts `stringified_splitting` into a list of strings, where each string represents the splitting of a direct factor of the original group. 
# These strings have been converted from list form into single strings with 'x' and '*' characters as appropriate.
def combine_strings(stringified_splitting, is_direct = True):
    for i, factor in enumerate(stringified_splitting):
        # If a factor is a list of non-lists, then replace the factor in `stringified_splitting` with its string representation.
        if isinstance(factor, list):
            for subfactor in factor:
                # If we find a list in factor's subfactors, then feed factor back into `combine_strings`.
                if isinstance(subfactor, list):
                    combine_strings(factor, not is_direct)
                    break
            # If `stringified splitting` is a direct splitting, then `factor` is a free splitting.
            if is_direct:
                # Collect together F_m terms and Z terms.
                free_rank = sum(int(subfactor[2:]) for subfactor in factor if subfactor.startswith('F_')) + len([subfactor for subfactor in factor if subfactor == 'Z'])
                factor = [f'F_{k}' for k in [free_rank] if free_rank > 1] + ['Z' for k in [free_rank] if free_rank == 1] + [subfactor for subfactor in factor if not subfactor.startswith('F_') and not subfactor == 'Z']
                # Collect together like terms and express in simplified notation.
                factor_count = defaultdict(int)
                for subfactor in factor:
                    factor_count[subfactor] += 1
                factor = [f'*^{factor_count[subfactor]}(' + subfactor + ')' for subfactor in factor_count if factor_count[subfactor] > 1] + [subfactor for subfactor in factor_count if factor_count[subfactor] == 1]
                if len(factor) > 1:
                    stringified_splitting[i] = '(' + ' * '.join(factor) + ')'
                else:
                    stringified_splitting[i] = factor[0]
            else:
                # Collect together like terms and express as powers.
                factor_count = defaultdict(int)
                for subfactor in factor:
                    factor_count[subfactor] += 1
                factor = [subfactor + f'^{factor_count[subfactor]}' for subfactor in factor_count if factor_count[subfactor] > 1 and len(subfactor) == 1] + ['(' + subfactor + ')' + f'^{factor_count[subfactor]}' for subfactor in factor_count if factor_count[subfactor] > 1 and len(subfactor) > 1] + [subfactor for subfactor in factor_count if factor_count[subfactor] == 1]
                if len(factor) > 1:
                    stringified_splitting[i] = '(' + ' x '.join(factor) + ')'
                else:
                    stringified_splitting[i] = factor[0]


#############
# Main code #
#############

def main():
    print('Checking gbg_data.txt...')
    file = open('gbg_data.txt')
    file_as_list = file.readlines()
    file.close()
    for i, line in enumerate(file_as_list):
        file_as_list[i] = line.strip()
    non_empty_lines = [line for line in file_as_list if line != '']
    dotted_line = non_empty_lines.index('-' * 50)
    gbg_data = non_empty_lines[dotted_line + 1:]

    # If data has not been entered in `gbg_data.txt`, prompt user to enter it manually.
    if len(gbg_data) == 3:
        (num_particles, adj_matrix, initial_config) = enter_data_manually()
    # Otherwise, get data from `gbg_data.txt` and verify it is formatted correctly.
    else: 
        (num_particles, adj_matrix, initial_config) = get_data_from_file(gbg_data)

    gbg = GraphBraidGroup(Graph(adj_matrix), num_particles, initial_config)
    if not gbg.is_reduced():
        print('''WARNING: In order to perform computations for B_n(\Gamma), the graph \Gamma must satisfy the following conditions:
1. All cycles must have length at least n+1.
2. All paths between vertices of degree not equal to 2 must have length at least n-1.
At least one of these conditions is not satisfied by your graph.
If you choose to continue, any results obtained will only be true for the reduced braid group RB_n(\Gamma).
Do you wish to continue? (y/n)''')
        while True:
            response = input().lower()
            if response == 'n':
                while True:
                    exit_sequence = input('Please amend your data then run the script again. Press [Enter] to exit.')
                    if not exit_sequence:
                        sys.exit()
            elif response == 'y':
                for comp in gbg.num_initial_particles_per_component:
                    if len(comp[0]) < gbg.num_initial_particles_per_component[comp]:
                        raise VertexException
                break

    print('Verified successfully.\nSearching for free splittings...')

    if gbg.is_trivial():
        if gbg.is_reduced():
            print(f'B_{num_particles} = 1')
            start_exit_sequence()
        else:
            print(f'RB_{num_particles} = 1')
            start_exit_sequence()

    # Returns a nested list, where the odd levels of the list correspond to direct factors and the even levels correspond to free factors.
    splitting = gbg.get_splitting()

    print('Search complete.')

    # length == 1 means it does not split as a direct product.
    if len(splitting) == 1:
        if not isinstance(splitting[0], str):
            if not isinstance(splitting[0], list):
                print('No splittings found.')
                start_exit_sequence()
            elif len(splitting[0]) == 1 and not isinstance(splitting[0][0], str):
                print('No splittings found.')
                start_exit_sequence()

    print('Splitting found. Converting data to readable format...')

    # Makes a fresh copy of `splitting.txt`, the file containing detailed splitting information (for the graph of group factors and the graph braid group factors).
    file = open('splitting.txt', 'w')
    file.write('')
    file.close()

    # Turns all factors in the splitting into strings and adds detailed data to `splitting.txt` where appropriate.
    stringify_factors(splitting, 1, 1)

    combine_strings(splitting)
    final_string = ' x '.join(splitting)

    # Prints the splitting and adds it to the beginning of `splitting.txt`
    with open('splitting.txt','r') as contents:
        save = contents.read()
    if gbg.is_reduced():
        with open('splitting.txt','w') as contents:
            contents.write(f'B_{num_particles} = ' + final_string + '\n\n')
        with open('splitting.txt','a') as contents:
            contents.write(save)
        print(f'B_{num_particles} = ' + final_string)
    else:
        with open('splitting.txt','w') as contents:
            contents.write(f'RB_{num_particles} = ' + final_string + '\n\n')
        with open('splitting.txt','a') as contents:
            contents.write(save)
        print(f'RB_{num_particles} = ' + final_string)

    if 'G_' in final_string or 'Gamma_' in final_string:
        print('In the above splitting, Gamma_i are graphs and G_i are fundamental groups of graphs of groups. More detailed data can be found in splitting.txt.')
    start_exit_sequence()

if __name__ == '__main__':
    main()