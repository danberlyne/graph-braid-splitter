#!/usr/bin/env python3
# splitter.py - Detects free splittings of graph braid groups
"""
A module that preprocesses users inputs, detects free splittings using the other modules, then outputs the results to the terminal.

Module contents: 
- `start_exit_sequence`: Prompts user to press Enter to exit the program.
- `enter_data_manually`: Handles manual entry of graph braid group data.
- `prompt_data_entry`: Prompts user to enter data in `gbg_data.txt`.
- `enter_particles`: Handles manual entry of particle data for the graph braid group.
- `enter_matrix`: Handles manual entry of the adjacency matrix for the graph braid group.
- `verify_matrix_formatting`: Verifies the entered adjacency matrix has valid dimensions and values.
- `verify_symmetric`: Verifies the entered adjacency matrix is symmetric.
- `enter_config`: Handles manual entry of initial configuration for the graph braid group.
- `verify_config_formatting`: Verifies the entered configuration has valid length and values.
- `get_data_from_file`: Parses data obtained from `gbg_data.txt` and extracts graph braid group information.
- `convert_particle_data`: Converts the line of `gbg_data.txt` defining the number of particles into the correct format.
- `is_list_format`: Detects whether the adjacency matrix in `gbg_data.txt` was entered as a list of lists.
- `convert_matrix_data_list`: If the adjacency matrix was entered as a list of lists, converts it into the correct format.
- `convert_matrix_data_lines`: If the adjacency matrix was entered line by line, converts it into the correct format.
- `convert_config_data`: Converts the line of `gbg_data.txt` defining the initial configuration into the correct format.
- `stringify_factors`: Converts factors of the detected free splitting into strings.
- `get_known_gbgs_from_file`: Reads and parses graph braid group data from `known_gbgs.txt`.
- `is_same`: Detects whether two adjacency matrices define the same graph.
- `combine_strings`: Converts the stringified splitting from list of lists format by inserting 'x' and '*' characters.
- `main`: Code that executes when running splitter from terminal.
"""

import sys
import ast
import itertools
from collections import defaultdict
from graph import Graph
from graph_of_groups import GraphOfGroups
from graph_braid_group import GraphBraidGroup
from typing import Union

# Recursive types used for free splittings.
NestedFactorList = list[Union[GraphOfGroups, 'GraphBraidGroup', 'NestedFactorList']]
StringifiedNestedFactorList = list[Union[str, 'StringifiedNestedFactorList']]

MatrixDimensionException = Exception('Adjacency matrix must be square.')
MatrixFormatException = Exception('Adjacency matrix must either only contain positive integers and spaces or use the list of lists format.')
MatrixSymmetryException = Exception('Adjacency matrix must be symmetric.')
GraphFormatException = Exception('Graph cannot contain loops or multi-edges.')
ParticleFormatException = Exception('Number of particles must be a positive integer.')
ConfigDimensionException = Exception('Number of integers in initial configuration must equal number of particles.')
ConfigFormatException = Exception('Initial configuration must only contain spaces and integers between 1 and the number of rows in the adjacency matrix.')
VertexException = Exception('Each connected component of the graph must have at least as many vertices as particles starting in that component.')

def start_exit_sequence() -> None:
    """Prompts the user to press Enter to exit the program."""
    while True:
        exit_sequence = input('Press [Enter] to exit.')
        if not exit_sequence:
            sys.exit()

def enter_data_manually() -> tuple[int, list[list[int]], list[int] | None]:
    """Prompts user to enter data manually if it has not been entered in `gbg_data.txt`."""
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

def prompt_data_entry() -> None:
    """Prompts user to enter data in `gbg_data.txt` and exits the program."""
    while True:
        exit_sequence = input('Please enter your data in gbg_data.txt then run the script again. Press [Enter] to exit.')
        if not exit_sequence:
            sys.exit()

def enter_particles() -> int:
    """Prompts manual entry of number of particles of the graph braid group."""
    while True:
        particle_response = input('Please enter the number of particles/strands in the graph braid group:\n')
        if particle_response.isdecimal():
            if int(particle_response) > 0:
                return int(particle_response)
            else:
                print('The number of particles must be a positive integer.')    
        else:
            print('The number of particles must be a positive integer.')

def enter_matrix() -> list[list[int]]:
    """Prompts manual entry of adjacency matrix of the graph braid group, one row at a time."""
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
        
def verify_matrix_formatting(initial_matrix_response: list[str], current_matrix_response: list[str], row_num: int) -> bool:
    """
    Returns `True` if the given row of the adjacency matrix has been entered in the correct format.
    
    `initial_matrix_response` should be the first row of the matrix that was entered, as a list of strings.
    `current_matrix_response` should be the row of the entered matrix that is currently being checked, as a list of strings.
    `row_num` should be the row number (starting from 1) of the row of the entered matrix that is currently being checked, as an integer.

    Returns `True` if: 
    - `current_matrix_response` contains only values of '1' and '0', plus spaces;
    - the length of `current_matrix_response` matches that of `initial_matrix_response`;
    - the value of `current_matrix_response` that lies on the diagonal of the matrix is '0'.
    """
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

def verify_symmetric(matrix_response: list[list[str]], current_matrix_response: list[str], row_num: int) -> bool:
    """
    Returns `True` if the given row of the adjacency matrix maintains symmetry of the overall matrix.
    
    `matrix_response` should be the rows of the matrix that have already been checked, as a list of lists of strings.
    `current_matrix_response` should be the row of the entered matrix that is currently being checked, as a list of strings.
    `row_num` should be the row number (starting from 1) of the row of the entered matrix that is currently being checked.

    Returns `True` if the values of `current_matrix_response` match the values of `matrix_response` that lie in column `row_num`.
    """
    for j in range(row_num - 1):
        if current_matrix_response[j] != matrix_response[j][row_num - 1]:
            print('Matrix must be symmetric.')
            return False
    return True
        
def enter_config(num_particles: int, adj_matrix: list[list[int]]) -> list[int] | None:
    """
    Prompts manual user entry of initial configuration of the graph braid group and returns the configuration if entered correctly.
    
    `num_particles` should be the number of particles of the graph braid group.
    `adj_matrix` should be the adjacency matrix of the graph braid group.

    If the user entered the configuration in the correct format, converts the data to a list of integers, indexing from 0, and return the list.
    """
    if Graph(adj_matrix).get_num_connected_components() == 1:
        return None
    else:
        while True:
            config_response = input(f'Please enter the initial configuration of the {num_particles} particles on the vertices of the graph. Enter this data in the form of {num_particles} integers separated by single spaces, corresponding to the row numbers of the vertices in the adjacency matrix:\n').split()
            correct_formatting = verify_config_formatting(config_response, num_particles, adj_matrix)
            if correct_formatting:
                return [int(entry) - 1 for entry in config_response]
            
def verify_config_formatting(config_response: list[str], num_particles: int, adj_matrix: list[list[int]]) -> bool:
    """
    Returns `True` if the initial configuration has been entered in the correct format.
    
    `config_response` should be the configuration that was entered, as a list of strings.
    `num_particles` should be the number of particles of the graph braid group, as an integer.
    `adj_matrix` should be the adjacency matrix of the graph braid group, as a list of lists of integers.

    Returns `True` if: 
    - `config_response` has length `num_particles`;
    - `config_response` contains only decimal string values between '1' and the number of vertices.
    """
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

def get_data_from_file(gbg_data: list[str]) -> tuple[int, list[list[int]], list[int] | None]:
    """
    Returns the number of particles, adjacency matrix, and initial configuration that were entered in `gbg_data.txt`.

    `gbg_data` should be a list of strings corresponding to the lines of `gbg_data.txt` below the dotted line in the file.

    Returns a 3-tuple consisting of the number of particles, adjacency matrix, and the initial configuration in the correct format.
    """
    num_particles = convert_particle_data(gbg_data)
    if is_list_format(gbg_data):
        adj_matrix = convert_matrix_data_list(gbg_data)
    else:
        adj_matrix = convert_matrix_data_lines(gbg_data)
    initial_config = convert_config_data(gbg_data, num_particles, adj_matrix)
    return (num_particles, adj_matrix, initial_config)

def convert_particle_data(gbg_data: list[str]) -> int:
    """
    Converts the number of particles found in `gbg_data.txt` to integer format and returns it.

    `gbg_data` should be a list of strings corresponding to the lines of `gbg_data.txt` below the dotted line in the file.
    """
    particle_data = gbg_data[0].lstrip('Number of particles:').strip()
    if particle_data.isdecimal():
        if int(particle_data) > 0:
            return int(particle_data)
        else:
            raise ParticleFormatException
    else: 
        raise ParticleFormatException

def is_list_format(gbg_data: list[str])-> bool:
    """
    Returns `True` if the adjacency matrix found in `gbg_data.txt` was entered in list of lists format.

    `gbg_data` should be a list of strings corresponding to the lines of `gbg_data.txt` below the dotted line in the file.
    """
    if gbg_data[2].replace(' ', '')[:2] == '[[':
        return True
    else:
        return False

def convert_matrix_data_list(gbg_data: list[str]) -> list[list[int]]:
    """
    Converts the adjacency matrix found in `gbg_data.txt` to a list of lists of integers and returns it.

    `gbg_data` should be a list of strings corresponding to the lines of `gbg_data.txt` below the dotted line in the file.
    The adjacency matrix in `gbg_data` should be in list of lists format (but as a string).
    """
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

def convert_matrix_data_lines(gbg_data: list[str]) -> list[list[int]]:
    """
    Converts the adjacency matrix found in `gbg_data.txt` to a list of lists of integers and returns it.

    `gbg_data` should be a list of strings corresponding to the lines of `gbg_data.txt` below the dotted line in the file.
    The adjacency matrix in `gbg_data` should be in the format of lines of integer characters separated by spaces.
    """
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

def convert_config_data(gbg_data: list[str], num_particles: int, adj_matrix: list[list[int]]) -> list[int] | None:
    """
    Converts the initial configuration found in `gbg_data.txt` (if it exists) to a list of integers and returns it.

    `gbg_data` should be a list of strings corresponding to the lines of `gbg_data.txt` below the dotted line in the file.
    `num_particles` should be the number of particles of the graph braid group, as an integer.
    `adj_matrix` should be the adjacency matrix of the graph braid group, as a list of lists of integers.

    If an initial configuration was entered in the correct format, converts it to a list of integers (reindexed from 0) and returns it.
    If no initial configuration was entered, returns `None`.
    """
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

def stringify_factors(splitting: NestedFactorList, braid_factor_counter: int, gog_factor_counter: int) -> None:
    """
    Recursively replaces all factors in a splitting with string versions and adds detailed data to `splitting.txt` where appropriate.

    `splitting` should be a list of lists of... of lists of factors of the graph braid group.
    The levels of list should alternate between direct factors and free factors.
    The lists should contain only `GraphBraidGroup` objects, `GraphOfGroups` objects, strings, or lists.
    `braid_factor_counter` is an integer that keeps track of how many `GraphBraidGroup` factors have been stringified so far.
    `gog_factor_counter` is an integer that keeps track of how many `GraphOfGroups` factors have been stringified so far.

    Modifies `splitting` in place, replacing all `GraphBraidGroup` and `GraphOfGroups` objects with string versions.
    When an object is replaced with a string version, further data on the object (e.g. adjacency matrix) is added to `splitting.txt`.
    Returns `None`.
    """
    known_gbgs = get_known_gbgs_from_file()
    for i, factor in enumerate(splitting):
        if isinstance(factor, GraphBraidGroup):
            if factor.is_reduced():
                essential_graph = factor.graph.make_essential()
                essential_adj_matrix_hashable = tuple(tuple(row) for row in essential_graph.adj_matrix)
                if (essential_adj_matrix_hashable, factor.num_particles) in known_gbgs: 
                    splitting[i] = known_gbgs[(essential_adj_matrix_hashable, factor.num_particles)]
                else:
                    with open('splitting.txt', 'r') as file:
                        splitting_data = file.readlines()
                    previous_matrices = [ast.literal_eval(line.partition('adjacency matrix: ')[2].strip()) for line in splitting_data if line.startswith('Ga')]
                    is_new_graph = True
                    for previous_matrix in previous_matrices:
                        if is_same(essential_graph.adj_matrix, previous_matrix):
                            splitting[i] = f'B_{factor.num_particles}(Gamma_{previous_matrices.index(previous_matrix) + 1})'
                            is_new_graph = False
                            break
                    if is_new_graph:
                        splitting[i] = f'B_{factor.num_particles}(Gamma_{braid_factor_counter})'
                        with open('splitting.txt', 'a') as file:
                            file.write(f'Gamma_{braid_factor_counter} adjacency matrix: {essential_graph.adj_matrix} \n')
                        braid_factor_counter += 1
            else:
                with open('splitting.txt', 'r') as file:
                    splitting_data = file.readlines()
                previous_matrices = [ast.literal_eval(line.partition('adjacency matrix: ')[2].strip()) for line in splitting_data if line.startswith('Ga')]
                is_new_graph = True
                for previous_matrix in previous_matrices:
                    if is_same(factor.adj_matrix, previous_matrix):
                        splitting[i] = f'RB_{factor.num_particles}(Gamma_{previous_matrices.index(previous_matrix) + 1})'
                        is_new_graph = False
                        break
                if is_new_graph:
                    splitting[i] = f'RB_{factor.num_particles}(Gamma_{braid_factor_counter})'
                    with open('splitting.txt', 'a') as file:
                        file.write(f'Gamma_{braid_factor_counter} adjacency matrix: {factor.adj_matrix} \n')
                    braid_factor_counter += 1
        elif isinstance(factor, GraphOfGroups):
            splitting[i] = f'G_{gog_factor_counter}'
            with open('splitting.txt', 'a') as file:
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
            gog_factor_counter += 1      
        elif isinstance(factor, list):
            stringify_factors(factor, braid_factor_counter, gog_factor_counter)

def get_known_gbgs_from_file() -> dict[tuple[tuple[tuple[int, ...], ...], int], str]:
    """
    Gets list of known graph braid groups from `known_gbgs.txt` and returns them as a dictionary.

    In the returned dictionary:
    - keys are (`adj_matrix`, `num_particles`) 2-tuples;
    - values are the braid groups as strings.
    """
    with open('known_gbgs.txt') as file:
        known_gbg_data = [line.strip() for line in file.readlines()]
    # Dictionary where keys are (adj_matrix, num_particles) tuples and values are the braid groups as strings.
    known_gbgs = {(ast.literal_eval(known_gbg_data[3*i + 1]), int(known_gbg_data[3*i + 2])): known_gbg_data[3*i] for i in range(int(len(known_gbg_data) / 3))}
    # Adds all possible permutations of the vertex labels of the graph to the dictionary.
    known_gbgs_unordered = {braid_data: known_gbgs[braid_data] for braid_data in known_gbgs}
    for (matrix, particles) in known_gbgs:
        for perm in itertools.permutations(range(len(matrix)), len(matrix)):
            permuted_matrix = tuple(tuple(matrix[perm[i]][perm[j]] for i in range(len(matrix))) for j in range(len(matrix)))
            known_gbgs_unordered.setdefault((permuted_matrix, particles), known_gbgs[(matrix, particles)])
    return known_gbgs_unordered

def is_same(adj_matrix_1: list[list[int]] | tuple[tuple[int, ...], ...], adj_matrix_2: list[list[int]] | tuple[tuple[int, ...], ...]) -> bool:
    """
    Takes as input two adjacency matrices and returns `True` if they define the same graph.

    `adj_matrix_1` and `adj_matrix_2` should be lists of lists of integers or tuples of tuples of integers.

    Checks if any permutation of the rows and columns of `adj_matrix_2` makes it equal to `adj_matrix_1`.

    Tip: To check if two `Graph` objects are homeomorphic, first use `make_essential` on both and then use `is_same` on their adjacency matrices.
    """
    matrix_2_permutations = [[[adj_matrix_2[perm[i]][perm[j]] for i in range(len(adj_matrix_2))] for j in range(len(adj_matrix_2))] for perm in itertools.permutations(range(len(adj_matrix_2)), len(adj_matrix_2))]
    if adj_matrix_1 in matrix_2_permutations:
        return True
    else:
        return False

def combine_strings(stringified_splitting: StringifiedNestedFactorList, is_direct: bool = True) -> list[str]:
    """
    Converts a stringified splitting into a list of strings and returns it, where each string represents the splitting of a direct factor of the original group. 
    
    `stringified_splitting` should be a list of lists of lists of... of lists of strings, representing direct/free factors of the group.
    The first list level represents the direct factors of the group.
    The second level represents the free factors of those direct factors. 
    This continues in this way, alternating between free and direct factors as we proceed deeper into the lists.

    The returned list has a single level, with all of its entries being strings.
    The strings in the returned list therefore represent the direct factors of the group.
    These strings have been converted from list form into single strings with 'x' and '*' characters as appropriate.
    """
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

def main() -> None:
    """Implements command line interface."""
    print('Checking gbg_data.txt...')
    with open('gbg_data.txt') as file:
        file_as_list = file.readlines()
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
    with open('splitting.txt', 'w') as file:
        file.write('')

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