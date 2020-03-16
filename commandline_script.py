from networkx import all_shortest_paths
from Class_wrDCJ_Node import Node

from Class_extremities_and_adjacencies import Extremities_and_adjacencies
import New_Network_wrDCJ

import time
import argparse
import sys
t0 = time.time()

def run(args):
    genomeA_file = args.source_genome
    genomeB_file = args.target_genome
    weight_ratios_file = args.ratios
    stdoutOrigin = sys.stdout
    sys.stdout = open(args.output_file, 'w')
    #outfile = open(args.output_file, 'w')
    with open(genomeA_file) as csv:
        line = [element.strip('\n').split(',') for element in csv]
    genomeA = []

    for element in line:
        element = list(map(int, element))
        genomeA.append(element)

    with open(genomeB_file) as csv:
        line = [element.strip('\n').split(',') for element in csv]
    genomeB = []

    for element in line:
        element = list(map(int, element))
        genomeB.append(element)

    with open(weight_ratios_file) as csv:
        line = [element.strip('\n').split(',') for element in csv]
    weight_ratios = []

    for element in line:
        element = list(map(int, element))
        weight_ratios.append(element)

    get_adjacencies = Extremities_and_adjacencies()
    adjacencies_genomeA = get_adjacencies.adjacencies_ordered_and_sorted(genomeA)


    adjacencies_genomeB = get_adjacencies.adjacencies_ordered_and_sorted(genomeB)

    #Create start and target node
    start_node = Node(adjacencies_genomeA)
    target_node = Node(adjacencies_genomeB)

    hash_table = {}
    hash_key_start = hash(str(start_node.state))
    hash_key_target = hash(str(target_node.state))
    hash_table.update({hash_key_start:start_node})
    hash_table.update({hash_key_target:target_node})

    #finding rearrangement weights
    max_number = max(weight_ratios[0])
    weights = []
    for element in weight_ratios[0]:
        if element == 0:
            weights.append(max_number^2)
        else:
            weights.append(max_number/element)

    New_Network_wrDCJ.build_hash_table(start_node, hash_table, adjacencies_genomeB, weights)

    network = New_Network_wrDCJ.build_network(hash_table)

    shortest_paths = (list(all_shortest_paths(network, start_node, target_node, weight='weight')))

    j = 1
    tot_b_trl = 0
    tot_u_trl = 0
    tot_inv = 0
    tot_trp1 = 0
    tot_trp2 = 0
    tot_fus = 0
    tot_fis = 0

    Paths_state = []
    Paths_state_weight = []
    # print(shortest_paths[0][4].children_weights[2])
    for path in shortest_paths:
        path_state = []
        path_state_weight = []

        i = 0
        while i < len(path):
            current = path[i]
            if i == 0:
                operation_type = 'none, this is the source genome'
                operation_weight = 'N/A'
                operation = 'N/A'
            else:
                x = path[i - 1].children.index(current)

                operation_type = path[i - 1].children_operations[x][1]
                operation_weight = path[i - 1].children_weights[x]
                operation = path[i - 1].children_operations[x][0]

            adjacencies = current.state
            genome = get_adjacencies.adjacencies_to_genome(adjacencies)
            path_state_weight.append((genome, ((operation_type, operation), operation_weight)))

            path_state.append((genome, (operation_type, operation)))

            i += 1
        Paths_state.append((path_state))
        Paths_state_weight.append(path_state_weight)

    for path in shortest_paths:

        i = 0
        b_trl = 0
        u_trl = 0
        inv = 0
        trp1 = 0
        trp2 = 0
        fus = 0
        fis = 0
        while i < len(path):

            current = path[i]
            if i == 0:
                pass
            else:
                x = path[i - 1].children.index(current)
                operation_type = path[i - 1].children_operations[x][1]
                if operation_type == 'b_trl':
                    b_trl += 1
                elif operation_type == 'u_trl':
                    u_trl += 1
                elif operation_type == 'inv':
                    inv += 1
                elif operation_type == 'trp1':
                    trp1 += 1
                elif operation_type == 'trp2':
                    trp2 += 1
                elif operation_type == 'fus':
                    fus += 1
                elif operation_type == 'fis':
                    fis += 1
            i += 1

        tot_b_trl += b_trl
        tot_u_trl += u_trl
        tot_inv += inv
        tot_trp1 += trp1
        tot_trp2 += trp2
        tot_fus += fus
        tot_fis += fis
        j += 1


    print('############################################################################################################')
    print()
    print('Source Genome: ', genomeA)
    print('Target Genome: ', genomeB)
    print()
    print('Number of most parsimonious solutions: ', len(shortest_paths))
    print()
    print('Average number of each operation per solution:')
    print('Inversions: ', int(tot_inv/len(shortest_paths)), '  Transpositions type 1: ', int(tot_trp1/len(shortest_paths)), '  Transpositions type 2: ', int(tot_trp2/len(shortest_paths)), '  Balanced translocations: ', int(tot_b_trl/len(shortest_paths)), '  Unbalanced translocations: ', int(tot_u_trl/len(shortest_paths)),
          '  Fusions: ', int(tot_fus/len(shortest_paths)),
          '  Fissions: ', int(tot_fis/len(shortest_paths)))
    print()
    print()
    print('Solutions: ')
    print()
    path_counter = 1
    for path in Paths_state:
        print('Solution number ', path_counter)
        for genome in path:
            print(genome)
        path_counter+=1
        print()
    print()
    print('############################################################################################################')



    ###############################
    # JUST FOR TESTING

    solution = [([[1, 2, 3, 4, 15], [-8, -7, 6, -5, -14, -13, -12], [9, 11],
                  [-20, -19, -18, -17, -16, -32, 10, -31, -30, -29, -28, -27], [21, 22, 23, 24, 25, 26], [-33],
                  [34, 35, 36, 37, 38, 39, 40]], ('none, this is the source genome', 'N/A')), (
                [[1, 2, 3, 4, 15], [-8, -7, -6, -5, -14, -13, -12], [9, 11],
                 [-20, -19, -18, -17, -16, -32, 10, -31, -30, -29, -28, -27], [21, 22, 23, 24, 25, 26], [-33],
                 [34, 35, 36, 37, 38, 39, 40]], ('inv', (((5.5, 6.5), (6, 7)), ((5.5, 6), (6.5, 7))))), (
                [[1, 2, 3, 4, 15], [-8, -7, -6, -5, -14, -13, -12], [9, 11], [16, 17, 18, 19, 20],
                 [21, 22, 23, 24, 25, 26], [27, 28, 29, 30, 31, -10, 32, 33], [34, 35, 36, 37, 38, 39, 40]],
                ('u_trl', (((16, 32.5), 33), ((32.5, 33), 16)))), (
                [[1, 2, 3, 4, 5, 6, 7, 8], [9, 11], [12, 13, 14, 15], [16, 17, 18, 19, 20], [21, 22, 23, 24, 25, 26],
                 [27, 28, 29, 30, 31, -10, 32, 33], [34, 35, 36, 37, 38, 39, 40]],
                ('b_trl', (((4.5, 15), (5, 14.5)), ((4.5, 5), (14.5, 15))))), (
                [[1, 2, 3, 4, 5, 6, 7, 8], [9, 11], [12, 13, 14, 15], [16, 17, 18, 19, 20], [21, 22, 23, 24, 25, 26],
                 [27, 28, 29, 30, 31, 32, 33], [34, 35, 36, 37, 38, 39, 40], ['o', 10]],
                ('trp0', (((10, 32), (10.5, 31.5)), ((10, 10.5), (31.5, 32))))), (
                [[1, 2, 3, 4, 5, 6, 7, 8], [9, 10, 11], [12, 13, 14, 15], [16, 17, 18, 19, 20],
                 [21, 22, 23, 24, 25, 26], [27, 28, 29, 30, 31, 32, 33], [34, 35, 36, 37, 38, 39, 40]],
                ('trp1', (((9.5, 11), (10, 10.5)), ((9.5, 10), (10.5, 11)))))]

    paths_operations = []
    for element in Paths_state:
        path_operations = [y for (x, y) in element]

        paths_operations.append(path_operations)

    solution_operations = [d for (c, d) in solution]

    path_types = []
    sol_types = [a for a, b in solution_operations]
    for element in paths_operations:
        types = [c for c, d in element]
        path_types.append(types)

    indexes = []
    counter = 0
    for element in path_types:
        if element == sol_types:
            indexer = path_types.index(element)
            counter += 1
            indexes.append(indexer)

    # for element in indexes:
    #     for x in Paths_state[element]:
    #         print(x)
    #     print()
    # print('*****')
    # for element in solution:
    #     print(element)

    print('sol len: ', len(solution))
    print('shortest path len: ',len(shortest_paths[0]))
    print('counter', counter)

    print('And the answer is... ', solution_operations in paths_operations)

    print('Source genome: ',genomeA)
    print('Target genome: ', genomeB)
    print()
    print('Solution: ', solution)
    ##########################################################################################################
    sys.stdout.close()
    sys.stdout=stdoutOrigin
def main():
    parser=argparse.ArgumentParser(description='A program that outputs all the optimal set of rearrangment operations that can descripe the evolution of one genome into another')
    parser.add_argument("-genB", help="this is the set of genes representing the target genome", dest='target_genome', required=True)
    parser.add_argument("-genA", help="this is the set of genes representing the source genome",
                        dest='source_genome', required=True,)
    parser.add_argument("-list_of_rearrangement_ratios", help='the ratios in which each rearrangement is expected to occur in the order inversions:transpositions:balanced translocations:unbalanced translocations:fissions:fusions', dest='ratios', required=True)
    parser.add_argument("-output", help="the name of the output file that will contain the set of rearrangements", dest='output_file', required=True)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main()