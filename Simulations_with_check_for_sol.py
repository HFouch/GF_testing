from networkx import all_shortest_paths
from Class_wrDCJ_Node import Node
from Class_extremities_and_adjacencies import Extremities_and_adjacencies
import New_Network_wrDCJ
import GenomeEvolve
import time
t0 = time.time()

number_of_simulations = 10000

results = []
genomeB = [[1,2,3,4,5,6,7,8,9,10], [11,12,13,14,15,16,17,18,19,20], [21, 22,23,24,25,26,27,28,29,30], [31,32,33,34,35,36,37,38,39,40], [41,42,43,44,45,46,47,48,49,50]]

weight_ratios = [[1,1,1,1,1,1]]

get_adjacencies = Extremities_and_adjacencies()
adjacencies_genomeB = get_adjacencies.adjacencies_ordered_and_sorted(genomeB)

number_of_solutions_found = 0

for i in range(0, number_of_simulations):
    print('simultation: ', i)
    genomeB_copy = genomeB[:]
    evolution_simulation = GenomeEvolve.get_evolved_genome_and_solution(genomeB_copy)
    genomeA = evolution_simulation[1]
    sorting_scenario = evolution_simulation[2]

    adjacencies_genomeA = get_adjacencies.adjacencies_ordered_and_sorted(genomeA)

    # Create start and target node
    start_node = Node(adjacencies_genomeA)
    target_node = Node(adjacencies_genomeB)

    hash_table = {}
    hash_key_start = hash(str(start_node.state))
    hash_key_target = hash(str(target_node.state))
    hash_table.update({hash_key_start: start_node})
    hash_table.update({hash_key_target: target_node})

    max_number = max(weight_ratios[0])
    weights = []

    for element in weight_ratios[0]:
        if element == 0:
            weights.append(max_number ^ 2)
        else:
            weights.append(max_number / element)

    New_Network_wrDCJ.build_hash_table(start_node, hash_table, adjacencies_genomeB, weights)

    network = New_Network_wrDCJ.build_network(hash_table)

    shortest_paths = (list(all_shortest_paths(network, start_node, target_node, weight='weight')))

    Paths_state = []
    Paths_state_weight = []
    # print(shortest_paths[0][4].children_weights[2])

    j = 1
    tot_b_trl = 0
    tot_u_trl = 0
    tot_inv = 0
    tot_trp1 = 0
    tot_trp2 = 0
    tot_fus = 0
    tot_fis = 0
    ave_b_trl = 0
    ave_u_trl = 0
    ave_inv = 0
    ave_trp1 = 0
    ave_trp2 = 0
    ave_fus = 0
    ave_fis = 0



    for path in shortest_paths:
        path_state = []
        #path_state_weight = []

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
                operation_type = 'none, this is the source genome'
                #operation_weight = 'N/A'
                operation = 'N/A'
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


                ###
                x = path[i - 1].children.index(current)

                operation_type = path[i - 1].children_operations[x][1]
                #operation_weight = path[i - 1].children_weights[x]
                operation = path[i - 1].children_operations[x][0]

                ###


            tot_b_trl += b_trl
            tot_u_trl += u_trl
            tot_inv += inv
            tot_trp1 += trp1
            tot_trp2 += trp2
            tot_fus += fus
            tot_fis += fis
            j += 1


            ###
            adjacencies = current.state
            genome = get_adjacencies.adjacencies_to_genome(adjacencies)
            #path_state_weight.append((genome, ((operation_type, operation), operation_weight)))

            path_state.append((genome, (operation_type, operation)))
            ###

            i += 1



        Paths_state.append((path_state))
        #Paths_state_weight.append(path_state_weight)


    ave_b_trl = tot_b_trl / len(shortest_paths)
    ave_u_trl = tot_u_trl / len(shortest_paths)
    ave_inv = tot_inv / len(shortest_paths)
    ave_trp1 = tot_trp1 / len(shortest_paths)
    ave_trp2 = tot_trp2 / len(shortest_paths)
    ave_fus = tot_fus / len(shortest_paths)
    ave_fis = tot_fis / len(shortest_paths)



    #number_of_operations = ave_b_trl + ave_fis + ave_fus + ave_inv + ave_trp1 + (ave_trp2 * 2) + ave_u_trl
    number_of_operations = len(shortest_paths[0])-2
    results.append(
        [number_of_operations, len(shortest_paths), ave_inv, ave_trp1, ave_trp2, ave_b_trl, ave_u_trl, ave_fus,
         ave_fis])

    network.clear()

    paths_operations = []
    for element in Paths_state:
        path_operations = [y for (x, y) in element]

        paths_operations.append(path_operations)

    solution_operations = [d for (c, d) in sorting_scenario]

    if solution_operations in paths_operations:
        number_of_solutions_found+=1





simulation_average = []
average_number_of_operations = 0
average_number_of_paths = 0
average_number_of_inv = 0
average_number_of_trp1 = 0
average_number_of_trp2 = 0
average_number_of_b_trl = 0
average_number_of_u_trl = 0
average_number_of_fus = 0
average_number_of_fis = 0
for element in results:
    average_number_of_operations += element[0]
    average_number_of_paths += element[1]
    average_number_of_inv += element[2]
    average_number_of_trp1 += element[3]
    average_number_of_trp2 += element[4]
    average_number_of_b_trl += element[5]
    average_number_of_u_trl += element[6]
    average_number_of_fus += element[7]
    average_number_of_fis += element[8]

average_number_of_operations = average_number_of_operations / number_of_simulations
average_number_of_paths = average_number_of_paths / number_of_simulations
average_number_of_inv = average_number_of_inv / number_of_simulations
average_number_of_trp1 = average_number_of_trp1 / number_of_simulations
average_number_of_trp2 = average_number_of_trp2 / number_of_simulations
average_number_of_b_trl = average_number_of_b_trl / number_of_simulations
average_number_of_u_trl = average_number_of_u_trl / number_of_simulations
average_number_of_fus = average_number_of_fus / number_of_simulations
average_number_of_fis = average_number_of_fis / number_of_simulations

ave_number_per_op = [average_number_of_inv, average_number_of_trp1, average_number_of_b_trl, average_number_of_u_trl, average_number_of_fis, average_number_of_fus, average_number_of_trp2]
num_per_op_ratios = []
tot_ave_op_num = 0
for element in ave_number_per_op:
    tot_ave_op_num+= element

for element in ave_number_per_op:
    ratio = element/tot_ave_op_num*average_number_of_operations
    num_per_op_ratios.append(ratio)



t1 = time.time()
print('##########################################################################################################')
print()
print('time: ', t1-t0)
print()
print('average number of solutions: ', average_number_of_paths)
print()
print('average solution length: ' ,average_number_of_operations)
print()
print('percentage time true solution was found: ',(number_of_solutions_found/number_of_simulations)*100)

print()
print('ratios as: inv:trp1:b_trl:u_trl:fis:fus:trp2')
print(num_per_op_ratios)


print('##########################################################################################################')
