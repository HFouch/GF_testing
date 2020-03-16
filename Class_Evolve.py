import random
from Class_extremities_and_adjacencies import Extremities_and_adjacencies

class Node_evolve:

    def __init__(self, state=None):
        self.state = state
        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]

    def find_next_extremity(self, current, next_extremity):

        if current[0] == next_extremity:
            if current[1] % 1 == 0:
                next = current[1] + 0.5
            else:
                next = current[1] - 0.5
        else:
            if current[0] % 1 == 0:
                next = current[0] + 0.5
            else:
                next = current[0] - 0.5
        return next

    def find_next_adjacency(self, next_extremity, chromosome, not_telomeres):
        for element in not_telomeres:
            if element[0] == next_extremity or element[1] == next_extremity:
                current = element
                chromosome.append(current)
                not_telomeres.remove(current)
                next_extremity = Node_evolve.find_next_extremity(self, current, next_extremity)
                return next_extremity, chromosome, not_telomeres

        return [next_extremity]

    def find_adjacency_cycle(self, next_extremity, chromosome, not_telomeres):

        next_adjacency = Node_evolve.find_next_adjacency(self, next_extremity, chromosome, not_telomeres)

        while len(next_adjacency) != 1:
            next_extremity = next_adjacency[0]
            next_adjacency = Node_evolve.find_next_adjacency(self, next_extremity, chromosome, not_telomeres)
        else:
            next_extremity = next_adjacency[0]

            return next_extremity, chromosome, not_telomeres

    def find_chromosomes(self, adjacencies):

        telomeres = [element for element in adjacencies if type(element) is not tuple]
        not_telomeres = [element for element in adjacencies if type(element) is tuple]
        linear_chromosomes = []
        circular_chromosomes = []
        chromosome = []
        i = 0

        # find linear chromosomes
        while len(telomeres) > 0:

            i += 1
            current = telomeres[0]

            telomeres.remove(current)
            chromosome.append(current)

            if current % 1 == 0:
                next_extremity = current + 0.5
            else:
                next_extremity = current - 0.5

            # if single gene chromosome
            if next_extremity in telomeres:
                current = next_extremity

                telomeres.remove(current)
                chromosome.append(current)
                linear_chromosomes.append(chromosome)
                chromosome = []

            # else find adjacency cycle
            else:
                adjacency_cycle = Node_evolve.find_adjacency_cycle(self, next_extremity, chromosome, not_telomeres)
                next_extremity = adjacency_cycle[0]

                if next_extremity in telomeres:
                    current = next_extremity
                    telomeres.remove(current)
                    chromosome.append(current)
                    linear_chromosomes.append(chromosome)
                    chromosome = []

        # find circular chromosomes
        while len(not_telomeres) > 0:
            current = not_telomeres[0]
            not_telomeres.remove(current)
            chromosome.append(current)

            # find next extremity:
            if current[0] % 1 == 0:
                next_extremity = current[0] + 0.5
            else:
                next_extremity = current[0] - 0.5

            # if it is a single gene chromosome:
            if next_extremity == current[1]:
                circular_chromosomes.append(chromosome)
                chromosome = []

            # go find adjacency cycle
            else:
                adjacency_cycle = Node_evolve.find_adjacency_cycle(self, next_extremity, chromosome, not_telomeres)
                next_extremity = adjacency_cycle[0]

                # if at end of circular chromosome
                if next_extremity == chromosome[0][1]:
                    circular_chromosomes.append(chromosome)
                    chromosome = []

        return linear_chromosomes, circular_chromosomes

    def evolve(self, num_inv=0, num_fus=0, num_fis=0, num_b_trl=0, num_u_trl=0, num_trp1=0, num_trp2=0):

        scenario = []
        total = num_b_trl+num_fis+num_fus+num_inv+num_trp1+num_trp2+num_u_trl
        inv_n = num_inv
        fus_n = num_fus
        fis_n = num_fis
        b_trl_n = num_b_trl
        u_trl_n = num_u_trl
        trp1_n = num_trp1
        trp2_n = num_trp2

        operations = ['inv' , 'trp1', 'trp2', 'b_trl', 'u_trl', 'fis', 'fus']

        while total > 0:
            operation = random.choice(operations)

            if operation == 'inv':
                if inv_n > 0:
                    execution = self.inversion()

                    scenario.append(execution)

                    inv_n -=1
                    total -=1
                else:
                    pass

            elif operation == 'trp1':
                if trp1_n > 0:
                    execution = self.transposition1()
                    for element in execution:
                        scenario.append(element)
                    trp1_n-=1
                    total-=1
                else:
                    pass

            elif operation == 'trp2':
                if trp2_n > 0:
                    self.transposition2()
                    trp2_n-=1
                    total-=1
                else:
                    pass

            elif operation == 'b_trl':
                if b_trl_n > 0:
                    execution = self.balanced_translocation()
                    scenario.append(execution)
                    b_trl_n-=1
                    total-=1
                else:
                    pass

            elif operation == 'u_trl':
                if u_trl_n >0:
                    execution =self.unbalanced_translocation()
                    scenario.append(execution)
                    u_trl_n-=1
                    total-=1
                else:
                    pass

            elif operation=='fis':
                if fis_n >0 :
                    execution = self.fission()
                    scenario.append(execution)
                    fis_n-=1
                    total-=1
                else:
                    pass

            elif operation=='fus':
                if fus_n>0:
                    execution = self.fusion()
                    scenario.append(execution)
                    fus_n-=1
                    total-=1
                else:
                    pass

        return scenario


    def inversion(self):

        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]
        chrm_num = random.randint(0, len(self.linear_chromosomes)-1)

        # if it is a single gene or two gene chromosome, choose another one
        while len(self.linear_chromosomes[chrm_num]) < 4:
            chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

        adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]

        adj1_num = random.randint(0, len(adjacencies)-1)
        adj2_num = random.randint(0, len(adjacencies)-1)

        while adj1_num == adj2_num:
            adj2_num = random.randint(0, len(adjacencies)-1)

        adj1 = adjacencies[adj1_num]
        adj2 = adjacencies[adj2_num]

        #order extremities in the adjacency tuple
        if adj1[0] < adj2[0]:
            new_adj1 = (adj1[0], adj2[0])
        else:
            new_adj1 = (adj2[0], adj1[0])

        if adj1[1] < adj2[1]:
            new_adj2 = (adj1[1], adj2[1])
        else:
            new_adj2 = (adj2[1], adj1[1])


        #perfrom operation
        self.state.remove(adj1)
        self.state.remove(adj2)
        self.state.append(new_adj1)
        self.state.append(new_adj2)

        if adj1[0] < adj2[0]:
            operationA = (adj1, adj2)
        else:
            operationA = (adj2, adj1)

        if new_adj1[0] < new_adj2[0]:
            operationB = (new_adj1, new_adj2)
        else:
            operationB = (new_adj2, new_adj1)

        operation = ((operationA, operationB), 'inv')


        intermediate = self.state[:]

        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]

        while len(self.circular_chromosomes) > 0:

            self.state.append(adj1)
            self.state.append(adj2)
            self.state.remove(new_adj1)
            self.state.remove(new_adj2)

            get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
            self.linear_chromosomes = get_chromosomes[0]
            self.circular_chromosomes = get_chromosomes[1]

            chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)
            while len(self.linear_chromosomes[chrm_num]) < 4:
                chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)
            adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]

            adj1_num = random.randint(0, len(adjacencies) - 1)
            adj2_num = random.randint(0, len(adjacencies) - 1)

            while adj1_num == adj2_num:
                adj2_num = random.randint(0, len(adjacencies) - 1)

            adj1 = adjacencies[adj1_num]
            adj2 = adjacencies[adj2_num]

            #order extremities in the adjacency tuple
            if adj1[0] < adj2[0]:
                new_adj1 = (adj1[0], adj2[0])
            else:
                new_adj1 = (adj2[0], adj1[0])

            if adj1[1] < adj2[1]:
                new_adj2 = (adj1[1], adj2[1])
            else:
                new_adj2 = (adj2[1], adj1[1])

            # perfrom operation
            self.state.remove(adj1)
            self.state.remove(adj2)
            self.state.append(new_adj1)
            self.state.append(new_adj2)

            get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
            self.linear_chromosomes = get_chromosomes[0]
            self.circular_chromosomes = get_chromosomes[1]

            if adj1[0] < adj2[0]:
                operationA = (adj1, adj2)
            else:
                operationA = (adj2, adj1)

            if new_adj1[0] < new_adj2[0]:
                operationB = (new_adj1, new_adj2)
            else:
                operationB = (new_adj2, new_adj1)

            operation = ((operationA, operationB),'inv')
            intermediate = self.state[:]

        return (operation, intermediate)



    def transposition1(self):

        # exsision and circularization
        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]
        chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

        #if it is a single gene or two gene chromosome, choose another one
        while len(self.linear_chromosomes[chrm_num]) < 4:
            chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

        adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]

        adj1_num = random.randint(0, len(adjacencies) - 1)
        adj2_num = random.randint(0, len(adjacencies) - 1)

        while adj1_num == adj2_num:
            adj2_num = random.randint(0, len(adjacencies) - 1)

        adj1 = adjacencies[adj1_num]
        adj2 = adjacencies[adj2_num]




        # order extremities in the adjacency tuple
        if adj1[0] < adj2[1]:
            new_adj1 = (adj1[0], adj2[1])
        else:
            new_adj1 = (adj2[1], adj1[0])

        if adj1[1] < adj2[0]:
            new_adj2 = (adj1[1], adj2[0])
        else:
            new_adj2 = (adj2[0], adj1[1])

        # perfrom operation
        self.state.remove(adj1)
        self.state.remove(adj2)
        self.state.append(new_adj1)
        self.state.append(new_adj2)

        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]

        if adj1[0] < adj2[0]:
            operationA = (adj1, adj2)
        else:
            operationA = (adj2, adj1)

        if new_adj1[0] < new_adj2[0]:
            operationB = (new_adj1, new_adj2)
        else:
            operationB = (new_adj2, new_adj1)

        operation = ((operationA, operationB), 'trp0')


        intermediate = self.state[:]

        if len(self.circular_chromosomes) == 0:
            self.state.append(adj1)
            self.state.append(adj2)
            self.state.remove(new_adj1)
            self.state.remove(new_adj2)

            #order extremities in the adjacency tuple
            if adj1[0] < adj2[0]:
                new_adj1 = (adj1[0], adj2[0])
            else:
                new_adj1 = (adj2[0], adj1[0])

            if adj1[1] < adj2[1]:
                new_adj2 = (adj1[1], adj2[1])
            else:
                new_adj2 = (adj2[1], adj1[1])

            #perform operation
            self.state.remove(adj1)
            self.state.remove(adj2)
            self.state.append(new_adj1)
            self.state.append(new_adj2)

            get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
            self.linear_chromosomes = get_chromosomes[0]
            self.circular_chromosomes = get_chromosomes[1]

            if adj1[0] < adj2[0]:
                operationA = (adj1, adj2)
            else:
                operationA = (adj2, adj1)

            if new_adj1[0] < new_adj2[0]:
                operationB = (new_adj1, new_adj2)
            else:
                operationB = (new_adj2, new_adj1)

            operation = ((operationA, operationB), 'trp0')


            intermediate = self.state[:]



        if new_adj1 in self.circular_chromosomes[0]:
            join = new_adj1
            excision = new_adj2
        else:
            join = new_adj2
            excision = new_adj1

        operation1 = operation
        intermediate1 = intermediate

        #decircularization and reinsertion
        adj1 = join
        chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

        #if it is a single gene chromosome, choose another chromosome
        while len(self.linear_chromosomes[chrm_num]) == 2:
            chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

        if len(self.linear_chromosomes[chrm_num]) == 3:
            adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]
            adj2 = adjacencies[0]

        else:
            adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]
            adj2_num = random.randint(0, len(adjacencies) - 1)
            adj2 = adjacencies[adj2_num]

        while adj2 == excision:
            chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

            # if it is a single gene chromosome, choose another chromosome
            while len(self.linear_chromosomes[chrm_num]) == 2:
                chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

            if len(self.linear_chromosomes[chrm_num]) == 3:
                adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]
                adj2 = adjacencies[0]

            else:
                adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]
                adj2_num = random.randint(0, len(adjacencies) - 1)
                adj2 = adjacencies[adj2_num]

        #order extremities in adjacency tuple
        if adj1[0] < adj2[1]:
            new_adj1 = (adj1[0], adj2[1])
        else:
            new_adj1 = (adj2[1], adj1[0])

        if adj1[1] < adj2[0]:
            new_adj2 = (adj1[1], adj2[0])
        else:
            new_adj2 = (adj2[0], adj1[1])

        # perfrom operation
        self.state.remove(adj1)
        self.state.remove(adj2)
        self.state.append(new_adj1)
        self.state.append(new_adj2)

        if adj1[0] < adj2[0]:
            operationA = (adj1, adj2)
        else:
            operationA = (adj2, adj1)

        if new_adj1[0] < new_adj2[0]:
            operationB = (new_adj1, new_adj2)
        else:
            operationB = (new_adj2, new_adj1)

        operation2 = ((operationA, operationB), 'trp1')
        intermediate2 = self.state[:]

        return [(operation1, intermediate1), (operation2, intermediate2)]




    def transposition2(self):

        # exsision and circularization
        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]
        chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

        # if it is a single gene or two gene chromosome, choose another one
        while len(self.linear_chromosomes[chrm_num]) < 6:
            chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)


        adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]

        adj1_num = random.randint(0, len(adjacencies) - 1)
        adj2_num = random.randint(0, len(adjacencies) - 1)

        # To ensure single gene circular chromosomes are not formed because then there is only one adj so cut and not != join
        while adj1_num == adj2_num or (int(adjacencies[adj2_num][0])) == int(adjacencies[adj1_num][0]) or (
        int(adjacencies[adj2_num][0])) == int(adjacencies[adj1_num][1]) or (int(adjacencies[adj2_num][1])) == int(
                adjacencies[adj1_num][0]) or (int(adjacencies[adj2_num][1])) == int(adjacencies[adj1_num][1]):
            adj2_num = random.randint(0, len(adjacencies) - 1)

        adj1 = adjacencies[adj1_num]
        adj2 = adjacencies[adj2_num]

        # order extremities in adjacency tuple
        if adj1[0] < adj2[1]:
            new_adj1 = (adj1[0], adj2[1])
        else:
            new_adj1 = (adj2[1], adj1[0])

        if adj1[1] < adj2[0]:
            new_adj2 = (adj1[1], adj2[0])
        else:
            new_adj2 = (adj2[0], adj1[1])

        # perfrom operation
        self.state.remove(adj1)
        self.state.remove(adj2)
        self.state.append(new_adj1)
        self.state.append(new_adj2)

        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]


        # if the above does not result in the formation of a cicular chromosome join the adjacencies differently
        # NTS you can actually take the above out and just start with this while loop...
        while len(self.circular_chromosomes) == 0 :
            self.state.append(adj1)
            self.state.append(adj2)
            self.state.remove(new_adj1)
            self.state.remove(new_adj2)

            if adj1[0] < adj2[0]:
                new_adj1 = (adj1[0], adj2[0])
            else:
                new_adj1 = (adj2[0], adj1[0])

            if adj1[1] < adj2[1]:
                new_adj2 = (adj1[1], adj2[1])
            else:
                new_adj2 = (adj2[1], adj1[1])

            self.state.remove(adj1)
            self.state.remove(adj2)
            self.state.append(new_adj1)
            self.state.append(new_adj2)

            get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
            self.linear_chromosomes = get_chromosomes[0]
            self.circular_chromosomes = get_chromosomes[1]

        if new_adj1 in self.circular_chromosomes[0]:
            join = new_adj1
            excision = new_adj2
        else:
            join = new_adj2
            excision = new_adj1

        # decircularization and reinsertion
        adj1_num = random.randint(0, len(self.circular_chromosomes[0])-1)
        adj1 = self.circular_chromosomes[0][adj1_num]

        while adj1 == join:
            adj1_num = random.randint(0, len(self.circular_chromosomes[0])-1)
            adj1 = self.circular_chromosomes[0][adj1_num]

        chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

        while len(self.linear_chromosomes[chrm_num]) == 2:
            chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

        if len(self.linear_chromosomes[chrm_num]) == 3:
            adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]
            adj2 = adjacencies[0]

        else:
            adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]
            adj2_num = random.randint(0, len(adjacencies) - 1)
            adj2 = adjacencies[adj2_num]

        while adj2 == excision:
            chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

            # if it is a single gene chromosome, choose another chromosome
            while len(self.linear_chromosomes[chrm_num]) == 2:
                chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)

            if len(self.linear_chromosomes[chrm_num]) == 3:
                adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]
                adj2 = adjacencies[0]

            else:
                adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]
                adj2_num = random.randint(0, len(adjacencies) - 1)
                adj2 = adjacencies[adj2_num]

        #oreder extremities in the adjacency tuple
        if adj1[0] < adj2[1]:
            new_adj1 = (adj1[0], adj2[1])
        else:
            new_adj1 = (adj2[1], adj1[0])

        if adj1[1] < adj2[0]:
            new_adj2 = (adj1[1], adj2[0])
        else:
            new_adj2 = (adj2[0], adj1[1])

        # perfrom operation
        self.state.remove(adj1)
        self.state.remove(adj2)
        self.state.append(new_adj1)
        self.state.append(new_adj2)

    def balanced_translocation(self):

        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]
        chrm_num_1 = random.randint(0, len(self.linear_chromosomes) - 1)

        while len(self.linear_chromosomes[chrm_num_1]) < 3:
            chrm_num_1 = random.randint(0, len(self.linear_chromosomes) - 1)

        chrm_num_2 = random.randint(0, len(self.linear_chromosomes) - 1)
        while chrm_num_1==chrm_num_2 or len(self.linear_chromosomes[chrm_num_2]) < 3:
            chrm_num_2 = random.randint(0, len(self.linear_chromosomes) - 1)

        adjacencies_chrm1 = [element for element in self.linear_chromosomes[chrm_num_1] if type(element) is tuple]
        adjacencies_chrm2 = [element for element in self.linear_chromosomes[chrm_num_2] if type(element) is tuple]
        adj1_num = random.randint(0, len(adjacencies_chrm1) - 1)
        adj2_num = random.randint(0, len(adjacencies_chrm2) - 1)

        adj1 = adjacencies_chrm1[adj1_num]
        adj2 = adjacencies_chrm2[adj2_num]

        # order extremities in the adjacency tuple
        if adj1[0] < adj2[1]:
            new_adj1 = (adj1[0], adj2[1])
        else:
            new_adj1 = (adj2[1], adj1[0])

        if adj1[1] < adj2[0]:
            new_adj2 = (adj1[1], adj2[0])
        else:
            new_adj2 = (adj2[0], adj1[1])



        # perfrom operation
        self.state.remove(adj1)
        self.state.remove(adj2)
        self.state.append(new_adj1)
        self.state.append(new_adj2)

        if adj1[0] < adj2[0]:
            operationA = (adj1, adj2)
        else:
            operationA = (adj2, adj1)

        if new_adj1[0] < new_adj2[0]:
            operationB = (new_adj1, new_adj2)
        else:
            operationB = (new_adj2, new_adj1)

        operation = ((operationA, operationB), 'b_trl')
        intermediate = self.state[:]

        return (operation, intermediate)

    def unbalanced_translocation(self):

        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]

        chrm_num_1 = random.randint(0, len(self.linear_chromosomes) - 1)

        # chromosome 1 should not be a single gene chromosome otherwise the operation would amount to a fusion
        while len(self.linear_chromosomes[chrm_num_1]) < 3:
            chrm_num_1 = random.randint(0, len(self.linear_chromosomes) - 1)

        chrm_num_2 = random.randint(0, len(self.linear_chromosomes) - 1)
        while chrm_num_1==chrm_num_2:
            chrm_num_2 = random.randint(0, len(self.linear_chromosomes) - 1)


        adjacencies_chrm1 = [element for element in self.linear_chromosomes[chrm_num_1] if type(element) is tuple]
        telomeres_chrm2 = [element for element in self.linear_chromosomes[chrm_num_2] if type(element) is not tuple]
        adj1_num = random.randint(0, len(adjacencies_chrm1) - 1)
        adj2_num = random.randint(0, len(telomeres_chrm2) - 1)

        adj1 = adjacencies_chrm1[adj1_num]
        telo2 = telomeres_chrm2[adj2_num]

        # order extremities in the adjacency tuple
        if telo2 < adj1[0]:
            new_adj1 = (telo2, adj1[0])
        else:
            new_adj1 = (adj1[0], telo2)

        # perfrom operation
        self.state.remove(adj1)
        self.state.remove(telo2)
        self.state.append(new_adj1)
        self.state.append(adj1[1])

        operation = (((adj1, telo2), (new_adj1, adj1[1])), 'u_trl')
        intermediate = self.state[:]

        return (operation, intermediate)

    def fission(self):

        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]
        chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)
        adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]

        while len(adjacencies) == 0:
            chrm_num = random.randint(0, len(self.linear_chromosomes) - 1)
            adjacencies = [element for element in self.linear_chromosomes[chrm_num] if type(element) is tuple]

        adj_num = random.randint(0, len(adjacencies) - 1)

        adj = adjacencies[adj_num]

        self.state.remove(adj)
        self.state.append(adj[0])
        self.state.append(adj[1])

        operation = ((adj, adj[0], adj[1]), 'fis')
        intermediate = self.state[:]
        return (operation, intermediate)


    def fusion(self):

        get_chromosomes = Node_evolve.find_chromosomes(self, self.state)
        self.linear_chromosomes = get_chromosomes[0]
        self.circular_chromosomes = get_chromosomes[1]
        chrm_num_1 = random.randint(0, len(self.linear_chromosomes) - 1)
        chrm_num_2 = random.randint(0, len(self.linear_chromosomes) - 1)
        while chrm_num_1 == chrm_num_2:
            chrm_num_2 = random.randint(0, len(self.linear_chromosomes) - 1)

        telomeres_chrm1 = [element for element in self.linear_chromosomes[chrm_num_1] if type(element) is not tuple]
        telomeres_chrm2 = [element for element in self.linear_chromosomes[chrm_num_2] if type(element) is not tuple]
        telo1_num = random.randint(0, len(telomeres_chrm1) - 1)
        telo2_num = random.randint(0, len(telomeres_chrm2) - 1)

        telo1 = telomeres_chrm1[telo1_num]
        telo2 =telomeres_chrm2[telo2_num]

        if telo1 < telo2:
            new_adj = (telo1, telo2)
        else:
            new_adj = (telo2, telo1)

        # perfrom operation
        self.state.remove(telo1)
        self.state.remove(telo2)
        self.state.append(new_adj)

        operation = ((telo1, telo2, new_adj), 'fus')
        intermediate = self.state[:]

        return (operation, intermediate)







