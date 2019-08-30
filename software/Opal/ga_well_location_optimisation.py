# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#Authors: P. Salinas
import operator

from opal_classes import *
from deap import tools, base, creator, cma
from importance_map_tools import *
import libspud
from exodus2gmsh import convertExodusII2MSH

#Read Opal options, for the time being just so we can call the same subroutines
# read in the opal options
opal_options, fwd_options, nirom_options = get_opal_options()
# get xml extension (mpml or flml)
xml_extension = get_xml_extension(opal_options.input_file)
# load the options for the forward model from xml file
libspud.load_options(opal_options.input_file)
# Production file
prod_file = libspud.get_option('/simulation_name') + "_outfluxes.csv"
# Final time to extract from the csv file
final_time = libspud.get_option('/timestepping/finish_time')
#Retrieve extra variables required to run the optimisation
MIN = opal_options.ga_variables[0].min_limit
MAX = opal_options.ga_variables[0].max_limit
spatial_precision = int((abs(MIN) + abs(MAX)) * opal_options.precision)
##Global variable to be used for gradient convergence
previous_convergence = [0.0, 0.0]
##Geothermal if we have temperature field
geothermal = libspud.have_option('material_phase['+str(0)+']/scalar_field::Temperature')
##Already visited studied
explored_locations = []

## TODO: 1) FIND OPTIMAL COMBINATION OF ALGORITHMS, FOR EXAMPLE DIFFERENT CROSSOVER?
## TODO: 1.5) IMPLEMENT SOME SORT OF HILL CLIMBER/GRADIENT DESCENT METHOD FOR CROSSOVER?


#Call the creators in the root of the program to ensure it works with scoop and in parallel
if  opal_options.ga_Minimise:
    creator.create("FitnessMin", base.Fitness, weights=(-1,))
    creator.create("Individual", list, fitness=creator.FitnessMin)
else:
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

def analyse_results(instance):
    cwd = os.getcwd()
    foldername = get_folder_name(instance)
    try:
        timeList, prod, inject, prod_temp, inject_temp = get_productions_and_injections(prod_file, geothermal, foldername, cwd)
    except:
        os.chdir(cwd)
        print "WARNING: Failed to load production files."
        timeList = 0.
        prod =0.
        inject=0.
        prod_temp=0.
        inject_temp =0.


    walltime= 0
    try:
        walltime = get_walltime(foldername, cwd)
    except:
        os.chdir(cwd)
        print "WARNING: Failed to obtain walltime"

    return timeList, prod, inject, prod_temp, inject_temp, walltime


def get_folder_name(instance):
    #To try to compact as much as possible the length of the folder while being unique
    #we convert the number representation to base64 from base10
    def get_digit(d):
        if 0 <= d <= 9:
            # 0 - 9
            c = 48 + d
        elif 10 <= d <= 35:
            # A - Z
            c = 55 + d
        elif 36 <= d <= 61:
            # a - z
            c = 61 + d
        elif d == 62:
            # -
            c = 45
        elif d == 63:
            # +
            c = 43
        return chr(c)

    def encode(n):
        out = []
        while n:
            n, r = n // 64, n % 64
            out.append(get_digit(r))
        while len(out) < 6:
            out.append('0')
        return ''.join(out)

    #To concatenate all the numbers we need them to be positive (no minus sign between them)
    ensure_positve = 0
    for var in opal_options.ga_variables:
        ensure_positve = max(abs(var.min_limit), ensure_positve)
    instance_bak = convert_instance(instance)
    foldername =  ''
    for val in instance_bak:
        foldername +=  str(val+ensure_positve)
    #We use base 64 representation of the numbers to try to reduce the length of the folder
    foldername = opal_options.input_file[:-5]+"_"+ encode(int(foldername))

    return foldername

def modify_experiment(input_file, foldername, instance, cwd):

    def modify_input_file(INPUT_FILE):
        # Read in the file
        with open(INPUT_FILE, 'r') as file:
            filedata = file.read()
        # Convert to physical values
        instance_bak = convert_instance(instance)
        ##Substitute a given pattern by the instance value
        for i in range(opal_options.ga_locations_to_study):
            for j in range(len(opal_options.ga_variables)):
                normaliser = opal_options.ga_variables[j].normaliser

                filedata = filedata.replace(opal_options.ga_variables[j].variable_pattern + str(i + 1),
                             str( float(instance_bak[j + len(opal_options.ga_variables) * i])/float(normaliser) ))

        # Overwrite file
        with open(INPUT_FILE, 'w') as file:
            file.write(filedata)


    ##Get into the folder
    os.chdir(foldername)
    #Modify trelis input file if requested
    if len(input_file)>1:
        modify_input_file(input_file)
        #Run trelis to create new exodusII file given the instance
        string = opal_options.trelis_path + " -nographics -nojournal -batch " + opal_options.trelis_input_file
        os.system(string)
        #Now convert the output .e file into .msh
        meshfile = libspud.get_option('/geometry/mesh::CoordinateMesh/from_file/file_name')
        convertExodusII2MSH(meshfile)
        #Decompose the mesh is required
        if opal_options.MPI_runs > 1: os.system("fldecomp -n " + str(opal_options.MPI_runs) + " " + meshfile)
        # Now convert the outp
    #Modify the mpml file if requested
    if opal_options.optimise_input: modify_input_file(opal_options.input_file)


    ##Return to original path
    os.chdir(cwd)
    return

def get_walltime(foldername, cwd):
    from fluidity_tools import stat_parser as stat
    walltime =1e50
    ##Get into the folder
    os.chdir(foldername)
    output_name = libspud.get_option('/simulation_name')
    walltime = stat('./' + output_name + '.stat')["ElapsedWallTime"]["value"][-1]
    ##Return to original path
    os.chdir(cwd)
    return walltime

def get_productions_and_injections(prod_file, geothermal, foldername, cwd):
    import csv
    nPhases = libspud.option_count('/material_phase')
    ##Get into the folder
    os.chdir(foldername)
    String_id = "-S"
    String_prod = "- Volume rate"

    #Replace spaces by commas to make it consistent
    opal_options.producer_ids.replace(' ', ',')
    opal_options.injector_ids.replace(' ', ',')
    #Create a list with the integers
    producer_list = [int(e) if e.isdigit() else e for e in opal_options.producer_ids.split(',')]
    injector_list = [int(e) if e.isdigit() else e for e in opal_options.injector_ids.split(',')]

    with open(prod_file, 'rb') as csvfile:
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        counter = -1 #To account for the header
        for row in datareader:
            try:
                counter+=1
            except:
                continue

    timeList = np.zeros(counter)
    prod = np.zeros((nPhases, counter))
    inject = np.zeros((nPhases, counter))
    prod_temp = np.zeros((nPhases, counter))
    inject_temp = np.zeros((nPhases, counter))
    counter = 0
    with open(prod_file, 'rb') as csvfile:
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        prod_cols =[]
        inject_cols = []

        header = True
        previousT = 0.0
        for row in datareader:
            try:
                if header:
                    #Find the positions of interest of the producers
                    for ids in producer_list:
                        for i in range(len(row)):
                            if (String_id+str(ids)+String_prod) in row[i]:
                                prod_cols.append(i)
                    #Find the positions of interest of the injectors
                    for ids in injector_list:
                        for i in range(len(row)):
                            if (String_id+str(ids)+String_prod) in row[i]:
                                inject_cols.append(i)
                    header = False
                    continue

                phase_prod = 0
                phase_inject = 0
                timeList[counter] = float(row[0])
                # Calculate deltaT from absolute times
                deltaT = float(row[0]) - previousT
                previousT = float(row[0])

                for i in range(len(row)):
                    # Update production
                    for j in range(len(prod_cols)):
                        if i == prod_cols[j]:
                            if geothermal:
                                #For geothermal I need to do temperature * production * rho * Cp
                                #CURRENTLY NOT INCLUDED Cp nor RHO
                                #First temperature, so we can calculate the production in this time level
                                prod_temp[phase_prod, counter] = abs(float(row[i]) * deltaT * float(row[i+nPhases*2]))
                            prod[phase_prod, counter] = abs(float(row[i]) * deltaT)
                            phase_prod += 1
                            if phase_prod == nPhases: phase_prod = 0

                    # Update Injection
                    for j in range(len(inject_cols)):
                        if i == inject_cols[j]:
                            if geothermal:
                                #For geothermal I need to do temperature * production * rho * Cp
                                #CURRENTLY NOT INCLUDED Cp nor RHO
                                #First temperature, so we can calculate the production in this time level
                                inject_temp[phase_inject, counter] = abs(float(row[i]) * deltaT * float(row[i+nPhases*2]))

                            inject[phase_inject, counter] = abs(float(row[i]) * deltaT)
                            phase_inject += 1
                            if phase_inject == nPhases: phase_inject = 0
                counter +=1
            except:
                continue
    ##Return to original path
    os.chdir(cwd)
    return timeList, prod, inject, prod_temp, inject_temp

#Function that evaluates the functional
#A penalty function can easily added here by providing a bad value if a requirement is not fulfilled
def fitness(instance):
    timeList, prod, inject, prod_temp, inject_temp, walltime = analyse_results(instance)  # Only analyse
    if len(opal_options.ga_fitness_functional) > 1:
        try:
            exec(opal_options.ga_fitness_functional)
        except:
            print "#############################################################################"
            print "ERROR evaluating the user functional. Please check the functional introduced."
            print "#############################################################################"
            exit()
    else:
        if geothermal:
            #Need to add rho_cp_to the mix
            val = np.sum(prod_temp)
        else:
            val = np.sum(prod)
    return val,

def eaAlgorithm_by_steps(pop, toolbox, CXPB, MUTPB, NGEN, halloffame):

    # Evaluate the entire population
    # Run the simulations in parallel first
    # for i in range(0, len(pop), opal_options.number_processors):
    #     run_in_batches(pop[i:i + opal_options.number_processors])
    counter = 0
    for i in range(0, len(pop), opal_options.number_processors):
        counter = +run_in_batches(pop[counter:])
        if counter >= len(pop): break
    # Evaluate
    fitnesses = map(toolbox.evaluate, pop)
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    for g in range(NGEN):
        if opal_options.ga_CMA:
            #If CMA then generate new population like this
            offspring = toolbox.generate()
        else: #Otherwise mate and mutation
            # Select the next generation individuals
            if opal_options.ga_evol_algorithm == 1:
                offspring = toolbox.select(pop, len(pop))
            else:
                offspring = toolbox.select(pop, opal_options.ga_lambda)
            # Clone the selected individuals
            offspring = map(toolbox.clone, offspring)
            #Now we are applying a VarAnd method since we do, mate AND mutation
            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < CXPB:
                    try: #By default the internal method, but sometimes it has failed
                        toolbox.mate(child1, child2)
                    except:#if it fails, use our version that is prepared for single variables
                        child1, child2 = One_point_crossover(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values
            for mutant in offspring:
                #Deciding if a mutation must be performed is decided here
                if random.random() < MUTPB:
                    toolbox.mutate(mutant)
                    del mutant.fitness.values
        # Evaluate the individuals with an invalid fitness
        # Run the simulations in parallel first
        counter = 0
        for i in range(0, len(offspring), opal_options.number_processors):
            counter = +run_in_batches(offspring[counter:])
            if counter >= len(offspring): break


            # run_in_batches(offspring[i:i + opal_options.number_processors])
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        #Update hall of fame
        halloffame.update(pop)
        best_so_far = halloffame[0]
        best_distance_so_far = fitness(best_so_far)
        # Check current state of convergence
        terminate(pop,best_distance_so_far, g)
        #Print best result so far:
        best_so_far_bak = convert_instance(best_so_far)
        print "Best simulation at generation " + str(g+1) + ": ", best_so_far_bak, "; Produces: ", best_distance_so_far[0]
        #Update as well best hall of fame
        create_output_file(halloffame)
        #create_best_folders(halloffame)

        # The population is entirely replaced by the offspring
        if opal_options.ga_evol_algorithm == 2:
            pop[:] = toolbox.select(offspring, opal_options.mu)
        elif opal_options.ga_evol_algorithm == 3:
            pop[:] = toolbox.select(pop + offspring, opal_options.mu)
        else:
            pop[:] = offspring

        if opal_options.ga_CMA: toolbox.update(pop)  # Update strategy


def One_point_crossover(ind1, ind2):
    size = min(len(ind1), len(ind2))
    if size == 1:
        point =1
    else:
        point = random.randint(1, size - 1) #if size == 1 this gives serious problems
    ind1[point:], ind2[point:] = ind2[point:], ind1[point:]
    return ind1, ind2


def run_in_batches(sub_pop):
    from subprocess import Popen
    from copy import deepcopy

    #For the time being, do not rerun simulations
    sub_population = deepcopy(sub_pop)
    pop_to_run = []
    #Check that the locations haven't been checked already
    # for instance in sub_population:
    counter = 0
    for instance in sub_population:
        foldername = get_folder_name(instance)
        #We check that a simulation has been run if there is a production file already in place
        if os.path.isfile(foldername+"/"+prod_file):
            explored_locations.append(instance)

        if instance in explored_locations:
            try:
                pop_to_run.remove(instance)
            except:
                pass
        else:
            counter += 1
            pop_to_run.append(instance)
            explored_locations.append(instance)
        #Once we fill up the total number of CPUs, we continue
        if len(pop_to_run) == opal_options.number_processors: break


    #Create folders
    folders = []
    input_file = opal_options.trelis_input_file
    # if libspud.have_option('/opal_operation/ga_optimisation/Trelis_integration'):
    #     input_file = opal_options.trelis_input_file#opal_options.input_file[:-5]
    cwd = os.getcwd()
    #Generate the sub-folders and create all the input files
    for instance in pop_to_run:
        cwd = os.getcwd()
        have_wells = True
        foldername = get_folder_name(instance)
        folders.append(foldername)
        existed = create_folder_to_run(opal_options, cwd, have_wells, foldername)
        if existed:#When rerunning simulations with folders already in place
            sub_population.remove(instance)
            continue
        modify_experiment(input_file, foldername, instance, cwd)

    if len(pop_to_run)>=1 and len(folders)>0:
        # Prepare the commands to run the batch of simulations
        commands = []
        # Run the simulations in parallel
        for k in range(len(folders)):
            if opal_options.MPI_runs > 1:
                string = "mpirun -n "+str(opal_options.MPI_runs)+ " " + opal_options.executable + " " + opal_options.input_file
            else:
                string = opal_options.executable + " " + opal_options.input_file
            commands.append('cd ' + folders[k] + ' && ' + string)


        processes = [Popen(cmd, shell=True) for cmd in commands]
        # wait for completion
        for p in processes: p.wait()
        # Cleanup the list
        del commands[:]
    #Just in case return to the root folder
    os.chdir(cwd)

    return counter

#Mutate method for integer
def mutate(instance, mutpb):

    if len(opal_options.mutation_method)>1:
        #Execute python code from the user
        exec(opal_options.mutation_method)
    else:
        #Wether to do the mutation or not is decided outside
        if random.random() <= 1.0:
            index = random.randint(0, len(instance) - 1)
            instance[index] += space_search_random(MIN, MAX)


    return instance,

def space_search_random(MINval, MAXval):

    val = (random.randint(MINval,MAXval)/spatial_precision ) * spatial_precision
    return val

def checkBounds(MIN, MAX):
    #Ensure that the results are bounded
    # If many locations are studied and a clear distance is specified, then ensure that
    # this is satisfied by the locations
    Nvar = len(opal_options.ga_variables)
    Nwells = opal_options.ga_locations_to_study
    def decorator(func):
        def wrapper(*args, **kargs):

            offspring = func(*args, **kargs)
            for child in offspring:
                #Ensure that wells are not too close
                # Perform this by "pushing" wells that are too close
                if opal_options.mind_the_gap > 0:
                    # First two variables of a location to study are the X and Y locations
                    for i in range(Nwells - 1):
                        Xorig = np.asanyarray(child[i * Nwells:i * Nwells + 2])
                        for j in range(Nwells):
                            if j == i: continue  # Ignore itself!
                            Xother = np.asanyarray(child[j * Nwells:j * Nwells + 2])
                            dist = (np.dot(Xorig - Xother, Xorig - Xother)) ** 0.5
                            if dist < opal_options.mind_the_gap:
                                if abs(dist) < spatial_precision:  # Move the node away
                                    Xother[0] = Xother[0] + (-1) ** random.randrange(2) * opal_options.mind_the_gap
                                    Xother[1] = Xother[1] + (-1) ** random.randrange(2) * opal_options.mind_the_gap
                                else:
                                    # Move node to new position
                                    Xother = (Xother - Xorig) * opal_options.mind_the_gap / dist + Xother
                                # Ensure that the node is in the space of search
                                child[j * Nwells] = (np.int(Xother[0]) / spatial_precision) * spatial_precision
                                child[j * Nwells + 1] = (np.int(Xother[1]) / spatial_precision) * spatial_precision
                for i in range(len(child)):
                    if child[i] > MAX:
                        child[i] = MAX
                    elif child[i] < MIN:
                        child[i] = MIN
            return offspring
        return wrapper
    return decorator



def convert_instance(instance):
    j =0
    instance_back = copy.deepcopy(instance)
    for i in range(len(instance)):
        #if i == 0: continue #Ignore first because it is the reference
        instance_back[i] = linear_converter(instance_back[i], MIN, MAX,
                         opal_options.ga_variables[j].min_limit, opal_options.ga_variables[j].max_limit)
        j += 1
        #Restart j so it iterates over values per location
        if j%(len(instance_back)/ opal_options.ga_locations_to_study)==0: j = 0
    return instance_back
def convert_instance_back(instance):
    j = 0
    instance_back = copy.deepcopy(instance)
    for i in range(len(instance)):
        #if i == 0: continue #Ignore first because it is the reference
        instance_back[i] = linear_converter(instance_back[i], opal_options.ga_variables[j].min_limit,
                        opal_options.ga_variables[j].max_limit, MIN, MAX)
        j += 1
        #Restart j so it iterates over values per location
        if j%(len(instance_back)/ opal_options.ga_locations_to_study)==0: j = 0
    return instance_back
def linear_converter(old_var, old_min, old_max, new_min, new_max):
    if new_max == old_max and new_min == old_min:
        val = old_var
    else:
        val = int((float(old_var - old_min) / float(old_max - old_min) * float(new_max - new_min) + new_min))
    return val

# Drop the first element,
# which will be replace by our initial guess.
def set_initial_guess(population, init_guess_string):
    #Convert string from diamond into list with integers
    initial_guess = [int(e) if e.isdigit() else e for e in init_guess_string.split(',')]
    if len(initial_guess) != len(opal_options.ga_variables)*opal_options.ga_locations_to_study:
        print "ERROR: The initial guess must be of the size of locations to study times number of variables"
        print "Initial guess ignored."
        return
    population.pop()
    #Convert to internal space
    initial_guess = convert_instance_back(initial_guess)
    guess_ind = creator.Individual(initial_guess)
    #Make sure it is in the correct space

    population.insert(0, guess_ind)

####################CHANGE CONVERGENCE CRITERIA######################
#Finds the fittest
def get_best_result(population):

    if isinstance(population[0], list):
        fitness_values = list(map(fitness, population))
        if opal_options.ga_Minimise:
            index = fitness_values.index(min(fitness_values))
        else:
            index = fitness_values.index(max(fitness_values))
        return population[index]
    else:
        if opal_options.ga_Minimise:
            return min(population, key=operator.attrgetter('fitness'))
        else:
            return max(population, key=operator.attrgetter('fitness'))



def terminate(population, best_distance_so_far, generation):
    global previous_convergence
    reslt =fitness(get_best_result(population))
    if opal_options.ga_Minimise:
        result = min(reslt[0],best_distance_so_far)
    else:
        result = max(reslt[0],best_distance_so_far)
    if opal_options.ga_gradient_convergence > 0.:
        if abs(2.*result - sum(previous_convergence))/abs(result) < opal_options.ga_gradient_convergence:
            if generation>1: raise StopIteration
        else:
            previous_convergence[1] = previous_convergence[0]
            previous_convergence[0] = result

    if opal_options.ga_absolute_convergence > 0.:
        if opal_options.ga_Minimise and result <= opal_options.ga_absolute_convergence:
            raise StopIteration
        elif not opal_options.ga_Minimise and result >= opal_options.ga_absolute_convergence:
            raise StopIteration
    return False

#I presume it evaluates some sort of residual
def distance_from_best_result(population):
    result = get_best_result(population)
    return fitness(result)[0]

#####################################################################


#Outputs the best
def output(best_instance):
    best_instance_bak = convert_instance(best_instance)
    print 'Best result:', best_instance_bak
    distance = fitness(best_instance)
    if distance > 0.1:
        print "The best result produces: ", abs(distance[0])

#This functions creates an output file in csv format that includes the information of the best results
#Untested
def create_output_file(halloffame):
    import csv

    HallOfFame_file = "Best_results_"+opal_options.output_filename+".csv"

    with open(HallOfFame_file, 'w') as csvfile:
        datarwriter = csv.writer(csvfile, delimiter=',', quotechar='|')
        #First write header
        writeList = ["Ranking"]
        writeList.append("Fitness evaluation")
        for k in range(opal_options.ga_locations_to_study):
            for variable in opal_options.ga_variables:
                writeList.append("Location " + str(k+1)+ ": " + variable.name)

        writeList.append("Folder name")
        datarwriter.writerow(writeList)
        #Now proceed to write the data
        i = 0
        for instance in halloffame:
            i += 1
            instance_bak = convert_instance(instance)
            writeList = [i] + [fitness(instance)[0]] + instance_bak + [get_folder_name(instance)]
            datarwriter.writerow(writeList)

    return

def create_best_folders(halloffame):


    #Create a folder with the best result
    try:
        foldername = get_folder_name(halloffame[0])
    except:
        foldername = ""
    #New name for the folder
    newfoldername = "Best_result"
    os.system("cp -r "+ foldername + " "+ newfoldername)

    vtufile = libspud.get_option('/simulation_name')+"_1.vtu"
    #Now we copy the initial vtu of all the results in the hall of fame, to be able to easily explore
    # the different well locations as well as the production .csv files
    newfoldername = "Best_well_configurations"
    os.system("mkdir " + newfoldername)
    for k in range(len(halloffame)):
        foldername = get_folder_name(halloffame[k])
        #Copy vtu file in order
        os.system("cp "+ foldername + "/" +vtufile + " " + newfoldername+ "/" +"best_result_ranked_" + str(k+1)+".vtu")
        #Copy production .csv file
        os.system("cp "+ foldername + "/" +prod_file + " " + newfoldername+ "/" +"best_result_ranked_" + str(k+1)+".csv")

    return

#Initialise the genetic algorithm
def setup():

    toolbox = base.Toolbox()
    toolbox.register("attribute", space_search_random, MIN, MAX)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
                     toolbox.attribute, n=len(opal_options.ga_variables)*opal_options.ga_locations_to_study)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxOnePoint)
    toolbox.register("mutate", mutate, mutpb=opal_options.ga_mutation_prob)
    toolbox.decorate("mate", checkBounds(MIN, MAX))
    toolbox.decorate("mutate", checkBounds(MIN, MAX))


    if opal_options.ga_selection_method == 1:
        toolbox.register("select", tools.selBest)
    elif opal_options.ga_selection_method == 2:
        toolbox.register("select", tools.selNSGA2)
    else:
        opal_options.toolbox.register("select", tools.selSPEA2)

    toolbox.register("evaluate", fitness)

    #To parallelise using SCOOP
    #toolbox.register("map", futures.map)

    if opal_options.ga_CMA:
        #Covariance Matrix Adaptation Evolution Strategy (CMA-ES) [Hansen2001]
        strategy = cma.Strategy(centroid=opal_options.ga_centroid, sigma=opal_options.ga_sigma)
        toolbox.register("generate", strategy.generate, creator.Individual)
        toolbox.register("update", strategy.update)


    return toolbox


def main():


    #The population of a generation has to be a multiple of the number of cpus used,  not doing this is un-optimal
    parProcess = opal_options.number_processors
    opal_options.ga_population_generation = int(ceil(float(opal_options.ga_population_generation)/float(parProcess))  * parProcess)


    if not os.path.isfile(opal_options.executable):
        print "#############################################################################"
        print "ERROR: IC-FERST executble not found in the given path."
        print "#############################################################################"
        exit()




    toolbox = setup()
    #Create an initial population
    population = toolbox.population(n=opal_options.ga_population_generation)
    #Specify that an initial seed value
    if opal_options.ga_initial_guess:
        set_initial_guess(population, opal_options.ga_initial_guess)
    stats = tools.Statistics()
    stats.register("best_instance_of_population", get_best_result)
    stats.register("distance", distance_from_best_result)
    stats.register("terminate", terminate)
    halloffame = tools.HallOfFame(min(opal_options.ga_hall_of_fame,
                                    opal_options.ga_population_generation* opal_options.ga_max_generations-1))
    try:

        eaAlgorithm_by_steps(population, toolbox, opal_options.ga_breeding_prob, opal_options.ga_mutation_prob,
                                       opal_options.ga_max_generations, halloffame)

    except StopIteration:
       pass
    finally:
        #Create list with best results
        create_output_file(halloffame)
        create_best_folders(halloffame)
        best_instance = halloffame[0]
        output(best_instance)
        return best_instance


# if __name__ == '__main__':
#     main()
