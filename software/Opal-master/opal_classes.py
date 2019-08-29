#!/usr/bin/env python

#    Copyright (C) 2017 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Prof. C Pain
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    amcgsoftware@imperial.ac.uk
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation,
#    version 2.1 of the License.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA

#Authors: C. Heaney, C.C. Pain, P. Salinas, and others.....

import sys
import libspud
import datetime

class Data_Assimilation():
    def __init__(self):
        self.fwd_input_file=""

class Opal_input():
    def __init__(self):
        self.output_filename = ""                       # fwd
        self.executable = ""                            # fwd
        self.input_file = ""                            # fwd
        self.avoid_perturbing_near_wells  = False       # opal
        self.opal_operation = ""                        # opal
        self.random_seed = -1                           # sensitivity
        self.time_window=-1.0                           # sensitivity
        self.pass_down_ginit = False                    # sensitivity
        self.functional = Functional()                  # func
        self.da_gamma = 0.0                             # da-soa
        self.da_gamma_step = 0.0                        # da-soa
        self.da_lambda = 0.0                            # da-soa
        self.da_lambda_step = 0.0                       # da-soa
        self.da_nu = 0.0                                # da-soa
        self.da_nits_CG = 1                             # da-soa
        self.da_nits_nonlinear = 1                      # da-soa
        self.da_nits_inner = 1                          # da-soa
        self.da_error_tolerance = 0.1                   # da-soa
        self.da_nits_precon = 1                         # da-soa
        self.da_w_relax = 0.0                           # da-soa
        self.da_inverse_toler = 0.0                     # da-soa
        self.data_assimilation = Data_Assimilation()
        self.Field_To_study_list=[]
        ###   GA optimisation    ###
        self.ga_max_generations = 0
        self.ga_locations_to_study = 0
        self.ga_initial_guess = ""
        self.ga_gradient_convergence = -1
        self.ga_absolute_convergence = -1
        self.ga_population_generation = 0
        self.ga_breeding_prob = 0.0
        self.ga_mutation_prob = 0.0
        self.ga_fitness_functional = ""
        self.ga_Minimise = True
        self.ga_evol_algorithm = 1
        self.ga_selection_method = 1
        self.ga_mu = 0
        self.ga_lambda = 0
        self.ga_CMA = False
        self.ga_centroid = 0
        self.ga_sigma = 0
        self.ga_variables = []
        self.ga_hall_of_fame = 10
        self.number_processors = 1
        self.MPI_runs = 1
        self.precision = 0.01
        self.mutation_method = ""
        self.generation_method = ""
        self.mind_the_gap = -1
        self.trelis_path = ""
        self.trelis_input_file = ""
        self.producer_ids = ""
        self.injector_ids = ""
        self.optimise_input = False
class Functional_field():
    def __init__(self):
        self.name = ""                                  # func field
        self.phase = -1                                 # func field
        self.field_type = ""                            # func field
        self.mesh_type = ""                             # func field
        self.indices_of_interest = -1                   # func field
        self.DGified = ""                               # func field
        self.reference_value = -1                       # func field


class Functional():
    def __init__(self):
        self.location_of_interest=0.0, 0.0, 0.0         # func
        self.time = ""                                  # func
        self.type = ""                                  # func
        self.Functional_field = []                      # func
        self.square_integrand = False

class Field_To_study():
    def __init__(self):
        self.field_type=""                              # sensitivity
        self.perturbations_to_do=0                      # sensitivity
        self.parallel_processors=1                      # sensitivity
        self.initial_condition=True                     # sensitivity
        self.boundary_condition=False                   # sensitivity
        self.source=False                               # sensitivity
        self.gram_schmidt=False                         # sensitivity
        self.nSmooth=0                                  # sensitivity
        self.use_G = False                              # sensitivity
        self.threshold_use_G = 3                        # sensitivity
        self.sigma=0.                                   # sensitivity
        self.phases_to_perturb=0                        # sensitivity
        self.diagonal_tweak=False                       # sensitivity
        self.diagonal_tweak_method=''                   # sensitivity
        self.diagonal_tweak_value=0.                    # sensitivity
        self.name=''                                    # sensitivity
        self.pycode=''                                  # sensitivity

class Nirom_input():
    def __init__(self):
        self.snapshots = Snapshots()
        self.compression = Compression()
        self.training = Train_NN()
        self.prediction = Prediction()

class Forward_model_options():
    def __init__(self):
        self.output_filename = ""                       # fwd
        self.executable = ""                            # fwd
        self.input_file = ""                            # fwd
        self.path_to_results = ""                       # fwd
        self.results_filebase = ""                      # fwd
        self.nTime = 0                                  # fwd
        self.results_extension=''                       # fwd
        self.ndim = 0
	#self.compression = Compression()

class Snapshots():
    def __init__(self):
        self.create=False
        self.location=False
        self.values=[]
        self.nParameters=1
        self.nComponents=1
        self.nNodes = 0

class Compression():
    def __init__(self):
        self.svd=False
        self.eigh=False
        self.DD_eigh = False
        self.autoencoder=False
        self.num_sub_base_2 = 1
        self.epochs=[]
        self.batch_size=[]
        self.neurons_each_layer=[]
        self.nlatent=[]
        self.cumulative_tol=[]
        self.field=[]
        self.nNodes=[]
        self.nPOD=[]
        self.basis_functions=[]
        self.s_values=[]
        self.svd_autoencoder=False
        self.whichd = []
        self.local_list = []
        self.order = []
        self.num_nodes_in_sub = []
        self.olde2new = []
        self.dd_snapshot = []
        self.dd_basis_functions = []
        #self.dd_s_values = []

    def write_sing_values(self):
        if self.DD_eigh is False:
            f= open('singular_values.dat',"w+")
            f.write('# index, s_values, normalised s_values, cumulative energy \n' )
            for k in range(len(self.s_values)):
                f.write('# field: %s\n' % self.field[k])
                total = 0.0
                s_values = self.s_values[k]
                for i in range(len(s_values)):
                    total = total + s_values[i]*s_values[i]

                running_total = 0.0
                for i in range(len(s_values)):
                    running_total = running_total + s_values[i]*s_values[i]
                    f.write ('%d %g %g %18.10g \n' % (i, s_values[i], s_values[i]/s_values[0], running_total/total) )
            f.close()
        elif self.DD_eigh:
            for j in range(len(self.order)):
                f= open('singular_values_'+str(j)+'.dat',"w+")
                f.write('# index, s_values, normalised s_values, cumulative energy \n' )
                for k in range(self.s_values.shape[1]):
                    f.write('# field: %s\n' % self.field[k])
                    total = 0.0
                    s_values = self.s_values[j][k]
                    for i in range(s_values.shape[0]):
                        total = total + s_values[i]*s_values[i]

                    running_total = 0.0
                    for i in range(s_values.shape[0]):
                        running_total = running_total + s_values[i]*s_values[i]
                        f.write ('%d %g %g %18.10g \n' % (i, s_values[i], s_values[i]/s_values[0], running_total/total) )
                f.close()


class Train_NN():
    def __init__(self):
        self.GPR = False
        self.DD_GPR = False
        self.LSTM = False
        self.DD_LSTM = False
        self.GPR_scaling = []
        self.GPR_RBF_length_scale = 1.0
        self.GPR_RBF_length_bounds = [1e-5, 1e5]
        self.GPR_constant_value = 1.0
        self.GPR_constant_bounds = [1e-5, 1e5]

class Prediction():
    def __init__(self):
        self.nTime = 0

class Timings():
    def __init__(self):
        self.generate_snapshots = -1.0
        self.read_in_snapshots = -1.0
        self.compression = -1.0
        self.training = -1.0
        self.replication = -1.0
        self.prediction = -1.0

    def print_values(self,svd_eigh):
        now = datetime.datetime.now()
        print ('nirom timings (in seconds), recorded  on', now.strftime("%d-%m-%Y %H:%M:%S"))
        if self.generate_snapshots > 0:
            print ('snapshot generation took: ',self.generate_snapshots)
        if self.read_in_snapshots > 0:
            print ('reading in snapshots took:',self.read_in_snapshots)
        if self.compression > 0:
            print ('compression (' + svd_eigh + ') took: ',self.compression)
        if self.training > 0:
            print ('training the NN took:     ',self.training)
        if self.prediction > 0:
            print ('prediction with NN took:  ',self.prediction)

    def write_values(self,filename, svd_eigh):
        now = datetime.datetime.now()
        f= open(filename,"a+")
        f.write('nirom timings (in seconds), recorded  on %s\n' % now.strftime("%d-%m-%Y %H:%M:%S"))
        if self.generate_snapshots > 0:
            f.write ('snapshot generation took:  %g\n' % self.generate_snapshots)
        if self.read_in_snapshots > 0:
            f.write ('reading in snapshots took: %g\n' % self.read_in_snapshots)
        if self.compression > 0:
            string = 'compression (' + svd_eigh + ') took:   %g\n'
            f.write (string % self.compression)
        if self.training > 0:
            f.write ('training the NN took:      %g\n' % self.training)
        if self.prediction > 0:
            f.write ('predicting with the NN took: %g\n' % self.prediction)
        f.close()


class ga_variable():
    def __init__(self):
        self.name=""
        self.min_limit = 0
        self.max_limit = 0
        self.variable_pattern = ""
        self.normaliser = 1.0


def read_in_importance_map_options(path,opal_options):

    opal_options.output_filename = libspud.get_option(path + 'Output_filename')
    opal_options.executable = libspud.get_option(path + 'Executable')
    opal_options.input_file = libspud.get_option(path + 'Input_file')

    if (libspud.have_option(path + 'Time_window')):
        opal_options.time_window = libspud.get_option(path + 'Time_window')
        if (libspud.have_option(path + 'Time_window/pass_down_ginit')):
            opal_options.pass_down_ginit = True
    else:
        opal_options.time_window = 0
        opal_options.pass_down_ginit = False

    if libspud.have_option(path + 'random_seed'):
        opal_options.random_seed = 200
    else:
        opal_options.random_seed = -1

    opal_options.functional.location_of_interest = libspud.get_option(path + 'functional/location_of_interest')
    if libspud.have_option(path + 'functional/all_time'):
        print "has option all time"
        opal_options.functional.time = "all_time"
    else:
        opal_options.functional.time = "end_time"
    if libspud.have_option(path + 'functional/type/standard'):
        opal_options.functional.type = "standard"
    elif libspud.have_option(path + 'functional/type/geothermal_over_time'):
        opal_options.functional.type = "geothermal_over_time"
    nffields = libspud.option_count(path + 'functional/Functional_field')
    for i in range(nffields):
        temp_field = Functional_field()
        fieldpath = path + 'functional/Functional_field['+str(i)+']'
        temp_field.name = libspud.get_option(fieldpath+'/name')
        if libspud.have_option(fieldpath+'/phase'):
            temp_field.phase = libspud.get_option(fieldpath+'/phase')
        if (libspud.have_option(fieldpath+'/field_type/scalar_field')):
            temp_field.field_type = "scalar_field"
        elif (libspud.have_option(fieldpath+'/field_type/vector_field')):
            temp_field.field_type = "vector_field"
        elif (libspud.have_option(fieldpath+'/field_type/tensor_field')):
            temp_field.field_type = "tensor_field"
        if (libspud.have_option(fieldpath+'/mesh_type/element')):
            temp_field.mesh_type = "element"
        elif (libspud.have_option(fieldpath+'/mesh_type/node')):
            temp_field.mesh_type = "node"

        opal_options.functional.Functional_field.append(temp_field)

    if libspud.have_option(path + 'functional/square_integrand'):
        opal_options.functional.square_integrand = True
    else:
        opal_options.functional.square_integrand = False

    nfields = libspud.option_count(path + 'Field_to_study')
    for i in range(nfields):
        temp_field = Field_To_study()
        fieldpath = path + 'Field_to_study['+str(i)+']'
        temp_field.name=libspud.get_option(fieldpath+'/name')
        if (libspud.have_option(fieldpath+'/field_type/scalar_field')):
            temp_field.field_type = "scalar_field"
        elif (libspud.have_option(fieldpath+'/field_type/vector_field')):
            temp_field.field_type = "vector_field"
        elif (libspud.have_option(fieldpath+'/field_type/tensor_field')):
            temp_field.field_type = "tensor_field"
        else:
            if (libspud.have_option(fieldpath+'/field_type/scalar_field')):
                temp_field.field_type = "porous_scalar"
            elif (libspud.have_option(fieldpath+'/field_type/vector_field')):
                temp_field.field_type = "porous_vector"
            elif (libspud.have_option(fieldpath+'/field_type/tensor_field')):
                temp_field.field_type = "porous_tensor"

        #temp_field.location_of_interest = libspud.get_option(fieldpath+'/location_of_interest')
        temp_field.perturbations_to_do = libspud.get_option(fieldpath+'/perturbations_to_do')
        temp_field.parallel_processors = 1 #Serial by default
        if (libspud.have_option(fieldpath+'/perturbations_to_do/parallel_processors')):
            temp_field.parallel_processors = libspud.get_option(fieldpath+'/perturbations_to_do/parallel_processors')
            if libspud.have_option(path + '/Number_processors/MPI_runs'):
                opal_options.MPI_runs = libspud.get_option(path + '/Number_processors/MPI_runs')
                opal_options.number_processors /= opal_options.MPI_runs
        temp_field.initial_condition = libspud.have_option(fieldpath+'/Initial_condition')
        temp_field.boundary_condition = libspud.have_option(fieldpath+'/Boundary_condition')
        temp_field.source = libspud.have_option(fieldpath+'/source')
        if (libspud.have_option(fieldpath+'/Sigma')):
            temp_field.sigma = libspud.get_option(fieldpath+'/Sigma')
        else:
            temp_field.sigma = 0.01
        if (libspud.have_option(fieldpath+'/python_function')):
            temp_field.pycode = libspud.get_option(fieldpath+'/python_function')
        temp_field.phases_to_perturb = libspud.get_option(fieldpath+'/Phases_to_perturb')
        if libspud.have_option(fieldpath+'/stabilise_through_diagonal'):
            temp_field.diagonal_tweak = True
            temp_field.diagonal_tweak_value = libspud.get_option(fieldpath+'/stabilise_through_diagonal/epsilon')
            if libspud.have_option(fieldpath+'/stabilise_through_diagonal/constrain_evalue'):
                temp_field.diagonal_tweak_method = 'constrain_evalue'
            elif libspud.have_option(fieldpath+'/stabilise_through_diagonal/add_to_diag'):
                temp_field.diagonal_tweak_method = 'add_to_diag'
            else:
                print "error, option not expected" , libspud.get_option(fieldpath+'/stabilise_through_diagonal/')
        else:
            # can probably remove these three lines now __init__ is working....
            temp_field.diagonal_tweak = False
            temp_field.diagonal_tweak_method = ''
            temp_field.diagonal_tweak_value = 0.0
        if libspud.have_option(fieldpath+'/Gram-Schmidt'):
            temp_field.gram_schmidt = True
        else:
            temp_field.gram_schmidt = False
        if libspud.have_option(fieldpath+'/smoothing'):
            temp_field.nSmooth = libspud.get_option(fieldpath+'/smoothing')
        else:
            temp_field.nSmooth = 0
        if libspud.have_option(fieldpath+'/use_G'):
            temp_field.use_G = True
            temp_field.threshold_use_G =  libspud.get_option(fieldpath+'/use_G')
        else:
            temp_field.use_G = False
            temp_field.threshold_use_G =  -1
        opal_options.Field_To_study_list.append(temp_field)

    return opal_options

def read_in_nirom_options(path, opal_options, nirom_options, fwd_options):

    if libspud.have_option('opal_operation/nirom/snapshots_create'):
        nirom_options.snapshots.create = True
        fwd_options.exectuable = libspud.get_option('opal_operation/nirom/snapshots_create/executable')
        fwd_options.path_to_results = libspud.get_option('opal_operation/nirom/snapshots_create/path')
        fwd_options.input_file = libspud.get_option('opal_operation/nirom/snapshots_create/input_file')
    elif libspud.have_option('opal_operation/nirom/snapshots_location'):
        nirom_options.snapshots.location = True
        fwd_options.path_to_results = libspud.get_option('opal_operation/nirom/snapshots_location/path')
        fwd_options.input_file = libspud.get_option('opal_operation/nirom/snapshots_location/input_file')
        if libspud.have_option('opal_operation/nirom/snapshots_location/output_filebase'):
           fwd_options.results_filebase = libspud.get_option('opal_operation/nirom/snapshots_location/output_filebase')
        if libspud.have_option('opal_operation/nirom/snapshots_location/extension'):
            fwd_options.results_extension = libspud.get_option('opal_operation/nirom/snapshots_location/extension')
        if libspud.have_option('opal_operation/nirom/snapshots_location/nParameters'):
            nirom_options.snapshots.nParameters = libspud.get_option('opal_operation/nirom/snapshots_location/nParameters')

    if libspud.have_option('/opal_operation/nirom/compression'):
        if libspud.have_option('/opal_operation/nirom/compression/svd_type'):
            if libspud.have_option('/opal_operation/nirom/compression/svd_type/svd_method'):
                nirom_options.compression.svd = True
            elif libspud.have_option('/opal_operation/nirom/compression/svd_type/eigh_method'):
                nirom_options.compression.eigh = True
            elif libspud.have_option('/opal_operation/nirom/compression/svd_type/DD_eigh_method'):
                nirom_options.compression.DD_eigh = True
                nirom_options.compression.num_sub_base_2 = libspud.get_option('/opal_operation/nirom/compression/svd_type/DD_eigh_method/num_sub_base_2')
            nFields = libspud.option_count('/opal_operation/nirom/compression/svd_type/field_name')
            for ii in range(nFields):
                my_string = libspud.get_option('/opal_operation/nirom/compression/svd_type/field_name[' + str(ii) + ']/name')
                nirom_options.compression.field.append(my_string)
                value_ct = 0.0
                value_nPOD = -1
                if libspud.have_option('/opal_operation/nirom/compression/svd_type/field_name['+str(ii)+']/cumulative_tol'):
                    value_ct = libspud.get_option('/opal_operation/nirom/compression/svd_type/field_name['+str(ii)+']/cumulative_tol')
                elif libspud.have_option('/opal_operation/nirom/compression/svd_type/field_name['+str(ii)+']/nPOD'):
                    value_nPOD = libspud.get_option('/opal_operation/nirom/compression/svd_type/field_name['+str(ii)+']/nPOD')
                nirom_options.compression.cumulative_tol.append(value_ct)
                nirom_options.compression.nPOD.append(value_nPOD)

            # svd_autoencoder option should go here
            if libspud.have_option('/opal_operation/nirom/compression/svd_type/svd_autoencoder'):
                nirom_options.compression.svd_autoencoder = True
                nirom_options.compression.epochs = libspud.get_option('/opal_operation/nirom/compression/svd_type/svd_autoencoder/number_epochs')
                nirom_options.compression.batch_size = libspud.get_option('/opal_operation/nirom/compression/svd_type/svd_autoencoder/batch_size')
                nirom_options.compression.neurons_each_layer = libspud.get_option('/opal_operation/nirom/compression/svd_type/svd_autoencoder/neurons_each_layer')

        elif libspud.have_option('/opal_operation/nirom/compression/autoencoder_type'):
            nirom_options.compression.autoencoder = True

            nFields = libspud.option_count('/opal_operation/nirom/compression/autoencoder_type/field_name')
            for ii in range(nFields):
                my_string = libspud.get_option('/opal_operation/nirom/compression/autoencoder_type/field_name['+str(ii)+']/name')
                nirom_options.compression.field.append(my_string)
                batch_size = libspud.get_option('/opal_operation/nirom/compression/autoencoder_type/field_name['+str(ii)+']/batch_size')
                nirom_options.compression.batch_size.append(batch_size)
                number_epochs = libspud.get_option('/opal_operation/nirom/compression/autoencoder_type/field_name['+str(ii)+']/number_epochs')
                nirom_options.compression.epochs.append(number_epochs)
                neurons_each_layer = libspud.get_option('/opal_operation/nirom/compression/autoencoder_type/field_name['+str(ii)+']/neurons_each_layer')
                nirom_options.compression.neurons_each_layer.append(neurons_each_layer) # neurons per layer

    if libspud.have_option('/opal_operation/nirom/training'):
        if libspud.have_option('/opal_operation/nirom/training/GPR'):
            nirom_options.training.GPR = True
            nirom_options.training.GPR_scaling = libspud.get_option('/opal_operation/nirom/training/GPR/scaling_bounds')
            nirom_options.training.GPR_RBF_length_scale = libspud.get_option('/opal_operation/nirom/training/GPR/RBF_length_scale')
            nirom_options.training.GPR_RBF_length_bounds = libspud.get_option('/opal_operation/nirom/training/GPR/RBF_length_scale_bounds')
            nirom_options.training.GPR_constant_value = libspud.get_option('/opal_operation/nirom/training/GPR/constant_value')
            nirom_options.training.GPR_constant_bounds = libspud.get_option('/opal_operation/nirom/training/GPR/constant_bounds')
        if libspud.have_option('/opal_operation/nirom/training/LSTM'):
            nirom_options.training.LSTM = True
        if libspud.have_option('/opal_operation/nirom/training/DD_LSTM'):
            nirom_options.training.DD_LSTM = True
            #print "LSTM not currently implemented ... aborting"
            #sys.exit(0)
        if libspud.have_option('/opal_operation/nirom/training/DD_GPR'):
            nirom_options.training.DD_GPR = True

    if libspud.have_option('opal_operation/nirom/prediction'):
        nirom_options.prediction.nTime = libspud.get_option('opal_operation/nirom/prediction/nTime')

    return  opal_options, nirom_options, fwd_options

def read_in_soda_options(path,opal_options):

    opal_options.da_gamma = libspud.get_option(path + 'gamma')
    opal_options.da_gamma_step = libspud.get_option(path + 'gamma_step')
    opal_options.da_lambda = libspud.get_option(path+'lambda')
    opal_options.da_lambda_step = libspud.get_option(path+'lambda_step')
    opal_options.da_nu = libspud.get_option(path+'nu')
    opal_options.da_nits_CG = libspud.get_option(path+'nits_CG')
    opal_options.da_nits_nonlinear = libspud.get_option(path+'nits_nonlinear')
    opal_options.da_nits_inner = libspud.get_option(path+'nits_inner')
    opal_options.da_error_tolerance = libspud.get_option(path+'error_tolerance')
    opal_options.da_nits_precon = libspud.get_option(path+'nits_precon')
    opal_options.da_w_relax = libspud.get_option(path+'w_relax')
    opal_options.da_inverse_toler = libspud.get_option(path+'inverse_toler')

    return opal_options

def get_opal_options():
    "Initialises the structure for the options."
    opal_options = Opal_input()
    nirom_options = Nirom_input()
    fwd_options = Forward_model_options()

    #Read the input data from the user
    opal_input_file = str(sys.argv[1])
    #Load the .opal file
    libspud.load_options(opal_input_file)
    ##Create Opal Variable
#    opal_options = Opal_input()
    ##Start reading options
    print "reading opal input file", opal_input_file
    if  libspud.have_option('/avoid_perturbing_near_wells'):
        opal_options.avoid_perturbing_near_wells = True

    if (libspud.have_option('/opal_operation/importance_map')):

        opal_options.opal_operation = 'importance_map'
        im_path = '/opal_operation/importance_map/'
        opal_options = read_in_importance_map_options(im_path, opal_options)

    elif (libspud.have_option('/opal_operation/data_assimilation')):

        opal_options.opal_operation = 'da'
        opal_options.data_assimilation.fwd_input_file = libspud.get_option('/opal_operation/data_assimilation/fwd_input_file')
        #print "opal_options.data_assimilation.fwd_input_file", opal_options.data_assimilation.fwd_input_file

    elif (libspud.have_option('/opal_operation/nirom')):

        opal_options.opal_operation = 'nirom'
        nirom_path = 'opal_operation/nirom/'
        opal_options, nirom_options, fwd_options = read_in_nirom_options(nirom_path, opal_options, nirom_options, fwd_options)

    elif (libspud.have_option('/opal_operation/second_order_da')):

        opal_options.opal_operation = 'second_order_da'
        da_path = '/opal_operation/second_order_da/'
        opal_options = read_in_soda_options(da_path,opal_options)

    elif (libspud.have_option('/opal_operation/ga_optimisation')):
        opal_options.opal_operation = 'ga_optimisation'
        path = '/opal_operation/ga_optimisation'

        # previously these three options were at the top level
        opal_options.output_filename = libspud.get_option(path + '/Output_filename')
        opal_options.executable = libspud.get_option(path + '/Executable')
        opal_options.input_file = libspud.get_option(path + '/Input_file')

        opal_options.ga_max_generations = libspud.get_option(path+'/convergence_settings/Maximum_generations')

        opal_options.ga_locations_to_study = libspud.get_option(path+'/Locations_to_study/value')
        if libspud.have_option(path+'/Locations_to_study/Initial_guess'):
            opal_options.ga_initial_guess = libspud.get_option(path+'/Locations_to_study/Initial_guess')
        if libspud.have_option(path + '/Locations_to_study/Mind_the_gap'):
            opal_options.mind_the_gap = libspud.get_option(path + '/Locations_to_study/Mind_the_gap')
        if libspud.have_option(path + '/Trelis_integration'):
            opal_options.trelis_path = libspud.get_option(path + '/Trelis_integration/trelis_path')
            opal_options.trelis_input_file = libspud.get_option(path + '/Trelis_integration/trelis_input_file')
        opal_options.optimise_input = libspud.have_option(path + '/optimise_input')
        if libspud.have_option(path+'/convergence_settings/Gradient_convergence'):
            opal_options.ga_gradient_convergence = libspud.get_option(path+'/convergence_settings/Gradient_convergence')
        if libspud.have_option(path+'/convergence_settings/Absolute_convergence'):
            opal_options.ga_absolute_convergence = libspud.get_option(path+'/convergence_settings/Absolute_convergence')
        if libspud.have_option(path+'/Number_processors'):
            opal_options.number_processors = libspud.get_option(path+'/Number_processors')
            if libspud.have_option(path + '/Number_processors/MPI_runs'):
                opal_options.MPI_runs = libspud.get_option(path + '/Number_processors/MPI_runs')
                opal_options.number_processors /= opal_options.MPI_runs
        opal_options.ga_population_generation = libspud.get_option(path+'/Population_generation')
        opal_options.ga_breeding_prob = libspud.get_option(path+'/Breeding_probability/value')
        if libspud.have_option(path+'/Breeding_probability/Generation_method'):
            opal_options.generation_method = libspud.get_option(path+'/Breeding_probability/Generation_method')
        opal_options.ga_mutation_prob = libspud.get_option(path+'/Mutation_probability/value')
        if libspud.have_option(path+'/Mutation_probability/Mutation_method'):
            opal_options.mutation_method = libspud.get_option(path+'/Mutation_probability/Mutation_method')
        opal_options.precision = libspud.get_option(path + '/Precision')
        if libspud.have_option(path+'/Fitness/python_function'):
            opal_options.ga_fitness_functional = libspud.get_option(path+'/Fitness/python_function')
        opal_options.ga_Minimise = libspud.have_option(path+'/Fitness/Optimisation/Minimise')
        opal_options.producer_ids = libspud.get_option(path+'/Fitness/producer_ids')
        opal_options.injector_ids = libspud.get_option(path+'/Fitness/injector_ids')
        if libspud.have_option(path+'/GA_methods/Evolutionary_algorithm/eaSimple'):
            opal_options.ga_evol_algorithm = 1
        elif libspud.have_option(path+'/GA_methods/Evolutionary_algorithm/eaMuCommaLambda'):
            opal_options.ga_evol_algorithm = 2
            opal_options.ga_mu = libspud.get_option(path+'/GA_methods/Evolutionary_algorithm/eaMuCommaLambda/Mu')
            opal_options.ga_lambda = libspud.get_option(path+'/GA_methods/Evolutionary_algorithm/eaMuCommaLambda/Lambda')
        elif libspud.have_option(path+'/GA_methods/Evolutionary_algorithm/eaMuPlusLambda'):
            opal_options.ga_evol_algorithm = 3
            opal_options.ga_mu = libspud.get_option(path+'/GA_methods/Evolutionary_algorithm/eaMuPlusLambda/Mu')
            opal_options.ga_lambda = libspud.get_option(path+'/GA_methods/Evolutionary_algorithm/eaMuPlusLambda/Lambda')
        if libspud.have_option(path+'/GA_methods/Selection_method/selBest'):
            opal_options.ga_selection_method = 1
        elif libspud.have_option(path+'/GA_methods/Selection_method/selNSGA2'):
            opal_options.ga_selection_method = 2
        elif libspud.have_option(path+'/GA_methods/Selection_method/selSPEA2'):
            opal_options.ga_selection_method = 3

        opal_options.ga_CMA = False
        if libspud.have_option(path+'/GA_methods/Use_CMA'):
            opal_options.ga_CMA = True
            opal_options.ga_centroid = libspud.get_option(path+'/GA_methods/Use_CMA/centroid')
            opal_options.ga_sigma = libspud.get_option(path+'/GA_methods/Use_CMA/sigma')
        if libspud.have_option(path+'Hall_of_fame'):
            opal_options.ga_hall_of_fame = libspud.get_option(path+'Hall_of_fame')
        nfields = libspud.option_count(path+'/Variables')
        for i in range(nfields):
            temp_field = ga_variable()
            fieldpath = path+'/Variables[' + str(i) + ']'
            temp_field.name = libspud.get_option(fieldpath+'/name')
            temp_field.min_limit = libspud.get_option(fieldpath+'/Min_limit')
            temp_field.max_limit = libspud.get_option(fieldpath+'/Max_limit')
            temp_field.variable_pattern = libspud.get_option(fieldpath+'/Variable_pattern')
            temp_field.normaliser = 1.0
            if libspud.have_option(fieldpath+'/normaliser'):
                temp_field.normaliser = libspud.get_option(fieldpath+'/normaliser')
            ##Now append to the variables list
            opal_options.ga_variables.append(temp_field)

    libspud.clear_options()

    return opal_options, fwd_options, nirom_options
