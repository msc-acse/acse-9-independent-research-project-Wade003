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

#Authors: C. Heaney, C.C. Pain and P. Salinas

import copy
import vtk
import vtktools
import sys
import os
import subprocess
from subprocess import Popen
import matplotlib.pyplot as plt
import numpy as np
from math import *
from scipy import interpolate
from scipy.interpolate import interp1d
from vtk.util.numpy_support import vtk_to_numpy
from importance_map_calc  import *
from importance_map_tools import *
import time
import random
import libspud


def generate_importance_map(opal_options):
    random_seed = opal_options.random_seed

    # Number of fields that are going to be studied
    nfields = len(opal_options.Field_To_study_list)

    if nfields>1:
        print "this might not work as to modify vol fraction libspud is used, and to modify Velocity"
        print " readlines, .insert is used"

    cwd = os.getcwd()

    # get xml extension (mpml or flml)
    xml_extension = get_xml_extension(opal_options.input_file)

    # load the options for the forward model from xml file
    libspud.load_options(opal_options.input_file)

    # count the number of material_phases
    nPhases = 0
    if libspud.have_option('/material_phase[0]'): 
        nPhases = libspud.option_count('/material_phase')

#    vol_frac = np.zeros(nPhases)
#    for iPhase in range(nPhases): 
#        vol_frac[iPhase] = libspud.get_option('/material_phase['+str(iPhase)+']/scalar_field::PhaseVolumeFraction/prognostic/initial_condition/constant') 

    phases_to_perturb = opal_options.Field_To_study_list[0].phases_to_perturb
    nPhases_to_perturb = len(phases_to_perturb)
    nPhases_orig = nPhases
    have_wells = False
    if libspud.have_option('/porous_media/wells_and_pipes'):
        nPhases = nPhases/2
        nWells = nPhases   
        have_wells = True

    # modify vol fracs so that perturbations will keep within the saturation curve
    # do this for all fields so don't have to worry about taking abs() of random perturbations???
    # SHOULD DO THIS FOR EVERY FIELD?
    modified_vol_frac = np.zeros(nPhases_to_perturb)

    if nPhases>1:
        modified_vol_frac = modify_volume_fractions(opal_options.Field_To_study_list[0].phases_to_perturb,opal_options.Field_To_study_list)
    
    # get the dimension of the problem geometry 
    nDim = libspud.get_option('/geometry[0]/dimension/')

    # create a copy of the mpml file 
    # (in meld, an mpml file generated in this way looks different to the original mpml file  )
    # new_xml = opal_options.input_file + '-copy'
    # libspud.write_options(new_xml)

    fwd_output_name =  libspud.get_option('/simulation_name')
    meshname_base =  libspud.get_option('/geometry[0]/mesh[0]/from_file/file_name')

    # for flml should get rid of this if set by user
    if xml_extension == 'mpml':
        # modify the file to include <disable_dump_at_end/>
        f = open(opal_options.input_file, "r")
        contents = f.readlines()
        f.close()
 
        disable = False
        for iline,line in enumerate(contents):

            if "<disable_dump_at_end" in line:
                disable = True

            if disable:
                break

            if "</io>" in line and not disable:
                contents.insert(iline, "    <disable_dump_at_end/>\n")
                break
 
        if not disable:
            f = open(opal_options.input_file, "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()

    elif xml_extension == 'flml':
        # modify the file to NOT include <disable_dump_at_end/>
        f = open(opal_options.input_file, "r")
        contents = f.readlines()
        f.close()
 
        for iline,line in enumerate(contents):
           if "<disable_dump_at_end" in line:
                del contents[iline]
                break
        f = open(opal_options.input_file, "w")
        contents = "".join(contents)
        f.write(contents)
        f.close()



    # run the code with checkpointing at t0 to have the original (unperturbed) initial 
    # conditions in vtu form and to have the unperturbed solution to calculate the 
    # change in solution field due to the perturbation
    exit = False
    try:
        for field in opal_options.Field_To_study_list:

            if field.name == 'Tracer' or field.name == 'Temperature' or field.name == 'PhaseVolumeFraction' or field.name == 'PhaseVolumeFractionSource': 
                if field.initial_condition:      
                    print "generating unperturbed ic with checkpointing"
                    nTime, n_in_window, nwindow = create_unperturbed_initial_conditions_vtu_and_run(opal_options,cwd,field.name,modified_vol_frac,xml_extension,have_wells,phases_to_perturb)

                elif field.source:
                    print "create unperturbed source" # modified_vol_frac = np.zeros(nPhases)
                    nTime, n_in_window, nwindow = create_unperturbed_initial_conditions_vtu_and_run(opal_options,cwd,field.name,np.zeros(nPhases_to_perturb),xml_extension,have_wells,phases_to_perturb)
            elif field.name == 'Velocity':
                print "not checked for vector fields ... aborting"
                exit = True
            else:
                print "field not recognised ... aborting"
                exit = True
    except:
        pass

    if exit:
        sys.exit(0)

    # does the code DGify the field of interest in the output files? (the ic is taken from the checkpoint files on 
    # the original mesh for that field, the output might be projected onto a DG mesh for "convenience".
    output_name =  libspud.get_option('/simulation_name')
    DGified, nNodes = have_cg_fields_been_dgified(opal_options.Field_To_study_list, output_name, phases_to_perturb, nPhases, nPhases_orig, 0) # iwindow = 0...
    # might have to change this to initial file of a window   

    libspud.clear_options()

    # file to store max value of G_initial
    create_max_G_file = True
    filename = 'max_G_t0.dat'
    max_G_t0 = open(filename, "w")
    max_G_t0.write( "# number of ensembles, phase number, max value of G" )
    max_G_t0.write( " \n" )
      

    time_fwd = 0.0
    time_imp_map = 0.0

    # write condition number to file
    filename = 'condition_number.dat'
    cond = open(filename, "w")
    cond.write( "# time index, condition number before add_to_diag and after" )
    cond.write( " \n" )

    # this code will only work for one field (iedally each field would have its own Importance map)
    directories = []
    #For free we make nPert a multiple of the number of processors used
    parProcess = opal_options.Field_To_study_list[0].parallel_processors
    nPert = opal_options.Field_To_study_list[0].perturbations_to_do
    nPert = int(ceil(float(nPert)/float(parProcess))  * parProcess)
    DF = np.zeros((nPert*nPhases_to_perturb,nwindow))
    # the above assumes these things don't depend on each field
    # atw all time windows ?
    G_initial_atw = np.ones((nNodes,nPhases_to_perturb*nwindow))
    importance_map_array = np.zeros((nfields,nNodes,nTime)) # won't work for time windows of nPhases>1
    Ms_array = np.zeros((nfields,nNodes,nPert,nTime)) #only works for nPhases_to_perturb = 1
    Pm = np.zeros((nfields,nNodes,nPert)) #matrix of perturbations for all fields
    Gm= np.zeros((nfields,nNodes,nPert))    
    
    if nPhases_to_perturb > 1:
        print "only coded for nPhases_to_perturb = 1  in Ms"
        exit(0)

    path = 'unperturbed/'
    functional = set_up_functional(opal_options,nPhases_orig,nPhases,path,fwd_output_name)
    
    for iwindow in reversed(range(nwindow)): # go backwards through time windows starting from last

        beginTimeWindow = iwindow*(n_in_window-1)
        endTimeWindow = (iwindow+1)*(n_in_window-1)

        if iwindow == nwindow - 1:
            endTimeWindow = (iwindow+1)*(n_in_window-1)

        # initialise G_initial_atw with G_initial_atw from future time window
        if iwindow < nwindow-1 and opal_options.pass_down_ginit: 
            print "passing down G initial"
            G_initial_atw[:,nPhases_to_perturb*iwindow:nPhases_to_perturb*iwindow+nPhases_to_perturb] = G_initial_atw[:,nPhases_to_perturb*(iwindow+1):nPhases_to_perturb*(iwindow+1)+nPhases_to_perturb] 

        for ifield,field in enumerate(opal_options.Field_To_study_list):
            diag_tweak = field.diagonal_tweak
            diag_method = field.diagonal_tweak_method
            diag_value = field.diagonal_tweak_value
            gram_schmidt = field.gram_schmidt
            nSmooth = field.nSmooth
            parProcess = field.parallel_processors
            use_G = field.use_G
            phases_to_perturb = field.phases_to_perturb
            nPhases_to_perturb = len(field.phases_to_perturb)
                    
            Initial_cond_matrix = np.zeros((nNodes,nPhases_to_perturb*nPert)) 
            M_initial = np.zeros((nNodes,nPhases_to_perturb*nPert)) # ALTER FOR VECTOR FIELDS
            M_final = np.zeros((nNodes,nPhases_to_perturb*nPert)) # ALTER FOR VECTOR FIELDS

            time_fwd_0 = time.time()

            # get begin and end time levels from the checkpoint mpml
            path = 'unperturbed/'
            iPhase = 0
            xml_checkpoint_file = get_name_of_checkpoint_xml_file(field.field_type,path,fwd_output_name,iwindow,xml_extension)
            tlevel_begin, tlevel_end = get_time_interval_from_xml(xml_checkpoint_file,opal_options.time_window)
            print "TIME ", tlevel_begin, tlevel_end
                    
            if iwindow == nwindow-1: # work out the lengths of fields
            # find the nearest node to the location of interest - need to use the mesh which is associated with the field of interest
                vtk_checkpoint_file = get_name_of_checkpoint_file(iPhase,nPhases_orig,field.field_type,path,fwd_output_name,iwindow)
                # this should be taken out the loop as there is only one loc of interest for all the fields  
                # done elsewhere 
                #location_of_interest = opal_options.functional.location_of_interest
                #node_of_interest = find_node_nearest_to_point(location_of_interest,vtk_checkpoint_file)
                    

            # create each ensemble/perturbation
            # iPert will run from 0,1,2,...,nPert-1 
            # jPert will run from 0,N,2N,... where N=parProcess
            jPert = 0
            while jPert < nPert:
            ##    for iPert in range (nPert):
                for iPert in range(jPert,jPert+parProcess):
                    strings_for_xml = []

                    #Create directory and copy msh file to directory 
                    string = create_directory('Window_' + str(iwindow) + '_Ensemble_'+field.name+'_'+str(iPert))
                    directories.append(string)
                    idir = ((nwindow-1)-iwindow)*nPert + iPert
                    string = 'cp ' + meshname_base + '.msh ./' + directories[idir]
                    os.system(string)
                    if have_wells:
#                        for file in os.listdir("../"):
                        for file in os.listdir("./"):
                            foundfile = False
                            if file.endswith(".jou"):
                                foundfile = True
                            if foundfile:
                                string = 'cp *.jou ./' + directories[idir]
                                os.system(string)
                            foundfile = False
                            if file.endswith(".bdf"):
                                foundfile = True
                            if foundfile:  
                                string = 'cp *.bdf ./' + directories[idir]
                                os.system(string)
                        
                    if len(field.pycode)>0:
                    # perturb initial condition by using python code supplied by the user - this option needs some looking at since recent changes
                        for iphase,phase in enumerate(phases_to_perturb):
                            # OPENS the Advection.mpml file 
                            # this will need a bit of work for it to run properly
                            perturb_IC_from_python_function_and_create_vtu(fwd_output_name+'_' + str(iwindow) + '.vtu', 'phase'+str(iphase+1)+'::'+field.name, 'Icond_'+field.name+str(iphase)+'.vtu', field.pycode)      
		 #                   perturb_IC_from_python_function_and_create_vtu(fwd_output_name+'_0.vtu', 'phase'+str(iphase)+'::'+field.name, 'Icond_'+field.name+str(iphase)+'.vtu', field.pycode)
                    else:
                        if field.name == 'Tracer' or field.name == 'Temperature' or field.name == 'PhaseVolumeFraction' or field.name == 'PhaseVolumeFractionSource':
                            # perturb the initial condition from a vtu file and write to a vtu file (three stages)
                            if field.initial_condition or field.source:
                                 #(1/3) create the perturbed vtu files
                                 path = 'unperturbed/'
                                 idir = ((nwindow-1)-iwindow)*nPert + iPert
                                 Initial_cond_matrix, random_seed = perturb_initial_condition(path, fwd_output_name, phases_to_perturb, nPhases_to_perturb, nPhases_orig, field.name, field.field_type, directories[idir] + '/perturbations.vtu', field.sigma, Initial_cond_matrix, iPert, gram_schmidt, nSmooth, nPert, modified_vol_frac, G_initial_atw[:, nPhases_to_perturb*iwindow:nPhases_to_perturb*iwindow + nPhases_to_perturb],use_G, random_seed,iwindow,opal_options.avoid_perturbing_near_wells)

                                 # Initial_cond_matrix, random_seed = perturb_initial_condition(path, fwd_output_name, nPhases, field.name, field.field_type, directories[iPert] + '/perturbations.vtu', field.sigma,Initial_cond_matrix,iPert,gram_schmidt,nSmooth,nPert,modified_vol_frac,G_initial,use_G,random_seed)

                                 # the following option needs updating for gs and smoothing
                                 #perturb_initial_condition_to_one_or_zero(path, fwd_output_name, nPhases, field.name, field.field_type, directories[iPert] + '/ic.vtu', iPert) 

                                 #(2/3) modify mpml file to read in IC from vtu - prepare the string - couldn't change this through libspud for some reason
                                 strings_for_xml = write_string_ic_from_file(opal_options.input_file, phases_to_perturb, nPhases_orig, field.name, strings_for_xml)

                            else:
                                print 'can only modify initial condition for velocity, temperature or pvf'
                                sys.exit(0)

                        else:
                            print "can only modify PhaseVolumeFraction, PhaseVolumeFractionSource, Velocity or Temperature"
                            sys.exit(0)


                    # modify the current_time and final_time when using windows
                    libspud.set_option('/timestepping/current_time',tlevel_begin)
                    libspud.set_option('/timestepping/finish_time', tlevel_end)                    
                    

                    # write the mpml options to file
                    new_xml = opal_options.input_file[:-5] +'_'+ field.name+'_'+str(iPert)+'.mpml' 
                    slash = '/'
                    idir = ((nwindow-1)-iwindow)*nPert + iPert
                    full_path = [cwd, directories[idir], new_xml] 
                    location = slash.join(full_path)
                    libspud.write_options(location)
                    #(3/3) modify mpml file to read in IC from vtu - write strings to mpml
                    if len(strings_for_xml)>0:
                        # bug in libspud so need to write these particular options
                        write_strings_to_xml(iPert, nPert, iwindow, nwindow, directories, new_xml, phases_to_perturb, field, strings_for_xml)

                #Prepare the commands to run the batch of simulations
                commands = []
                  
                for k in range(parProcess):
                    input_file = opal_options.input_file[:-5] +'_'+ field.name+'_'+str(jPert+k)+'.mpml'
                    string = opal_options.executable + " " + input_file
                    idir = ((nwindow-1)-iwindow)*nPert + jPert+k
                    commands.append('cd '+directories[idir] + ' && ' +string)
                    
                processes = [Popen(cmd, shell=True) for cmd in commands]
                # wait for completion
                for p in processes: p.wait()
                #Cleanup the list
                del commands[:]

                # modify time level of files from "an index local to the window to the global index
                if iwindow > 0:
                    for k in range(parProcess):
                        input_file = opal_options.input_file[:-5] +'_'+ field.name+'_'+str(jPert+k)+'.mpml'
                        string = opal_options.executable + " " + input_file
                        idir = ((nwindow-1)-iwindow)*nPert + jPert + k
                        os.chdir(directories[idir])
                        for i in reversed(range(n_in_window)):
                            ii = iwindow * (n_in_window-1) + i
                            string = 'mv ' + fwd_output_name + '_' + str(i) + '.vtu ' + fwd_output_name +  '_' + str(ii) + '.vtu '
                            os.system(string)
                        os.chdir(cwd)
               
                # if we are using G to modify the initial conditions we calculate G based on the ensemble solutions at iTime = 0 
                # as they become available
                if use_G or create_max_G_file:

                    for iPhase,phase in enumerate(phases_to_perturb):

                        if nPhases_orig == 1:
                            field_name = field.name
                        elif nPhases_orig > 1:
                            field_name = "phase" +str(iPhase+1)+ "::" + field.name
                        else:
                            print 'Not expecting nPhases =',nPhases
                            sys.exit(0)

                        # store solutions in M_initial and in DF
                        for iPert in range(jPert,jPert+parProcess):
                            jTime = iwindow * (n_in_window - 1) # initial level time in a window
                            kTime = (iwindow + 1 ) * (n_in_window - 1) # final time level in a window
                            # store differences in solution at the first time step
                            idir = ((nwindow-1)-iwindow)*nPert + iPert
                            ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(jTime) + ".vtu"      
                            solution_field, dummy        = get_field_from_vtk(ensemble, field_name)        		    
                            solution_field_unpert, dummy = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(jTime) + '.vtu', field_name)
	                
                            # if CG solutions have been projected onto a DG mesh, undo this
                            if DGified[ifield]:
                                solution_field = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, solution_field)
                                solution_field_unpert = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, solution_field_unpert)


                            M_initial[:,nPert*iPhase + iPert] = solution_field - solution_field_unpert  

                            # DF and M at the final time step
#                            idir =  iPert # the final time level is in the "first to be calculated / last chronologically" time window
#                            ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(nTime-1) + ".vtu"
#                            print "nTime", nTime, "n_in_window", n_in_window
#                            print "ensemble", ensemble
#                            print "field name", field_name
#                            solution_field        = get_field_from_vtk(ensemble, field_name)
#                            solution_field_unpert = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(nTime-1) + '.vtu', field_name)
#                            if DGified[ifield]:
#                                solution_field = transform_dg2cg(nPhases, field.field_type, fwd_output_name, iwindow, solution_field)
#                                solution_field_unpert = transform_dg2cg(nPhases, field.field_type, fwd_output_name, iwindow, solution_field)

                            if iwindow == nwindow-1: 

                                idir =  iPert # the final time level is in the "first to be calculated / last chronologically" time window
                                DF[iPhase*nPert + iPert,iwindow] = calculate_DF_explicitly(opal_options,directories,idir, fwd_output_name, beginTimeWindow+1,endTimeWindow+1,iwindow,n_in_window,functional,nPhases_orig, nPhases,iPhase) 
                            else: 
 
                                # store M_final (must move this - only done if use_G on, should always be done)
                                idir = ((nwindow-1)-iwindow)*nPert + iPert
                                ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(kTime) + ".vtu"      
                                solution_field, dummy        = get_field_from_vtk(ensemble, field_name)        		    
                                solution_field_unpert, dummy = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(kTime) + '.vtu', field_name)	                
                                # if CG solutions have been projected onto a DG mesh, undo this
                                if DGified[ifield]:
                                    solution_field = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, solution_field)
                                    solution_field_unpert = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, solution_field_unpert)

                                M_final[:,nPert*iPhase + iPert] = solution_field - solution_field_unpert 

                                # Does this DF need modifying?
                                DF[nPert*iPhase : nPert*iPhase + iPert+1,iwindow] = DF[iPhase*nPert + iPert,iwindow] + np.dot(np.transpose(M_final[:, nPert*iPhase : nPert*iPhase + iPert+1]) , G_initial_atw[:,nPhases_to_perturb*(iwindow+1) + iPhase]) #+ calculate_DF_explicitly(opal_options,directories,idir, fwd_output_name, beginTimeWindow+1,endTimeWindow+1,iwindow,n_in_window,functional,nPhases_orig, nPhases,iPhase) 

                        jjPert = jPert + parProcess
                        #print "DF[iPhase*nPert:iPhase*nPert+jjPert,iwindow]", DF[iPhase*nPert:iPhase*nPert+jjPert,iwindow]
                        if jPert > field.threshold_use_G:
                            G_initial_atw[:,nPhases_to_perturb*iwindow + iPhase], dummy, dummy = calculate_G_non_chaotic(M_initial[:,nPert*iPhase:nPert*iPhase+jjPert], DF[iPhase*nPert:iPhase*nPert+jjPert,iwindow], diag_tweak, diag_method, diag_value)
                        
                        # overwrites until end of time window / all ensembles have been run
                        
                        max_G_t0.write( "{:10.4f} {:10.4f} {:22.16e}\n".format(jjPert, iPhase, np.max(G_initial_atw[:,nPhases_to_perturb*iwindow + iPhase]) ) ) 

                #Advance the number of processes
                #print "finished ensemble number", iPert+1
                jPert +=parProcess

                # -----------------------------------------------------------------------------------------------------------------

            print "Finished running simulations.\n"
            libspud.clear_options()
            time_fwd = time_fwd + time.time() - time_fwd_0
                    
            #Allocate matrix that contains the information
            #C_matrix = np.zeros((vtk_field.shape[0], len(data_name), numb_files, field.perturbations_to_do))
            #Nodes, fields of interest, time-levels, number of ensembles
                    
            #for ensemble in range(field.perturbations_to_do):
            #    for field in range(len(data_name)):
            #        for i in range(numb_files):
            #            filename = filename_base+'_'+str(i)+'.vtu'
            #            C_matrix[:,field, i, ensemble] = get_field_from_vtk(filename, str(data_name[field]))
                    
            # -------------------------------------------------------------------------------------------------
            # Calculate importance map
                    
            M = np.zeros((nNodes, nPert)) 

            time_imp_map_0 = time.time() 

            print "Calculating importance map and writing to vtu \n"
            
            # create a vtk file from any results file which has the geometry information - find one with the mesh associated wiith the field of interest....?
            clean_vtk = get_clean_vtk_file(vtk_checkpoint_file)
            imp_map_filename = 'Importance_map'
            #it seems there is a limit of 3 dimensions... wonderful
            #    C_all = np.zeros(int(nNodes), int(nfields*nPhases), int(nTime*field.perturbations_to_do))

            #---------------------------------------------------BELOW WILL ONLY WORK FOR SCALAR FIELDS------
            # Calculate DF at final time step
# need to do this if (use_G or max_G0_file)  =  False
# Claire: check if we need the contents of this if - we have calculated DF already? 
#            if not (create_max_G_file or use_G):
#                if iwindow == nwindow-1: 
#                    iTime = nTime # want to calculate DF at the final step
#                    #prepend phase to field_name if nPhases > 1
#                    for iPhase,phase in enumerate(phases_to_perturb):
#                        field_name = prepend_phase_to_field_name(iPhase,nPhases_orig,field.name)
#                        mDim = get_num_of_field_components(field.field_type,nDim)
#
#		                # order of loops must be iDim then iPert
#                        for iDim in range(mDim):
#                            for iPert in range(nPert):
#                                idir = ((nwindow-1)-iwindow)*nPert + iPert
#
#                                DF[iPhase*nPert + iPert,iwindow] = calculate_DF_explicitly(opal_options,directories,idir,ensemble, fwd_output_name, nTime, iwindow,n_in_window,functional,nPhases_orig, nPhases,iPhase)

############################Inserted this for the second order code       	
#To calculate F at the final time step
            #F = DF[iPhase*nPert + iPert,iwindow] + solution_field_unpert[node_of_interest]

#To calculate F as the value of half of the square of the unperturbed solution, F = 0.5(psi^2) for the unperturbed values
#This works for the value of the functional evaluated at the end of time            
            # 3D ify point as paraview does this
            point = functional.location_of_interest
            if point.shape[1]==2:
            #print "3D ifying point"
                point_3d = np.zeros((1,3))
                point_3d[0,0:2] = point
                point = point_3d            

            path_to_unperturbed  = './unperturbed/' + fwd_output_name + '_' + str(nTime-1) + '.vtu'
            vtu_data = vtktools.vtu(path_to_unperturbed)
            probe = VTU_Probe(vtu_data.ugrid, point)
            f_name = prepend_phase_to_field_name(iPhase,nPhases_orig,field.name)
            solution_at_point_unpert = probe.GetField(f_name)

            if opal_options.functional.time == "all_time" and opal_options.opal_operation == 'second_order_da':
                print "ERROR! Not yet coded for calculating the functional for all time"
                sys.exit()

            if opal_options.functional.square_integrand == True:
                F = 0.5*((solution_at_point_unpert)**2)   ######## Add value to the functional. Remove for a general functional. F will not always be the unperturbed solution field value. F should be a vector of size number of perturbations
            else:
                F = solution_at_point_unpert ######## Add value to the functional. Remove for a general functional. F will not always be the unperturbed solution field value. F should be a vector of size number of perturbations

      
#To calculate Gm, which is the mapping from dm to dm_s. Gm = ((CT_C)_inv)_CT, where C is tha matrix of the perturbations, CT is C transpose and inv refers to inverse.
            if opal_options.opal_operation == 'second_order_da':
                Pm[ifield,:,:] = Initial_cond_matrix #matrix of the perturbations
                C = Pm[ifield,:,:] # C is the 2D array of pertubations for each field
                CT_C_inv = np.linalg.inv(np.matmul(np.transpose(C),C))
                CT_C_inv_CT = np.matmul(CT_C_inv,(np.transpose(C)))
                Gm[ifield, :,:] = np.transpose(CT_C_inv_CT)
###############Inserted for the second order code

            #---------------------------------------------------ABOVE WILL ONLY WORK FOR SCALAR FIELDS------
            print "time levels range in current window", iwindow*(n_in_window-1),(iwindow+1)*(n_in_window-1)
            print "loop                               ", iwindow*(n_in_window-1),(iwindow+1)*(n_in_window-1)
            for iTime in range(iwindow*(n_in_window-1),(iwindow+1)*(n_in_window-1)):
		    #  for iTime in range(0,nTime): # include/omit initial condition from loop?
                #if not (create_max_G_file or use_G):
                    DF = np.zeros((nPert*nPhases_to_perturb,nwindow))
                    endTimeWindow = (iwindow+1)*(n_in_window-1)
                    beginTimeWindow = iwindow*(n_in_window-1)
                    #print "iTime in main", iTime
                    #print "beginTimeWindow", beginTimeWindow, "endTimeWindow", endTimeWindow

## remove this code - rewritten below                    
#                # don't always need to recalculate this, but for safety we are recalculating...
#                #if not (create_max_G_file or use_G):
#                    if iwindow == nwindow-1: 

#                        endTimeWindow = (iwindow+1)*(n_in_window-1)

#                        #prepend phase to field_name if nPhases > 1
#                        for iPhase,phase in enumerate(phases_to_perturb):
#                            field_name = prepend_phase_to_field_name(iPhase,nPhases_orig,field.name)
#                            mDim = get_num_of_field_components(field.field_type,nDim)

#                            # order of loops must be iDim then iPert
#                            for iDim in range(mDim):
#                                for iPert in range(nPert):
#                                    idir = ((nwindow-1)-iwindow)*nPert + iPert 
##                                    if iPert == 0: 
##                                        print "calling DF for times ", iTime+1, endTimeWindow-1, "incl"
#                                    DF[iPhase*nPert + iPert,iwindow] = calculate_DF_explicitly(opal_options,directories,idir, fwd_output_name, iTime+1, endTimeWindow+1, iwindow,n_in_window,functional,nPhases_orig, nPhases,iPhase)
##                                    #print "check: iTime", iTime, "DF", DF[iPhase*nPert + iPert,iwindow]

#                    else:
#                        kTime = endTimeWindow 

#                        #prepend phase to field_name if nPhases > 1
#                        for iPhase,phase in enumerate(phases_to_perturb):
#                            field_name = prepend_phase_to_field_name(iPhase,nPhases_orig,field.name)
#                            mDim = get_num_of_field_components(field.field_type,nDim)

#                            # order of loops must be iDim then iPert
#                            for iDim in range(mDim):
#                                for iPert in range(nPert):
# 
#                                # store M_final (must move this - only done if use_G on, should always be done)
#                                    idir = ((nwindow-1)-iwindow)*nPert + iPert
#                                    ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(kTime) + ".vtu" 
#                                    #print "ensemble", ensemble     
#                                    solution_field, dummy        = get_field_from_vtk(ensemble, field_name)        		    
#                                    solution_field_unpert, dummy = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(kTime) + '.vtu', field_name)	                
#                                    # if CG solutions have been projected onto a DG mesh, undo this
#                                    if DGified[ifield]:
#                                        solution_field = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, solution_field)
#                                        solution_field_unpert = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, solution_field_unpert)

#                                    M_final[:,nPert*iPhase + iPert] = solution_field - solution_field_unpert 


#                                    #DF[nPert*iPhase : nPert*iPhase + iPert+1,iwindow] = np.dot(np.transpose(M_final[:, nPert*iPhase : nPert*iPhase + iPert+1]) , G_initial_atw[:,nPhases_to_perturb*(iwindow+1) + iPhase]) + calculate_DF_explicitly(opal_options,directories,idir, fwd_output_name, iTime+1, endTimeWindow+1, iwindow,n_in_window,functional,nPhases_orig, nPhases,iPhase)
#                                    DF[nPert*iPhase + iPert,iwindow] = np.dot(np.transpose(M_final[:, nPert*iPhase + iPert]) , G_initial_atw[:,nPhases_to_perturb*(iwindow+1) + iPhase]) + calculate_DF_explicitly(opal_options,directories,idir, fwd_output_name, iTime+1, endTimeWindow+1, iwindow,n_in_window,functional,nPhases_orig, nPhases,iPhase)

                                  
## new version
                # don't always need to recalculate this, but for safety we are recalculating...
                #if not (create_max_G_file or use_G):
                    #if iwindow == nwindow-1: 
                    #prepend phase to field_name if nPhases > 1
                    for iPhase,phase in enumerate(phases_to_perturb):
                        field_name = prepend_phase_to_field_name(iPhase,nPhases_orig,field.name)
                        mDim = get_num_of_field_components(field.field_type,nDim)

                        # order of loops must be iDim then iPert
                        for iDim in range(mDim):
                            for iPert in range(nPert):
                                idir = ((nwindow-1)-iwindow)*nPert + iPert 
#                                if iPert == 0: 
#                                    print "calling DF for times ", iTime+1, endTimeWindow-1, "incl"
                                #print "calling DF for final TW, iPert", iPert

                                #print "DF index", iPhase*nPert + iPert,iwindow 
                                DF[iPhase*nPert + iPert,iwindow] = calculate_DF_explicitly(opal_options,directories,idir, fwd_output_name, iTime+1, endTimeWindow+1, iwindow,n_in_window,functional,nPhases_orig, nPhases,iPhase)
#                                #print "check: iTime", iTime, "DF", DF[iPhase*nPert + iPert,iwindow]

                    if iwindow < nwindow - 1:
                        kTime = endTimeWindow 

                        #prepend phase to field_name if nPhases > 1
                        for iPhase,phase in enumerate(phases_to_perturb):
                            field_name = prepend_phase_to_field_name(iPhase,nPhases_orig,field.name)
                            mDim = get_num_of_field_components(field.field_type,nDim)

                            # order of loops must be iDim then iPert
                            for iDim in range(mDim):
                                for iPert in range(nPert):
 
                                # store M_final (must move this - only done if use_G on, should always be done)
                                    idir = ((nwindow-1)-iwindow)*nPert + iPert
                                    ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(kTime) + ".vtu" 
                                    #print "ensemble", ensemble     
                                    solution_field, dummy        = get_field_from_vtk(ensemble, field_name)        		    
                                    solution_field_unpert, dummy = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(kTime) + '.vtu', field_name)	                
                                    # if CG solutions have been projected onto a DG mesh, undo this
                                    if DGified[ifield]:
                                        solution_field = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, solution_field)
                                        solution_field_unpert = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, solution_field_unpert)

                                    M_final[:,nPert*iPhase + iPert] = solution_field - solution_field_unpert 


                                    print  "writing to row ",nPert*iPhase + iPert,"column",iwindow  
                                    DF[nPert*iPhase + iPert,iwindow] = DF[nPert*iPhase + iPert,iwindow] + np.dot(np.transpose(M_final[:, nPert*iPhase + iPert]) , G_initial_atw[:,nPhases_to_perturb*(iwindow+1) + iPhase])  
## end new version                                


		    # prepare vtu file
		    new_vtu = vtktools.vtu()
		    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
		    new_vtu.filename = imp_map_filename + '_' + str(iTime) + '.vtu'
		    # loop over fields of interest
                    for iPhase,phase in enumerate(phases_to_perturb):
		        field_name = prepend_phase_to_field_name(iPhase,nPhases_orig,field.name)
                        mDim = get_num_of_field_components(field.field_type,nDim)
                        # order of loops must be iDim then iPert
                        for iDim in range(mDim):
                            for iPert in range(nPert):
                                idir = ((nwindow-1)-iwindow)*nPert + iPert
                                ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(iTime) + ".vtu"

                                if mDim==1: 
                                    # use this method instead? 
                                    #vtu_data = vtktools.vtu(ensemble)
                                    #probe = VTU_Probe(vtu_data.ugrid, points)
                                    #solution_field = probe.GetField(f_name)

                                    #path_to_unperturbed  = './unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu'
                                    #vtu_data = vtktools.vtu(path_to_unperturbed)
                                    #probe = VTU_Probe(vtu_data.ugrid, points)
                                    #solution_field_unpert = probe.GetField(f_name)

                                    #delta_solution = solution_field - solution_field_unpert
		                    #M[:,iPert] = delta_solution

                                    solution_field, dummy        = get_field_from_vtk(ensemble, field_name)
                                    solution_field_unpert, dummy = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu', field_name)
                                    delta_solution = solution_field - solution_field_unpert
                                    if DGified[ifield]:
                                        delta_solution = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, delta_solution)
                                    M[:,iPert] = delta_solution
		                else:
                                    solution_field, dummy        = get_field_from_vtk(ensemble, field_name)
                                    solution_field_unpert, dummy = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu', field_name)
                                    delta_solution = solution_field[:,iDim] - solution_field_unpert[:,iDim]
                                    if DGified[ifield]:
                                        delta_solution = transform_dg2cg(nPhases_orig, field.field_type, fwd_output_name, iwindow, delta_solution)
                                    M[:,iPert] = delta_solution

		                #Start storing everything into the big C matrix
		                #C[:,counter_field*nPhases + iPhase, iTime, iPert] = get_field_from_vtk(ensemble, field)                
		                #This has to move to an external .py file so it is compiled with .pyc
                            G, cond_before, cond_after = calculate_G_non_chaotic(M,DF[iPhase*nPert:(iPhase+1)*nPert,iwindow],diag_tweak, diag_method, diag_value)


                            if mDim>1:  
                                new_vtu.AddScalarField(field_name + '_' + str(iDim+1),G) 
                            else:
                                new_vtu.AddScalarField(field_name,G)
                            # write condition number to file
                            cond.write( "{:10.4f} {:10.4f} {:22.16e} {:22.16e}\n".format(iwindow, iTime, cond_before,  cond_after) )                 
                    new_vtu.Write() 

                    importance_map_array[ifield,:,iTime] = G
                    Ms_array[ifield,:,:,iTime] = M



            create_pvd_file(imp_map_filename,nTime) # for different fields of interest this will be overwritten

            time_imp_map = time_imp_map + time.time() - time_imp_map_0 

    cond.close()
    max_G_t0.close()

    print "Time (s) to calculate importance map:", time_imp_map
    print "Time (s) for forward simulations:", time_fwd
    path = 'unperturbed/'
    iwindow = 0
    iPhase = 0
    field = opal_options.Field_To_study_list[0]
    vtk_checkpoint_file = get_name_of_checkpoint_file(iPhase,nPhases_orig,field.field_type,path,fwd_output_name,iwindow)
    location_of_interest = opal_options.functional.location_of_interest
    node_of_interest = find_node_nearest_to_point(location_of_interest,vtk_checkpoint_file)

    return importance_map_array, Ms_array, F, Gm, Pm, node_of_interest
