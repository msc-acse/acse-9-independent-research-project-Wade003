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
import matplotlib.pyplot as plt
import numpy as np
from math import *
from scipy import interpolate
from scipy.interpolate import interp1d
from vtk.util.numpy_support import vtk_to_numpy
from vtktools import VTU_Probe
import time
import random
from opal_classes import Functional, Functional_field
import libspud

#this is a method to extract a field from a vtk file and also stores it as a numpy array
#also it seems that the coordinates can be extracted from this, not yet fully working
def get_field_from_vtk(filename, fieldname):
    "Gets a field from a vtu/pvtu file."

    # load a vtk file as input
    if filename.endswith('pvtu'):
        reader = vtk.vtkXMLPUnstructuredGridReader()
    else:   
        reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()

    #Grab a scalar from the vtk file
    my_vtk_array = reader.GetOutput().GetPointData().GetArray(fieldname)
    mesh_type = "node"
    if my_vtk_array is None:
        mesh_type = "element"
        my_vtk_array = reader.GetOutput().GetCellData().GetArray(fieldname) #self.ugrid.GetCellData().GetArray(fieldname)
        if my_vtk_array is None:
          raise Exception("ERROR: couldn't find point or cell field data with name "+fieldname+".")

    my_numpy_array = vtk_to_numpy(my_vtk_array)
    return my_numpy_array, mesh_type

def get_node_coords_from_vtk(filename):
    "Gets nodal coordinates from a vtu file (N by nDim)"
    vtu_data = vtktools.vtu(filename)
    coords = vtu_data.GetLocations()
    return coords
    
def write_vtk_from_vtk(template_file, fieldname, output_name, np_field):
    """This function reads in a vtk file. It removes all the fields to finally
       add the field np_field named as fieldname in the output file output_name"""
    #Clean the input vtu file
    clean_vtk = get_clean_vtk_file(template_file)
    # prepare vtu file
    new_vtu = vtktools.vtu()
    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
    new_vtu.filename = output_name
    #if field_type == 'scalar':
    #    new_vtu.AddScalarField(fieldname,np_field)
    #elif field_type == 'vector':
    #    new_vtu.AddVectorField(fieldname,np_field)
    new_vtu.AddField(fieldname,np_field)
    new_vtu.Write()   
    

def get_clean_vtk_file(filename):
    "Removes fields and arrays from a vtk file, leaving the coordinates/connectivity information."
    vtu_data = vtktools.vtu(filename)
    clean_vtu = vtktools.vtu()
    clean_vtu.ugrid.DeepCopy(vtu_data.ugrid)
    fieldNames = clean_vtu.GetFieldNames()
# remove all fields and arrays from this vtu
    for field in fieldNames:
        clean_vtu.RemoveField(field)
        fieldNames = clean_vtu.GetFieldNames()
        vtkdata=clean_vtu.ugrid.GetCellData()
        arrayNames = [vtkdata.GetArrayName(i) for i in range(vtkdata.GetNumberOfArrays())]
    for array in arrayNames:
        vtkdata.RemoveArray(array)
    return clean_vtu

def get_filename_and_number_files(path):
    "Counts the vtu files in a directory and finds the filename base."
    numb_files = 0
    
    for files in os.listdir(path):
        if files.endswith(".vtu"):
            pos = files.rfind('_')
            pos2 = files.rfind('.')
            if pos != -1: # to ignore the perturbations.vtu file where "_" isn't found
                filename = files[:pos]
                numb_files = max(numb_files, int(files[pos+1:pos2]))
    numb_files = numb_files + 1 # as there will be a Test_0.vtu file
    return (filename, numb_files)

def get_mpml_filename(path):
    "Finds the first mpml file in a given directory and returns its filename base."
    mpml_file_found = False
    for files in os.listdir(path):
        if files.endswith(".mpml"):
            pos = files.rfind('.')
            filename = files[:pos]
            mpml_file_found = True
            break
    if not mpml_file_found:
        print ".mpml file not found ... exiting"
        sys.exit(0)
    return filename

def get_xml_extension(input_file):
    if input_file.endswith('flml'):
        xml_extension = 'flml'
    elif input_file.endswith('mpml'):
        xml_extension = 'mpml'
    else:
        print "extension of this file not compatible with OPAL (yet):", input_file
        sys.exit(0)
    return xml_extension

def get_volume_fraction(phases_to_perturb):
    # returns the volume fractions from the original mpml
    nPhases = len(phases_to_perturb)
    vol_frac = np.zeros(nPhases)
    for iPhase in range(nPhases): 
        vol_frac[iPhase] = libspud.get_option('/material_phase['+str(iPhase)+']/scalar_field::PhaseVolumeFraction/prognostic/initial_condition/constant')   
    #assert abs(sum(vol_frac) - 1.0) < 1e-10,  'Volume fractions in original mpml do not sum to 1.0!' 
    return vol_frac

# modify vol fracs so that perturbations will keep within the saturation curve
def modify_volume_fractions(phases_to_perturb, field_list):

    vol_frac = np.zeros(len(phases_to_perturb))

    modified_vol_frac = np.zeros(len(phases_to_perturb))
    for field in field_list:
        if field.name == "PhaseVolumeFraction":
            vol_frac = get_volume_fraction(phases_to_perturb)
            modified_vol_frac = copy.deepcopy(vol_frac)
            imin = np.argmin(vol_frac)
            imax = np.argmax(vol_frac)
            
            min_value = vol_frac[imin]
            for i in range(len(vol_frac)):
                if abs(vol_frac[i] - min_value) < 1e-6:
                    modified_vol_frac[i] = modified_vol_frac[i] + field.sigma 

            max_value = vol_frac[imax]
            for i in range(len(vol_frac)):
                if abs(vol_frac[i] - max_value) < 1e-6:
                    modified_vol_frac[i] = modified_vol_frac[i] - field.sigma

    return  modified_vol_frac

def create_directories(npert):
    "Creates npert directories called ensemble_number."
    directories = []
    for i in range(npert):
        iPert = i + 1
        path = "ensemble_" + str(iPert)
        if os.path.exists(path):
            print "folder already exists", path, " exiting..."
            sys.exit(0) 
        os.mkdir( path, 0775 );
        directories.append(path)
    return directories 

def create_directory(path):
    "Creates the given directory "
    os.system('mkdir ' + path)
    #If I use anything different from os.system it complains if the folder already exists
    #subprocess.Popen('mkdir ' + path, shell=True).wait()
        
    return path

def set_up_functional(opal_options,nPhases_orig,nPhases,path,fwd_output_name):
    '''Fills the functional object with info from opal_options and by calculating some other options'''
    functional = Functional()
    functional.time = opal_options.functional.time
    point = opal_options.functional.location_of_interest
    location_of_interest = np.zeros((1,len(point)))
    for i in range(len(point)):
        location_of_interest[0,i] = point[i]
    functional.location_of_interest = location_of_interest
    functional.type = opal_options.functional.type
    for i in range (len(opal_options.functional.Functional_field)):
        temp_field = Functional_field()
        temp_field.name = opal_options.functional.Functional_field[i].name
        temp_field.field_type = opal_options.functional.Functional_field[i].field_type
        temp_field.phase = opal_options.functional.Functional_field[i].phase
        if temp_field.phase == -1:
            iphase = 0
        else:
            iphase = temp_field.phase - 1
        iwindow = 0 
        temp_field.DGified, dummy = has_cg_field_been_dgified(temp_field, fwd_output_name, iphase, nPhases, nPhases_orig, iwindow)

        functional.Functional_field.append(temp_field)

    #vtk_checkpoint_file = get_name_of_checkpoint_file(0,nPhases_orig,"scalar_field",path,fwd_output_name,iwindow)  
    #location_of_interest = opal_options.functional.location_of_interest
    #node_of_interest = locate_point_of_interest_in_mesh(location_of_interest,vtk_checkpoint_file)
    #print "node_of_interest", node_of_interest
    #sys.exit(0)

    return functional

def find_node_nearest_to_point(location_of_interest,vtk_template_path):
    """ Finds the node nearest to a given point. The nodal locations 
        are found from a vtu file."""

    p = location_of_interest
    coords = get_node_coords_from_vtk(vtk_template_path)

    # the vtu seems always to be 3D, p might be 2D or 3D so find the dimension from p
    nDim = len(p)

    # get some idea of the maximum distance between nodes
    xn= max(coords[:,0])
    yn= max(coords[:,1])
    zn= max(coords[:,2])

    x0= min(coords[:,0])
    y0= min(coords[:,1])
    z0= min(coords[:,2])

    min_dist = max([abs(xn-x0), abs(yn-y0), abs(zn-z0)])

    # check location_of_interest is inside domain
    outside = False
    # 2D case
    if p[0]<x0 or p[0]>xn:
        outside = True
    if p[1]<y0 or p[1]>yn:
        outside = True
    # 3D case 
    if nDim==3:
        if p[2]<z0 or p[2]>zn:
            outside = True
    if outside:
        print "location of interest is outside the computational domain" 
        print p
        sys.exit(0) 
    # find node nearest to the location of interest (p)
    jCoord = -1
    for iCoord in range(coords.shape[0]):
        x = coords[iCoord,0:nDim]
        distance = np.linalg.norm(x-p)
        if distance < min_dist:
            min_dist = distance
            jCoord = iCoord
    if jCoord == -1:
        print "could not find a node near the location of interest"
        sys.exit(0) 

    #print "location of interest", coords[jCoord,0:nDim]

    return jCoord
   
def calculate_DF_explicitly(opal_options,directories,idir,fwd_output_name,jTime,nTime,iwindow,n_in_window,functional,nPhases_orig,nPhases,iPhase):
    #print "calculating DF explicitly"

    #print "nTime", nTime
    #print range(jTime,nTime)

    # 3D ify point as paraview does this
    point = functional.location_of_interest
    if point.shape[1]==2:
        #print "3D ifying point"
        point_3d = np.zeros((1,3))
        point_3d[0,0:2] = point
        point = point_3d

    if functional.type == "standard":

        field = functional.Functional_field[0]

        f_name = prepend_phase_to_field_name(iPhase,nPhases_orig,field.name) #"phase1::"+field.name # HARDWIRED
        f_type = field.field_type


        if opal_options.functional.time == "end_time":

            ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(nTime-1) + ".vtu"
            vtu_data = vtktools.vtu(ensemble)
            probe = VTU_Probe(vtu_data.ugrid, point)
            solution_at_point = probe.GetField(f_name)

            path_to_unperturbed  = './unperturbed/' + fwd_output_name + '_' + str(nTime-1) + '.vtu'
            vtu_data = vtktools.vtu(path_to_unperturbed)
            probe = VTU_Probe(vtu_data.ugrid, point)
            solution_at_point_unpert = probe.GetField(f_name)

             # calculating the functional as the square of the integrand
            if opal_options.functional.square_integrand == True:
                DF = 0.5*((solution_at_point)**2 - (solution_at_point_unpert)**2) ######## Add value to the functional. Remove for a general functional. F will not always be the unperturbed solution field value. F should be a vector of size number of perturbations
	    else:
                DF = ((solution_at_point)-(solution_at_point_unpert)) ######## Add value to the functional. Remove for a general functional. F will not always be the unperturbed solution field value. F should be a vector of size number of perturbations


#            solution_field, dummy        = get_field_from_vtk(ensemble, f_name)
#            solution_field_unpert, dummy = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(nTime-1) + '.vtu', f_name)
#
#            if DGified:
#                solution_field = transform_dg2cg(nPhases_orig, f_type, fwd_output_name, iwindow, solution_field)
#                solution_field_unpert = transform_dg2cg(nPhases_orig, f_type, fwd_output_name, iwindow, solution_field_unpert)

#            DF = (solution_field[indices]-solution_field_unpert[indices])
        elif opal_options.functional.time == "all_time":

            # get delta_t
            libspud.load_options(opal_options.input_file)
            dt = libspud.get_option('/timestepping/timestep')
            libspud.clear_options()

            DF = 0.0
            print "calculate_DF_explicitly ...summing...from", jTime, "to", nTime-1 , "inclusive"
            for iTime in range(jTime,nTime):
                #print "iTime in intergration", iTime
                # mid-point time integration
                delta_t = dt 
                #if iTime == 0 or iTime == nTime-1:
                #    delta_t = 0.5*dt

                ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(iTime) + ".vtu"
                vtu_data = vtktools.vtu(ensemble)
                #print "ens:", ensemble
                #print "field names", vtu_data.GetFieldNames()
                probe = VTU_Probe(vtu_data.ugrid, point)
                solution_at_point = probe.GetField(f_name)

                path_to_unperturbed  = './unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu'
                vtu_data = vtktools.vtu(path_to_unperturbed)
                probe = VTU_Probe(vtu_data.ugrid, point)
                solution_at_point_unpert = probe.GetField(f_name)
 
                #print "DF contribution", delta_t * (solution_at_point - solution_at_point_unpert)
                
                # calculating the functional as the square of the integrand
                if opal_options.functional.square_integrand == True:
                    value = 0.5*((solution_at_point)**2 - (solution_at_point_unpert)**2) ######## Add value to the functional. Remove for a general functional. F will not always be the unperturbed solution field value. F should be a vector of size number of perturbations
	        else:
                    value = ((solution_at_point)-(solution_at_point_unpert)) ######## Add value to the functional. Remove for a general functional. F will not always be the unperturbed solution field value. F should be a vector of size number of perturbations

                DF = DF + delta_t * value

#    if opal_options.functional.type == "geothermal_over_time":
#        #  functional based on the product of velocity and temperature

#        # get delta_t
#        libspud.load_options(opal_options.input_file)
#        dt = libspud.get_option('/timestepping/timestep')
#        libspud.clear_options()

#        temperature_field = get_functional_field("emperature",functional)
#        velocity_field  = get_functional_field("elocity",functional)

#        DF = 0.0
#        #for iTime in range(nTime):
#        for iTime in range(iwindow * (n_in_window-1)+1, iwindow * (n_in_window-1)+n_in_window ):
#            print "iTime", iTime
#            # mid-point time integration
#            delta_t = dt 
#            if iTime == 0 or iTime == nTime-1:
#                delta_t = 0.5*dt 

#            ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(iTime) + ".vtu"
#            vtu_data = vtktools.vtu(ensemble)
#            probe = VTU_Probe(vtu_data.ugrid, point)
#            field_name = prepend_phase_to_field_name(temperature_field.phase-1,nPhases_orig,temperature_field.name)
#            print "ensemble", ensemble, "field, ", field_name
#            temperature = probe.GetField(field_name)
#            print "temp", temperature

#            vtu_data = vtktools.vtu(ensemble)
#            probe = VTU_Probe(vtu_data.ugrid, point)
#            field_name = prepend_phase_to_field_name(velocity_field.phase-1,nPhases_orig,velocity_field.name)
 #           print "ensemble", ensemble, "field, ", field_name
 #           velocity = probe.GetField(field_name)

#            filename = './unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu'
#            vtu_data = vtktools.vtu(filename)
#            probe = VTU_Probe(vtu_data.ugrid, point)
#            field_name = prepend_phase_to_field_name(temperature_field.phase-1,nPhases_orig,temperature_field.name)
 #           print "ensemble", filename, "field, ", field_name
#            temperature_unpert = probe.GetField(field_name)
#            print "temp_unpert", temperature_unpert

#            filename = './unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu'
#            vtu_data = vtktools.vtu(filename)
#            probe = VTU_Probe(vtu_data.ugrid, point)
 #           field_name = prepend_phase_to_field_name(velocity_field.phase-1,nPhases_orig,velocity_field.name)
#            print "ensemble", filename, "field, ", field_name
#            velocity_unpert = probe.GetField(field_name)


##            temperature        = get_field_from_vtk(ensemble, temperature_field.name)
##            temperature_unpert = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu', temperature_field.name)

##            if temperature_field.DGified:
##                print "size temperature before dgcgtrasform", temperature_field.shape
##                temperature = transform_dg2cg(nPhases_orig, temperature_field.field_type, fwd_output_name, iwindow, temperature)
##                print "size temperature after dgcgtrasform", temperature_field.shape
##                temperature_unpert = transform_dg2cg(nPhases_orig, temperature_field.field_type, fwd_output_name, iwindow, temperature_unpert)[:,1]

##            velocity        = get_field_from_vtk(ensemble, velocity_field.name)[:,1]
##            velocity_unpert = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu', velocity_field.name)
##            print "size velocity before dgcgtransform", velocity_field.shape

##            if velocity_field.DGified:
##                velocity = transform_dg2cg(nPhases_orig, velocity_field.field_type, fwd_output_name, iwindow, velocity)
##                velocity_unpert = transform_dg2cg(nPhases_orig, velocity_field.field_type, fwd_output_name, iwindow, velocity_unpert)

#            print "np.linalg.norm(velocity)", np.linalg.norm(velocity)

#            DF = DF + delta_t * ( np.linalg.norm(velocity) * temperature  - np.linalg.norm(velocity_unpert) * temperature_unpert  ) # temperature reference field cancels


    if opal_options.functional.type == "geothermal_over_time":
        #  functional based on the product of velocity and temperature

        # get delta_t
        libspud.load_options(opal_options.input_file)
        dt = libspud.get_option('/timestepping/timestep')
        libspud.clear_options()

        temperature_field = get_functional_field("emperature",functional)
        velocity_field  = get_functional_field("elocity",functional)

        DF = 0.0
        #for iTime in range(nTime):
        for iTime in range(iwindow * (n_in_window-1)+1, iwindow * (n_in_window-1)+n_in_window ):
            print "iTime", iTime
            # mid-point time integration
            delta_t = dt 
            if iTime == 0 or iTime == nTime-1:
                delta_t = 0.5*dt 

            ensemble = './'+directories[idir]+'/' + fwd_output_name +'_'+ str(iTime) + ".vtu"
            vtu_data = vtktools.vtu(ensemble)
            field_name = prepend_phase_to_field_name(temperature_field.phase-1,nPhases_orig,temperature_field.name)
            temperature_values = vtu_data.GetField(field_name)
            temperature = np.max(temperature_values)
            print "temp", temperature

            vtu_data = vtktools.vtu(ensemble)
            field_name = prepend_phase_to_field_name(velocity_field.phase-1,nPhases_orig,velocity_field.name)
            velocity_values = vtu_data.GetField(field_name)
            velocity = np.max(velocity_values[:,2])
            print "velocity", velocity

            filename = './unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu'
            vtu_data = vtktools.vtu(filename)
            print "filename", filename
            print "temperature_field.phase", temperature_field.phase  
            field_name = prepend_phase_to_field_name(temperature_field.phase-1,nPhases_orig,temperature_field.name)
            temperature_unpert_values = vtu_data.GetField(field_name)
            temperature_unpert = np.max(temperature_unpert_values)
            print "temp_unpert", temperature_unpert

            filename = './unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu'
            vtu_data = vtktools.vtu(filename)
            field_name = prepend_phase_to_field_name(velocity_field.phase-1,nPhases_orig,velocity_field.name)
            velocity_unpert_values = vtu_data.GetField(field_name)
            velocity_unpert = np.max(velocity_unpert_values[:,2])
            print "velocity_unpert", velocity_unpert


#            temperature        = get_field_from_vtk(ensemble, temperature_field.name)
#            temperature_unpert = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu', temperature_field.name)

#            if temperature_field.DGified:
#                print "size temperature before dgcgtrasform", temperature_field.shape
#                temperature = transform_dg2cg(nPhases_orig, temperature_field.field_type, fwd_output_name, iwindow, temperature)
#                print "size temperature after dgcgtrasform", temperature_field.shape
#                temperature_unpert = transform_dg2cg(nPhases_orig, temperature_field.field_type, fwd_output_name, iwindow, temperature_unpert)[:,1]

#            velocity        = get_field_from_vtk(ensemble, velocity_field.name)[:,1]
#            velocity_unpert = get_field_from_vtk('./unperturbed/' + fwd_output_name + '_' + str(iTime) + '.vtu', velocity_field.name)
#            print "size velocity before dgcgtransform", velocity_field.shape

#            if velocity_field.DGified:
#                velocity = transform_dg2cg(nPhases_orig, velocity_field.field_type, fwd_output_name, iwindow, velocity)
#                velocity_unpert = transform_dg2cg(nPhases_orig, velocity_field.field_type, fwd_output_name, iwindow, velocity_unpert)

#            print "np.linalg.norm(velocity)", np.linalg.norm(velocity)

            DF = DF + delta_t * ( velocity * (temperature - 310)  - velocity_unpert * (temperature_unpert - 310)  ) # temperature reference field cancels
 
        
    return DF        


def get_functional_field(string,f):

    found = False
    for i in range(len(f.Functional_field)):
        if string in f.Functional_field[i].name:
            field = f.Functional_field[i]
            found = True
            break
    if not found: 
        print "expecting ...",string," field but none found in the list of functional fields"
        sys.exit(0)
    return field

def perturb_IC_from_python_function_and_create_vtu(base_vtu, Field_name, output_path, pycode):
    """This function reads the code written in the diamond interface and runs it to create 
       a vtu file to be used as initial condition.
       From diamond the user has access to the field and the coordinates
       Extract field and coordinates"""
    field, dummy = get_field_from_vtk(base_vtu, Field_name)
    coordinates = get_node_coords_from_vtk(base_vtu)  
    #Run code from the user
    exec (pycode)
    #Create the input file
    write_vtk_from_vtk(base_vtu, Field_name, output_path, field)


def locate_point_of_interest_in_mesh(location_of_interest,vtk_template_path):#,field):
    """ Finds the node nearest to a given point. The nodal locations 
        are found from a vtu file."""

    p = location_of_interest
    coords = get_node_coords_from_vtk(vtk_template_path)

    # the vtu seems always to be 3D, p might be 2D or 3D so find the dimension from p
    nDim = len(p)

    # check location of interest is inside the computational domain
    outside, min_dist = is_point_outside_domain(p,coords,nDim)
    if outside:
        print "location of interest is outside the computational domain" 
        print p
        sys.exit(0) 

    # find node nearest to the location of interest 
    jCoord = find_node_nearest_to_location_of_interest(p,coords,nDim,min_dist)
    if jCoord == -1:
        print "could not find a node near the location of interest"
        sys.exit(0)

#    # find out if the field is node- or cell-based
#    dummy,mesh_type = get_field_from_vtk(vtk_template_path,field.name)

#    if mesh_type == "element":
#        print "vtk_template_path",vtk_template_path
#        print "point", p  
#        tol = 1e-10
#        vtu_data = vtktools.vtu(vtk_template_path)
#        element_list = []
#        elements = vtu_data.GetPointCells(jCoord)
#        for i in range(len(elements)):
#            element_coords = coords[vtu_data.GetCellPoints(elements[i]),:]
#            print "elements[i]", elements[i]
#            print "elements_coords", element_coords
#            print "elements_coords", element_coords.shape
#            bary_coords = get_barycentric_coords(p,element_coords)
#            print "bary coords", bary_coords 
#            keep = True
#            for j in range(len(bary_coords)): 
#                if bary_coords[j]<-tol or bary_coords[j]-1>tol:
#                    keep = False 
#            if keep:
#                element_list.append(elements[i]) 

#    #    loop over elements linked to node
#    #    get barycentric coords
#    #    if not <0 and not >1 choose that element
#        print "element list", element_list
#        sys.exit(0)

    return jCoord#, mesh_type

def get_barycentric_coords(p,coords):
    N = coords.shape[0]
    M = coords.shape[1]
    Tmatrix = np.zeros((N-1,M))
    r = p - coords[N-1,:]
    for i in range(N-1):
       Tmatrix[i,:] = coords[i,:] - coords[N-1,:]

    barycentric = np.zeros(M+1)
    barycentric[0:M] = np.dot(r,np.linalg.inv(Tmatrix))
    barycentric[M] = 1-sum(np.dot(r,np.linalg.inv(Tmatrix)))   

    return barycentric

def is_point_outside_domain(p,coords,nDim):

    # get some idea of the maximum distance between nodes
    xn= max(coords[:,0])
    yn= max(coords[:,1])
    zn= max(coords[:,2])

    x0= min(coords[:,0])
    y0= min(coords[:,1])
    z0= min(coords[:,2])

    min_dist = max([abs(xn-x0), abs(yn-y0), abs(zn-z0)])

    # check location_of_interest is inside domain
    outside = False
    # 2D case
    if p[0]<x0 or p[0]>xn:
        outside = True
    if p[1]<y0 or p[1]>yn:
        outside = True
    # 3D case 
    if nDim==3:
        if p[2]<z0 or p[2]>zn:
            outside = True

    return outside, min_dist

def find_node_nearest_to_location_of_interest(p,coords,nDim,min_dist):
    
    jCoord = -1
    for iCoord in range(coords.shape[0]):
        x = coords[iCoord,0:nDim]
        distance = np.linalg.norm(x-p)
        if distance < min_dist:
            min_dist = distance
            jCoord = iCoord
    if jCoord == -1:
        print "could not find a node near the location of interest"
        sys.exit(0) 

    #print "location of interest", coords[jCoord,0:nDim]

    return jCoord

def get_name_of_checkpoint_xml_file(field_type,path,test_name,iwindow,xml_ext):

    xml_file = path + test_name + '_' + str(iwindow) + '_checkpoint.' + xml_ext   

    return xml_file


def get_name_of_checkpoint_file(iPhase,nPhases_orig,field_type,path,test_name,iwindow):
    # the field_type determines which checkpoint file (Pressure mesh or Velocity mesh) we need
    if field_type[0:6] == 'scalar':
        meshtype = 'Pressure'
    elif field_type[0:6] == 'vector':
        meshtype = 'Velocity'
    else:
        "not expecting this field type", field_type
        sys.exit(0)

    if nPhases_orig == 1:
        vtu_file = path + test_name + '_' + meshtype + 'Mesh_' + str(iwindow) + '_checkpoint.vtu'   
#        vtu_file = path + test_name + '_' + meshtype + 'Mesh_0_checkpoint.vtu'      
    elif nPhases_orig > 1:
        vtu_file = path + test_name + '_phase' + str(iPhase+1) + '_' + meshtype + 'Mesh_' + str(iwindow)  + '_checkpoint.vtu' 
#        vtu_file = path + test_name + '_phase' + str(iPhase+1) + '_' + meshtype + 'Mesh_0_checkpoint.vtu' 
        #'../unperturbed/QuickTest_DG_phase1_VelocityMesh_0_checkpoint.vtu'
    else:
        print 'nPhases_orig=', nPhases_orig,'not yet coded'
        sys.exit(0)
    return vtu_file


#def perturb_initial_condition_to_one_or_zero(path, test_name, nPhases, fieldname, fieldtype, output_path, iPert):

#    # (could check field of interest is in this checkpoint vtu first)
#    # get a checkpoint vtu of any phase
#    iPhase = 0 
#    checkpoint_vtu = get_name_of_checkpoint_file(iPhase,nPhases_orig,fieldtype,path,test_name,iwindow)

#    #Get a clean vtu file
#    clean_vtk = get_clean_vtk_file(checkpoint_vtu)
#    new_vtu = vtktools.vtu()
#    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
#    new_vtu.filename = output_path

#    for iPhase in range(nPhases):
#        checkpoint_vtu = get_name_of_checkpoint_file(iPhase,nPhases_orig,fieldtype,path,test_name,iwindow)
#        #full_fieldname = 'phase' + str(iPhase+1) + '::' + fieldname
#        Initial_cond = get_field_from_vtk(checkpoint_vtu, fieldname)
#        #Introduce the desired perturbation

#        nDoF = Initial_cond.shape[0]
#        if Initial_cond.size == nDoF:
#            nDim = 1
#        else:
#            nDim = Initial_cond.shape[1]

#        if nDim == 1:
#            Initial_cond = np.zeros((nDoF))
#            Initial_cond[iPert] = 1.0
#        else: 
#            Initial_cond = np.zeros((nDoF,nDim))
#            for iDim in range(nDim):
#                Initial_cond[iPert,iDim] = 1.0

#        if nPhases == 1:
#            string = fieldname
#        elif nPhases>1:
#            string = 'phase'+str(iPhase+1)+'::' + fieldname

#        new_vtu.AddField(string,Initial_cond) # [:,0:2]

#    #Create the input file
#    #write_vtk_from_vtk(template_vtu, fieldname, output_path, Initial_cond)
#    
#    # prepare vtu file
#    #if field_type == 'scalar':
#    #    new_vtu.AddScalarField(fieldname,np_field)
#    #elif field_type == 'vector':
#    #    new_vtu.AddVectorField(fieldname,np_field)
#    new_vtu.Write()

def smooth(ic, filename):

    old_field = copy.deepcopy(ic)
    coordinates = get_node_coords_from_vtk(filename)  
    vtu_data = vtktools.vtu(filename)
    nCoord = coordinates.shape[0]
    for k in range(nCoord):
        # get indices of connecting nodes and node itself
        connecting_nodes = vtu_data.GetPointPoints(k)
        coord = coordinates[k,:]
        N = connecting_nodes.shape[0] - 1.0
        # assign weights to the solution at each connecting node
        # the central node (k) has weight N-1 and the surrounding N-1 nodes have weight 1.
        total = sum(old_field[connecting_nodes]) - old_field[k] + N * old_field[k]
        ic[k] = total / (2.0 * N)
    #old_field = ic # if doing more than one sweep
#  vtu_data.AddScalarField('the smooth temp. '+str(j+1),new_field)
    return ic

def perturb_initial_condition(path, test_name, phases_to_perturb, nPhases, nPhases_orig, fieldname, fieldtype, output_path, sigma, ic_matrix, iPert, gs, nSmooth, nPert, vol_frac, G_init, use_G,random_seed,iwindow,avoid_perturbing_near_wells):

    # (could check field of interest is in this checkpoint vtu first)
    # get a checkpoint vtu of any phase
    iPhase = 0
    checkpoint_vtu = get_name_of_checkpoint_file(iPhase,nPhases_orig,fieldtype,path,test_name,iwindow)

    #Get a clean vtu file
    clean_vtk = get_clean_vtk_file(checkpoint_vtu)
    new_vtu = vtktools.vtu()
    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
    new_vtu.filename = output_path

    perturb_directions = [[] for iPhase in range(len(phases_to_perturb))]
    #if fieldname == "PhaseVolumeFraction":
    #    perturb_directions = determine_type_of_python_function_for_vf(nPhases,vol_frac)
 
    # if random seed is positive, use this option to generate the same sequences of random numbers
    # useful for comparing different sets of opal simulations
    if random_seed>0:
        np.random.seed(random_seed)
        
    for iPhase,phase in enumerate(phases_to_perturb):
        perturb_dir = perturb_directions[iPhase]
        
        checkpoint_vtu = get_name_of_checkpoint_file(iPhase,nPhases_orig,fieldtype,path,test_name,iwindow)
        Initial_cond_unperturbed,dummy = get_field_from_vtk(checkpoint_vtu, fieldname)
        perturbation = copy.deepcopy(Initial_cond_unperturbed) # or set to zeros, same size as Initial_cond_unperturbed

        coords = get_node_coords_from_vtk(checkpoint_vtu)
        function_avoiding_wells = np.ones_like(perturbation)  
  
        #Introduce the desired perturbation

        nDoF = Initial_cond_unperturbed.shape[0]
        if Initial_cond_unperturbed.size == nDoF:
            nDim = 1
        else:
            nDim = Initial_cond_unperturbed.shape[1]

        # find min and max values - probably should do this independently for each dimension
        #maxval = np.max(Initial_cond)
        #minval = np.min(Initial_cond)

        if nDim>1:
            print "Not yet coded"
            sys.exit(0)
#            for x in range(nDoF):
#                mu = Initial_cond[x,0:nDim]
#                delta_mu = np.random.normal(mu, sigma, nDim) - mu

#                if mu + delta_mu - maxval > 1e-06:
#                    Initial_cond[x,0:nDim] = mu - delta_mu
#                elif minval - (mu + delta_mu) > 1e-06:
#                    Initial_cond[x,0:nDim] = mu - delta_mu
#                else:
#                    Initial_cond[x,0:nDim] = mu + delta_mu    
        elif nDim==1:
            for x in range(nDoF):
                    # this code is BAD
#                    mu = Initial_cond[x]
#                    #delta_mu = np.random.normal(mu, sigma, nDim) - mu
#                    delta_mu = np.random.uniform(-0.5, 0.5, nDim) 
#                    # temperature field is really a concentration field which we constrain to be in [0,1]

#                    if fieldname == 'Temperature': # adjusts pointwise
#                        if abs(mu) < 1e-6:
#                            perturb_dir_orig = copy.deepcopy(perturb_dir)
#                            perturb_dir = "up"
#                        elif abs(mu-1.0) < 1e-6:
#                            perturb_dir_orig = copy.deepcopy(perturb_dir)
#                            perturb_dir = "down"

#                    if perturb_dir == "up":
#                        Initial_cond[x] = abs(delta_mu)
#                    elif perturb_dir == "perturb_dir":
#                        Initial_cond[x] = - abs(delta_mu)
#                    else:
#                        Initial_cond[x] = delta_mu    

#                    if fieldname == 'Temperature': # reset 
#                        if abs(mu) < 1e-6  or abs(mu-1.0) < 1e-6:
#                            perturb_dir = perturb_dir_orig

                    mu = Initial_cond_unperturbed[x]
                    delta_mu = np.random.uniform(-0.5, 0.5, nDim) 
                    perturbation[x] = delta_mu 
                    
                    # hardwiring
                    if avoid_perturbing_near_wells:
                        in_well1 = np.abs(coords[x,0]-250)<50 and np.abs(coords[x,1]-650)<50 and coords[x,2]>-350  
                        #if in_well1:
                        #    print "in well 1 ", coords[x,:]
                        in_well2 = np.abs(coords[x,0]-800)<50 and np.abs(coords[x,1]-450)<50 and coords[x,2]>-350  
                        #if in_well2:
                        #    print "in well 2 ", coords[x,:]
                    
                        if in_well1 or in_well2:
                            function_avoiding_wells[x] = 0.   

        else:
            print 'nDim ', nDim, ' unexpected'
            sys.exit(0)

        #print "min max perturbation begin", np.min(perturbation), np.max(perturbation)

        # smoothing
        for s in range(nSmooth):
            perturbation = smooth(perturbation, checkpoint_vtu)


        # function avoiding perturbing near wells
        if avoid_perturbing_near_wells:
            perturbation = np.multiply(perturbation,function_avoiding_wells) 

        # multiply initial condition by G
        if use_G:
            G_init[:,iPhase] = abs(G_init[:,iPhase]) 
            perturbation = np.multiply(perturbation,G_init[:,iPhase]) / np.max(G_init[:,iPhase])

        # orthogonalise Initial_cond against entries in ic_matrix
        if gs:
            for col in range(iPert):    
                perturbation = perturbation - np.dot(perturbation,ic_matrix[:,iPhase*nPert + col]) * ic_matrix[:,iPhase*nPert + col] / np.dot(ic_matrix[:,iPhase*nPert + col],ic_matrix[:,iPhase*nPert + col])
            #ic_matrix[:,iPhase*nPert + iPert] = Initial_cond 

        # rescale perturbation
        b = np.max(perturbation)
        a = np.min(perturbation) 

        #print "min max perturbation 2:", np.min(perturbation), np.max(perturbation)
        ab = sigma / (2.0*(b-a))
        perturbation = ab*perturbation
        #print "a,b,sigma",a,b,sigma
        #print "min max perturbation 3:", np.min(perturbation), np.max(perturbation)
        
        if gs:
            ic_matrix[:,iPhase*nPert + iPert] = perturbation


        if nPhases_orig == 1:
            string = fieldname
        elif nPhases_orig > 1:
            string = 'phase'+str(iPhase+1)+'::' + fieldname

        Initial_cond = perturbation + Initial_cond_unperturbed  
        new_vtu.AddField(string,Initial_cond) 

        #print "min max perturbation end", np.min(perturbation), np.max(perturbation)
   
    # overwrite perturbation for one phase so vol fracs sum to 1
    if fieldname == 'PhaseVolumeFraction':
        # might only work if vol_frac is ordered (increasing) - think more about this
        if nPhases_orig == 2:
            phase_sum_to_one = 1
        elif nPhases_orig >2:
             phase_sum_to_one = nPhases - 2
        else:
             print "error - there should be more than one phase and volume fraction..."
             print "nPhases_orig=", nPhases_orig
             exit(0)        

        sum_vf = np.zeros((nDoF,1))
        for iPhase in range(nPhases_orig):
            if iPhase == phase_sum_to_one:
                break
            if nPhases_orig == 1:
                string = fieldname
            elif nPhases_orig > 1:
                string = 'phase'+str(iPhase+1)+'::' + fieldname
            sum_vf = sum_vf + new_vtu.GetField(string)

        Initial_cond = 1.0-sum_vf

        if nPhases_orig == 1:
            string = fieldname
        elif nPhases_orig > 1:
            string = 'phase'+str(phase_sum_to_one+1)+'::' + fieldname   
        new_vtu.AddField(string,Initial_cond)


    new_vtu.Write() 

    if random_seed>0:
        random_seed = random_seed + 1

    return ic_matrix, random_seed



def get_time_interval_from_xml(xml_checkpoint_file,time_window):

    libspud.load_options(xml_checkpoint_file)
    
    t0 = libspud.get_option('/timestepping/current_time')
    if time_window < 0:
        tN = libspud.get_option('/timestepping/finish_time')
    else:
        tN = t0 + time_window

    libspud.clear_options()

    return t0, tN

def create_unperturbed_initial_conditions_vtu_and_run(opal_options,cwd,field_name,modified_vol_frac,xml_ext,have_wells,phases_to_perturb):
    """Runs the unperturbed problem with checkpointing to produce 
       the unpertubed initial conditions in a vtu file. These are 
       now easy to modify / smooth etc. Also we get the unperturbed 
       solution. """

    libspud.load_options(opal_options.input_file)
    fwd_output_name =  libspud.get_option('/simulation_name')
    meshname_base =  libspud.get_option('/geometry[0]/mesh[0]/from_file/file_name')
    create_directory('unperturbed')
    new_xml =  fwd_output_name +'_unperturbed.' + xml_ext
    t0 = libspud.get_option('/timestepping/current_time')
    dt = libspud.get_option('/timestepping/timestep')
    t_final = libspud.get_option('/timestepping/finish_time')

    if opal_options.time_window == 0:
         opal_options.time_window = t_final - t0# ensuring one time window

    nwindow = np.rint( (t_final - t0) / opal_options.time_window ) 
    nwindow = nwindow.astype(np.int)
    nwindow = max(nwindow, 1)#At least 1 time-window covering the whole time

    n_in_window = np.rint( opal_options.time_window / dt ) + 1
    n_in_window = n_in_window.astype(np.int)

    nTime = np.rint( (t_final - t0) / dt ) + 1
    nTime = nTime.astype(np.int)
  
    #Big checkpoint period by default, so we only get one at the beginning
    checkpoint_period = 100000
    #Check if we are using time windows. If that is the case then ensure that dump_period is used instead of
    #dump_period_in_timesteps
    if (opal_options.time_window > 0):
        #if (libspud.have_option('io/dump_period_in_timesteps')):
        #    #libspud.copy_option('io/dump_period_in_timesteps/', 'io/dump_period/')
        #    libspud.delete_option('io/dump_period_in_timesteps')
        #    try:
        #        libspud.add_option('io/dump_period/constant/')
        #    except:
        #        pass

        #Ensure that there is a checkpoint at the beginning of every time-window
        checkpoint_period = n_in_window - 1#<= This ensures that at the beginning of each time window there will be a checkpoint, which is what we want if we have time windows - it would be nice to remove the checkpoints created at the end of time....  
        try:
            libspud.set_option('io/dump_period_in_timesteps/constant', 1)
        except:
            pass

    if libspud.have_option('io/checkpointing/checkpoint_period'):
        libspud.delete_option('io/checkpointing/checkpoint_period')
    try:
        libspud.set_option('io/checkpointing/checkpoint_period_in_dumps',checkpoint_period) 
    except:
        pass

################################################################
# do for all phases_to_perturb???????????????????????????????

    if field_name == 'PhaseVolumeFractionSource':
        # add a fix to force fluidity to add PhaseVolumeFractionSource to the checkpoint vtu
        #print libspud.get_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/scalar_field::Source/prescribed/value/constant')
        try:
            print libspud.add_option('/material_phase[0]/scalar_field::PhaseVolumeFraction/prognostic/scalar_field::Source/prescribed/consistent_interpolation')
        except:
            pass 

    # increase the minimum volume fraction and decrease the maximum by the perturbation size in order to
    # perturb (+/-) and apply smooth and orthogonalisation more easily (ie without having to worry about 
    # going outside the definition of the saturation curve)
    if field_name == 'PhaseVolumeFraction':
        for jPhase in range(len(phases_to_perturb)):
            iPhase = phases_to_perturb[jPhase]
            try: 
                libspud.set_option('/material_phase['+str(jPhase)+']/scalar_field::PhaseVolumeFraction/prognostic/initial_condition/constant',modified_vol_frac[iPhase]) 
            except:
                 pass 
    os.chdir('unperturbed')
    libspud.write_options(new_xml)

    # modify the file to include <checkpoint_at_start/>
    f = open(new_xml, "r")
    contents = f.readlines()
    f.close()
 
    for iline,line in enumerate(contents):
        if "</checkpointing>" in line:
            string = "            <checkpoint_at_start/>\n"
            contents.insert(iline, string)
            break
    
    f = open(new_xml, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

    # run the unperturbed, forward model
    string = 'cp ../' + meshname_base + '.msh ./'
    os.system(string)

    if opal_options.opal_operation == 'second_order_da':
        string = 'cp ../initial_condition.vtu ./'
        os.system(string)

    if have_wells:
        for file in os.listdir("../"):
            foundfile = False
            if file.endswith(".jou"):
                foundfile = True
            if foundfile:
                string = 'cp ../*.jou ./'
                os.system(string)
            foundfile = False
            if file.endswith(".bdf"):
                foundfile = True
            if foundfile:  
                string = 'cp ../*.bdf ./'
                os.system(string)

    string = opal_options.executable + " " + new_xml
    print "running ",string , "\n"
    os.system(string)
    os.chdir(cwd)

    return nTime, n_in_window, nwindow

def write_strings_to_xml(iPert, nPert, iwindow, nwindow, directories, new_xml,phases_to_perturb, field, strings_for_xml):
    field_name = field.name 
    field_type = field.field_type
    idir = ((nwindow-1)-iwindow)*nPert + iPert
    f = open(directories[iPert] + '/' + new_xml, "r")
    contents = f.readlines()
    f.close()

    field_string = '<' + field_type + ' name="' + field_name + '"'

    # needs amending if not ALL phases are to be modified
    if field.initial_condition:
        search_tree = ['<material_phase name=', field_string, '<initial_condition', '</initial_condition>'] 
        i = 0
        kPhase = -1 # change if not ALL phases to be modified
        for iline,line in enumerate(contents):
            # identify a material_phase
            if search_tree[0] in line:
                kPhase = kPhase + 1
                i = 1  
            if search_tree[i] in line and i<2:
                # locate "<material_phase" (i=0) and "<scalar_field name=NAME" (i=1)
                i = i + 1
            if search_tree[i] in line and i==2:
                # locates <initial condition>
                i = i + 1
                begin_ic = iline
            if search_tree[i] in line and i==3:
                # locates </initial condition>
                i = 0
                end_ic = iline

                if kPhase+1 in phases_to_perturb:
                   del contents[begin_ic+1:end_ic]
                   contents.insert(begin_ic+1,strings_for_xml[kPhase])

# search tree for source#material phase , scalar field PhaseVolumeFraction, prognostic, scalar_field source, prescribed, value(whole mesh) , constant
    if field.source:
        search_tree = ['<material_phase name=', '<scalar_field name="PhaseVolumeFraction', '<prognostic', '<scalar_field name="Source', '<prescribed', '<value name=','</value>'] 
        i = 0
        kPhase = -1 # change if not ALL phases to be modified
        for iline,line in enumerate(contents):
            # identify a material_phase
            if search_tree[0] in line:
                kPhase = kPhase + 1
                i = 1   
            if search_tree[i] in line and i<5:
                # locate "<material_phase" (i=0) and "<scalar_field name=NAME" (i=1)
                i = i + 1
            if search_tree[i] in line and i==5:
                # locates <value name=
                i = i + 1
                begin_ic = iline
            if search_tree[i] in line and i==6:
                # locates </value>
                i = 0
                end_ic = iline

                if kPhase+1 in phases_to_perturb:
                   del contents[begin_ic+1:end_ic]
                   contents.insert(begin_ic+1,strings_for_xml[kPhase])


    #print "writing to", directories[idir] + '/' + new_xml
    f = open(directories[idir] + '/' + new_xml,  "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

    return

def get_nNodes_from_gmsh_file(filename):
    f = open(filename)
    f.readline() # '$MeshFormat\n'
    f.readline() # '2.2 0 8\n'
    f.readline() # '$EndMeshFormat\n'
    f.readline() # '$Nodes\n'
    nNodes = int(f.readline()) # '8\n'
    f.close()
    return nNodes

def write_string_ic_from_file(input_file,phases_to_perturb,nPhases,field_name,strings):
    libspud.load_options(input_file)
    for iPhase,phase in enumerate(phases_to_perturb):
        s1 = '                    <from_file file_name="perturbations.vtu">\n'
        s2 = '                        <format name="vtu">\n'
        s3 = '                            <field_name>\n' 
        if nPhases == 1: 
            s4 = '                                <string_value lines="1">'+field_name+'</string_value>\n'
        elif nPhases>1:
            s4 = '                                <string_value lines="1">phase'+str(iPhase+1)+'::'+field_name+'</string_value>\n'
        else:
            print 'nPhases=',nPhases,'unexpected'
            sys.exit(0)
        s5 = '                            </field_name>\n' 
        s6 = '                        </format>\n'
        s7 = '                    </from_file>\n'
        string = "%s%s%s%s%s%s%s" % (s1,s2,s3,s4,s5,s6,s7)
        strings.append(string)
# old way of trying to modify the mpml file to read ic from vtu
#                    for iPhase in range(nPhases):
#                        try:
#                            libspud.delete_option('/material_phase['+str(iPhase)+']/vector_field::Velocity/prognostic/initial_condition::WholeMesh/constant')
#                            libspud.set_option_attribute('/material_phase['+str(iPhase)+']/vector_field::Velocity/prognostic/initial_condition::WholeMesh/from_file/file_name', 'perturbations.vtu')
#                            libspud.set_option('/material_phase['+str(iPhase)+']/vector_field::Velocity/prognostic/initial_condition::WholeMesh/from_file/format/field_name', 'phase'+str(iPhase+1)+'::Velocity')
#
#                        except libspud.SpudNewKeyWarning:
#                            pass
#                        except libspud.SpudKeyError:
#                            pass
    return strings


def create_pvd_file(filename,nTime):
    pvd_filename = filename + '.pvd'  
    pvd_file = open(pvd_filename,"w")
    pvd_file.write( '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">' )
    pvd_file.write( " \n" )
    pvd_file.write( '<Collection>' )
    pvd_file.write( " \n" )

    for iTime in range(0,nTime):
        pvd_file.write( '<DataSet part="0" timestep="' + str(float(iTime)) + ' " file="' + filename + '_' + str(iTime) + '.vtu"/>' ) 
        pvd_file.write( " \n" )

    pvd_file.write( '</Collection>' )
    pvd_file.write( " \n" )
    pvd_file.write( '</VTKFile>' )
    pvd_file.write( " \n" )
    pvd_file.close()
    return

def have_cg_fields_been_dgified(fields_list, test_name, phases_to_perturb, nPhases, nPhases_orig, iwindow):
    """ CG fields can be projected onto DG fields for visulaisation in fluidity. Detect whether
        this has happened for each field of interest by comparing a checkpoint vtu file with a 
        standard vtu"""

    nfields = len(fields_list)
    nPhases_to_perturb = len(phases_to_perturb)
    DGified = np.zeros(nfields*nPhases_to_perturb, dtype=bool)
    for ifield,field in enumerate(fields_list):
        for iPhase,phase in enumerate(phases_to_perturb):

            path = 'unperturbed/'
            filename_input = get_name_of_checkpoint_file(iPhase,nPhases_orig,field.field_type,path,test_name,iwindow)

            vtu_data_input = vtktools.vtu(filename_input)
            nNodes_input = vtu_data_input.ugrid.GetNumberOfPoints() 

            vtu_data_output = vtktools.vtu('unperturbed/' + test_name + '_0.vtu') 
            nNodes_output = vtu_data_output.ugrid.GetNumberOfPoints() 

            if nNodes_input < nNodes_output:
                DGified[iPhase + nPhases_to_perturb*ifield] = True
    return DGified, nNodes_input

def has_cg_field_been_dgified(field, test_name, phase, nPhases, nPhases_orig, iwindow):
    """ CG fields can be projected onto DG fields for visulaisation in fluidity. Detect whether
        this has happened for each field of interest by comparing a checkpoint vtu file with a 
        standard vtu"""

    DGified = False
    path = 'unperturbed/'
    filename_input = get_name_of_checkpoint_file(phase,nPhases_orig,field.field_type,path,test_name,iwindow)

    vtu_data_input = vtktools.vtu(filename_input)
    nNodes_input = vtu_data_input.ugrid.GetNumberOfPoints() 

    vtu_data_output = vtktools.vtu('unperturbed/' + test_name + '_0.vtu') 
    nNodes_output = vtu_data_output.ugrid.GetNumberOfPoints() 

    if nNodes_input < nNodes_output:
        DGified = True

    return DGified, nNodes_input


def transform_dg2cg(nPhases_orig, fieldtype, test_name, iwindow, dg_field):
    """ Transform  a solution field from a DG mesh to its original CG mesh"""

#   eg. two_well_test_phase1_PressureMesh_1_checkpoint.vtu
    path = 'unperturbed/'
    iPhase = 0
    filename_input = get_name_of_checkpoint_file(iPhase,nPhases_orig,fieldtype,path,test_name,iwindow)

    vtu_data_input = vtktools.vtu(filename_input)
    nNodes_input = vtu_data_input.ugrid.GetNumberOfPoints()

    vtu_data_output = vtktools.vtu('unperturbed/' + test_name + '_0.vtu') 
    nNodes_output = vtu_data_output.ugrid.GetNumberOfPoints() 

    nElements = vtu_data_output.ugrid.GetNumberOfCells()

    cg_field = np.zeros((nNodes_input))
    cg_count = np.zeros((nNodes_input), dtype=np.int)

    for iElement in range(nElements):
        cg_node_numbers = vtu_data_input.GetCellPoints(iElement)
        dg_node_numbers = vtu_data_output.GetCellPoints(iElement)
        for i in range(len(cg_node_numbers)):
            nod_cg = cg_node_numbers[i]
            nod_dg = dg_node_numbers[i]
            cg_field[nod_cg] = cg_field[nod_cg] + dg_field[nod_dg]
            cg_count[nod_cg] = cg_count[nod_cg] + 1
      
    cg_field = np.divide(cg_field,cg_count)

    return cg_field

def prepend_phase_to_field_name(iPhase,nPhases_orig,field_name):

    if nPhases_orig == 1:
        phase_field_name = field_name
    elif nPhases_orig > 1:
        phase_field_name = "phase" +str(iPhase+1)+ "::" + field_name
    else:
        print 'Not expecting nPhases =',nPhases
        sys.exit(0)

    return phase_field_name

def get_num_of_field_components(field_type,nDim):
    if field_type == 'scalar_field':
        mDim = 1
    elif field_type == 'vector_field':
        mDim = nDim
    else:
        print "field type not yet coded - sorry"
        sys.exit(0) 
    return mDim


def create_folder_to_run(opal_options, cwd, have_wells, foldername):
    """Copy the necessary files to run the problem inside the given name folder"""
    existed = os.path.isdir(foldername)
    if existed: return existed
    libspud.load_options(opal_options.input_file)
    meshname_base = libspud.get_option('/geometry[0]/mesh[0]/from_file/file_name')
    create_directory(foldername)
    ################################################################


    os.chdir(foldername)

    # run the forward model
    string = 'cp ../' + meshname_base + '.msh ./'
    os.system(string)



    # string = 'cp ../themal_porous_template.vtu ./'##ONLY FOR THE CURRENT TEST CASE
    # os.system(string)##ONLY FOR THE CURRENT TEST CASE
    string = 'cp ../' + opal_options.input_file + ' ./'
    os.system(string)


    if have_wells:
        string = 'cp ../*.jou ./'
        os.system(string)
        string = 'cp ../*.bdf ./'
        os.system(string)
        #Copy as well step files folder
        string = 'cp -r ../step_files ./'
        os.system(string)
    os.chdir(cwd)

    return existed

def run_forward_model(opal_options, cwd, foldername):
    """Runs the forward model problem inside the given name folder"""
    #Go into de folder
    os.chdir(foldername)
    string = opal_options.executable + " " + opal_options.input_file
    print
    "running ", string, "\n"
    os.system(string)
    os.chdir(cwd)

    return



