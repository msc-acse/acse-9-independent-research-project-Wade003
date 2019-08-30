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
#import libspud

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





def get_time_interval_from_xml(xml_checkpoint_file,time_window):

    libspud.load_options(xml_checkpoint_file)
    
    t0 = libspud.get_option('/timestepping/current_time')
    if time_window < 0:
        tN = libspud.get_option('/timestepping/finish_time')
    else:
        tN = t0 + time_window

    libspud.clear_options()

    return t0, tN


