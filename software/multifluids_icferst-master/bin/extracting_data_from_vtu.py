#!/usr/bin/env python

import vtk, vtktools
import numpy as np
import sys,os

from importance_map_reduced import get_node_coords_from_vtk, get_field_from_vtk

# this file has results from a fluid flowing in a pipe (cylinder)
filename = "slug_214.vtu"
# open slug_214.vtu in Paraview and you'll see that this is one of the possible fields
fieldname = "Component1::ComponentMassFractionPhase1" 

# find the mesh coordinates from the results file called slug_214.vtu
coords = get_node_coords_from_vtk(filename)
print "shape of coords array", coords.shape
print "the 415th node has coordinates", coords[414,:]

# get the volume_fraction field from the results file called slug_214.vtu 
volume_fraction = get_field_from_vtk(filename, fieldname)
print "shape of volume_fraction array", (volume_fraction[0].shape)
print "maximum value of volume fraction", np.max(volume_fraction[0])
print "minimum value of volume fraction", np.min(volume_fraction[0])

# find the bounding box of the cylinder
print "x in", np.min(coords[:,0]), np.max(coords[:,0]) 
print "y in", np.min(coords[:,1]), np.max(coords[:,1]) 
print "z in", np.min(coords[:,2]), np.max(coords[:,2]) 

# generate a random point within the cylinder
radius = 0.5 * (np.max(coords[:,0]) - np.min(coords[:,0])) * np.random.rand()
cos_theta = np.cos(2 * np.pi * np.random.rand())
xp = radius * cos_theta 
yp = radius * np.sqrt (1 - cos_theta**2)
h = np.random.rand()
point = np.array((h,xp,yp)) # the "height" of the cylinder is along the x-axis
point = np.atleast_2d(point) # point needs to be of dimension nPoints by nDim where nPoints>=1  
#print point

# find out what the value of volume fraction is at this point
vtu_data = vtktools.vtu(filename)
probe = vtktools.VTU_Probe(vtu_data.ugrid, point)
solution_at_point = probe.GetField(fieldname)
print "at point", point, "..."
print "...the volume fraction is ", solution_at_point



