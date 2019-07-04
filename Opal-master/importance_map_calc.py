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

import vtk
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from math import *
from scipy import interpolate
from scipy.interpolate import interp1d
from vtk.util.numpy_support import vtk_to_numpy


def calculate_G_non_chaotic(M,DF,diag_tweak, diag_method, diag_value):
    "Calculates an importance map at a particular time for a particular field"
    MTM = np.dot(np.transpose(M),M)
    condition_number_before = np.linalg.cond(MTM)
    nPert = M.shape[1] 

    if diag_tweak:
        if diag_method == 'add_to_diag':
            for iPert in range(nPert):
                MTM[iPert,iPert] = MTM[iPert,iPert] + diag_value
        elif diag_method == 'constrain_evalue':
            D,V = np.linalg.eig(MTM)
            for i in range(D.shape[0]):
                if D[i]<diag_value:
                    D[i] = diag_value
            MTM = np.dot(V,np.dot(np.diag(D),np.transpose(V)))

    condition_number_after = np.linalg.cond(MTM)
    MTM_INV = np.linalg.inv(MTM)
    G = np.dot(M, np.dot(MTM_INV,DF) )
    return G, condition_number_before, condition_number_after 

def G_CALC(C, weight_nod, subtract_mean_C, add_epsilon):
#! This sub calculates the value of G at every time level of 
#! the time window AND also G_SUM. It does this by assuming the flow is chaotic and 
#! thus there are natual fluctuations that can be used to calculate 
#! sensitivities with. 
#! If SUBTRACT_MEAN_C then subtract out the mean from the concentration field 
#! as this then follows the perturbation theory a bit more may be - try both.
#! ADD_EPSILON = the no added onto diagonal diagonal for conditioning of the Mel Petros sudo inverse...
#! if -ve then use the defaut value. 
#! C is the concentration fields
#! weight_nod contains the weights of the nodes used to calculate 
#! the functional F=sum(WEIGHT_NOD,C(:,ICON,ITIME))
#! NONODS is the no of nodes/CVs in the mesh
#! NCON is the no of concentration fields solved for each time level
#! NTIM is the no of time levels. 
#! NWINDOW is the no of time levels in the time window used to calculate 
#! the senitivies (importance map) so there will be NWINDOW sensitivties e.g. 40
    nonods = C.shape[0]
    nfields = C.shape[1]
    ntime = C.shape[2]
    nsembles = C.shape[3]
    C_mean = np.zeros(nonods, nfields)#Nodes, fields of interest, time-levels, number of ensembles
    C_win = np.zeros(nonods, nsembles)
    DF = np.zeros(nsembles)
    if (subtract_mean_C):
        for ifield in range(nfields):
            for itime in range(ntime):
                C_mean[:, ifield] += C[:, ifield, itime, 0]  #for the time being I do not consider ensembles
            C_mean[:,ifield] = C_mean[:, ifield]/float(ntime)
            
    G_sum = np.zeros(nonods)
    for isemble in range(nsemble):
        i_ensemb = 0
        G = np.zeros(nonods, nsemble)
        for itime in range(nwindow, ntime):
            for ifield in range(nfield):
                C_win[:, i_ensemb] = C[:,ifield,itime,isemble] - C_mean[:,ifield]
                DF[isemble] = np.sum(weight_nod[:]*(C[:,ifield,itime,isemble])-C_mean[:,icon])
                #SHOULDN@T THESE TWO EQUATIONS ABOVE BE ADDITIVE?
        G[:,isemble] = G_window(C_win, DF, add_epsilon)
        if (isemble!=nsemble): 
            G_sum[:] += G[:,isemble]
    return G_sum
            
def g_window(C_win, DF, add_epsilon):
# This sub calculates the value of g at every time level of 
# the time window.    
    
    
    nonods = C_win.shape[0]
    nsembles = C_win.shape[1]
    
    for i in range(nsemble): 
        for j in range(nsemble):
            MTM[i,j] = np.sum(C_win[:,i]*C_win[:,j])
#Add a small no onto diagonal for conditioning of the Mel Petros sudo inverse...    
    for i in range(nsemble):
        MTM[i,i]+=add_epsilon
    
    MTM_INV = np.linalg.inv(MTM)
    G_short[:]=np.matmul(MTM_INV,DF)
    G = np.zeros(nonods)
    #Calculate g: G=K*DF with K=C_WIN^T * MTM_INV
    for i in range(nsemble):
        G+=C_win[:,i]*G_short[i]

    return G


