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
import csv
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
from importance_map import * # generate_importance_map
from second_order_da_tools import *
import time
import random
#sys.path.insert(1, '/usr/local/lib/python2.7/dist-packages/')
import libspud


#SUBROUTINE NONLINEAR_INVERSE(M_CONTROLS, F, NU, LAMBDA,GAMMA, X,NDIM, E, G, Mm, Ms, W, DELTA_T, HPRECON, FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,E,NTIME,ERROR_TOLERANCE,CONVERGING_WELL,NITS_CG,NITS_PRECON, W_RELAX, NITS_NONLINEAR, NITS_INNER)
# This subroutine solve for the controls M_CONTROLS using a Hessian eqn H*DM=-G approach.
#
# Variables:
# Gm - is the mapping from m_s to m
# Ms - the mapping of from m_s to \Psi the solution variables
# G - is the gradient
# E - no of ensambles
# HPRECON is the regularization matrix which is also used as a preconditioner.
# NTIME - no of time steps stored in the approximate Hessian calculation.
# DELTA_T - is the normalized sclaing factor for each time level stored.
# CONVERGING_WELL - =.true. if the preconditioner is 'good' and system is converging well.
# E - no of ensambles used to calculate the gradient.
# NSCALAR - no of scalars solved for in the inversion method.
# NONODS - no of FEM or CV nodes.
# FINDRM0,MIDM0,COLM0,NCOLM0 - sparcity partern of the matrix.
# X - the nodal coordinates.
# NDIM - the number of dimensions.
#
# Options for DA (to be specified in diamond):
# GAMMA - LAMBDA annealing parameter e.g. 0.95
# LAMBDA - the magnitude of the reularization.
# NU -  the magnutude of the step length damping term.
# NITS_NONLINEAR - no of non-linear iteratons to solve the invere problem
# NITS_INNER - no of inner iterations used to adjust the step length damping without running more forward models.
# if NITS_INNER = 1 then dont do any inner iterations.
# ERROR_TOLERANCE - max difference in the controls before assuming convergence.
# NITS_PRECON - no of iterations used for the SSOR preconditioner.
# W_RELAX - relaxation coefficient for the SSOR preconditioner.
# INVERSE_TOLER - inverse toerance - maximum update.
# npert_set - no of perturbations used to calculate the gradient.
#
#            INTEGER NONODS,NSCALAR,NCOLM0,E,NTIME,NITS_CG,NITS_PRECON,NDIM
#            INTEGER FINDRM0(NONODS+1),MIDM0(NONODS),COLM0(NCOLM0)
#            REAL F
#            REAL X(NDIM,NONODS)
#            REAL HPRECON(NSCALAR,NCOLM0)
#            REAL ERROR_TOLERANCE, W_RELAX(NSCALAR), LAMBDA(NSCALAR),NU(NSCALAR),GAMMA(NSCALAR)
#            REAL G(NSCALAR,NONODS), Gm(NSCALAR,NONODS,E), Ms(NSCALAR,NONODS,E,NTIME), W(NSCALAR,NONODS,NTIME), DELTA_T(NTIME)
#            LOGICAL FINISHED
# Local variables...
#            LOGICAL CONVERGING_WELL
#            INTEGER ITS
#            REAL FOLD
#            REAL ALLOCATABLE K(:,:),ML(:,:),DM(:,:),M_CONTROLS_OLD(:,:)

#            ALLOCATE(K(NSCALAR,NCOLM0),ML(NSCALAR,NCOLM0),DM(NSCALAR,NONODS),M_CONTROLS_OLD(NSCALAR,NONODS))


def second_order_da(opal_options):

    print "\nChris's baby\n"

    print "gamma", opal_options.da_gamma
    print "gamma_step", opal_options.da_gamma_step
    print "lambda", opal_options.da_lambda
    print "lambda_step", opal_options.da_lambda_step
    print "nu", opal_options.da_nu
    print "nits_CG",    opal_options.da_nits_CG
    print "nits_nonlinear",    opal_options.da_nits_nonlinear
    print "nits_inner",    opal_options.da_nits_inner
    print "error_tolerance", opal_options.da_error_tolerance
    print "nits_precon",    opal_options.da_nits_precon
    print "w_relax",    opal_options.da_w_relax
    print "inverse_toler",    opal_options.da_inverse_toler


    NSCALAR = len(opal_options.Field_To_study_list)
    E = opal_options.Field_To_study_list[0].perturbations_to_do * NSCALAR #nsembles = C.shape[3], C= concentration fields
    NONODS = get_nNodes_from_gmsh_file(meshname_base + '.msh') ###check if there's a more efficient way
    #M_CONTROLS_OLD = np.zeros((NSCALAR, NONODS)) #this is never used in the code
    M_CONTROLS = np.zeros((NSCALAR, NONODS))
    M_CONTROLS_old = np.zeros((NSCALAR, NONODS))
    DM = np.zeros((NSCALAR, NONODS))
    TOTALLY_CONVERGED = False
    ALPHA_TOLER = 0.01
    NITS_CG = opal_options.da_nits_CG
    NITS_PRECON = opal_options.da_nits_precon
    NITS_NONLINEAR = opal_options.da_nits_nonlinear
    NITS_INNER =  opal_options.da_nits_inner
    LAMBDA = opal_options.da_lambda
    LAMBDA_STEP = opal_options.da_lambda_step
    GAMMA = opal_options.da_gamma
    GAMMA_STEP = opal_options.da_gamma_step
    NU = opal_options.da_nu
    W_RELAX = opal_options.da_w_relax
    ERROR_TOLERANCE = opal_options.da_error_tolerance
    INVERSE_TOLER = opal_options.da_inverse_toler
    F = 1.E+19 # initialize to ensure it gets smaller
    #F = np.zeros(E)

    for ITS in range(NITS_NONLINEAR):
    # Calculate gradient G & functional F & Ms and Mm ...

        #if ITS == NITS_NONLINEAR-1 or TOTALLY_CONVERGED == True:
        #    npert = 0 # Just run the forward model on the last update so we have access to that.
        #else:
        #   npert = opal_options.Field_To_study_list[0].perturbations_to_do  # not sure if this block of code is needed, but can be modified to imply that the frward model is run just once for every non-linear iteration

        G, Ms, F, Gm, Perturbations, node_i = generate_importance_map(opal_options)
        #W should be passed down from the Importance map code and into the MULMAT_HESSIAN and SOLVER_CG_INVERSE routines as defined in Chris's original codes
        #Get the initial conditions from the unperturbed results

        Initial_cond_unperturbed,dummy = get_field_from_vtk("unperturbed/Advection_PressureMesh_0_checkpoint.vtu", "Temperature")
        M_CONTROLS[0,:] = Initial_cond_unperturbed
        print "Initial Conditions", M_CONTROLS

        if ITS == 0:
            G_old=G; Ms_old=Ms; F_old=F; Gm_old=Gm; M_CONTROLS_old=M_CONTROLS
            print F
        #to get the coooordintaes and information for the connecting nodes
            vtu_data, coordinates = get_coordinates("unperturbed/Advection_PressureMesh_0_checkpoint.vtu")

        else:

            accept_forward = F < F_old
            print "F", F , "F_old", F_old, "accept_forward", accept_forward

            if accept_forward:
                NU = NU*1.0
                G_old=G; Ms_old=Ms; F_old=F; Gm_old=Gm; M_CONTROLS_old=M_CONTROLS
            else:
                NU = NU*10.0
                G=G_old; Ms=Ms_old; F=F_old; Gm=Gm_old; M_CONTROLS=M_CONTROLS_old

        print "G", G
        print "accepted Initial Conditions", M_CONTROLS

        if ITS == 0:
            FINDRM0,COLM0,MIDM0,NCOLM0 = map_matrix(vtu_data,coordinates, NONODS)
            HPRECON = np.zeros((NSCALAR,NCOLM0))
            K,ML = GET_K_ML(NDIM, FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,coordinates) #NSCALAR is not needed for this function

        #Obtain the gradient that is regularized and genralized for with NSCALARS > 1, and for the time level for the inversion
        G = G_REGULARIZED(G,NONODS,FINDRM0,MIDM0,COLM0,NCOLM0,NSCALAR,LAMBDA,K,M_CONTROLS)

        if ITS != NITS_NONLINEAR or TOTALLY_CONVERGED == False: #we only ran the forward model on the last iteration.
        # DM is the model update...
        # Solve H*DM=-G for model update DM

            ZERO_GUESS = False ###use ZER0_GUESS = False when things start working
            #if F > F_old:
            #    ALPHA_TOLER = ALPHA_TOLER + 0.01 # This will adjust NU so it increases the step length damping
            #    ZERO_GUESS = True

            for ITS_INNER in range(NITS_INNER): # INNER ITERATIONS THAT DONT NEED A FUNCTIONAL EVALUATION.
        # Calculate the regularization and preconditioning matrix HPRECON
                for COUNT in range(NCOLM0):  # form new preconditioner.
                    HPRECON[:,COUNT] = (LAMBDA + LAMBDA_STEP)*K[COUNT] + NU*ML[COUNT]

                DM,CONVERGING_WELL = SOLVER_CG_INVERSE(opal_options,G, Gm, Ms, DELTA_T, HPRECON,FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,NSCALAR,E,NTIME,ERROR_TOLERANCE,NITS_CG,NITS_PRECON, W_RELAX,ZERO_GUESS,ALPHA_TOLER,node_i,DM)
                #ZERO_GUESS = False # use when things start working

                print ITS_INNER, "NU", NU, "\n", "LAMBDA_STEP", LAMBDA_STEP
                print "CONVERGING_WELL", CONVERGING_WELL

                if CONVERGING_WELL:
                    NU = NU*0.5
                else:
                    NU = NU*2.0

                print "DM", DM

            LAMBDA_STEP = GAMMA_STEP*LAMBDA_STEP # anneal down the regularization coefficient that is included only in the preconditioner for the step-length damping
            LAMBDA = GAMMA*LAMBDA # anneal down the regularization coefficient that is included in the preconditioner and the gradient
            M_CONTROLS = M_CONTROLS + DM

            if np.max(abs(DM)) < INVERSE_TOLER: ##### use absolute values or not???
                TOTALLY_CONVERGED = True

        else:
            break

       #Re shape and read M_CONTROLS to vtu files
        M_CONTROLS = np.reshape(M_CONTROLS, (NONODS,NSCALAR)) # Changed the shape of the output from rows to columns for storing in csv and vtu files
        write_vtk_from_vtk("unperturbed/Advection_PressureMesh_0_checkpoint.vtu", "M_CONTROLS", "M_CONTROLS_" + str(ITS) + ".vtu", M_CONTROLS) #To output M_CONTROLS to a vtu file
        write_vtk_from_vtk("unperturbed/Advection_PressureMesh_0_checkpoint.vtu", "Temperature", "initial_condition.vtu", M_CONTROLS) #To output M_CONTROLS to a vtu file called initial_conditions.vtu. This is used for running the forward model in the next iteration.

       #To write to a csv file
        with open('output.csv', 'w') as csvfile:
            csvwriter = csv.writer(csvfile)
            for col1 in M_CONTROLS:
                csvwriter.writerow(col1)

        #Reshape M_CONTROLS to correspond to size for subsequent non-linear iterations. M_CONTROLS was initially reshaped for storing in vtu files
        M_CONTROLS = np.reshape(M_CONTROLS, (NSCALAR,NONODS))
        print "M_controls", "ITS", ITS, M_CONTROLS

    return M_CONTROLS
