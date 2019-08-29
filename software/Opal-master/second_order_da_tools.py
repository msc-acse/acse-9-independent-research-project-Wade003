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
from importance_map_calc  import *
from importance_map_tools import *
from importance_map import *
import time
import random
#sys.path.insert(1, '/home/zdt16/MultiFluids_Dev/libspud/')
import libspud

#To load the vtu files and read the coordinates of the nodes
#vtu_data = vtktools.vtu("Input_Importance_map_1.vtu") ###Hardcoded
#coordinates = vtu_data.GetLocations()


#To get the mpml file
path = os.getcwd()
mpmlfile = get_mpml_filename(path)


libspud.load_options(mpmlfile + '.mpml')
meshname_base =  libspud.get_option('/geometry[0]/mesh[0]/from_file/file_name')
NDIM = libspud.get_option('/geometry/dimension')
Dt_n = libspud.get_option('/timestepping/timestep') #timestep
tlevel_begin = libspud.get_option('/timestepping/current_time') #start time
tlevel_end = libspud.get_option('/timestepping/finish_time') #end time
NSCALAR = libspud.option_count('/Field_to_study')
DELTA_T = libspud.get_option('/io/dump_period_in_timesteps/constant') ###change to point to checkpint file in unperturbed
libspud.clear_options()


DT = tlevel_end - tlevel_begin  #total simulation period
NTIME = int(DT/Dt_n) + 1  # number of timelevels



def get_coordinates(vtu_filename):
    vtu_data = vtktools.vtu(vtu_filename)
    coordinates = vtu_data.GetLocations()
    return vtu_data,coordinates



#def map_matrix(coordinates, NONODS):
def map_matrix(vtu_data,coordinates,NONODS):
    #This is the matrix for the compressed row storage (CSR) format
    #NCOLM0 gives the number of number of non-zero elements
    #COLM0 is the matrix of the non-zero elements (connecting nodes for each node)
    #FINDRM0 is gives the number of non-zero elements at the start of each row
    NCOLM0 = 0
    n = 0
    # get indices of connecting nodes and node itself
    for k in range(NONODS):
        connecting_nodes = vtu_data.GetPointPoints(k)
        NCOLM0 = NCOLM0 + connecting_nodes.shape[0]

    COLM0 = np.zeros(NCOLM0, dtype = "int32") #data type is specified as an integer
    FINDRM0=np.zeros(NONODS+1,dtype='int32')
    MIDM0=np.zeros(NONODS,dtype='int32')

    for k in range(NONODS):
        FINDRM0[k] = n
        connecting_nodes = vtu_data.GetPointPoints(k)
        connecting_nodes.sort()
        length_row = connecting_nodes.shape[0]
        COLM0[FINDRM0[k]:FINDRM0[k] + length_row] = connecting_nodes ###### 1  needs to be added as connecting nodes being displayed have values with 1 unit less.
        n = n + length_row
    FINDRM0[NONODS] = n

    for k in range(NONODS):
        for count in range(FINDRM0[k],FINDRM0[k+1]):
            if k == COLM0[count]:
                MIDM0[k] = count

    return FINDRM0,COLM0, MIDM0, NCOLM0



def GETVOL(V,NCOL,NDIM,nodi):
    #Calculate the approximate volume associated with a node VOLUME from a given set of vectors V across the node.
    #INTEGER NCOL,NDIM
    #REAL V(NDIM,NCOL),VOLUME
    #Local variables
    #INTEGER IDIM,JDIM
    #NCOL = Number of neighbours/connecting_nodes(FINDRM0[k+1] - FINDRM0[k]) to the node

    if NCOL <= 2:
        print 'NCOL =', NCOL,'SOMETHING HAS GONE WRONG'
        exit(0)
    VT = np.matrix.transpose(V)
    VV = np.matmul(V,VT)
    VV_det= np.linalg.det(VV) #Determinant of a matrix
    volume = np.sqrt(VV_det)*float(NDIM)/float(NCOL-1)
    return volume



def GET_K_ML(NDIM, FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,coordinates):
    #Calculate the regularization matrices K (the diffusion matrix) and ML the step length damping matrix or diagonal mass matrix.
    #INTEGER NONODS,NCOLM0,NDIM
    #INTEGER FINDRM0(NONODS+1),MIDM0(NONODS),COLM0(NCOLM0)
    #REAL X(NDIM,NONODS),K(NCOLM0),ML(NCOLM0)
    #REAL PARAMETER(TOLER=1.E-14)
    #Local variables...
    #REAL VEC(NDIM),VOLUME,MAGNITUDE
    #INTEGER COUNT,COUNT2, NODI,NODJ,NODJ2
    #REAL ALLOCATABLE V(:,:),ML_DIAG(:)
    #ALLOCATE(V(NDIM,NONODS),ML_DIAG(NONODS))
    TOLER = 1.E-14
    K = np.zeros((NCOLM0))
    ML = np.zeros((NCOLM0))
    ML_DIAG = np.zeros((NONODS))
    VEC = np.zeros((NDIM))
    X = np.zeros((NDIM,NONODS))
    Method_1 = True


    for i in range(NONODS):
        X[:,i] = coordinates[i,:][0:NDIM]

    for nodi in range(NONODS):
        num = 0
        V = np.zeros((NDIM,FINDRM0[nodi+1]-FINDRM0[nodi]))
        for count in range(FINDRM0[nodi], FINDRM0[nodi+1]):
            nodj = COLM0[count]
            VEC[:] = X[:,nodj]-X[:,nodi]
            magnitude = np.sqrt(sum(VEC*VEC))
            K[count] = -1./max(magnitude*magnitude,TOLER)
            V[:, num] = VEC[:]
            num += 1
        volume = GETVOL(V,FINDRM0[nodi+1]-FINDRM0[nodi],NDIM,nodi)
        K[FINDRM0[nodi]:FINDRM0[nodi+1]] = volume*K[FINDRM0[nodi]:FINDRM0[nodi+1]]
        K[MIDM0[nodi]] = 0.0
        ML[MIDM0[nodi]] = volume
        ML_DIAG[nodi] = volume

    if Method_1:
        #To compute the old diffusion
        #Scale K by the mean of the volume of nodes I and J
        for nodi in range(NONODS):
            for count in range(FINDRM0[nodi],FINDRM0[nodi+1]):
                nodj=COLM0[count]
                K[count] = K[count]*0.5*(ML_DIAG[nodi]+ML_DIAG[nodj]) #K is symmetric by construction.
            K[MIDM0[nodi]] = 0.0 #To make the each row of K sum to zero
            K[MIDM0[nodi]] = -sum(K[FINDRM0[nodi]:FINDRM0[nodi+1]])

    else:
        # To compute the new diffusion matrix by pre and post multiplication with the square root of the mass matrix
        sqrt_ML_DIAG = np.sqrt(ML_DIAG)
        for nodi in range(NONODS):
            for count in range(FINDRM0[nodi],FINDRM0[nodi+1]):
                K[count] = -sqrt_ML_DIAG[nodi]*sqrt_ML_DIAG[COLM0[count]]

            K[MIDM0[nodi]] = 0.0 #To make the each row of K sum to zero
            K[MIDM0[nodi]] = -sum(K[FINDRM0[nodi]:FINDRM0[nodi+1]])

    return K,ML



def MULMAT_HESSIAN(opal_options,Gm, Ms, DELTA_T, HPRECON,FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,NSCALAR,E,NTIME, P,node_i):
    #this sub performs matrix vector multiplication involving the approximate Hessian: HP = H*P
    # using the second Hessian method H = M^T W M + \Lambda K + \no ML
    #            INTEGER NONODS,NSCALAR,NCOLM0,E,NTIME
    #            INTEGER FINDRM0(NONODS+1),MIDM0(NONODS),COLM0(NCOLM0)
    #            REAL Z(NSCALAR,NONODS),R(NSCALAR,NONODS)
    #            REAL Gm(NSCALAR,NONODS,E), Ms(NSCALAR,NONODS,E,NTIME), W(NSCALAR,NONODS,NTIME), DELTA_T(NTIME)
    #            REAL HP(NSCALAR,NONODS),P(NSCALAR,NONODS)
    # Local variables...
    #            REAL ALLOCATABLE V1(:),V2(:,:),V3(:,:),V4(:,:)
    ###            INTEGER N,I,J
    # using the second Hessian method H = M^T W M + \Lambda K + \no ML
    #            INTEGER NONODS,NSCALAR,NCOLM0,E,NTIME
    #            INTEGER FINDRM0(NONODS+1),MIDM0(NONODS),COLM0(NCOLM0)
    #            REAL Z(NSCALAR,NONODS),R(NSCALAR,NONODS)
    #            REAL Gm(NSCALAR,NONODS,E), Ms(NSCALAR,NONODS,E,NTIME), W(NSCALAR,NONODS,NTIME), DELTA_T(NTIME)
    #            REAL HP(NSCALAR,NONODS),P(NSCALAR,NONODS)
    # Local variables...
    #            REAL ALLOCATABLE V1(:),V2(:,:),V3(:,:),V4(:,:)
    #            INTEGER N,I,J

    V1 = np.zeros((NSCALAR,E))
    V2 = np.zeros((NSCALAR,NONODS))
    V3 = np.zeros((NSCALAR,NONODS))
    V4 = np.zeros((NSCALAR,E))
    W = np.zeros((NSCALAR,NONODS,NTIME))
    HP = np.zeros((NSCALAR,NONODS))


    for nod in range(NONODS):
        for count in range(FINDRM0[nod],FINDRM0[nod+1]):
            HP[:,nod] = HP[:,nod] + HPRECON[:,count]*P[:,COLM0[count]]


    #To compute the weight matrix for the functional evaluated at all time
    if opal_options.functional.time == "all_time":
        W[:,node_i,:] = 1#Dt_n/DT
    elif opal_options.functional.time == "end_time":
        W[:,node_i,NTIME-1] = 1#Dt_n/DT

    for N in range(NTIME):

        for I in range(NONODS):
            for J in range(E):
                V1[:,J] = V1[:,J] + Ms[:,I,J,N]*P[:,I]

        for I in range(NONODS):
            for J in range(E):
                V2[:,I] = V2[:,I] + Gm[:,I,J]*V1[:,J]

        for I in range(NONODS):
                V3[:,I] = W[:,I,N]*V2[:,I]

        for I in range(NONODS):
            for J in range(E):
                V4[:,J] = V4[:,J] + Gm[:,I,J]*V3[:,I]

        for I in range(NONODS):
            for J in range(E):
                HP[:,I] = HP[:,I] + 1.0*DELTA_T*Ms[:,I,J,N]*V4[:,J]

    return HP



def PRECON_HESSIAN(HPRECON,FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,NSCALAR,R,NITS_PRECON, W_RELAX):
    #Solve for a preconditioner (HPRECON*Z=R) using SSOR on matrix HPRECON and solve for Z.
    #
    # Variables:
    # FINDRM0,MIDM0,COLM0,NCOLM0 contains the matrix sparcity
    # W_RELAX is the relaxation coefficient for SSOR preconditioning.
    # NITS_PRECON is the number of SSOR preconditioning iterations.
    # NSCALAR is the number of scalars.
    # NONODS is the no of FEM or CV nodes.
    #            INTEGER NONODS,NSCALAR,NCOLM0,NITS_PRECON
    #            INTEGER FINDRM0(NONODS+1),MIDM0(NONODS),COLM0(NCOLM0)
    #            REAL Z(NSCALAR,NONODS),R(NSCALAR,NONODS)
    #            REAL W_RELAX(NSCALAR)
    # Local variables...
    V = np.zeros((NSCALAR,NONODS))
    Z = np.zeros((NSCALAR,NONODS))
    ISTART = 0
    RSUM = np.zeros(NSCALAR)
    Fast=True

    for ITS_PRECON in range(NITS_PRECON):

        if Fast == True:
            for nod in range(0,NONODS): # forward sweep
                IMID = MIDM0[nod]
                RSUM[:] = - HPRECON[:,IMID]*Z[:,nod]
                for count in range(FINDRM0[nod],FINDRM0[nod+1]):
                    RSUM[:] = RSUM[:] + HPRECON[:,count]*Z[:,COLM0[count]]
                Z[:,nod] = ((1.0 - W_RELAX) * Z[:,nod]) + (W_RELAX*(R[:,nod] - RSUM[:] - V[:,nod]))/HPRECON[:,IMID]

            for nod in reversed(range(NONODS)): # forward sweep
                IMID = MIDM0[nod]
                RSUM[:] = - HPRECON[:,IMID]*Z[:,nod]
                for count in range(FINDRM0[nod],FINDRM0[nod+1]):
                    RSUM[:] = RSUM[:] + HPRECON[:,count]*Z[:,COLM0[count]]
                Z[:,nod] = ((1.0 - W_RELAX) * Z[:,nod]) + (W_RELAX*(R[:,nod] - RSUM[:] - V[:,nod]))/HPRECON[:,IMID]

        else:

            for nod in range(ISTART,NONODS): # forward sweep
                IMID = MIDM0[nod]
                RSUM[:]= 0.0
                for count in range(FINDRM0[nod],IMID):
                    RSUM[:] = RSUM[:] + HPRECON[:,count]*Z[:,COLM0[count]]
                Z[:,nod] = ((1.0 - W_RELAX) * Z[:,nod]) + (W_RELAX*(R[:,nod] - RSUM[:] - V[:,nod]))/HPRECON[:,IMID]
                V[:,nod] = RSUM[:]

                V[:,0] = 0 #######change to np.zeros((NSCALAR,NONODS))? Make sure we overwrite the results from the backwards sweep

            for nod in reversed(range(NONODS)): # backward sweep (have already just done nod NONODS)
                RSUM[:] = 0.0
                IMID = MIDM0[nod]
                for count in range(IMID+1,FINDRM0[nod+1]):
                    RSUM[:]=RSUM[:] + HPRECON[:,count]*Z[:,COLM0[count]]
                Z[:,nod] = ((1.0 - W_RELAX)* Z[:,nod]) + (W_RELAX*(R[:,nod] - RSUM[:] -V[:,nod]))/HPRECON[:,IMID]
                V[:,nod] = RSUM[:]

                V[:,NONODS-1] = 0 # Make sure we overwrite the results from the forward sweep

            ISTART = 0 ########should this be ISTART + 1???  ISTART is used so we dont try to repeat the soln on the node 1 that we have just done
    return Z



def G_REGULARIZED(G,NONODS,FINDRM0,MIDM0,COLM0,NCOLM0,NSCALAR,LAMBDA,K, M_CONTROLS):

    #This subroutine adds the regularization term to the generalized gradient. G = G - Lambda*K*M_CONTROLS
    #VSUM = np.zeros(NSCALAR) is an allocatable

    VSUM = np.zeros(NSCALAR)
    for N in range(NTIME):
        for nod in range(0,NONODS):
            VSUM[:] = 0.0
            for count in range(FINDRM0[nod],FINDRM0[nod+1]):
                VSUM[:] = VSUM[:] + LAMBDA*K[count]*M_CONTROLS[:,COLM0[count]]
            G[:,nod,N] = G[:,nod,N] + VSUM[:]

    G_inv = np.zeros((NSCALAR,NONODS))
    for I in range(NSCALAR):
        G_inv[:,:] = G[I,:,0]  #instead of 0, t_inv can be used subsequently
        G = G_inv

    return G



def SOLVER_CG_INVERSE(opal_options,G, Gm, Ms, DELTA_T, HPRECON,FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,NSCALAR,E,NTIME,ERROR_TOLERANCE,NITS_CG,NITS_PRECON, W_RELAX,ZERO_GUESS,ALPHA_TOLER,node_i,DM):
    # This subroutine solve the Hessian eqn H*DM=-G in which G is the gradient and H the approximate Hessian and solves for the
    # correction of the controls DM.
    #
    # Variables:
    # Gm - is the mapping from m_s to m
    # Ms - the mapping of from m_s to \Psi the solution variables
    # G - is the gradient
    # HPRECON is the regularization matrix which is also used as a preconditioner.
    # NTIME - no of time steps stored in the approximate Hessian calculation.
    # DELTA_T - is the normalized sclaing factor for each time level stored.
    # ERROR_TOLERANCE - max difference in the controls before assuming convergence.
    # CONVERGING_WELL - =.true. if the preconditioner is 'good' and system is converging well.
    # NITS_PRECON - no of iterations used for the SSOR preconditioner.
    # W_RELAX - relaxation coefficient for the SSOR preconditioner.
    # E - no of ensambles used to calculate the gradient.
    # NSCALAR - no of scalars solved for in the inversion method.
    # NONODS - no of FEM or CV nodes.
    # FINDRM0,MIDM0,COLM0,NCOLM0 - sparcity partern of the matrix.
    # CONVERGING_WELL IF ALPHA > ALPHA_TOLER on final iteration.
    #            INTEGER NONODS,NSCALAR,NCOLM0,E,NTIME,NITS_CG,NITS_PRECON
    #            INTEGER FINDRM0(NONODS+1),MIDM0(NONODS),COLM0(NCOLM0)
    #            REAL HPRECON(NSCALAR,NCOLM0)
    #            REAL ERROR_TOLERANCE, W_RELAX
    #            REAL G(NSCALAR,NONODS), Gm(NSCALAR,NONODS,E), Ms(NSCALAR,NONODS,E,NTIME), W(NSCALAR,NONODS,NTIME), DELTA_T(NTIME)
    #            REAL DM(NSCALAR,NONODS)
    #            LOGICAL CONVERGING_WELL
    # Local variables...
    #            REAL ALLOCATABLE P(:,:),HP(:,:),Z(:,:),R(:,:),ROLD(:,:),ZOLD(:,:)
    #            REAL ALPHA,BETA,ERROR
    #            INTEGER ITS


    P = np.zeros((NSCALAR,NONODS))
    Z = np.zeros((NSCALAR,NONODS))
    R = np.zeros((NSCALAR,NONODS))
    ROLD = np.zeros((NSCALAR,NONODS))
    ZOLD = np.zeros((NSCALAR,NONODS))
    #DM = np.zeros((NSCALAR,NONODS))

    # for initialization:
    if (ZERO_GUESS): #Zero initial guess
        DM = np.zeros((NSCALAR,NONODS))
        HP = np.zeros((NSCALAR,NONODS))
        CONVERGING_WELL = False
    else:
        HP = MULMAT_HESSIAN(opal_options,Gm, Ms, DELTA_T, HPRECON, FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,NSCALAR,E,NTIME,DM,node_i)
    R = -G - HP

    # Solve preconditioned system (HPRECON*Z=R) for Z:
    Z = PRECON_HESSIAN(HPRECON,FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,NSCALAR,R,NITS_PRECON, W_RELAX)
    P = Z


    for ITS_CG in range(NITS_CG):
    # for each iteration:
        HP = MULMAT_HESSIAN(opal_options,Gm, Ms, DELTA_T, HPRECON,FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,NSCALAR,E,NTIME,P,node_i)
        ALPHA = np.sum(Z*R)/np.sum(P*HP)
        ROLD = R
        R = R - ALPHA * HP
        ERROR = ALPHA*(np.max(P))
        DM = DM + ALPHA*P
        if ERROR < ERROR_TOLERANCE:
            break
    # Solve preconditioned system (HPRECON*Z=R)for Z:
        ZOLD = Z
        Z = PRECON_HESSIAN(HPRECON,FINDRM0,MIDM0,COLM0,NCOLM0,NONODS,NSCALAR,R,NITS_PRECON, W_RELAX)
        BETA = np.sum(Z*R)/np.sum(ZOLD*ROLD) #######changed sum to np.inner as sum returns an array instead of a scalar
        P = Z + (BETA*P)

    print "ALPHA", ALPHA

    CONVERGING_WELL = ALPHA > ALPHA_TOLER # ALPHA says how good the preconditioner is

    return DM,CONVERGING_WELL,
