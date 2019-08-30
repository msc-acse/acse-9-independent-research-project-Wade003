#!/usr/bin/env python

#    Copyright (C) 2017 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
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

#Authors: Claire Heaney, Pablo Salinas, Dunhui Xiao, Christopher Pain, Wade Song

import numpy as np
import sys
import os
import vtk, vtktools
from importance_map_tools import *
#from vtk.util.numpy_support import vtk_to_numpy
import time
#import matplotlib.pyplot as plt
#sys.path.insert(1, '/usr/local/lib/python2.7/dist-packages/')
import libspud
#import copy
import fnmatch
import f2py

from sklearn.preprocessing import MinMaxScaler
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from sklearn import model_selection
#import pickle # not recommended for keras
from sklearn.externals import joblib

from autoencoder_basic import autoencoder
os.environ['KERAS_BACKEND'] = 'tensorflow'
from keras.models import load_model

from sklearn import __version__
os.system('echo sklearn ' + str(__version__) + ' >> ml_lib_versions')


def get_nirom_settings(fwd_options, nirom_options):

    # find the results filename so the snapshots can be read in
    if fwd_options.results_filebase == '':
        input_xml = fwd_options.input_file
        print "fwd_options.path_to_results", fwd_options.path_to_results
        print "input_xml", input_xml
        libspud.load_options(fwd_options.path_to_results + "/" + input_xml)
        fwd_options.results_filebase = libspud.get_option('/simulation_name')
        print "ndim before",fwd_options.ndim
        fwd_options.ndim = libspud.get_option('/geometry/dimension')
        print "ndim after",fwd_options.ndim

        libspud.clear_options()

    if fwd_options.results_extension == '':
    # work out the results_extension and get the number of (p)vtu files
        files = os.listdir(fwd_options.path_to_results)
        nvtu = len(fnmatch.filter(files, '*.vtu'))
        npvtu = len(fnmatch.filter(files, '*.pvtu'))
        if nvtu>0 and npvtu>0:
            print "both vtu and pvtu files are present"
            print "don't know which to use ... aborting"
            sys.exit(0)
        if nvtu > 0:
            fwd_options.results_extension = "vtu"
            nTime = nvtu
        elif npvtu > 0:
            fwd_options.results_extension = "pvtu"
            nTime = npvtu#
    else:
        files = os.listdir(fwd_options.path_to_results)
        if fwd_options.results_extension == 'pvtu':
            npvtu = len(fnmatch.filter(files, '*.pvtu'))
            nTime = npvtu
        elif fwd_options.results_extension == 'vtu':
            nvtu = len(fnmatch.filter(files, '*.vtu'))
            nTime = nvtu

    fwd_options.nTime = nTime/nirom_options.snapshots.nParameters # this is the number of snapshots not the nTime


    # superfluos now as done above
    # find number of dumps (i.e. nTime) by looking in snapshots folder
    #ext = fwd_options.results_extension
    #iCount = 0
    #for file in files:
    #    if file.startswith(fwd_options.results_filebase) and file.endswith('.vtu'):
    #       iCount = iCount + 1
    #    elif file.startswith(fwd_options.results_filebase) and file.endswith('.pvtu'):
    #       iCount = iCount + 1
    #nTime = iCount
    #fwd_options.nTime = nTime

    # finding nNodes from a snapshot
    ext = fwd_options.results_extension
    # make more general than this
    #filename = fwd_options.path_to_results + '/' + fwd_options.results_filebase + '_0.' + ext
    files = os.listdir(fwd_options.path_to_results)
    for file in files:
        if file.startswith(fwd_options.results_filebase) and file.endswith(ext):
            filename = fwd_options.path_to_results + "/" + file
            break

    vtu_data = vtktools.vtu(filename)
    for ii in range(len(nirom_options.compression.field)):
        my_field = vtu_data.GetField(nirom_options.compression.field[ii])
        nNodes = len(my_field) # use shape[0]?
        nirom_options.compression.nNodes.append(nNodes)

    return fwd_options, nirom_options

def read_in_snapshots(fwd_options, nirom_options):

    nTime = fwd_options.nTime # nSNapshots really
    nFields = len(nirom_options.compression.field)

    filebase = fwd_options.results_filebase
    print "reading in snapshots from ", filebase, "..."

    nParameters = nirom_options.snapshots.nParameters

    t0 = time.time()

    ext = fwd_options.results_extension

    # it might be more efficient to open one vtu file and then get all the required fields...?
    for iField in range(nFields):
      nNodes = nirom_options.compression.nNodes[iField]

      for p in range(nParameters):


        for iTime in range(nTime):
            if nParameters == 1:
                file_name = fwd_options.path_to_results + '/' + filebase + '_' + str(iTime) + '.' + ext
            elif nParameters>1:
                file_name = fwd_options.path_to_results + '/' + filebase + '_param_' + str(p) + "_time_" + str(iTime) + '.' + ext
            values,dummy = get_field_from_vtk(file_name, nirom_options.compression.field[iField])
            #print "getting values from  ", file_name

            # is there a better way to do this?
            if iTime==0 and p==0:
                if values.size != values.shape[0]:
                    nComponents = fwd_options.ndim #values.shape[1]
                    field_values = values
                else:
                    nComponents = 1
                    field_values = np.atleast_2d(values).T

                #nComp = nirom_options.compression.nComp
                snapshots = np.zeros((nComponents*nNodes,nTime*nParameters))

            if nComponents > 1:
                field_values = values
            else:
                nComponents = 1
                field_values = np.atleast_2d(values).T

            #print "saving data from ", file_name
            for j in range(nComponents):
                #print "filling column", iTime+p*nTime
                snapshots[j*nNodes:(j+1)*nNodes,iTime+p*nTime] = field_values[:,j]

      nirom_options.snapshots.values.append(snapshots)



#    # it might be more efficient to open one vtu file and then get all the required fields...?
#    files = os.listdir(fwd_options.path_to_results)
#    for iField in range(nFields):
#        nNodes = nirom_options.compression.nNodes[iField]
#        iCount = 0
#        #for iTime in range(nTime):
#        for file in (files):
#            new_data_found = False
#            if file.startswith(fwd_options.results_filebase) and file.endswith(ext):
#
#                iCount = iCount + 1
#                file_name = fwd_options.path_to_results + "/" + file
#                #file_name = fwd_options.path_to_results + '/' + filebase + '_' + str(iTime) + '.' + ext
#                values,dummy = get_field_from_vtk(file_name, nirom_options.compression.field[iField])
#                print "getting values ffrom  ", file_name
#                new_data_found = True
#
#            # is there a better way to do this?
#            if iCount==1:
#                if values.size != values.shape[0]:
#                    nComponents = values.shape[1]
#                    field_values = values
#                else:
#                    nComponents = 1
#                    field_values = np.atleast_2d(values).T
#
#                #nComp = nirom_options.compression.nComp
#                snapshots = np.zeros((nComponents*nNodes,nTime))
#
#            if new_data_found:
#                if nComponents > 1:
#                    field_values = values
#                else:
#                    nComponents = 1
#                    field_values = np.atleast_2d(values).T
#
#                print "saving data from ", file_name
#                for j in range(nComponents):
#                    snapshots[j*nNodes:(j+1)*nNodes,iCount-1] = field_values[:,j]
#
#        nirom_options.snapshots.values.append(snapshots)

    for i in range(len(nirom_options.snapshots.values)):
        print "snapshots from field", i, "have shape", nirom_options.snapshots.values[i].shape
        #for j in range(nirom_options.snapshots.values[i].shape[1]):
        #    print "col", j, "has min, max, sum", np.min(nirom_options.snapshots.values[i][:,j]), np.max(nirom_options.snapshots.values[i][:,j])

    t1 = time.time()

    return nirom_options, t1-t0

def compress_snapshots(fwd_options, nirom_options):

    # truncate SVD of snapshots and get basis functions and singular values

    plot_sing_values = False

    if nirom_options.compression.svd:

        t = 0
        for i in range(len(nirom_options.compression.field)):

            # SVD truncation - percentage of information captured or number
            cumulative_tol = nirom_options.compression.cumulative_tol[i]

            t0 = time.time()
            snapshots = nirom_options.snapshots.values[i]
            Usvd, s_values, Vsvd = np.linalg.svd(snapshots)
            t1 = time.time()
            t = t + t1-t0
            #print "singular vals"
            #for ii in range(len(s_values)):
            #    print s_values[ii]

            cumulative_info = np.zeros(len(s_values))
            cumulative_info[0] = s_values[0]*s_values[0]
            for j in range(1,len(s_values)):
                cumulative_info[j] = cumulative_info[j-1] + s_values[j]*s_values[j]

            cumulative_info = cumulative_info / cumulative_info[-1]

            if nirom_options.compression.nPOD[i] == -1:
                # SVD truncation - percentage of information captured or number
                cumulative_tol = nirom_options.compression.cumulative_tol[i]
                nPOD = sum(cumulative_info <= cumulative_tol) #tolerance
                nirom_options.compression.nPOD.append(nPOD)
            else:
                nPOD = nirom_options.compression.nPOD[i]

            print "retaining", nPOD, "basis functions of a possible", len(s_values)

            # repeated code - turn into function
            #nNodes = nirom_options.compression.nNodes[i]
            #if basis_functions.shape[0] != nNodes:
            #    nComponents = int(basis_functions.shape[0] / nNodes)
            #else:
            #    nComponents = 1
            #basis_functions = np.zeros((nNodes*nComponents,nPOD))
            basis_functions = Usvd[:,:nPOD]

            nirom_options.compression.basis_functions.append(basis_functions)
            nirom_options.compression.s_values.append(s_values)

    elif nirom_options.compression.eigh:

        t = 0
        for i in range(len(nirom_options.compression.field)):

            snapshots = nirom_options.snapshots.values[i]
            print "len(nirom_options.compression.field)", len(nirom_options.compression.field)
            print "snapshots.shape", snapshots.shape


            t0 = time.time()
            STS = np.dot(np.transpose(snapshots),snapshots)
            t1 = time.time()
            t = t + t1 - t0

            eigvalues, v = np.linalg.eigh(STS)
            eigvalues =  eigvalues[::-1]
            s_values = np.sqrt(eigvalues)
            #print "singular vals"
            #for ii in range(STS.shape[0]):
            #    print s_values[ii]
            #print "e vals"
            #for ii in range(STS.shape[0]):
            #    print eigvalues[ii]

            cumulative_info = np.zeros(len(eigvalues))
            for j in range(len(eigvalues)):
                if j==0:
                    cumulative_info[j] = eigvalues[j]
                else:
                    cumulative_info[j] = cumulative_info[j-1] + eigvalues[j]

            cumulative_info = cumulative_info / cumulative_info[-1]

            if nirom_options.compression.nPOD[i] == -1:
                # SVD truncation - percentage of information captured or number
                cumulative_tol = nirom_options.compression.cumulative_tol[i]
                nPOD = sum(cumulative_info <= cumulative_tol) #tolerance
                nirom_options.compression.nPOD[i] = nPOD
            else:
                nPOD = nirom_options.compression.nPOD[i]

            print "retaining", nPOD, "basis functions of a possible", len(eigvalues)

            nAll = len(eigvalues)
            # repeated code - turn into function
            nNodes = nirom_options.compression.nNodes[i]
            if snapshots.shape[0] != nNodes:
                nComponents = int(snapshots.shape[0] / nNodes)
            else:
                nComponents = 1
            basis_functions = np.zeros((nNodes*nComponents,nPOD))
            for j in reversed(range(nAll-nPOD,nAll)):
                Av = np.dot(snapshots,v[:,j])
                basis_functions[:,nAll-j-1] = Av/np.linalg.norm(Av)

            nirom_options.compression.basis_functions.append(basis_functions)
            nirom_options.compression.s_values.append(s_values)
    elif nirom_options.compression.DD_eigh:
        print("DD_eigh")
        t = 0             # timing
        for i in range(len(nirom_options.compression.field)):
            snapshots = nirom_options.snapshots.values[i]
            print "len(nirom_options.compression.field)", len(nirom_options.compression.field)
            print "snapshots.shape", snapshots.shape
            # Do DD here
            nsplt = nirom_options.compression.num_sub_base_2

            nNodes = nirom_options.compression.nNodes[i]
            if snapshots.shape[0] != nNodes:
                nComponents = int(snapshots.shape[0] / nNodes)
            else:
                nComponents = 1
            #np.set_printoptions(threshold = sys.maxsize)
            nirom_options.compression.whichd, ncola, cola, fina = f2py.DD(nNodes, nsplt,fwd_options, nirom_options)    # DD's crutial data
            print("whichd", nirom_options.compression.whichd)

            ext = fwd_options.results_extension
            filename = fwd_options.path_to_results + '/' + fwd_options.results_filebase + '_0.' + ext
            clean_vtk = get_clean_vtk_file(filename)
            new_vtu = vtktools.vtu()
            new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
            new_vtu.filename = "whichd.vtu"
            new_vtu.AddField("whichd",nirom_options.compression.whichd)
            new_vtu.Write()


            nirom_options.compression.local_list = f2py.get_local_list(nirom_options.compression.whichd, fina, cola, nsplt, nNodes)
            print("local_list", nirom_options.compression.local_list)
            nirom_options.compression.order = f2py.reorde(nirom_options.compression.local_list, int(pow(2, nsplt)))
            n_sub = len(nirom_options.compression.order)
            num_nodes_in_sub = np.zeros(n_sub, dtype = int)
            for sub in range(n_sub):
                num_nodes_in_sub[sub] = nirom_options.compression.whichd.tolist().count(sub)  # count node number in each sudomain
            nirom_options.compression.num_nodes_in_sub = num_nodes_in_sub
            print("num_nodes_in_sub", num_nodes_in_sub)
            t0 = time.time()
            dd_snapshots = []                                  # num_domains*timestpes(201)*nComponents
            for subdomain in range(int(pow(2, nsplt))):        #nirom_options.compression.order
                dd_snapshots.append([])
                for time_step in range(snapshots.shape[1]):
                    dd_snapshots[subdomain].append([])

            print("dd_snapshots", len(dd_snapshots[0]))
            old2new = np.zeros(nNodes, dtype = int)
            be_record = np.zeros(n_sub)
            for node in range(nNodes):
                belong_to = nirom_options.compression.whichd[node]
                for subs in range(n_sub):
                    if subs < belong_to:                                 # so the place for smaller index subdomain shall be record
                        old2new[node] += num_nodes_in_sub[belong_to]     # add the number of nodes for smaller index subdomain
                    elif subs==belong_to:
                        old2new[node] += be_record[belong_to]            # add the number of nodes for itself subdomain has been recorded
                be_record[belong_to] += 1                                # one more node be record in 'belong_to' subdomain
                for component in range(nComponents):
                    for time_step in range(snapshots.shape[1]):         # get shnapshots data for each subdomain
                        dd_snapshots[belong_to][time_step].append(snapshots[component*nNodes + node][time_step])
            new_dd_snapshots = np.ones_like(np.array(dd_snapshots))
            for sub in range(n_sub):
                for time_step in range(snapshots.shape[1]):
                    for node in range(num_nodes_in_sub[sub]):
                        for component in range(nComponents):
                            new_dd_snapshots[sub][time_step][node + component*num_nodes_in_sub[sub]] = dd_snapshots[sub][time_step][nComponents*node+component]
            dd_snapshots = np.array(new_dd_snapshots)
            nirom_options.compression.dd_snapshot = dd_snapshots
            print("old2new", old2new)
            nirom_options.compression.old2new = old2new
            print("snapshots.shape", snapshots.shape)
            print("dd_snapshots.shape", dd_snapshots.shape)
            #print("dd_snapshots[0][0]", len(dd_snapshots[0][0]))
            #print(dd_snapshots[0][1])
            #for i in range(int(pow(2, nsplt))):
            #    print("witchd", i, nirom_options.compression.whichd.tolist().count(i))


            for subdomain in range(int(pow(2, nsplt))):                           # do DD in each subdomain
                STS = np.dot(dd_snapshots[subdomain], np.transpose(dd_snapshots[subdomain]))        # matrix dot preparing for eigh
                eigvalues, v = np.linalg.eigh(STS)              # now we already got the eigvalues
                eigvalues =  eigvalues[::-1]
                s_values = np.sqrt(eigvalues)
                #print "singular vals"
                #for ii in range(STS.shape[0]):
                #    print s_values[ii]
                #print "e vals"
                #for ii in range(STS.shape[0]):
                #    print eigvalues[ii]
                cumulative_info = np.zeros(len(eigvalues))
                for j in range(len(eigvalues)):
                    if j==0:
                        cumulative_info[j] = eigvalues[j]
                    else:
                        cumulative_info[j] = cumulative_info[j-1] + eigvalues[j]
                cumulative_info = cumulative_info / cumulative_info[-1]
                if nirom_options.compression.nPOD[i] == -1:
                    # SVD truncation - percentage of information captured or number
                    cumulative_tol = nirom_options.compression.cumulative_tol[i]
                    nPOD = sum(cumulative_info <= cumulative_tol) #tolerance
                    nirom_options.compression.nPOD[i] = nPOD
                else:
                    nPOD = nirom_options.compression.nPOD[i]
                print "retaining", nPOD, "basis functions of a possible", len(eigvalues)
                nAll = len(eigvalues)
                basis_functions = np.zeros((dd_snapshots.shape[2],nPOD))
                #print("nAll", nAll)   201
                #print(np.transpose(dd_snapshots[1].shape))
                for j in reversed(range(nAll-nPOD,nAll)):
                    Av = np.dot(np.transpose(dd_snapshots[subdomain]),v[:,j])
                    basis_functions[:,nAll-j-1] = Av/np.linalg.norm(Av)
                print("basis_functions", basis_functions.shape)
                nirom_options.compression.basis_functions.append([basis_functions])
                nirom_options.compression.s_values.append([s_values])
                #print "orthogonal test"
                #print np.dot(basis_functions[:,2],basis_functions[:,3])
            #print("nirom_options.compression.basis_functions", nirom_options.compression.basis_functions)
            #print("nirom_options.compression.s_values", nirom_options.compression.s_values)
            nirom_options.compression.basis_functions = np.array(nirom_options.compression.basis_functions)
            nirom_options.compression.s_values = np.array(nirom_options.compression.s_values)
            print("nirom_options.compression.basis_functions", nirom_options.compression.basis_functions.shape)
            print("nirom_options.compression.s_values", nirom_options.compression.s_values.shape)
            t1 = time.time()
            t = t + t1 - t0
    else:
        print "error - compression method not recognised ... aborting"
        sys.exit(0)

    if plot_sing_values and nirom_options.compression.DD_eigh is False:
        plt.semilogy(s_values, 'rs') #ax.semilogx(t, np.exp(-t / 5.0))
        plt.grid(True)
        plt.ylabel('Singular values')
        plt.show()
        plt.close('all')

    #print "orthogonal test"
    #print np.dot(basis_functions[:,2],basis_functions[:,2])
    #print np.dot(basis_functions[:,2],basis_functions[:,3])

    # use write_vtk_from_vtk() to do this..... it needs modifying to be able to write multiple fields
    # prepare vtu file
    U_encoded = 0
    if nirom_options.compression.DD_eigh is False:
        ext =  fwd_options.results_extension
        files = os.listdir(fwd_options.path_to_results)
        for file in files:
            if file.startswith(fwd_options.results_filebase) and file.endswith(ext):
                filename = fwd_options.path_to_results + "/" + file
                break

        clean_vtk = get_clean_vtk_file(filename)
        new_vtu = vtktools.vtu()
        new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
        new_vtu.filename = "basis_functions.vtu"

        for i in range(len(nirom_options.compression.field)):
            basis_functions = nirom_options.compression.basis_functions[i]
            print "basis_functions.shape", basis_functions.shape
            nNodes = nirom_options.compression.nNodes[i]
            if basis_functions.shape[0] != nNodes:
                nComponents = int(basis_functions.shape[0] / nNodes)
                print "nComponents (integer hopefully)", nComponents
            else:
                nComponents = 1

            basis_function_field = np.zeros((nNodes,nComponents))
            nPOD = nirom_options.compression.nPOD[i]
            for k in range(nPOD):
                for j in range(nComponents):
                    basis_function_field[:,j] = basis_functions[j*nNodes:(j+1)*nNodes,k]
                new_vtu.AddField(nirom_options.compression.field[i] + '_' + str(k), basis_function_field)
        new_vtu.Write()

        print "finished writing basis functions"

        print "nirom_options.compression.svd_autoencoder ", nirom_options.compression.svd_autoencoder


        if nirom_options.compression.svd_autoencoder == True:

            ################# below copied from train_the_nirom
            nPOD_total = 0
            for ifield in range(len(nirom_options.compression.field)):
                nPOD_total = nPOD_total + nirom_options.compression.nPOD[ifield]

            nPOD_running_total = 0
            nTime = fwd_options.nTime

            pod_coeffs_all = np.zeros((nTime,nPOD_total))
            for ifield in range(len(nirom_options.compression.field)):

                snapshots = nirom_options.snapshots.values[ifield]
                basis_functions = nirom_options.compression.basis_functions[ifield]
                nPOD = nirom_options.compression.nPOD[ifield]
                nPOD_running_total = nPOD_running_total + nPOD
                # shape of pod_coeffs is nTime by number of POD coeffs
                pod_coeffs = np.dot(np.transpose(snapshots), basis_functions)
                pod_coeffs_all[:,nPOD_running_total-nPOD:nPOD_running_total] = pod_coeffs
            ################### above copied from train_the_nirom

            # training data for the autoencoder
            training_data = pod_coeffs_all #np.transpose(pod_coeffs_all)
            m = nPOD_total # 382 in this case # "number of snapshots"
            n = nTime      # number of variables (nodes, time levels etc)
            nr = nirom_options.compression.nlatent # 24 # 24 # number of recuded variables

            print "calling autoencoder...."

            U_predict, U_encoded = autoencoder(n,m,training_data,nr,nirom_options)

            print "norm U-train - U-Predict", np.linalg.norm(U_predict-training_data)

            print "n latent from options", nr
            print "training_data.shape", training_data.shape
            print " m by n", m, n

            print "training_data.shape ", training_data.shape
            print "U_predict.shape     ", U_predict.shape
            print "U_encoded.shape     ", U_encoded.shape
            print "pod_coeffs_all.shape", pod_coeffs_all.shape

    elif nirom_options.compression.DD_eigh:
        print("do this to save")
        print("nirom_options.compression.basis_functions", nirom_options.compression.basis_functions.shape)
        ext =  fwd_options.results_extension
        files = os.listdir(fwd_options.path_to_results)
        for file in files:
            if file.startswith(fwd_options.results_filebase) and file.endswith(ext):
                filename = fwd_options.path_to_results + "/" + file
                break
        clean_vtk = get_clean_vtk_file(filename)
        new_vtu = vtktools.vtu()
        print("new_vtu.ugrid.GetNumberOfCells()", new_vtu.ugrid.GetNumberOfCells())
        new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
        print("new_vtu.ugrid.GetNumberOfCells()", new_vtu.ugrid.GetNumberOfCells())
        new_vtu.filename = "basis_functions.vtu"
        basis_functions = nirom_options.compression.basis_functions
        for i in range(len(nirom_options.compression.field)):
            new_basis_function = []
            #for place_id in nirom_options.compression.order:
                #new_basis_function.append(basis_functions[place_id][i])
            #basis_functions_straight = basis_functions.reshape(basis_functions.shape[0] * basis_functions.shape[2], basis_functions.shape[3])

            for component in range(nComponents):
                for find_old_place in nirom_options.compression.old2new:
                    domain_belong_to = find_old_place/nirom_options.compression.num_nodes_in_sub[0]  # basis_functions[domain][field][data]
                    # for velocity, two that is u and v data, store seperately, in the old snapshots, it's seperated totally
                    # so first nNodes is for u and new nNodes is for v, u is the first component
                    new_basis_function.append(basis_functions[domain_belong_to][i][component*nirom_options.compression.num_nodes_in_sub[0] + find_old_place%nirom_options.compression.num_nodes_in_sub[0]])
                    #print(component*nirom_options.compression.num_nodes_in_sub[0] + find_old_place%nirom_options.compression.num_nodes_in_sub[0], component, domain_belong_to)
            #print("new_basis_function", np.array(new_basis_function).shape)
            print("new_basis_function", np.array(new_basis_function).shape)
            print(np.array(new_basis_function))
            new_basis_function = np.array(new_basis_function)
            #new_basis_function =  np.array(new_basis_function).reshape(basis_functions.shape[0] * basis_functions.shape[2], basis_functions.shape[3])
            nNodes = nirom_options.compression.nNodes[i]
            print("nNodes", nNodes)
            if new_basis_function.shape[0] != nNodes:
                nComponents = int(new_basis_function.shape[0] / nNodes)
                print "nComponents (integer hopefully)", nComponents
            else:
                nComponents = 1
            basis_function_field = np.zeros((nNodes,nComponents))
            print("basis_function_field", basis_function_field.shape)
            nPOD = nirom_options.compression.nPOD[i]
            for k in range(nPOD):
                for j in range(nComponents):
                    basis_function_field[:,j] = new_basis_function[j*nNodes:(j+1)*nNodes,k]
                print("nirom_options.compression.field[i]", nirom_options.compression.field[i])
                new_vtu.AddField(nirom_options.compression.field[i] + '_' + str(k), basis_function_field)
        new_vtu.Write()
    print "finished writing basis functions"

    return  nirom_options, t, U_encoded

def train_the_NIROM_GPR(nirom_options, fwd_options, U_encoded):

    print "training..."

    scaling_bounds = nirom_options.training.GPR_scaling
    RBF_length_scale = nirom_options.training.GPR_RBF_length_scale
    RBF_length_bounds = nirom_options.training.GPR_RBF_length_bounds
    constant_value = nirom_options.training.GPR_constant_value
    constant_bounds = nirom_options.training.GPR_constant_bounds


    if nirom_options.compression.svd_autoencoder == True:

        pod_coeffs_all = U_encoded #U_encoded.T
        nTime = fwd_options.nTime

    else:
        t = 0
        nPOD_total = 0

        for ifield in range(len(nirom_options.compression.field)):
            nPOD_total = nPOD_total + nirom_options.compression.nPOD[ifield]

        # 1. map snapshots to reduced space to give pod coefficients
        nPOD_running_total = 0
        nTime = fwd_options.nTime

        pod_coeffs_all = np.zeros((nTime*nirom_options.snapshots.nParameters,nPOD_total))
        for ifield in range(len(nirom_options.compression.field)):

            snapshots = nirom_options.snapshots.values[ifield]
            basis_functions = nirom_options.compression.basis_functions[ifield]
            nPOD = nirom_options.compression.nPOD[ifield]
            nPOD_running_total = nPOD_running_total + nPOD
            # shape of pod_coeffs is nTime by number of POD coeffs
            pod_coeffs = np.dot(np.transpose(snapshots), basis_functions)
            pod_coeffs_all[:,nPOD_running_total-nPOD:nPOD_running_total] = pod_coeffs

        #pod_coeffs_all_orig = copy.deepcopy(pod_coeffs_all)

    # 2. scale data
    print ("min max of pod coeffs before transform",np.min(pod_coeffs_all), np.max(pod_coeffs_all))
    scaling = MinMaxScaler(feature_range=scaling_bounds)
    #scaling = MinMaxScaler(feature_range=(0,10))
    pod_coeffs_all = scaling.fit_transform(pod_coeffs_all)
    print ("min max of pod coeffs after transform",np.min(pod_coeffs_all), np.max(pod_coeffs_all))

    # save the scaling settings for use in when predicting
    scaler_filename = "scaler.sav"
    joblib.dump(scaling, scaler_filename)

    # arrays for predicting the seen data (so not really prediction)
    pod_coeffs_all_prediction = np.zeros_like(pod_coeffs_all)
    pod_coeffs_all_prediction[0,:] = pod_coeffs_all[0,:]

    #for ifield in range(len(nirom_options.compression.field)):
    #    nPOD = nirom_options.compression.nPOD[ifield]
    #    nNodes = nirom_options.compression.nNodes[ifield]

    # 3. training (using all data at the moment)
    print ('training the GPR')
    t_train = 0
    t_replicate = 0
    nParameters = nirom_options.snapshots.nParameters
    trainX = np.zeros(( (nTime-1)*nParameters,nPOD_total ))
    #trainX = pod_coeffs_all[:nTime-1,:]
    for p in range(nParameters):
        trainX[p*(nTime-1):(p+1)*(nTime-1),:] = pod_coeffs_all[p*nTime:(p+1)*nTime-1,:] # n-1 by nPOD (same every iteration so could take out of loop)

    for i in range(pod_coeffs_all.shape[1]):

        t0 = time.time()

        #trainY = np.atleast_2d(pod_coeffs_all[1:,i]).T  # n-1 by 1 (different every iteration)
        trainY = np.zeros(( (nTime-1)*nParameters,1))
        for p in range(nParameters):
            #vector = np.atleast_2d(pod_coeffs_all[p*nTime+1:(p+1)*nTime,i]).T
            #print "shape output bits", vector.shape
            trainY[p*(nTime-1):(p+1)*(nTime-1),:] = np.atleast_2d(pod_coeffs_all[p*nTime+1:(p+1)*nTime,i]).T

        if i==0:
            print('shape trainX and Y', trainX.shape, trainY.shape)
            #print "training POD coeff", i+1
        #else:
        #    print "POD coeff", i+1

        kernel = C(constant_value, constant_bounds) * RBF(RBF_length_scale, RBF_length_bounds)
        #kernel = C(1.0, (1e-3, 1e3)) * RBF(100, (1e-2, 1e2))
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=100)

        # Fit to data using Maximum Likelihood Estimation of the parameters
        gp.fit(trainX, trainY)

        t1 = time.time()
        t_train = t_train + t1 - t0

        # gives a convergence warning
        #print "score of gp" , gp.score(trainX,trainY)

        # save the model
        joblib.dump(gp, 'gpr_model_' + str(i) + '.sav')

        if nParameters==1:
            # gives a warning - predicted variances < 0, setting to 0
            t0 = time.time()
            pod_pred_all, sigma = gp.predict(trainX, return_std=True)
            t1 = time.time()
            t_replicate = t_replicate + t1-t0

            # just replicating the snapshots not really predicting
            pod_coeffs_all_prediction[1:,i] = pod_pred_all[:,0]

    if nParameters==1:
        # inverse of the scaling
        pod_coeffs_all_prediction = scaling.inverse_transform(pod_coeffs_all_prediction)

        # map pod_coeffs_prediction to full space and write to file

        if nirom_options.compression.svd_autoencoder == True:

            decoder = load_model('decoder.h5')
            U_decoded = decoder.predict(pod_coeffs_all_prediction) ##### check shape of this
            print "U_decoded shape", U_decoded.shape
            pod_coeffs_all_prediction = U_decoded #U_decoded.T

        # modify write_vtk_from_vtk()
        filebase = fwd_options.results_filebase
        ext = fwd_options.results_extension
        files = os.listdir(fwd_options.path_to_results)
        for file in files:
            if file.startswith(filebase) and file.endswith(ext):
                filename = fwd_options.path_to_results + "/" + file
                break

        clean_vtk = get_clean_vtk_file(filename)
        nPOD_running_total = 0
        for ifield in range(len(nirom_options.compression.field)):
            nPOD = nirom_options.compression.nPOD[ifield]
            nPOD_running_total =  nPOD_running_total + nPOD
            nNodes = nirom_options.compression.nNodes[ifield]
            field = nirom_options.compression.field[ifield]
            basis_functions = nirom_options.compression.basis_functions[ifield]
            pod_coeffs_pred = pod_coeffs_all_prediction[:,nPOD_running_total-nPOD:nPOD_running_total]
            solution_pred = np.dot(basis_functions,np.transpose(pod_coeffs_pred))
            print "solution_pred", solution_pred.shape

            if solution_pred.shape[0] != nNodes:
                nComponents = int(solution_pred.shape[0] / nNodes)
                print "solution_pred.shape[0] / nNodes", solution_pred.shape[0] / nNodes
                print "nComponents (integer hopefully)", nComponents
            else:
                nComponents = 1


            for i in range(nTime):
                if ifield==0:
                    new_vtu = vtktools.vtu()
                    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
                    new_vtu.filename = "nirom_replication_" +str(i)+ ".vtu"
                else:
                    filename = "nirom_replication_" +str(i)+".vtu"
                    new_vtu = vtktools.vtu(filename)

                solution_pred_field = np.zeros((nNodes,nComponents))
                for j in range(nComponents):
                    solution_pred_field[:,j] = solution_pred[j*nNodes:(j+1)*nNodes,i]


                new_vtu.AddField(field,solution_pred_field)
                new_vtu.Write()


#    #for hyperparameter in kernel.hyperparameters:
#    #    print(hyperparameter)

#    #params = kernel.get_params()
#    #for key in sorted(params):
#    #    print("%s : %s" % (key, params[key]))

#    #print("theta", kernel.theta)
#    #print("bounds", kernel.bounds)


    return t_train, t_replicate

def train_the_NIROM_DD_GPR(nirom_options, fwd_options, U_encoded):

    print "DD_GPR"

    scaling_bounds = nirom_options.training.GPR_scaling
    RBF_length_scale = nirom_options.training.GPR_RBF_length_scale
    RBF_length_bounds = nirom_options.training.GPR_RBF_length_bounds
    constant_value = nirom_options.training.GPR_constant_value
    constant_bounds = nirom_options.training.GPR_constant_bounds


    if nirom_options.compression.svd_autoencoder == True:

        pod_coeffs_all = U_encoded #U_encoded.T
        nTime = fwd_options.nTime

    else:
        t = 0
        nPOD_total = 0

        for ifield in range(len(nirom_options.compression.field)):
            nPOD_total = nPOD_total + nirom_options.compression.nPOD[ifield]

        # 1. map snapshots to reduced space to give pod coefficients
        nPOD_running_total = 0
        nTime = fwd_options.nTime
        dd_pod_coeffs_all = []
        dd_snapshots = nirom_options.compression.dd_snapshot
        basis_functions = nirom_options.compression.basis_functions
        print("dd_snapshots", dd_snapshots.shape)
        print("basis_functions", basis_functions.shape)
        n_sub = len(nirom_options.compression.local_list)
        for sub in range(n_sub):
            nPOD_running_total = 0
            pod_coeffs_all = np.zeros((nTime*nirom_options.snapshots.nParameters,nPOD_total))
            for ifield in range(len(nirom_options.compression.field)):
                nPOD = nirom_options.compression.nPOD[ifield]
                nPOD_running_total = nPOD_running_total + nPOD
                # shape of pod_coeffs is nTime by number of POD coeffs
                #pod_coeffs = np.dot(np.transpose(snapshots), basis_functions)
                pod_coeffs = np.dot(dd_snapshots[sub], basis_functions[sub][ifield])
                pod_coeffs_all[:,nPOD_running_total-nPOD:nPOD_running_total] = pod_coeffs
            dd_pod_coeffs_all.append(pod_coeffs_all)
        dd_pod_coeffs_all = np.array(dd_pod_coeffs_all)
        print("dd_pod_coeffs_all", dd_pod_coeffs_all.shape)

    # 2. scale data
    #print ("min max of pod coeffs before transform",np.min(pod_coeffs_all), np.max(pod_coeffs_all))
    #scaling = MinMaxScaler(feature_range=scaling_bounds)
    #scaling = MinMaxScaler(feature_range=(0,10))
    #pod_coeffs_all = scaling.fit_transform(pod_coeffs_all)
    #print ("min max of pod coeffs after transform",np.min(pod_coeffs_all), np.max(pod_coeffs_all))

    # save the scaling settings for use in when predicting
    #scaler_filename = "scaler.sav"
    #joblib.dump(scaling, scaler_filename)

    # arrays for predicting the seen data (so not really prediction)
    pod_coeffs_all_prediction = np.zeros_like(dd_pod_coeffs_all)
    pod_coeffs_all_prediction[:,0,:] = dd_pod_coeffs_all[:,0,:]

    #for ifield in range(len(nirom_options.compression.field)):
    #    nPOD = nirom_options.compression.nPOD[ifield]
    #    nNodes = nirom_options.compression.nNodes[ifield]

    # 3. training (using all data at the moment)
    print ('training the GPR')
    t_train = 0
    t_replicate = 0
    nParameters = nirom_options.snapshots.nParameters
    trainX = np.zeros((n_sub, (nTime-1)*nParameters,nPOD_total ))
    #trainX = pod_coeffs_all[:nTime-1,:]
    for p in range(nParameters):
        trainX[:, p*(nTime-1):(p+1)*(nTime-1),:] = dd_pod_coeffs_all[:, p*nTime:(p+1)*nTime-1,:] # n-1 by nPOD (same every iteration so could take out of loop)

    for i in range(dd_pod_coeffs_all.shape[2]):

        t0 = time.time()

        #trainY = np.atleast_2d(pod_coeffs_all[1:,i]).T  # n-1 by 1 (different every iteration)
        trainY = np.zeros((n_sub, (nTime-1)*nParameters,1))
        for p in range(nParameters):
            #vector = np.atleast_2d(pod_coeffs_all[p*nTime+1:(p+1)*nTime,i]).T
            #print "shape output bits", vector.shape
            trainY[:, p*(nTime-1):(p+1)*(nTime-1),:] = np.atleast_3d(dd_pod_coeffs_all[:, p*nTime+1:(p+1)*nTime,i])

        if i==0:
            print('shape trainX and Y', trainX.shape, trainY.shape)
            #print "training POD coeff", i+1
        #else:
        #    print "POD coeff", i+1

        kernel = C(constant_value, constant_bounds) * RBF(RBF_length_scale, RBF_length_bounds)
        #kernel = C(1.0, (1e-3, 1e3)) * RBF(100, (1e-2, 1e2))
        gp_list = []
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=100)

        # Fit to data using Maximum Likelihood Estimation of the parameters
        for sub in range(n_sub):
            gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=100)
            #print("trainX[np.array(nirom_options.compression.local_list[sub])]", trainX[np.array(nirom_options.compression.local_list[sub])].shape)
            #print("trainY[sub]", trainY[sub].shape)
            train_X_sub = trainX[np.array(nirom_options.compression.local_list[sub])].transpose(1, 0, 2)
            gp.fit(train_X_sub.reshape(train_X_sub.shape[0],train_X_sub.shape[1]*train_X_sub.shape[2]), trainY[sub])
            gp_list.append(gp)

        t1 = time.time()
        t_train = t_train + t1 - t0

        # gives a convergence warning
        #print "score of gp" , gp.score(trainX,trainY)

        # save the model
        #joblib.dump(gp, 'gpr_model_' + str(i) + '.sav')

        if nParameters==1:
            # gives a warning - predicted variances < 0, setting to 0
            pod_pred_all_list = []
            t0 = time.time()
            for sub in range(n_sub):
                train_X_sub = trainX[np.array(nirom_options.compression.local_list[sub])].transpose(1, 0, 2)
                pod_pred_all, sigma = gp_list[sub].predict(train_X_sub.reshape(train_X_sub.shape[0],train_X_sub.shape[1]*train_X_sub.shape[2]), return_std=True)
                pod_pred_all_list.append(pod_pred_all)
            t1 = time.time()
            pod_pred_all_list = np.array(pod_pred_all_list)
            t_replicate = t_replicate + t1-t0

            # just replicating the snapshots not really predicting
            pod_coeffs_all_prediction[:, 1:,i] = pod_pred_all_list[:,:,0]

    if nParameters==1:
        # inverse of the scaling
        #pod_coeffs_all_prediction = scaling.inverse_transform(pod_coeffs_all_prediction)

        # map pod_coeffs_prediction to full space and write to file

        if nirom_options.compression.svd_autoencoder == True:

            decoder = load_model('decoder.h5')
            U_decoded = decoder.predict(pod_coeffs_all_prediction) ##### check shape of this
            print "U_decoded shape", U_decoded.shape
            pod_coeffs_all_prediction = U_decoded #U_decoded.T

        # modify write_vtk_from_vtk()
        filebase = fwd_options.results_filebase
        ext = fwd_options.results_extension
        files = os.listdir(fwd_options.path_to_results)
        for file in files:
            if file.startswith(filebase) and file.endswith(ext):
                filename = fwd_options.path_to_results + "/" + file
                break

        clean_vtk = get_clean_vtk_file(filename)
        nPOD_running_total = 0
        for ifield in range(len(nirom_options.compression.field)):
            nPOD = nirom_options.compression.nPOD[ifield]
            nPOD_running_total =  nPOD_running_total + nPOD
            nNodes = nirom_options.compression.nNodes[ifield]
            print("nNodes", nNodes)
            field = nirom_options.compression.field[ifield]

            solution_all = []
            for sub in range(n_sub):
                basis_functions = nirom_options.compression.basis_functions[sub][ifield]
                print("nirom_options.compression.basis_functions", nirom_options.compression.basis_functions.shape)
                print("basis_functions", basis_functions.shape)
                print("pod_coeffs_all_prediction", pod_coeffs_all_prediction.shape)
                pod_coeffs_pred = pod_coeffs_all_prediction[sub][:,nPOD_running_total-nPOD:nPOD_running_total]
                print("pod_coeffs_pred", pod_coeffs_pred.shape)
                solution_pred = np.dot(basis_functions,np.transpose(pod_coeffs_pred))
                solution_all.append(solution_pred)
                print "solution_pred", solution_pred.shape
            solution_all = np.array(solution_all)
            #solution_all = solution_all.reshape(solution_all.shape[0]*solution_all.shape[1], solution_all.shape[2])
            if solution_all.shape[0]*solution_all.shape[1] != nNodes:
                nComponents = int(solution_all.shape[0]*solution_all.shape[1] / nNodes)
                print "solution_pred.shape[0] / nNodes", solution_all.shape[0]*solution_all.shape[1]/ nNodes
                print "nComponents (integer hopefully)", nComponents
            else:
                nComponents = 1

            new_solution = []
            for component in range(nComponents):
                for find_old_place in nirom_options.compression.old2new:
                    domain_belong_to = find_old_place/nirom_options.compression.num_nodes_in_sub[0]
                    new_solution.append(solution_all[domain_belong_to][component*nirom_options.compression.num_nodes_in_sub[0] + find_old_place%nirom_options.compression.num_nodes_in_sub[0]])
            new_solution = np.array(new_solution)
            print("new_solution", new_solution.shape)

            for i in range(nTime-1):
                if ifield==0:
                    new_vtu = vtktools.vtu()
                    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
                    new_vtu.filename = "nirom_replication_" +str(i)+ ".vtu"
                else:
                    filename = "nirom_replication_" +str(i)+".vtu"
                    new_vtu = vtktools.vtu(filename)

                solution_pred_field = np.zeros((nNodes,nComponents))
                for j in range(nComponents):
                    solution_pred_field[:,j] = new_solution[j*nNodes:(j+1)*nNodes,i]


                new_vtu.AddField(field,solution_pred_field)
                new_vtu.Write()


    return t_train, t_replicate


# hard wired for svd-autoencoder and B&E test case
def predict_with_the_NIROM(nirom_options, fwd_options):

    nPOD = nirom_options.compression.nPOD[0]            # need to calculate these things!
    nPOD_total = nirom_options.compression.nPOD[0]
    nComponents = nirom_options.Snapshots.nComponents
    nNodes = 767559

    t0 = time.time()

    #pod_coeffs_initial_condition = np.zeros((1,nPOD_total))

    # get initial condition from which to start predicting
    # read in basis functions
    filename = "basis_functions.vtu"
    basis_functions_vtu = vtktools.vtu(filename)

    basis_functions = np.zeros((nComponents*nNodes,nPOD_total))
    for i in range(nPOD_total):
        field_name = 'Velocity_' + str(i)
        basis_functions_field = basis_functions_vtu.GetField(field_name)
        for j in range(nComponents):
            basis_functions[nNodes*j:nNodes*(j+1),i] = basis_functions_field[:,j]

# hard wired for svd-autoencoder and B&E test case
########################################################
    kTime = 419
    filebase = "LSBUv2"
    ext = "pvtu"
    # get final snapshot
    file_name = '../LSBU_results/' + filebase + '_' + str(kTime) + '.' + ext
    values,dummy = get_field_from_vtk(file_name, 'Velocity')
    final_snapshot = np.zeros((nNodes*nComponents))
    initial_condition_pod_coeffs = np.zeros((1,nPOD))
    for j in range(nComponents):
        final_snapshot[j*nNodes:(j+1)*nNodes] = values[:,j]

    result =  np.transpose(np.dot(np.transpose(basis_functions), final_snapshot))
    initial_condition_pod_coeffs[0,:] = result.T

    t1 = time.time()
    print "getting initial condition", t1-t0, "seconds"

######################################################################

    predicted_pod_coeffs = np.zeros((nirom_options.prediction.nTime,nPOD))
    # loop over the fields
    testX = np.zeros((1,nPOD))
    testX[0,0:nPOD] = initial_condition_pod_coeffs

    tload = 0.0
    tpred = 0.0
    # predict the next time step for each POD coeff
    print "starting prediction loop"
    for j in range(nirom_options.prediction.nTime):

        for i in range (nPOD_total):

            print "reading in previously trained model"
            # read in the previously trained model

            t0 = time.time()
            filename = 'gpr_model_' + str(i) + '.sav'
            gp_model = joblib.load(filename)
            t1 = time.time()
            tload = tload + t1 - t0

            t0 = time.time()
            testY, sigma = gp_model.predict(testX, return_std=True)
            t1 = time.time()
            tpred = tpred + (t1-t0)


            # prevent data from going out of bounds
            max_testY = 10.; min_testY = 0.0
            if testY > max_testY:
                testY = max_testY
            elif testY < min_testY:
                testY = min_testY


            predicted_pod_coeffs[j,i] = testY

        testX[0,0:nPOD] = predicted_pod_coeffs[j,:]


    # finished predicting for all POD coeffs

    t0 = time.time()
    # apply inverse scaling
    scaler_filename = "scaler.sav"
    scaling = joblib.load(scaler_filename)
    predicted_pod_coeffs = scaling.inverse_transform(predicted_pod_coeffs)

    # get solution and write to file
    predicted_solution = np.dot(basis_functions, np.transpose(predicted_pod_coeffs))
    clean_vtu = get_clean_vtk_file("basis_functions.vtu")

    new_vtu = vtktools.vtu()
    new_vtu.ugrid.DeepCopy(clean_vtu.ugrid)

    predicted_solution_field = np.zeros((nNodes, nComponents))
    for k in range(nirom_options.prediction.nTime):
        new_vtu.filename = "nirom_prediction_" +str(k)+ ".vtu"
        for j in range(nComponents):
            predicted_solution_field[:,j] = predicted_solution[nNodes*j:nNodes*(j+1),k]

        new_vtu.AddField('Velocity_NIROM',predicted_solution_field)
        new_vtu.Write()


    #result = gp_model.score(X_test, Y_test)
    #print(result)

    t1 = time.time()
    print "converting POD coeffs to solution and writing results", t1-t0, "seconds"

    return tload, tpred


# hard wired for GBNirom and advection case
def predict_with_the_NIROM_2D_advection_GBN(nirom_options, fwd_options):

    nPOD = 12           # need to calculate these things!
    nPOD_total = 12
    nComponents = 1
    nNodes = 600

    t0 = time.time()

    scaler_filename = "scaler.sav"
    scaling = joblib.load(scaler_filename)

    # 1. read in initial condition
    # pod_coeffs_last_time.dat
    # know nPOD_total

    #pod_coeffs_initial_condition = np.zeros((1,nPOD_total))

    # 1. read in basis functions
    filename = "basis_functions.vtu"
    basis_functions_vtu = vtktools.vtu(filename)

    basis_functions = np.zeros((nComponents*nNodes,nPOD_total))
    for i in range(nPOD_total):
        field_name = 'Temperature_' + str(i)
        basis_functions_field = basis_functions_vtu.GetField(field_name)
        for j in range(nComponents):
            basis_functions[nNodes*j:nNodes*(j+1),i] = basis_functions_field[:,j]

# hard wired for GBNirom and advection case
########################################################

    # 2. read in initial condition
    nTest = 100
    for iTest in range(nTest):
        file_name = '../unseen_tests/Window_0_Ensemble_Temperature_' + str(iTest) + "/Advection_0.vtu"
        print "file_name of initial conditions =========", file_name
        values,dummy = get_field_from_vtk(file_name, 'Temperature')
        #print "values", values.T
        initial_condition_from_vtu = np.zeros((nNodes*nComponents))
        initial_condition_pod_coeffs = np.zeros((1,nPOD))
        if nComponents==1:
            initial_condition_from_vtu[0:nNodes] = values
        else:
            for j in range(nComponents):
                initial_condition_from_vtu[j*nNodes:(j+1)*nNodes] = values[:,j]

        result =  np.transpose(np.dot(np.transpose(basis_functions), initial_condition_from_vtu))
        initial_condition_pod_coeffs[0,:] = result.T
        print "initial_condition_pod_coeffs", initial_condition_pod_coeffs
        #print "shape at least 2D", np.atleast_2d(initial_condition_pod_coeffs[0,:]).shape
        print "above scaled", scaling.transform(np.atleast_2d(initial_condition_pod_coeffs[0,:]))
        initial_condition_pod_coeffs[0,:] = scaling.transform(np.atleast_2d(initial_condition_pod_coeffs[0,:]))

######################################################################

        predicted_pod_coeffs = np.zeros((nirom_options.prediction.nTime,nPOD))
        # loop over the fields
        testX = np.zeros((1,nPOD))

        print "testX.shape", testX.shape
        print "initial_condition_pod_coeffs.shape", initial_condition_pod_coeffs.shape


        for j in range(nirom_options.prediction.nTime):

            if j==0:
                testX[0,0:nPOD] = initial_condition_pod_coeffs
            else:
                testX[0,0:nPOD] = predicted_pod_coeffs[j-1,:]

            for i in range (nPOD_total):     # loop over POD coefficients

                # 3. read in the previously trained model
                filename = 'gpr_model_' + str(i) + '.sav'
                gp_model = joblib.load(filename)

                testY, sigma = gp_model.predict(testX, return_std=True)

                max_testY = 10.; min_testY = 0.0
                if testY > max_testY:
                    testY = max_testY
                elif testY < min_testY:
                    testY = min_testY

                predicted_pod_coeffs[j,i] = testY

            testX[0,0:nPOD] = predicted_pod_coeffs[j,:]


        # finished predicting for all POD coeffs

        # 5. load scaling and apply inverse
        #scaler_filename = "scaler.sav"
        #scaling = joblib.load(scaler_filename)
        predicted_pod_coeffs = scaling.inverse_transform(predicted_pod_coeffs)

        # get solution and write to file
        predicted_solution = np.dot(basis_functions, np.transpose(predicted_pod_coeffs))
        #new_vtu = vtktools.vtu()
        #new_vtu.ugrid.DeepCopy(basis_functions_vtu.ugrid)
        clean_vtu = get_clean_vtk_file("basis_functions.vtu")

        print "copying ugrid from basis functions"
        new_vtu = vtktools.vtu()
        new_vtu.ugrid.DeepCopy(clean_vtu.ugrid)

        print "making test directory"
        string = "mkdir Test_" + str(iTest)
        os.system(string)

        predicted_solution_field = np.zeros((nNodes, nComponents))
        for k in range(nirom_options.prediction.nTime):
            new_vtu.filename = "Test_" + str(iTest) + "/" + "nirom_prediction_" +str(k)+ ".vtu"
            for j in range(nComponents):
                predicted_solution_field[:,j] = predicted_solution[nNodes*j:nNodes*(j+1),k]

            new_vtu.AddField('Temperature_NIROM',predicted_solution_field)
            new_vtu.Write()


    #result = gp_model.score(X_test, Y_test)
    #print(result)

    t1 = time.time()

    return t0, t1
