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

#Author: Wade Song

import numpy as np
import sys
import os
import vtk, vtktools
from importance_map_tools import *
import time
import libspud
import fnmatch
from sklearn.preprocessing import RobustScaler
from sklearn.externals import joblib

import LSTM_Global
import LSTM_DD

def train_the_NIROM_LSTM(fwd_options, nirom_options):
    print("fwd_options", fwd_options)
    print("nirom_options", nirom_options)
    print ("hello from Wade's LSTM")
    print("nirom_options.compression.svd_autoencoder", nirom_options.compression.svd_autoencoder)
    if nirom_options.compression.svd_autoencoder == True:
        #print("U_encoded",U_encoded)
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
            print("basis_functions", basis_functions.shape)
            print("pod_coeffs", pod_coeffs.shape)
            pod_coeffs_all[:,nPOD_running_total-nPOD:nPOD_running_total] = pod_coeffs
        #pod_coeffs_all_orig = copy.deepcopy(pod_coeffs_all)


    # 2. training (using 100 data at the moment)
    print ('training the LSTM')
    nParameters = nirom_options.snapshots.nParameters
    history_level = 3
    t_train = 0
    t_pred = 0
    if fwd_options.ndim == 3:
        history_level = 4
    t0 = time.time()
    lstmNN = LSTM_Global.LSTM_train(pod_coeffs_all, history_level)
    t1 = time.time()
    t_train = t_train + t1 - t0
    # 3. predicting (predict 1000 time steps)
    predict_times = 420
    t0 = time.time()
    pod_coeffs_all_prediction = np.zeros((predict_times, nPOD_total))
    pod_coeffs_all_prediction[:history_level] = pod_coeffs_all[:history_level]
    pod_coeffs_all_prediction[history_level:] = np.array(LSTM_Global.LSTM_predict(pod_coeffs_all, history_level, lstmNN, predict_times))
    t1 = time.time()
    t_pred = t_train + t1 - t0
    if nParameters==1:

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

            for i in range(predict_times - history_level):
                if ifield==0:
                    new_vtu = vtktools.vtu()
                    new_vtu.ugrid.DeepCopy(clean_vtk.ugrid)
                    new_vtu.filename = "nirom_replication_" +str(i+history_level)+ ".vtu"
                else:
                    filename = "nirom_replication_" +str(i+history_level)+".vtu"
                    new_vtu = vtktools.vtu(filename)

                solution_pred_field = np.zeros((nNodes,nComponents))
                for j in range(nComponents):
                    solution_pred_field[:,j] = solution_pred[j*nNodes:(j+1)*nNodes,i]


                new_vtu.AddField(field,solution_pred_field)
                new_vtu.Write()

    return t_train, t_pred


def train_the_NIROM_DD_LSTM(fwd_options, nirom_options):
    print ("hello from Wade's DD_LSTM")
    print("nirom_options.compression.svd_autoencoder", nirom_options.compression.svd_autoencoder)
    n_sub = len(nirom_options.compression.local_list)
    if nirom_options.compression.svd_autoencoder == True:
        pod_coeffs_all = U_encoded #U_encoded.T
        nTime = fwd_options.nTime

    else:
        t = 0
        nPOD_total = 0

        for ifield in range(len(nirom_options.compression.field)):
            nPOD_total = nPOD_total + nirom_options.compression.nPOD[ifield]

        # 1. map snapshots to reduced space to give pod coefficients
        nTime = fwd_options.nTime
        dd_pod_coeffs_all = []
        dd_snapshots = nirom_options.compression.dd_snapshot
        basis_functions = nirom_options.compression.basis_functions
        print("dd_snapshots", dd_snapshots.shape)
        print("basis_functions", basis_functions.shape)
        for sub in range(n_sub):
            nPOD_running_total = 0
            pod_coeffs_all = np.zeros((nTime*nirom_options.snapshots.nParameters,nPOD_total))
            for ifield in range(len(nirom_options.compression.field)):
                nPOD = nirom_options.compression.nPOD[ifield]
                nPOD_running_total = nPOD_running_total + nPOD
                # shape of pod_coeffs is nTime by number of POD coeffs
                #pod_coeffs = np.dot(np.transpose(snapshots), basis_functions)
                pod_coeffs = np.dot(dd_snapshots[sub], basis_functions[sub][ifield])
                print("pod_coeffs", pod_coeffs.shape)
                print("nPOD_running_total-nPOD", nPOD_running_total-nPOD)
                print("pod_coeffs_all", pod_coeffs_all.shape, pod_coeffs_all[:,nPOD_running_total-nPOD:nPOD_running_total].shape)
                pod_coeffs_all[:,nPOD_running_total-nPOD:nPOD_running_total] = pod_coeffs
            dd_pod_coeffs_all.append(pod_coeffs_all)
        dd_pod_coeffs_all = np.array(dd_pod_coeffs_all)
        print("dd_pod_coeffs_all", dd_pod_coeffs_all.shape)
        #pod_coeffs_all_orig = copy.deepcopy(pod_coeffs_all)


    # 3. training (using 100 data at the moment)
    print ('training the DD-LSTM')
    nParameters = nirom_options.snapshots.nParameters
    history_level = 3
    t_train = 0
    t_pred = 0
    t0 = time.time()
    lstmNN = LSTM_DD.LSTM_DD_global_train(dd_pod_coeffs_all, history_level, nirom_options.compression.local_list, nirom_options.compression.order)
    t1 = time.time()
    t_train = t_train + t1 - t0
    n_sub = dd_pod_coeffs_all.shape[0]
    # 3. predicting (predict 1000 time steps)
    predict_times = 420

    t0 = time.time()
    pod_coeffs_all_prediction = np.array(LSTM_DD.LSTM_DD_global_predict(dd_pod_coeffs_all, history_level, lstmNN, predict_times, nirom_options.compression.local_list, nirom_options.compression.order )).reshape(n_sub, predict_times-history_level, dd_pod_coeffs_all.shape[2])
    t1 = time.time()
    t_pred = t_pred + t1 - t0
    print("pod_coeffs_all_prediction", pod_coeffs_all_prediction.shape)
    if nParameters==1:

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

            for i in range(predict_times - history_level):
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

    return t_train, t_pred

def train_the_NIROM_DD_Global_LSTM(fwd_options, nirom_options):
    print ("hello from Wade's DD_Global_LSTM")
    print("nirom_options.compression.svd_autoencoder", nirom_options.compression.svd_autoencoder)
    n_sub = len(nirom_options.compression.local_list)
    if nirom_options.compression.svd_autoencoder == True:
        pod_coeffs_all = U_encoded #U_encoded.T
        nTime = fwd_options.nTime

    else:
        t = 0
        nPOD_total = 0

        for ifield in range(len(nirom_options.compression.field)):
            nPOD_total = nPOD_total + nirom_options.compression.nPOD[ifield]

        # 1. map snapshots to reduced space to give pod coefficients
        nTime = fwd_options.nTime
        dd_pod_coeffs_all = []
        dd_snapshots = nirom_options.compression.dd_snapshot
        basis_functions = nirom_options.compression.basis_functions
        print("dd_snapshots", dd_snapshots.shape)
        print("basis_functions", basis_functions.shape)
        for sub in range(n_sub):
            nPOD_running_total = 0
            pod_coeffs_all = np.zeros((nTime*nirom_options.snapshots.nParameters,nPOD_total))
            for ifield in range(len(nirom_options.compression.field)):
                nPOD = nirom_options.compression.nPOD[ifield]
                nPOD_running_total = nPOD_running_total + nPOD
                # shape of pod_coeffs is nTime by number of POD coeffs
                #pod_coeffs = np.dot(np.transpose(snapshots), basis_functions)
                pod_coeffs = np.dot(dd_snapshots[sub], basis_functions[sub][ifield])
                print("pod_coeffs", pod_coeffs.shape)
                print("nPOD_running_total-nPOD", nPOD_running_total-nPOD)
                print("pod_coeffs_all", pod_coeffs_all.shape, pod_coeffs_all[:,nPOD_running_total-nPOD:nPOD_running_total].shape)
                pod_coeffs_all[:,nPOD_running_total-nPOD:nPOD_running_total] = pod_coeffs
            dd_pod_coeffs_all.append(pod_coeffs_all)
        dd_pod_coeffs_all = np.array(dd_pod_coeffs_all)
        print("dd_pod_coeffs_all", dd_pod_coeffs_all.shape)
        #pod_coeffs_all_orig = copy.deepcopy(pod_coeffs_all)


    # 3. training (using 100 data at the moment)
    print ('training the DD-LSTM')
    nParameters = nirom_options.snapshots.nParameters
    t_train = 0
    t_pred = 0
    history_level = 3

    lstmNN = LSTM_DD.LSTM_DD_global_train(dd_pod_coeffs_all, history_level, nirom_options.compression.local_list, nirom_options.compression.order)


    # 3. predicting (predict 1000 time steps)
    predict_times = 200
    t0 = time.time()
    #pod_coeffs_all_prediction = np.array(LSTM_DD.LSTM_predict(dd_pod_coeffs_all, history_level, lstmNN, predict_times, nirom_options.compression.local_list, nirom_options.compression.order ))
    pod_coeffs_all_prediction = np.array(LSTM_DD.LSTM_DD_global_predict(dd_pod_coeffs_all, history_level, lstmNN, predict_times, nirom_options.compression.local_list, nirom_options.compression.order ))
    t1 = time.time()
    t_pred = t_pred + t1 - t0
    print("pod_coeffs_all_prediction", pod_coeffs_all_prediction.shape)
    if nParameters==1:

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

            for i in range(predict_times - history_level):
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

    return t_train, t_pred
