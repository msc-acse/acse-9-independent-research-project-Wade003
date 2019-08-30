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

#Authors: Claire Heaney, Pablo Salinas, Dunhui Xiao, Christopher Pain, Wade Song

import numpy as np
import sys
from nirom_tools import get_nirom_settings, read_in_snapshots, compress_snapshots, train_the_NIROM_GPR, train_the_NIROM_DD_GPR, predict_with_the_NIROM, predict_with_the_NIROM_2D_advection_GBN

from nirom_LSTM_tools import train_the_NIROM_LSTM, train_the_NIROM_DD_LSTM
from opal_classes import Timings, Compression


def nirom(fwd_options,nirom_options):

# temporarily checking options
#    print "svd            ", nirom_options.compression.svd
#    print "eigh           ", nirom_options.compression.eigh
#    print "autoencoder    ", nirom_options.compression.autoencoder
#    print "svd_autoencoder", nirom_options.compression.svd_autoencoder
#    print "fields         ", nirom_options.compression.field
#    print "epochs         ", nirom_options.compression.epochs
#    print "batch_size     ", nirom_options.compression.batch_size
#    print "nPod           ", nirom_options.compression.nPOD
#    print "cum_tol        ", nirom_options.compression.cumulative_tol
#    print "neurons        ", nirom_options.compression.neurons_each_layer

#    print "LSTM                ", nirom_options.training.LSTM
#    print "GPR                 ", nirom_options.training.GPR
#    print "GPR_scaling         ", nirom_options.training.GPR_scaling
#    print "RBF_length_scale    ", nirom_options.training.GPR_RBF_length_scale
#    print "RBF_length_bounds   ", nirom_options.training.GPR_RBF_length_bounds
#    print "RBF_constant_value  ", nirom_options.training.GPR_constant_value
#    print "RBF_constant_bounds ", nirom_options.training.GPR_constant_bounds

#    sys.exit()

    # for reading in snapshots / compressing and training
    if nirom_options.prediction.nTime == 0:

        timings = Timings()

        # generate snapshots if needs be ----------------------------------------------------
        if nirom_options.snapshots.create:
            print 'Sorry - you can only point at pre-generated snapshots at the moment'
            print 'Run the simulations and then point to where the snapshots are located'
            print 'exiting'
            sys.exit(0)

        # get some settings in order to read in the snapshots
        fwd_options, nirom_options  = get_nirom_settings(fwd_options, nirom_options)

        # (1) read in snapshots -----------------------------------------------------------------
        nirom_options, t = read_in_snapshots(fwd_options, nirom_options)

        timings.read_in_snapshots = t
        print "time for snapshots",t

        # (2) compression ----------------------------------------------------------------------
        nirom_options, t,  U_encoded = compress_snapshots(fwd_options, nirom_options)

        nirom_options.compression.write_sing_values()

        timings.compression = t

        # else read in basis functions

        # map snapshots to reduced space ready for the training -----------------------------
        # do this here and write to file

        # (3) training of the GPR or LSTM ---------------------------------------------------------------
        # (4) predicting
        if nirom_options.training.GPR :
            t_train, t_pred = train_the_NIROM_GPR(nirom_options, fwd_options, U_encoded)
            timings.training = t_train
            timings.replication = t_pred
        elif nirom_options.training.DD_GPR:
            t_train, t_pred = train_the_NIROM_DD_GPR(nirom_options, fwd_options, U_encoded)
            timings.training = t_train
            timings.replication = t_pred
        elif nirom_options.training.LSTM:
            t_train, t_pred = train_the_NIROM_LSTM(fwd_options, nirom_options)
            timings.training = t_train
            timings.replication = t_pred
        elif nirom_options.training.DD_LSTM:
            t_train, t_pred = train_the_NIROM_DD_LSTM(fwd_options, nirom_options)
            timings.training = t_train
            timings.replication = t_pred
        print "time for training", t_train
        print "time for prediction", t_pred
        if nirom_options.compression.svd:
            svd_eigh = 'svd'
        else:
            svd_eigh = 'eigh'
        timings.print_values(svd_eigh)
        timings.write_values('nirom_timings.dat',svd_eigh)

    # (4) for predicting forward in time from the final snapshot
    elif nirom_options.prediction.nTime > 0:

        # predicting unseen behaviour with the GPR ---------------------------------------------------------------

        #if nirom_options.compression.svd_autoencoder:
        #    # autoencoder / prediction for BE case
        #    t1, t2 = predict_with_the_NIROM(nirom_options, fwd_options)
        #else:
        #    # GB NIROM advection case
        #    t1, t2 = predict_with_the_NIROM_2D_advection_GBN(nirom_options, fwd_options)
        tload, tpred = predict_with_the_NIROM(nirom_options, fwd_options)

        print "time to load neural networks (.sav files)", tload
        print "time to predict", tpred
        timings.prediction = tpred


    return
