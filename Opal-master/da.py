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

#Authors: R. Arcucci, C. Heaney, C.C. Pain and P. Salinas

import sys
import os
import numpy as np
from importance_map_calc  import *
from importance_map_tools import *
import time
import libspud
from scipy.optimize import minimize


def run_forward_model_tk(k,executable,input_file):
    ''' run the forward model for one time step '''
    cwd = os.getcwd()
    os.chdir(str(k))
    
    string = executable + " " + input_file
    print "running ", string
    os.system(string)

    os.chdir(cwd)

    return

def modify_initial_conditions_at_tk(k,Vv): 
    ''' modify the initial conditions in the checkpoint vtu file'''
    cwd = os.getcwd()
    os.chdir(str(k))

    field = 'Tracer'
    vtu_file = '2d_canyon_PressureMesh_8_checkpoint.vtu'
    vtu_data = vtktools.vtu(vtu_file)
    tracer_values = vtu_data.GetScalarField(field)
    vtu_data.AddField('Tracer_Old',tracer_values)
    vtu_data.RemoveField('Tracer')
    # write Vv in the checkpoint to go one step forward, Vv at t_k
    vtu_data.AddField('Tracer',Vv)
    vtu_data.Write() 

    os.chdir(cwd)

    return field


def prepare_inputs_for_forward_model(k):
    ''' make directory to run the fowrard model in, copy across the input files and
        modify the time step in flml so the forward model will run just for one time step'''             
    # make a directory to run the code in and copy across input files
    cwd = os.getcwd()
    input_file_name = "fwd_model.flml" 
#             print "files in cwd"
#             for files in os.listdir(cwd):
#                print files
    if not os.path.isdir( str(k) ):
        print "attempting to mkdir" 
        os.mkdir(str(k))
        os.system('cp *.msh ' + str(k))
        os.system('cp *_' +str(k)+ '_checkpoint* ' + str(k))

        os.chdir(str(k))

        # modify the checkpoint file times
        for files in os.listdir('./'):
            # get the name of the checkpoint flml     
            if files.endswith(str(k)+"_checkpoint.flml"):
#               pos = files.rfind('.')
                checkpoint_file_name = files
#               print "checkpoint fname", checkpoint_file_name

                # load options from checkpoint file
                libspud.load_options(checkpoint_file_name)

                # change the name of the output 
                libspud.set_option('/simulation_name','2d_canyon')

                # change the final time so it runs from t_k to t_k+1 only
                t0 = libspud.get_option('/timestepping/current_time')
                dt = libspud.get_option('/timestepping/timestep')
                libspud.set_option('/timestepping/finish_time',t0+dt)
                # could check with vtu's that these are the correct times

                # rename input file
                libspud.write_options(input_file_name)
                libspud.clear_options()

    os.chdir(cwd)

    return input_file_name


def data_assimilation(opal_options):

   global Model_updated, iter_count

   # functions used within data_assimilation: J() and gradJ()
   def J(v):
	global Model_updated, iter_count

        iter_count = iter_count + 1

        vT = np.transpose(v)
        vTv = np.dot(vT,v)
        Vv = np.dot(V,v)
        HVv = np.dot(H,Vv)

        # we need the model results - check if these are already available, if not, run the forward model with Vv as the input
        if Model_updated:
             print "in J(): using pre-exiting forward model model solution"
             Model_updated = False
        else:
             print "in J(): updating the forward model solution"

             # prepare directory and input files for forward model
             input_file_name = prepare_inputs_for_forward_model(k)

             # modify initial condition of tracer
             field = modify_initial_conditions_at_tk(k,Vv)

             # run forward model
             run_forward_model_tk(k,opal_options.executable,input_file_name)

             Model_updated = True

        # retrieve the forward model results, MVv
        path_to_vtu_file = str(k) + '/2d_canyon_1.vtu'
        vtu_data = vtktools.vtu(path_to_vtu_file)
        MVv = vtu_data.GetScalarField(field)
#        print "J(): MVv[0:10]", MVv[0:10]

        ##MVv = np.dot(M,Vv)
        HMVv = np.dot(H,MVv)
        Jmis = np.subtract(HVv,d)
        JmisM = np.subtract(HMVv,d)
        invR = np.linalg.inv(R)
        JmisT = np.transpose(Jmis)
        JmisMT = np.transpose(JmisM)
        RJmis = np.dot(invR,JmisT)
        RJmisM = np.dot(invR,JmisMT)
        J1 = np.dot(Jmis,RJmis)
        JM1 = np.dot(JmisM,RJmisM)

        Jv = (vTv + J1 + JM1) / 2

        return Jv    
    
###############################################
#######         GRADIENT OF J          ########
###############################################


   def gradJ(v):
        global Model_updated

        Vv = np.dot(V,v)
        HVv = np.dot(H,Vv)

        # CODE COPIED FROM J() ###########################################
        # we need the model results - check if these are already available, 
        # if not, run the forward model with Vv as the input
        if Model_updated:
             print "in gradJ(): using pre-exiting forward model model solution"
             Model_updated = False
        else:
             print "in gradJ(): updating the forward model solution"

             # prepare directory and input files for forward model
             input_file_name = prepare_inputs_for_forward_model(k)

             # modify initial condition of tracer
             field = modify_initial_conditions_at_tk(k,Vv)

             # run forward model
             run_forward_model_tk(k,opal_options.executable,input_file_name)

             Model_updated = True

        # END OF CODE COPIED FROM J() ###########################################

        # MVv = np.dot(M,Vv)
        # retrieve the forward model results, MVv
        path_to_vtu_file = str(k) + '/2d_canyon_1.vtu'
        vtu_data = vtktools.vtu(path_to_vtu_file)
        MVv = vtu_data.GetScalarField('Tracer') #vtu_data.GetScalarField(field)
#        print "gradJ: MVv[0:10]", MVv[0:10]

        HMVv = np.dot(H,MVv)
        Jmis = np.subtract(HVv,d)
        JmisM = np.subtract(HMVv,d)
        invR = np.linalg.inv(R)
        RJmis = np.dot(invR,Jmis)
        RJmisM = np.dot(invR,JmisM)
        HT = np.transpose(H)
        g1 = np.dot(HT,RJmis)
        g1M = np.dot(HT,RJmisM)
        ##MT = ... MT(g1M) = from importance map t_k+1 , map at t_k
        VT = np.transpose(V)
        ##VTMT = np.dot(VT,MT)
        g2 = np.dot(VT,g1)


        ggJ = v + g2 ##+ VTMT

        return ggJ

    #print "exectuable which has been selected:", opal_options.executable
    
    # ......... read the input .........
    
###############################################
## TRUNCATION AND REGULARIZATION PARAMETERS ###
###############################################
   # inputs
   n = 852
   lam = 1  #REGULARIZATION PARAMETER
   m = 45  #TRUNCATION PARAMETER FROM buildV.py
   xB = np.ones(n)
   y = np.ones(n)

   k = 8 # time at which observation is known

###############################################
########  INTIAL RUN OF FLUIDITY #############
###############################################
   # put checkpointing on for file k
   print "name of fwd_input_file", opal_options.data_assimilation.fwd_input_file
   libspud.load_options('2d_canyon.flml')#(opal_options.data_assimilation.fwd_input_file)

   # don't need these currently
   if libspud.have_option('io/checkpointing/checkpoint_at_start'):
       libspud.delete_option('io/checkpointing/checkpoint_at_start')
   if libspud.have_option('io/checkpointing/checkpoint_at_end'):
       libspud.delete_option('io/checkpointing/checkpoint_at_end')
   if libspud.have_option('io/checkpointing/checkpoint_period_in_dumps'):
       libspud.set_option('io/checkpointing/checkpoint_period_in_dumps',k)
   else:
      print "checkpoint_period_in_dumps option missing from xml file"
      sys.exit(0)

   libspud.write_options(opal_options.data_assimilation.fwd_input_file)
   libspud.clear_options() 

   string = opal_options.executable + " " + opal_options.data_assimilation.fwd_input_file


   # run code which will checkpoint every "k" dumps at the moment....
   print string
   os.system(string)


###############################################
########  COVARIANCE MATRICES     #############
###############################################


   V = np.loadtxt('matrixVprec'+str(m)+'.txt', usecols=range(m))

   R = lam * 0.5 * np.identity(n)

   H = np.identity(n)


###############################################
####### FROM PHYSICAL TO CONTROL SPACE ########
###############################################

   x0 = np.ones(n)
   Vin = np.linalg.pinv(V)
   v0 = np.dot(Vin,x0)


###############################################
#######       COMPUTE THE MISFIT       ########
###############################################


   VT = np.transpose(V)
   HxB = np.dot(H,xB)
   # consider multiple observations later - just one for now
   d = np.subtract(y,HxB)


###############################################
#######    COMPUTE THE MINIMUM OF J    ########
###############################################

   t = time.time()
   iter_count = 0
   Model_updated = False

   res = minimize(J, v0, method='L-BFGS-B', jac=gradJ,
                options={'disp': True})

###############################################
####### FROM CONTROL TO PHYSICAL SPACE ########
###############################################


   vDA = np.array([])
   vDA = res.x
   deltaxDA = np.dot(V,vDA)
   xDA = xB + deltaxDA

   elapsed = time.time() - t

   print " iter_count", iter_count

   #return 

###############################################
####### PRECONDITIONED COST FUNCTION J ########
###############################################



   return


'''
from scipy.optimize import minimize

t = time.time()


res = minimize(J, v0, method='L-BFGS-B', jac=gradJ,
                options={'disp': True})

###############################################
####### FROM CONTROL TO PHYSICAL SPACE ########
###############################################


vDA = np.array([])
vDA = res.x
deltaxDA = np.dot(V,vDA)
xDA = xB + deltaxDA

elapsed = time.time() - t




#   ---- OUTPUT  ---- 


    return ......

   return 
'''



