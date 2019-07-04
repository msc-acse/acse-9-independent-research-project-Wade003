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

#Authors: R.Arcucci

def variational_data_assimilation(opal_options):
  
    #print "exectuable which has been selected:", opal_options.executable
    
    # ......... read the input .........
    
###############################################
## TRUNCATION AND REGULARIZATION PARAMETERS ###
###############################################

lam=1  #REGULARIZATION PARAMETER

m=145  #TRUNCATION PARAMETER FROM buildV.py

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
d = np.subtract(y,HxB)


###############################################
####### PRECONDITIONED COST FUNCTION J ########
###############################################


def J(v):

        vT = np.transpose(v)
        vTv = np.dot(vT,v)
        Vv = np.dot(V,v)
        HVv = np.dot(H,Vv)
        Jmis = np.subtract(HVv,d)
        invR = inv(R)
        JmisT = np.transpose(Jmis)
        RJmis = np.dot(invR,JmisT)
        J1 = np.dot(Jmis,RJmis)

        Jv = (vTv + J1) / 2

        return Jv

    
    
###############################################
#######         GRADIENT OF J          ########
###############################################


def gradJ(v):

        Vv = np.dot(V,v)
        HVv = np.dot(H,Vv)
        Jmis = np.subtract(HVv,d)
        invR = inv(R)
        RJmis = np.dot(invR,Jmis)
        HT = np.transpose(H)
        g1 = np.dot(HT,RJmis)
        VT = np.transpose(V)
        g2 = np.dot(VT,g1)

        ggJ = v + g2

        return ggJ



###############################################
#######    COMPUTE THE MINIMUM OF J    ########
###############################################


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





