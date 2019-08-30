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

#Authors: C. Heaney, C.C. Pain, P. Salinas, and others.....
from opal_classes import get_opal_options
from da import *
from importance_map import generate_importance_map


def call_optimisation():
    import ga_well_location_optimisation
    ga_well_location_optimisation.main()

def call_nirom():
    from nirom import nirom
    nirom(fwd_options, nirom_options)

def call_second_order_da():
    import  second_order_da
    second_order_da(opal_options)



# read in the opal options
opal_options, fwd_options,nirom_options = get_opal_options()

if opal_options.opal_operation == 'importance_map':
    generate_importance_map(opal_options)   
elif opal_options.opal_operation == 'da':
    data_assimilation(opal_options) 
elif opal_options.opal_operation == 'nirom':
    call_nirom()
elif opal_options.opal_operation == 'second_order_da':
    call_second_order_da()
elif opal_options.opal_operation == 'ga_optimisation':
    call_optimisation()


print "\nOPAL ending \n"

