import numpy as np
import sys
import f2py

nsplt = 1
nx = 100
ny = 1
np.set_printoptions(threshold = sys.maxsize)
whichd, ncola, cola, fina = f2py.DD(nx, ny, nsplt)    # DD's crutial data
print("whichd", whichd)
list_subs_all = f2py.get_subdomain_list(whichd, fina, cola, nsplt, nx*ny)
print("list_subs_all", list_subs_all)
order = f2py.reorde(list_subs_all, int(pow(2, nsplt)))
print("order", order)
