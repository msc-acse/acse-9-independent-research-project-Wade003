import decomp_new
import numpy as np 


#print(decomp.python_set_up_recbis.__doc__)
nx=500
ny=2000
n = nx*ny
nsplt = 7
#ncola, cola, fina = decomp_new.one_d_row_stor(n,1,n,(n*3)-2)
ncola, cola, fina = decomp_new.one_d_row_stor(nx,ny,n,(n*5))
split_levels = np.zeros((nsplt),dtype='int32')
split_levels[:] = 2

havwnod = 2

wnod = np.zeros(n)
wnod[:] = 1.0

havmat = 0
a = np.zeros(1)
exact = True 
iexact = 1
ii=1
na=0

print "OK up to here"

#witchd = decomp_new.python_set_up_recbis(splev,fina,cola, nsplt,ncola,n, havwnod,wnod,exact, havmat)

print "split_levels", split_levels
print "shape fina", fina.shape
print "shape cola", cola.shape
print "nsplit", nsplt
print "ncola", ncola
print "ii", ii
print decomp_new.python_set_up_recbis.__doc__
# we have to finish the sub clal with the variables used in the decleration within this sub.
witchd=decomp_new.python_set_up_recbis(split_levels,fina,cola, wnod,a, havwnod,havmat,iexact, nsplt,ncola,n,na)
#print(witchd)

for j in range(ny): #do j=1,ny
    print witchd[j*nx : (j+1)*nx]

