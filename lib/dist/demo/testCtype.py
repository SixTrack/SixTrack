from ctypes import *
double6 = c_double * 6
DOUBLE = c_double
PDOUBLE = POINTER(DOUBLE)
ptr = (DOUBLE*6*6)()
DBL5ARR = DOUBLE * 6
# An array of double* can be passed to your function as double**.
PDBL4ARR = PDOUBLE * 6

ptr = PDBL4ARR()
for i in range(6):
    # fill out each pointer with an array of doubles.
    ptr[i] = DBL5ARR()
    for j in range(6):
        ptr[i][j] = i + j  # just to initialize the actual doubles.
        print(i)
dim = c_int(6)

coord_c = double6(1,1,1,1,1,0)
dist = cdll.LoadLibrary("./libhello.so")
#dist.hello()
dist.printmatrix(dim, dim, ptr)
#dim_c = c_int(6)
#dist.addclosedorbit_(coord_c)
#mychar = c_char_p(b"wrrritiing")
#dist.printvector(mychar, dim_c, coord_c)
#for d in coord_c:
#	print(d)
