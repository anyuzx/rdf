import cython
import numpy as np
cimport numpy as np

cdef extern from "c_rdf.h":
	void c_rdfone(double* r_array, double* rdf_array, double* dim_array, int n, int natmm, int nHist, double rHistMax, double v)
	void c_rdftwo(double* r1_array, double* r2_array, double* rdf_array, double* dim_array, int n, int natmm1, int natmm2, int nHist, double rHistMax, double v)

@cython.boundscheck(False)
@cython.wraparound(False)

def rdfone(np.ndarray[double,ndim = 1,mode="c"] r_array not None, np.ndarray[double,ndim = 1,mode="c"] dim_array not None, int natmm, int nHist, double rHistMax, double v):
	cdef int n
	cdef np.ndarray[double,ndim = 1,mode="c"] rdf_array = np.ndarray((nHist,),dtype=np.float64)
	n = r_array.shape[0]/(3*natmm)

	c_rdfone(&r_array[0], &rdf_array[0], &dim_array[0], n, natmm, nHist, rHistMax, v)

	return rdf_array

def rdftwo(np.ndarray[double,ndim = 1,mode="c"] r1_array not None, np.ndarray[double,ndim = 1,mode="c"] r2_array not None, np.ndarray[double,ndim = 1,mode="c"] dim_array not None, int natmm1, int natmm2, int nHist, double rHistMax, double v):
	cdef int n
	cdef np.ndarray[double,ndim = 1,mode="c"] rdf_array = np.ndarray((nHist,),dtype=np.float64)
	n = r1_array.shape[0]/(3*natmm1)

	c_rdftwo(&r1_array[0], &r2_array[0], &rdf_array[0], &dim_array[0], n, natmm1, natmm2, nHist, rHistMax, v)

	return rdf_array