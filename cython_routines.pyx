import cython
cimport cython

import numpy as np
cimport numpy as np

dtypeInt = np.int
dtypeFloat = np.float

ctypedef np.int_t dtypeInt_t
ctypedef np.float_t dtypeFloat_t

@cython.boundscheck(False)
@cython.wraparound(False)

def compute_sim(np.ndarray[dtypeFloat_t,ndim=2] infmat, np.ndarray[dtypeFloat_t,ndim=1] h, np.ndarray[dtypeInt_t,ndim=1] indices, int p, int q):

    cdef np.ndarray[dtypeFloat_t,ndim=2] M = np.zeros((q,q), dtype=dtypeFloat)
    
    cdef unsigned int i, j
    
    for j in xrange(q):
        for i in xrange(q):
            M[i,j] = infmat[indices[i],indices[j]]*h[j]
            
    return M
    
def compute_sim_classic(np.ndarray[dtypeFloat_t,ndim=2] infmat, np.ndarray[dtypeFloat_t,ndim=1] h, np.ndarray[dtypeInt_t,ndim=1] indices, int p, int q):

    cdef np.ndarray[dtypeFloat_t,ndim=2] M = np.zeros((q,q), dtype=dtypeFloat)
    
    cdef unsigned int i, j
    
    for j in xrange(q):
        for i in xrange(j,q):
            M[i,j] = min(infmat[indices[i],indices[j]],infmat[indices[j],indices[i]])*max(h[i], h[j])
            M[j,i] = M[i,j]
            
    return M
