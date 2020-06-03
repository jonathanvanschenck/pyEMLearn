cimport numpy as cnp
import numpy as np
#cimport scipy.linalg.cython_lapack as clap
from libc.stdlib cimport malloc, free

ctypedef cnp.float64_t REAL_t

def inv(mat):
  cdef int k = mat.shape[0]
  cdef int info = 0

  cdef REAL_t* mat_pointer = <REAL_t*>cnp.PyArray_DATA(mat)
  cdef REAL_t* iden_pointer

  cdef int* piv_pointer = <int*>malloc(sizeof(int)*k)

  try:
    identity = np.eye(k)
    iden_pointer = <REAL_t*>cnp.PyArray_DATA(identity)

    clap.dgsev(&k, &k, mat_pointer, &k,
               piv_pointer, iden_pointer, &k, &info)

    # add info check?

    return identity

  finally:
    # Ensure the piv allocation gets freed
    free(piv_pointer)
