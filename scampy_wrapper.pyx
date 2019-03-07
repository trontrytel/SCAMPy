# wrapper from cdef (not accessible from Python) to def (accessible from Python)
# written to be able to use the thermodynamic functions from Pytest

cimport thermodynamic_functions as fun
include "parameters.pxi"

def exner_c(p0):
    return fun.exner_c(p0)

def pv_star(T):
    return fun.pv_star(T)

# assuming qt = qv + ql + qi
def qv_star_c(p0, qt, pv):
    return fun.qv_star_c(p0, qt, pv)
