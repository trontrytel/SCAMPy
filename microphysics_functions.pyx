import numpy as np
cimport numpy as np

from libc.math cimport fmax

from thermodynamic_functions cimport *

# instantly remove all cloud water that is above some threshold
# the threshold is specified as axcess saturation
cdef double acnv_rate(double ql, double qt, double sat_treshold, double T, double p0) nogil :

    cdef double psat = pv_star(T)
    cdef double qsat = qv_star_c(p0, qt, psat)

    return fmax(0.0, ql - sat_treshold * qsat)
