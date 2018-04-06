cdef double r2q(double r_, double qt) nogil
cdef double q2r(double q_, double qt) nogil
cdef double acnv_instant(double ql, double qt, double sat_treshold, double T, double p0) nogil
cdef double acnv_rate(double ql, double qt) nogil
cdef double accr_rate(double ql, double qr, double qt) nogil
cdef double evap_rate(double rho, double qv, double qr, double qt, double T, double p0) nogil
cdef double terminal_velocity(double rho, double rho0, double qr, double qt) nogil
