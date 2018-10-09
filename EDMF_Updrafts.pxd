cimport Grid
cimport ReferenceState
from Variables cimport GridMeanVariables
from NetCDFIO cimport NetCDFIO_Stats
from EDMF_Environment cimport EnvironmentVariables

cdef class UpdraftVariable:
    cdef:
        double [:,:] values
        double [:,:] new
        double [:,:] old
        double [:,:] tendencies
        double [:,:] flux
        double [:] bulkvalues
        str loc
        str kind
        str name
        str units
    cpdef set_bcs(self, Grid.Grid Gr)

cdef class UpdraftVariables:
    cdef:
        Grid.Grid Gr
        UpdraftVariable W
        UpdraftVariable Area
        UpdraftVariable QT
        UpdraftVariable QL
        UpdraftVariable H
        UpdraftVariable THL
        UpdraftVariable T
        UpdraftVariable B
        Py_ssize_t n_updrafts
        bint prognostic
        double updraft_fraction
        double [:] cloud_base
        double [:] cloud_top
        double [:] cloud_cover
    cpdef initialize(self, GridMeanVariables GMV)
    cpdef initialize_io(self, NetCDFIO_Stats Stats)
    cpdef io(self, NetCDFIO_Stats Stats)
    cpdef set_means(self, GridMeanVariables GMV)
    cpdef set_new_with_values(self)
    cpdef set_old_with_values(self)
    cpdef set_values_with_new(self)
    cpdef get_cloud_base_top_cover(self)

cdef class UpdraftRain:
    cdef:
        Grid.Grid Gr
        Py_ssize_t n_updrafts

        UpdraftVariable QR
        UpdraftVariable RainArea

        double puddle

        bint rain_model

    cpdef initialize_io(self, NetCDFIO_Stats Stats)
    cpdef io(self, NetCDFIO_Stats Stats)

    cpdef set_means(self, GridMeanVariables GMV)
    cpdef set_new_with_values(self)
    cpdef set_old_with_values(self)
    cpdef set_values_with_new(self)

cdef class UpdraftThermodynamics:
    cdef:
        double (*t_to_prog_fp)(double p0, double T,  double qt, double ql, double qi)   nogil
        double (*prog_to_t_fp)(double H, double pd, double pv, double qt ) nogil
        Grid.Grid Gr
        ReferenceState.ReferenceState Ref
        Py_ssize_t n_updraft

    #cpdef satadjust(self, UpdraftVariables UpdVar) TODO - not used?
    cpdef buoyancy(self, UpdraftVariables UpdVar, EnvironmentVariables EnvVar,
                   GridMeanVariables GMV, bint extrap)

cdef class UpdraftMicrophysics:
    cdef:
        Grid.Grid Gr
        ReferenceState.ReferenceState Ref
        Py_ssize_t n_updraft
        double [:,:] prec_source_h
        double [:,:] prec_source_qt
        double [:]  prec_source_h_tot
        double [:] prec_source_qt_tot
        double max_supersaturation

    cpdef compute_column_sources(self, UpdraftVariables UpdVar)
    cpdef update_column_UpdVar(self,   UpdraftVariables UpdVar)
    cpdef update_column_UpdRain(self,  UpdraftVariables UpdVar, UpdraftRain UpdRain)

    cpdef cleanup_column_UpdRain(self, UpdraftRain UpdRain, double eps)

    cdef void update_UpdVar(self, double *qt, double *ql, double *h, double *T,
                            double qr_src, double th_src, double qt_new, double ql_new, double T_new, double thl_new,
                            Py_ssize_t i, Py_ssize_t k) nogil

    cdef void update_UpdRain(self, double *upd_area, double *qr, double *rain_Area, double qr_new, Py_ssize_t i, Py_ssize_t k) nogil
