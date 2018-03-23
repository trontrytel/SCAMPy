from NetCDFIO cimport NetCDFIO_Stats
from Grid cimport  Grid
from ReferenceState cimport ReferenceState
from Variables cimport VariableDiagnostic, GridMeanVariables

cdef class EnvironmentVariable:
    cdef:
        double [:] values
        double [:] flux
        str loc
        str kind
        str name
        str units

cdef class EnvironmentVariables:
    cdef:
        EnvironmentVariable W
        EnvironmentVariable QT
        EnvironmentVariable QL
        EnvironmentVariable H
        EnvironmentVariable THL
        EnvironmentVariable T
        EnvironmentVariable B
        EnvironmentVariable TKE
        EnvironmentVariable Hvar
        EnvironmentVariable QTvar
        EnvironmentVariable HQTcov
        EnvironmentVariable CF
        Grid Gr
        bint use_tke
        bint use_scalar_var
        bint use_prescribed_scalar_var
        double prescribed_QTvar
        double prescribed_Hvar
        double prescribed_HQTcov

    cpdef initialize_io(self, NetCDFIO_Stats Stats )
    cpdef io(self, NetCDFIO_Stats Stats)

cdef class EnvironmentThermodynamics:
    cdef:
        Grid Gr
        ReferenceState Ref
        Py_ssize_t quadrature_order
        double (*t_to_prog_fp)(double p0, double T,  double qt, double ql, double qi)   nogil
        double (*prog_to_t_fp)(double H, double pd, double pv, double qt ) nogil
        void eos_update_SA_sgs(self, EnvironmentVariables EnvVar, VariableDiagnostic GMV_B)
        double [:] qt_dry
        double [:] th_dry
        double [:] t_cloudy
        double [:] qv_cloudy
        double [:] qt_cloudy
        double [:] th_cloudy
        double max_supersaturation

    cpdef satadjust(self, EnvironmentVariables EnvVar, GridMeanVariables GMV)
