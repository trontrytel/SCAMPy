cimport Grid
cimport ReferenceState
from Variables cimport GridMeanVariables
from EDMF_Environment cimport EnvironmentThermodynamics
from EDMF_Updrafts cimport UpdraftThermodynamics
from NetCDFIO cimport NetCDFIO_Stats
from TimeStepping cimport TimeStepping

cdef class PrecipVariable:
    cdef:
        str loc
        str kind
        str name
        str units

        double [:] values
        double [:] new
        double [:] flux

    cpdef set_bcs(self, Grid.Grid Gr)

cdef class PrecipVariables:
    cdef:
        str precip_model

        double mean_rwp
        double env_rwp
        double upd_rwp

        double mean_swp
        double env_swp
        double upd_swp

        Grid.Grid Gr

        PrecipVariable QR
        PrecipVariable RainArea
        PrecipVariable Upd_QR
        PrecipVariable Upd_RainArea
        PrecipVariable Env_QR
        PrecipVariable Env_RainArea

        PrecipVariable QS
        PrecipVariable SnowArea
        PrecipVariable Upd_QS
        PrecipVariable Upd_SnowArea
        PrecipVariable Env_QS
        PrecipVariable Env_SnowArea

    cpdef initialize_io(self, NetCDFIO_Stats Stats)
    cpdef io(self, NetCDFIO_Stats, ReferenceState.ReferenceState Ref)
    cpdef sum_subdomains_precip(self, UpdraftThermodynamics UpdThermo, EnvironmentThermodynamics EnvThermo)
    cpdef precip_diagnostics(self, ReferenceState.ReferenceState Ref)

cdef class PrecipPhysics:
    cdef :
        Grid.Grid Gr
        ReferenceState.ReferenceState Ref

        double [:] rain_evap_source_h
        double [:] rain_evap_source_qt

    cpdef solve_precip_fall(
        self, GridMeanVariables GMV, TimeStepping TS, PrecipVariable QP,
        PrecipVariable PrecipArea
    )

    cpdef solve_rain_evap(
        self, GridMeanVariables GMV, TimeStepping TS, PrecipVariable QR,
        PrecipVariable RainArea
    )
