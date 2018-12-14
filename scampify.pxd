cimport Grid
cimport ReferenceState
cimport TimeStepping
from Variables cimport GridMeanVariables

cimport EDMF_Environment
cimport EDMF_Updrafts
cimport EDMF_Rain

cdef class Scampify1d:
    cdef:
        Grid.Grid Gr
        ReferenceState.ReferenceState Ref
        TimeStepping.TimeStepping TS
        GridMeanVariables GMV

        EDMF_Rain.RainVariables rain_var
        EDMF_Updrafts.UpdraftVariables upd_var
        EDMF_Environment.EnvironmentVariables env_var

        EDMF_Updrafts.UpdraftThermodynamics upd_thr
        EDMF_Environment.EnvironmentThermodynamics env_thr
        EDMF_Rain.RainPhysics rain_phs

    cpdef do_updrafts(self)
