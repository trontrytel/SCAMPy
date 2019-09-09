cimport Grid
cimport ReferenceState
cimport TimeStepping

from Cases cimport CasesBase
from Variables cimport GridMeanVariables
from NetCDFIO cimport NetCDFIO_Stats

cimport EDMF_Environment
cimport EDMF_Updrafts
cimport EDMF_Rain

cdef class Scampify1d:
    cdef:
        int it
        int norm

        Grid.Grid Gr
        ReferenceState.ReferenceState Ref
        TimeStepping.TimeStepping TS
        GridMeanVariables GMV
        NetCDFIO_Stats Stats
        CasesBase Case

        EDMF_Rain.RainVariables rain_var
        EDMF_Updrafts.UpdraftVariables upd_var
        EDMF_Environment.EnvironmentVariables env_var

        EDMF_Updrafts.UpdraftThermodynamics upd_thr
        EDMF_Environment.EnvironmentThermodynamics env_thr
        EDMF_Rain.RainPhysics rain_phs

        dict r_dict
        dict p_dict
        dict v_dict

    cpdef do_updrafts(self)
    cpdef do_rain(self)
