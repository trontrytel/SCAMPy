cimport Grid
cimport ReferenceState
from Variables cimport GridMeanVariables
from NetCDFIO cimport NetCDFIO_Stats
from TimeStepping cimport TimeStepping

cdef class RainVariable:
    cdef:
        str loc
        str kind
        str name
        str units

        double [:] values
        double [:] new
        double [:] flux

    cpdef set_bcs(self, Grid.Grid Gr)

cdef class RainVariables:
    cdef:
        bint rain_model
        bint rain_const_area

        double puddle
        double rain_area_value
        double max_supersaturation

        Grid.Grid Gr

        RainVariable QR
        RainVariable RainArea
        RainVariable Upd_QR
        RainVariable Upd_RainArea
        RainVariable Env_QR
        RainVariable Env_RainArea

    cpdef initialize_io(self, NetCDFIO_Stats Stats)
    cpdef io(self, NetCDFIO_Stats)
    cpdef set_values_with_new(self)

cdef class RainPhysics:
    cdef :
        Grid.Grid Gr
        ReferenceState.ReferenceState Ref

    cpdef solve_rain_fall(self, GridMeanVariables GMV, TimeStepping TS,
                          RainVariable QR, RainVariable RainArea,
                          double rain_area_value)
