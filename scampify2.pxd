cimport Grid
cimport ReferenceState
cimport EDMF_Environment
cimport EDMF_Updrafts

cdef class Scampify1d:
    cdef:
        Grid.Grid Gr
        ReferenceState.ReferenceState ref_state

        EDMF_Environment.EnvironmentVariables env_var
        EDMF_Environment.EnvironmentThermodynamics env_thr

        EDMF_Updrafts.UpdraftVariables upd_var
        EDMF_Updrafts.UpdraftThermodynamics upd_thr
        EDMF_Updrafts.UpdraftMicrophysics upd_mcr

    cpdef do_updrafts(self)

