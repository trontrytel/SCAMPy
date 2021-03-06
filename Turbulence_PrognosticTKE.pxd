cimport EDMF_Updrafts
cimport EDMF_Environment
cimport EDMF_Rain

from Grid cimport Grid
from Variables cimport VariablePrognostic, VariableDiagnostic, GridMeanVariables
from Surface cimport  SurfaceBase
from ReferenceState cimport  ReferenceState
from Cases cimport CasesBase
from TimeStepping cimport  TimeStepping
from NetCDFIO cimport NetCDFIO_Stats
from turbulence_functions cimport *
from Turbulence cimport ParameterizationBase

cdef class EDMF_PrognosticTKE(ParameterizationBase):
    cdef:
        Py_ssize_t n_updrafts

        EDMF_Updrafts.UpdraftVariables UpdVar
        EDMF_Updrafts.UpdraftThermodynamics UpdThermo

        EDMF_Environment.EnvironmentVariables EnvVar
        EDMF_Environment.EnvironmentThermodynamics EnvThermo

        EDMF_Rain.RainVariables Rain
        EDMF_Rain.RainPhysics RainPhysics

        entr_struct (*entr_detr_fp) (entr_in_struct entr_in) nogil

        pressure_buoy_struct (*pressure_func_buoy) (pressure_in_struct press_in) nogil
        pressure_buoy_struct (*pressure_func_buoysin) (pressure_in_struct press_in) nogil
        pressure_drag_struct (*pressure_func_drag) (pressure_in_struct press_in) nogil

        bint use_const_plume_spacing
        bint similarity_diffusivity
        bint use_steady_updrafts
        bint calc_scalar_var
        bint calc_tke

        str asp_label
        bint drag_sign
        double surface_area
        double minimum_area
        double entrainment_factor
        double detrainment_factor
        double entrainment_sigma
        double entrainment_scale
        double updraft_mixing_frac
        double sorting_power
        double turbulent_entrainment_factor
        double vel_pressure_coeff # used by diagnostic plume option; now calc'ed from Tan et al 2018 coefficient set
        double vel_buoy_coeff # used by diagnostic plume option; now calc'ed from Tan et al 2018 coefficient set
        double pressure_buoy_coeff # Tan et al. 2018: coefficient alpha_b in Eq. 30
        double pressure_drag_coeff # Tan et al. 2018: coefficient alpha_d in Eq. 30
        double [:] pressure_plume_spacing # Tan et al. 2018: coefficient r_d in Eq. 30
        double pressure_normalmode_buoy_coeff1
        double pressure_normalmode_buoy_coeff2
        double pressure_normalmode_adv_coeff
        double pressure_normalmode_drag_coeff
        double dt_upd
        double constant_plume_spacing
        double aspect_ratio
        double [:,:] entr_sc
        double [:,:] detr_sc
        double [:,:] nh_pressure
        double [:,:] sorting_function
        double [:,:] nh_pressure_adv
        double [:,:] nh_pressure_drag
        double [:,:] nh_pressure_b
        double [:,:] asp_ratio
        double [:,:] b_coeff
        double [:,:] b_mix
        double [:,:] frac_turb_entr
        double [:,:] frac_turb_entr_full
        double [:,:] turb_entr_W
        double [:,:] turb_entr_H
        double [:,:] turb_entr_QT
        double [:,:] horizontal_KM
        double [:,:] horizontal_KH
        double [:] area_surface_bc
        double [:] h_surface_bc
        double [:] qt_surface_bc
        double [:] w_surface_bc
        double [:,:] m # mass flux
        double [:] massflux_h
        double [:] massflux_qt
        double [:] massflux_tke
        double [:] massflux_tendency_h
        double [:] massflux_tendency_qt
        double [:] diffusive_flux_h
        double [:] diffusive_flux_u
        double [:] diffusive_flux_v
        double [:] diffusive_flux_qt
        double [:] diffusive_tendency_h
        double [:] diffusive_tendency_qt
        double [:] mixing_length
        double [:] tke_buoy
        double [:] tke_dissipation
        double [:] tke_entr_gain
        double [:] tke_detr_loss
        double [:] tke_shear
        double [:] tke_pressure
        double [:] tke_transport
        double [:] tke_advection
        double max_area_factor
        double tke_ed_coeff
        double tke_diss_coeff
        double static_stab_coeff
        double lambda_stab

        double [:] Hvar_shear
        double [:] QTvar_shear
        double [:] Hvar_entr_gain
        double [:] QTvar_entr_gain
        double [:] Hvar_detr_loss
        double [:] QTvar_detr_loss
        double [:] Hvar_diss_coeff
        double [:] QTvar_diss_coeff
        double [:] HQTcov
        double [:] HQTcov_shear
        double [:] HQTcov_entr_gain
        double [:] HQTcov_detr_loss
        double [:] HQTcov_diss_coeff
        double [:] Hvar_dissipation
        double [:] QTvar_dissipation
        double [:] HQTcov_dissipation
        double [:] Hvar_rain
        double [:] QTvar_rain
        double [:] HQTcov_rain

        double [:] mls
        double [:] ml_ratio
        double [:] l_entdet
        double [:] b
        double [:] prandtl_nvec
        str mixing_scheme

    cpdef initialize(self, GridMeanVariables GMV)
    cpdef initialize_io(self, NetCDFIO_Stats Stats)
    cpdef io(self, NetCDFIO_Stats Stats)
    cpdef update(self,GridMeanVariables GMV, CasesBase Case, TimeStepping TS)
    cpdef compute_prognostic_updrafts(self, GridMeanVariables GMV, CasesBase Case, TimeStepping TS)
    cpdef compute_diagnostic_updrafts(self, GridMeanVariables GMV, CasesBase Case)
    cpdef update_inversion(self, GridMeanVariables GMV, option)
    cpdef compute_mixing_length(self, double obukhov_length, double ustar, GridMeanVariables GMV)
    cpdef compute_eddy_diffusivities_tke(self, GridMeanVariables GMV, CasesBase Case)
    cpdef compute_horizontal_eddy_diffusivities(self, GridMeanVariables GMV)
    cpdef reset_surface_covariance(self, GridMeanVariables GMV, CasesBase Case)
    cpdef compute_pressure_plume_spacing(self, GridMeanVariables GMV,  CasesBase Case)
    cpdef compute_nh_pressure(self)

    cpdef set_updraft_surface_bc(self, GridMeanVariables GMV, CasesBase Case)
    cpdef decompose_environment(self, GridMeanVariables GMV, whichvals)
    cpdef compute_turbulent_entrainment(self, GridMeanVariables GMV, CasesBase Case)
    cpdef compute_entrainment_detrainment(self, GridMeanVariables GMV, CasesBase Case)
    cpdef zero_area_fraction_cleanup(self, GridMeanVariables GMV)
    cpdef set_subdomain_bcs(self)
    cpdef solve_updraft_velocity_area(self)
    cpdef solve_updraft_scalars(self, GridMeanVariables GMV)
    cpdef update_GMV_MF(self, GridMeanVariables GMV, TimeStepping TS)
    cpdef update_GMV_ED(self, GridMeanVariables GMV, CasesBase Case, TimeStepping TS)
    cpdef compute_covariance(self, GridMeanVariables GMV, CasesBase Case, TimeStepping TS)

    cpdef initialize_covariance(self, GridMeanVariables GMV, CasesBase Case)
    cpdef cleanup_covariance(self, GridMeanVariables GMV)
    cpdef compute_tke_buoy(self, GridMeanVariables GMV)
    cpdef compute_tke_pressure(self)
    cdef void compute_covariance_dissipation(self, EDMF_Environment.EnvironmentVariable_2m Covar)
    cdef void compute_covariance_entr(self, EDMF_Environment.EnvironmentVariable_2m Covar, EDMF_Updrafts.UpdraftVariable UpdVar1,
                EDMF_Updrafts.UpdraftVariable UpdVar2, EDMF_Environment.EnvironmentVariable EnvVar1, EDMF_Environment.EnvironmentVariable EnvVar2,
                VariablePrognostic GmvVar1, VariablePrognostic GmvVar2)
    cdef void compute_covariance_detr(self, EDMF_Environment.EnvironmentVariable_2m Covar)
    cdef void compute_covariance_shear(self,GridMeanVariables GMV, EDMF_Environment.EnvironmentVariable_2m Covar,
                                       double *UpdVar1, double *UpdVar2, double *EnvVar1, double *EnvVar2)
    cpdef compute_covariance_rain(self, TimeStepping TS, GridMeanVariables GMV)
    cdef void compute_covariance_interdomain_src(self, EDMF_Updrafts.UpdraftVariable au, EDMF_Updrafts.UpdraftVariable phi_u, EDMF_Updrafts.UpdraftVariable psi_u,
                        EDMF_Environment.EnvironmentVariable phi_e,  EDMF_Environment.EnvironmentVariable psi_e, EDMF_Environment.EnvironmentVariable_2m covar_e)
    cdef void update_covariance_ED(self, GridMeanVariables GMV, CasesBase Case,TimeStepping TS, VariablePrognostic GmvVar1, VariablePrognostic GmvVar2,
            VariableDiagnostic GmvCovar, EDMF_Environment.EnvironmentVariable_2m Covar, EDMF_Environment.EnvironmentVariable  EnvVar1, EDMF_Environment.EnvironmentVariable  EnvVar2,
            EDMF_Updrafts.UpdraftVariable UpdVar1, EDMF_Updrafts.UpdraftVariable UpdVar2)
    cpdef compute_tke_transport(self)
    cpdef compute_tke_advection(self)
    cpdef update_GMV_diagnostics(self, GridMeanVariables GMV)
    cpdef double compute_zbl_qt_grad(self, GridMeanVariables GMV)
    cdef get_GMV_CoVar(self, EDMF_Updrafts.UpdraftVariable au,
                        EDMF_Updrafts.UpdraftVariable phi_u, EDMF_Updrafts.UpdraftVariable psi_u,
                        EDMF_Environment.EnvironmentVariable phi_e,  EDMF_Environment.EnvironmentVariable psi_e,
                        EDMF_Environment.EnvironmentVariable_2m covar_e,
                       double *gmv_phi, double *gmv_psi, double *gmv_covar)
    cdef get_env_covar_from_GMV(self, EDMF_Updrafts.UpdraftVariable au,
                                EDMF_Updrafts.UpdraftVariable phi_u, EDMF_Updrafts.UpdraftVariable psi_u,
                                EDMF_Environment.EnvironmentVariable phi_e, EDMF_Environment.EnvironmentVariable psi_e,
                                EDMF_Environment.EnvironmentVariable_2m covar_e,
                                double *gmv_phi, double *gmv_psi, double *gmv_covar)

    cdef void GMV_third_m(self, VariableDiagnostic Gmv_third_m, EDMF_Environment.EnvironmentVariable_2m env_covar,
                           EDMF_Environment.EnvironmentVariable  env_mean, EDMF_Updrafts.UpdraftVariable  upd_mean)