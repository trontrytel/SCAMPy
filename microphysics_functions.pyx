import numpy as np
cimport numpy as np
import sys
from libc.math cimport fmax, exp, sqrt, fmin
from thermodynamic_functions cimport *
include "parameters.pxi"

# cutoff microphysics parameter
cdef double max_supersaturation = 0.02


# CLIMA microphysics parameters

# assumed intercept parameters for rain and ice size distributions
# (assuming exponential distr)
cdef double MP_n_0_r = 16 * 1e6
cdef double MP_n_0_i = 1 * 1e7

# density of water and ice
cdef double rho_cloud_liq = 1e3
cdef double rho_cloud_ice = 900

# assumed relationship between particle radius r and mass m
# m(r) = a * r^b
cdef double a_rai = 4./3 * pi * rho_cloud_liq
cdef double b_rai = 3.
cdef double a_ice = 4./3 * pi * rho_cloud_ice
cdef double b_ice = 3.
cdef double a_sno = 2.5 * 1e-2 # Grabowski 1998
cdef double b_sno = 2.         # Grabowski 1998

# assumed timescales of the processes
cdef double tau_cond_evap_cloud = 10.
cdef double tau_autoconv_rain = 1e3
cdef double tau_autoconv_snow = 3 * 60

# drag coefficient for rain drops
cdef double C_drag = 0.55
cdef double terminal_velocity_single_drop_coeff(double rho) nogil:
    return sqrt(8./3 / C_drag * (rho_cloud_liq / rho - 1.))
# assumed relationship between particle radius r and terminal velocity v
# v(r) = c * r^d
cdef double c_rai(double rho) nogil:
    return terminal_velocity_single_drop_coeff(rho) * sqrt(g)
cdef double d_rai = 0.5
cdef double c_sno = 4. / pow(2, 0.25) # Grabowski 1998
cdef double d_sno = 0.25              # Grabowski 1998

# rain autoconversion threshold (TODO - test for smooth funstion instead)
cdef double q_liq_threshold = 5e-4

# assumed collision efficiency for rain
cdef double E_col = 0.8

# rain evaporation coefficients
cdef double nu_air = 1.6e-5
cdef double K_therm = 2.4e-2
cdef double D_vapor = 2.26e-5
cdef double a_vent = 1.5
cdef double b_vent = 0.53


cdef double precip_source_to_thetal(double p0, double T, double qp) nogil :
    """
    Source term for thetal because of qr or qs transitioning
    between the working fluid and rain (simple version to avoid exponents)
    """
    return latent_heat(T) * qp / exner_c(p0) / cpd

# TODO
#cdef double precip_source_to_thetal_detailed(double p0, double T, double qt, double ql, double qr) nogil :
#    """
#    Source term for thetal because of qr transitioning between the working fluid and rain
#    (more detailed version, but still ignoring dqt/dqr)
#    """
#    cdef double L = latent_heat(T)
#
#    old_source = L * qr / exner_c(p0) / cpd
#
#    new_source = old_source / (1.-qt) * exp(-L * ql / T / cpd / (1.-qt))
#
#    return new_source

cdef double acnv_instant(double ql, double qt, double T, double p0) nogil :
    """
    instantly convert all cloud water exceeding a threshold to rain water
    the threshold is specified as axcess saturation
    rain water is immediately removed from the domain
    """
    cdef double psat = pv_star(T)
    cdef double qsat = qv_star_c(p0, qt, psat)

    return fmax(0.0, ql - max_supersaturation * qsat)

# CLIMA microphysics rates

cdef double lambda_param(double a_, double b_, double MP_n_0,
                         double rho, double q_):
    """
    Slope parameter of the assumed exponential dstr.
    Input:
    - a_, b_ - m(r) relationship parameters
    - MP_n_0 - intercept parameter
    - air density
    - q_ - specific humidity
    """
    return pow(a_ * MP_n_0 * gamma(b_ + 1.) / rho / q_, 1. / (b_+1))

cdef double terminal_velocity(double q_rai, double rho) nogil:

    cdef double v_c = terminal_velocity_single_drop_coeff(rho)
    cdef double gamma_9_2 = 11.631728396567448

    cdef double term_vel = 0.

    if q_rai > 0:
        lambda_param = (8. * pi * rho_cloud_liq * MP_n_0 / rho / q_rai)**0.25
        term_vel = gamma_9_2 * v_c / 6. * sqrt(g / lambda_param)

    return term_vel


cdef double conv_q_vap_to_q_liq(double q_sat_liq, double q_liq) nogil:

  return (q_sat_liq - q_liq) / tau_cond_evap


cdef double conv_q_liq_to_q_rai_acnv(double q_liq) nogil:

  return fmax(0., q_liq - q_liq_threshold) / tau_acnv


cdef double conv_q_liq_to_q_rai_accr(double q_liq, double q_rai, double rho) nogil:

  cdef double v_c = terminal_velocity_single_drop_coeff(rho)
  cdef double gamma_7_2 = 3.3233509704478426

  cdef double accr_coeff = gamma_7_2 * 8.**(-7./8) * pi**(1./8) * v_c * E_col *\
                           (rho / rho_cloud_liq)**(7./8)

  return accr_coeff * MP_n_0**(1./8) * sqrt(g) * q_liq * q_rai**(7./8)


cdef double conv_q_rai_to_q_vap(double q_rai, double q_tot, double q_liq, double q_ice,
                                  double T, double p, double rho) nogil:

  cdef double L = latent_heat(T)
  cdef double gamma_11_4 = 1.6083594219855457
  cdef double v_c = terminal_velocity_single_drop_coeff(rho)
  cdef double N_Sc = nu_air / D_vapor

  cdef double av_param = sqrt(2. * pi) * a_vent * sqrt(rho / rho_cloud_liq)
  cdef double bv_param = 2.**(7./16) * gamma_11_4 * pi**(5./16) * b_vent *\
                         N_Sc**(1./3) * sqrt(v_c) * (rho / rho_cloud_liq)**(11./16)

  cdef double p_vs = pv_star(T)
  cdef double qv_sat = qv_star_c(p, q_tot, p_vs)
  cdef double q_v = q_tot - q_liq - q_ice
  cdef double S = q_v/qv_sat - 1

  cdef double G_param = 1. / (L / K_therm / T * (L / Rv / T - 1.) +\
                              Rv * T / D_vapor / p_vs\
                             )

  cdef double F_param = av_param * sqrt(q_rai) +\
                        bv_param * g**0.25 / MP_n_0**(3./16) /\
                          sqrt(nu_air) * q_rai**(11./16)

  return S * F_param * G_param * sqrt(MP_n_0) / rho


# TODO - recode to reuse the stuff from warm rain param
cdef double conv_q_ice_to_q_snw(double q_ice, double q_tot, double q_liq,
                                  double T, double p, double rho) nogil:

  cdef double r_0 = 125 * 1e-6

  cdef double lmbd = (8. * pi * rho_cloud_ice * MP_n_0_i / rho / q_ice)**0.25

  cdef double L = latent_heat(T)

  cdef double p_vs = pv_star(T)

  cdef double qv_sat = qv_star_c(p, q_tot, p_vs)
  cdef double q_v = q_tot - q_liq - q_ice
  cdef double S = q_v/qv_sat - 1

  cdef double G_param = 1. / (L / K_therm / T * (L / Rv / T - 1.) +\
                              Rv * T / D_vapor / p_vs\
                             )

  cdef double tmp = max(0.,\
                        4 * pi / rho * S * G_param * MP_n_0_i * exp(- lmbd * r_0) *\
                          (r_0**2 / 3. + (r_0 * lmbd + 1.) / lmbd**2))

  tmp *= 1e10

  #with gil:
  #    if q_ice > 0:
  #        print "qs_src = ", tmp

  return tmp


cdef mph_struct microphysics_precip_src(
                  str rain_model,
                  str snow_model,
                  double qt,
                  double ql,
                  double qr,
                  double qi,
                  double qs,
                  double area,
                  double T,
                  double p0,
                  double rho,
                  double dt) nogil:

    """
    compute the autoconversion and accretion rate for rain
    compute autoconversion rate for snow
    return:
      new values: qt, ql, qv, qi, thl, th, alpha
      rates: qr_src, qs_src, thl_precip_src
    """
    # TODO assumes no ice
    cdef mph_struct _ret
    _ret.qv    = qt - ql - qi
    _ret.thl   = t_to_thetali_c(p0, T, qt, ql, qi)
    _ret.th    = theta_c(p0, T)
    _ret.alpha = alpha_c(p0, T, qt, _ret.qv)

    #TODO - temporary way to handle different autoconversion rates
    # cython doesn't allow for string comparison without gil
    tmp_clima_acnv_flag = False
    tmp_cutoff_acnv_flag = False
    tmp_no_acnv_flag = False
    tmp_snow_acnv_flag = False
    with gil:
        if rain_model == 'clima_1m':
            tmp_clima_acnv_flag = True
        elif rain_model == 'cutoff':
            tmp_cutoff_acnv_flag = True
        elif rain_model == 'None':
            tmp_no_acnv_flag = True
        else:
            sys.exit('rain model not recognized')

    tmp_snow_acnv_flag = False
    with gil:
        if snow_model == 'snow_testing':
            tmp_snow_acnv_flag = True

    if area > 0.:
        if tmp_clima_acnv_flag:
            _ret.qr_src = fmin(ql,
                                  (conv_q_liq_to_q_rai_acnv(ql) +
                                   conv_q_liq_to_q_rai_accr(ql, qr, rho)) * dt
                              )

        if tmp_cutoff_acnv_flag:
            _ret.qr_src = fmin(ql, acnv_instant(ql, qt, T, p0))

        if tmp_no_acnv_flag:
            _ret.qr_src = 0.

        if tmp_snow_acnv_flag and qi > 0.:
            _ret.qs_src = fmin(qi, conv_q_ice_to_q_snw(qi, qt, ql, T, p0, rho))
        else:
            _ret.qs_src = 0.

        _ret.thl_precip_src = precip_source_to_thetal(p0, T, _ret.qr_src + _ret.qs_src)

    else:
        _ret.qr_src = 0.
        _ret.qs_src = 0.
        _ret.thl_precip_src = 0.

    _ret.qt = qt - _ret.qr_src - _ret.qs_src
    _ret.ql = ql - _ret.qr_src
    _ret.qi = qi - _ret.qs_src

    _ret.thl += _ret.thl_precip_src

    return _ret

cdef rain_struct rain_area(double source_area,  double source_qr,
                           double current_area, double current_qr ) nogil:
    """
    Source terams for rain and rain area
    assuming constant rain area fraction of 1
    """
    cdef rain_struct _ret

    if source_qr <= 0.:
        _ret.qr = current_qr
        _ret.ar = current_area
    else:
        _ret.qr = current_qr + source_area * source_qr
        _ret.ar = 1.

    # sketch of what to do for prognostic rain area fraction:

    #cdef double a_big, q_big, a_sml, q_sml
    #cdef double a_const = 0.2
    #cdef double eps     = 1e-5

    #if source_qr ==  0.:
    #    _ret.qr = current_qr
    #    _ret.ar = current_area
    #else:
    #    if current_area != 0.:
    #        if current_area >= source_area:
    #            a_big = current_area
    #            q_big = current_qr
    #            a_sml = source_area
    #            q_sml = source_qr
    #        else:
    #            a_sml = current_area
    #            q_sml = current_qr
    #            a_big = source_area
    #            q_big = source_qr

    #        _ret.qr = q_big + a_sml / a_big * q_sml
    #        _ret.ar = a_big

    #    else:
    #        _ret.qr = source_qr
    #        _ret.ar = source_area

    return _ret

