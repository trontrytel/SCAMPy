#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: initializedcheck=False
#cython: cdivision=True

import numpy as np
import sys
include "parameters.pxi"
import cython
from Grid cimport  Grid
from TimeStepping cimport TimeStepping
from ReferenceState cimport ReferenceState
from Variables cimport VariableDiagnostic, GridMeanVariables
from libc.math cimport fmax, fmin, sqrt, exp, erf
from thermodynamic_functions cimport  *
from microphysics_functions cimport *

cdef class EnvironmentVariable:
    def __init__(self, nz, loc, kind, name, units):
        self.values = np.zeros((nz,),dtype=np.double, order='c')
        self.flux = np.zeros((nz,),dtype=np.double, order='c')
        if loc != 'half' and loc != 'full':
            print('Invalid location setting for variable! Must be half or full')
        self.loc = loc
        if kind != 'scalar' and kind != 'velocity':
            print ('Invalid kind setting for variable! Must be scalar or velocity')
        self.kind = kind
        self.name = name
        self.units = units


cdef class EnvironmentVariables:
    def __init__(self,  namelist, Grid Gr  ):
        cdef Py_ssize_t nz = Gr.nzg
        self.Gr = Gr

        self.W = EnvironmentVariable(nz, 'full', 'velocity', 'w','m/s' )
        self.QT = EnvironmentVariable( nz, 'half', 'scalar', 'qt','kg/kg' )
        self.QL = EnvironmentVariable( nz, 'half', 'scalar', 'ql','kg/kg' )
        self.QR = EnvironmentVariable( nz, 'half', 'scalar', 'qr','kg/kg' )
        if namelist['thermodynamics']['thermal_variable'] == 'entropy':
            self.H = EnvironmentVariable( nz, 'half', 'scalar', 's','J/kg/K' )
        elif namelist['thermodynamics']['thermal_variable'] == 'thetal':
            self.H = EnvironmentVariable( nz, 'half', 'scalar', 'thetal','K' )
        self.THL = EnvironmentVariable(nz, 'half', 'scalar', 'thetal', 'K')
        self.T = EnvironmentVariable( nz, 'half', 'scalar', 'temperature','K' )
        self.B = EnvironmentVariable( nz, 'half', 'scalar', 'buoyancy','m^2/s^3' )
        self.CF = EnvironmentVariable(nz, 'half', 'scalar','cloud_fraction', '-')

        # TKE
        # TODO - kind of repeated from Variables.pyx logic
        if  namelist['turbulence']['scheme'] == 'EDMF_PrognosticTKE':
            self.use_tke = True
        else:
            self.use_tke = False

        try:
            self.use_scalar_var = namelist['turbulence']['EDMF_PrognosticTKE']['use_scalar_var']
        except:
            self.use_scalar_var = False
            print('Defaulting to non-calculation of scalar variances')

        try:
            self.EnvThermo_scheme = str(namelist['thermodynamics']['saturation'])
        except:
            self.EnvThermo_scheme = 'sa_mean'
            print('Defaulting to saturation adjustment with respect to environmental means')

        if self.use_tke:
            self.TKE = EnvironmentVariable( nz, 'half', 'scalar', 'tke','m^2/s^2' )

        if self.use_scalar_var:
            self.QTvar = EnvironmentVariable( nz, 'half', 'scalar', 'qt_var','kg^2/kg^2' )
            if namelist['thermodynamics']['thermal_variable'] == 'entropy':
                self.Hvar = EnvironmentVariable(nz, 'half', 'scalar', 's_var', '(J/kg/K)^2')
                self.HQTcov = EnvironmentVariable(nz, 'half', 'scalar', 's_qt_covar', '(J/kg/K)(kg/kg)' )
            elif namelist['thermodynamics']['thermal_variable'] == 'thetal':
                self.Hvar = EnvironmentVariable(nz, 'half', 'scalar', 'thetal_var', 'K^2')
                self.HQTcov = EnvironmentVariable(nz, 'half', 'scalar', 'thetal_qt_covar', 'K(kg/kg)' )
                if self.EnvThermo_scheme == 'sommeria_deardorff':
                    self.THVvar = EnvironmentVariable(nz, 'half', 'scalar', 'thetav_var', 'K^2' )

        #TODO  - most likely a temporary solution (unless it could be useful for testing)
        try:
            self.use_prescribed_scalar_var = namelist['turbulence']['sgs']['use_prescribed_scalar_var']
        except:
            self.use_prescribed_scalar_var = False
        if self.use_prescribed_scalar_var == True:
            self.prescribed_QTvar  = namelist['turbulence']['sgs']['prescribed_QTvar']
            self.prescribed_Hvar   = namelist['turbulence']['sgs']['prescribed_Hvar']
            self.prescribed_HQTcov = namelist['turbulence']['sgs']['prescribed_HQTcov']

        if (self.EnvThermo_scheme == 'sommeria_deardorff' or self.EnvThermo_scheme == 'sa_quadrature'):
            if (self.use_scalar_var == False and self.use_prescribed_scalar_var == False ):
                sys.exit('EDMF_Environment.pyx 96: scalar variance has to be specified for Sommeria Deardorff or quadrature saturation')

        return

    cpdef initialize_io(self, NetCDFIO_Stats Stats):
        Stats.add_profile('env_w')
        Stats.add_profile('env_qt')
        Stats.add_profile('env_ql')
        Stats.add_profile('env_qr')
        if self.H.name == 's':
            Stats.add_profile('env_s')
        else:
            Stats.add_profile('env_thetal')
        Stats.add_profile('env_temperature')
        if self.use_tke:
            Stats.add_profile('env_tke')
        if self.use_scalar_var:
            Stats.add_profile('env_Hvar')
            Stats.add_profile('env_QTvar')
            Stats.add_profile('env_HQTcov')
        if self.EnvThermo_scheme == 'sommeria_deardorff':
            Stats.add_profile('env_THVvar')
        return

    cpdef io(self, NetCDFIO_Stats Stats):
        Stats.write_profile('env_w', self.W.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
        Stats.write_profile('env_qt', self.QT.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
        Stats.write_profile('env_ql', self.QL.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
        Stats.write_profile('env_qr', self.QR.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
        if self.H.name == 's':
            Stats.write_profile('env_s', self.H.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
        else:
            Stats.write_profile('env_thetal', self.H.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])

        Stats.write_profile('env_temperature', self.T.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
        if self.use_tke:
            Stats.write_profile('env_tke', self.TKE.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
        if self.use_scalar_var:
            Stats.write_profile('env_Hvar', self.Hvar.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
            Stats.write_profile('env_QTvar', self.QTvar.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
            Stats.write_profile('env_HQTcov', self.HQTcov.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])
        if self.EnvThermo_scheme  == 'sommeria_deardorff':
            Stats.write_profile('env_THVvar', self.THVvar.values[self.Gr.gw:self.Gr.nzg-self.Gr.gw])

        #ToDo [suggested by CK for AJ ;]
        # Add output of environmental cloud fraction, cloud base, cloud top (while the latter can be gleaned from ql profiles
        # it is more convenient to simply have them in the stats files!
        # Add the same with respect to the grid mean
        return

cdef class EnvironmentThermodynamics:
    def __init__(self, namelist, paramlist, Grid Gr, ReferenceState Ref, EnvironmentVariables EnvVar):
        self.Gr = Gr
        self.Ref = Ref
        try:
            self.quadrature_order = namelist['condensation']['quadrature_order']
        except:
            self.quadrature_order = 5
        if EnvVar.H.name == 's':
            self.t_to_prog_fp = t_to_entropy_c
            self.prog_to_t_fp = eos_first_guess_entropy
        elif EnvVar.H.name == 'thetal':
            self.t_to_prog_fp = t_to_thetali_c
            self.prog_to_t_fp = eos_first_guess_thetal

        self.qt_dry = np.zeros(self.Gr.nzg, dtype=np.double, order='c')
        self.th_dry = np.zeros(self.Gr.nzg, dtype=np.double, order='c')

        self.t_cloudy  = np.zeros(self.Gr.nzg, dtype=np.double, order ='c')
        self.qv_cloudy = np.zeros(self.Gr.nzg, dtype=np.double, order ='c')
        self.qt_cloudy = np.zeros(self.Gr.nzg, dtype=np.double, order='c')
        self.th_cloudy = np.zeros(self.Gr.nzg, dtype=np.double, order='c')

        self.Hvar_rain_dt   = np.zeros(self.Gr.nzg, dtype=np.double, order='c')
        self.QTvar_rain_dt  = np.zeros(self.Gr.nzg, dtype=np.double, order='c')
        self.HQTcov_rain_dt = np.zeros(self.Gr.nzg, dtype=np.double, order='c')

        self.max_supersaturation = paramlist['turbulence']['updraft_microphysics']['max_supersaturation']

        return

    cdef void update_EnvVar(self, int k, EnvironmentVariables EnvVar, double T, double H, double qt, double ql, double qr, double alpha) nogil :

        EnvVar.T.values[k]   = T
        EnvVar.THL.values[k] = H
        EnvVar.H.values[k]   = H
        EnvVar.QT.values[k]  = qt
        EnvVar.QL.values[k]  = ql
        EnvVar.QR.values[k]  = qr
        EnvVar.B.values[k]   = buoyancy_c(self.Ref.alpha0_half[k], alpha)
        return

    cdef void update_cloud_dry(self, int k, EnvironmentVariables EnvVar, double T, double H, double qt, double ql, double qv) nogil :

        if ql > 0.0:
            EnvVar.CF.values[k] = 1.
            self.th_cloudy[k]   = H
            self.t_cloudy[k]    = T
            self.qt_cloudy[k]   = qt
            self.qv_cloudy[k]   = qv
        else:
            EnvVar.CF.values[k] = 0.
            self.th_dry[k]      = H
            self.qt_dry[k]      = qt
        return

    cdef void eos_update_SA_mean(self, EnvironmentVariables EnvVar, bint in_Env):

        cdef:
            Py_ssize_t k
            Py_ssize_t gw = self.Gr.gw

            eos_struct sa
            mph_struct mph

        if EnvVar.H.name != 'thetal':
            sys.exit('EDMF_Environment: rain source terms are defined for thetal as model variable')

        with nogil:
            for k in xrange(gw,self.Gr.nzg-gw):
                # condensation + autoconversion
                sa  = eos(self.t_to_prog_fp, self.prog_to_t_fp, self.Ref.p0_half[k], EnvVar.QT.values[k], EnvVar.H.values[k])
                mph = microphysics(sa.T, sa.ql, self.Ref.p0_half[k], EnvVar.QT.values[k], EnvVar.H.values[k], self.max_supersaturation, in_Env)

                self.update_EnvVar(   k, EnvVar, mph.T, mph.H, mph.qt, mph.ql, mph.qr, mph.alpha)
                self.update_cloud_dry(k, EnvVar, mph.T, mph.H, mph.qt, mph.ql, mph.qv)
        return

    cdef void eos_update_SA_sgs(self, EnvironmentVariables EnvVar, bint in_Env):
        a, w = np.polynomial.hermite.hermgauss(self.quadrature_order)

        #TODO - remember you output source terms multipierd by dt (bec. of instanteneous autoconcv)
        #TODO - read prescribed var/covar from file to compare with LES data
        #TODO - add tendencies for GMV H, QT and QR due to rain

        cdef:
            Py_ssize_t gw = self.Gr.gw
            Py_ssize_t k, m_q, m_h
            double [:] abscissas = a
            double [:] weights = w

            inner_env = {'ql':0.0, 'T':0.0, 'thl':0.0, 'alpha':0.0, 'cf':0.0, 'qr':0.0, 'qt_cld': 0.0, 'qt_dry':0.0, 'T_cld':0.0, 'T_dry':0.0}
            outer_env = {'ql':0.0, 'T':0.0, 'thl':0.0, 'alpha':0.0, 'cf':0.0, 'qr':0.0, 'qt_cld': 0.0, 'qt_dry':0.0, 'T_cld':0.0, 'T_dry':0.0}
            inner_src = {'SH_qt' : 0.0, 'Sqt_H': 0.0, 'SH_H':0.0, 'Sqt_qt':0.0 , 'Sqt':0.0, 'SH':0.0, 'H':0.0, 'qt':0.0}
            outer_src = {'SH_qt' : 0.0, 'Sqt_H': 0.0, 'SH_H':0.0, 'Sqt_qt':0.0 , 'Sqt':0.0, 'SH':0.0, 'H':0.0, 'qt':0.0}

            double h_hat, qt_hat, sd_h, sd_q, corr, mu_h_star, sigma_h_star, qt_var
            double sqpi_inv = 1.0/sqrt(pi)
            double sqrt2 = sqrt(2.0)
            double sd_q_lim
            eos_struct sa
            mph_struct mph

        if EnvVar.H.name != 'thetal':
            sys.exit('EDMF_Environment: rain source terms are only defined for thetal as model variable')

        # TODO - for testing (to be removed)
        if EnvVar.use_prescribed_scalar_var:
            for k in xrange(gw, self.Gr.nzg-gw):
                if k * self.Gr.dz <= 1500:
                    EnvVar.QTvar.values[k]  = EnvVar.prescribed_QTvar
                else:
                    EnvVar.QTvar.values[k]  = 0.
                if k * self.Gr.dz <= 1500 and k * self.Gr.dz > 500:
                    EnvVar.Hvar.values[k]   = EnvVar.prescribed_Hvar
                else:
                    EnvVar.Hvar.values[k]   = 0.
                if k * self.Gr.dz <= 1500 and k * self.Gr.dz > 200:
                    EnvVar.HQTcov.values[k] = EnvVar.prescribed_HQTcov
                else:
                    EnvVar.HQTcov.values[k] = 0.

        for k in xrange(gw, self.Gr.nzg-gw):
            if True: #EnvVar.QTvar.values[k] != 0.0 and EnvVar.Hvar.values[k] != 0.0 and EnvVar.HQTcov.values[k] != 0.0:
                sd_q = sqrt(EnvVar.QTvar.values[k])
                sd_h = sqrt(EnvVar.Hvar.values[k])
                corr = fmax(fmin(EnvVar.HQTcov.values[k]/fmax(sd_h*sd_q, 1e-13),1.0),-1.0)

                # limit sd_q to prevent negative qt_hat
                sd_q_lim = (1e-10 - EnvVar.QT.values[k])/(sqrt2 * abscissas[0])
                sd_q = fmin(sd_q, sd_q_lim)
                qt_var = sd_q * sd_q
                sigma_h_star = sqrt(fmax(1.0-corr*corr,0.0)) * sd_h

                for key in outer_env:
                    outer_env[key] = 0.0
                if in_Env:
                    for key in outer_src:
                        outer_src[key] = 0.0

                for m_q in xrange(self.quadrature_order):
                    qt_hat    = EnvVar.QT.values[k] + sqrt2 * sd_q * abscissas[m_q]
                    mu_h_star = EnvVar.H.values[k]  + sqrt2 * corr * sd_h * abscissas[m_q]

                    for key in outer_env:
                        inner_env[key] = 0.0
                    if in_Env:
                        for key in outer_src:
                            inner_src[key] = 0.0

                    for m_h in xrange(self.quadrature_order):
                        h_hat = sqrt2 * sigma_h_star * abscissas[m_h] + mu_h_star

                        # condensation + autoconversion
                        sa  = eos(self.t_to_prog_fp, self.prog_to_t_fp, self.Ref.p0_half[k], qt_hat, h_hat)
                        mph = microphysics(sa.T, sa.ql, self.Ref.p0_half[k], qt_hat, h_hat, self.max_supersaturation, in_Env)

                        # environmental variables
                        inner_env['ql']    += mph.ql    * weights[m_h] * sqpi_inv
                        inner_env['qr']    += mph.qr    * weights[m_h] * sqpi_inv
                        inner_env['T']     += mph.T     * weights[m_h] * sqpi_inv
                        inner_env['thl']   += mph.H     * weights[m_h] * sqpi_inv
                        inner_env['alpha'] += mph.alpha * weights[m_h] * sqpi_inv
                        # cloudy/dry categories for buoyancy in TKE
                        if mph.ql  > 0.0:
                            inner_env['cf']     +=          weights[m_h] * sqpi_inv
                            inner_env['qt_cld'] += mph.qt * weights[m_h] * sqpi_inv
                            inner_env['T_cld']  += mph.T  * weights[m_h] * sqpi_inv
                        else:
                            inner_env['qt_dry'] += mph.qt * weights[m_h] * sqpi_inv
                            inner_env['T_dry']  += mph.T  * weights[m_h] * sqpi_inv
                        # products for variance and covariance source terms
                        if in_Env:
                            inner_src['Sqt']    += -mph.qr                     * weights[m_h] * sqpi_inv
                            inner_src['SH']     +=  mph.thl_rain_src           * weights[m_h] * sqpi_inv
                            inner_src['Sqt_H']  += -mph.qr           * mph.thl * weights[m_h] * sqpi_inv
                            inner_src['Sqt_qt'] += -mph.qr           * mph.qt  * weights[m_h] * sqpi_inv
                            inner_src['SH_H']   +=  mph.thl_rain_src * mph.thl * weights[m_h] * sqpi_inv
                            inner_src['SH_qt']  +=  mph.thl_rain_src * mph.qt  * weights[m_h] * sqpi_inv

                    for key in outer_env:
                        outer_env[key] += inner_env[key] * weights[m_q] * sqpi_inv

                    if in_Env:
                        for key in outer_src:
                            outer_src[key] += inner_src[key] * weights[m_q] * sqpi_inv

                # update environmental variables
                self.update_EnvVar(k, EnvVar, outer_env['T'], outer_env['thl'],\
                                   outer_env['qt_cld']+outer_env['qt_dry'], outer_env['ql'],\
                                   outer_env['qr'], outer_env['alpha'])
                # update cloudy/dry variables for buoyancy in TKE
                EnvVar.CF.values[k]  = outer_env['cf']
                self.qt_dry[k]    = outer_env['qt_dry']
                self.th_dry[k]    = outer_env['T_dry'] / exner_c(self.Ref.p0_half[k])
                self.t_cloudy[k]  = outer_env['T_cld']
                self.qv_cloudy[k] = outer_env['qt_cld'] - outer_env['ql']
                self.qt_cloudy[k] = outer_env['qt_cld']
                self.th_cloudy[k] = outer_env['T_cld'] / exner_c(self.Ref.p0_half[k])

                if in_Env:
                    self.Hvar_rain_dt[k]   = outer_src['SH_H']   - outer_src['SH']  * EnvVar.H.values[k]
                    self.QTvar_rain_dt[k]  = outer_src['Sqt_qt'] - outer_src['Sqt'] * EnvVar.QT.values[k]
                    self.HQTcov_rain_dt[k] = outer_src['SH_qt']  - outer_src['SH']  * EnvVar.QT.values[k] + \
                                             outer_src['Sqt_H']  - outer_src['Sqt'] * EnvVar.H.values[k]

            else:
                # condensation + autoconversion
                sa  = eos(self.t_to_prog_fp, self.prog_to_t_fp, self.Ref.p0_half[k], EnvVar.QT.values[k], EnvVar.H.values[k])
                mph = microphysics(sa.T, sa.ql, self.Ref.p0_half[k], EnvVar.QT.values[k], EnvVar.H.values[k], self.max_supersaturation, in_Env)

                self.update_EnvVar(   k, EnvVar, mph.T, mph.H, mph.qt, mph.ql, mph.qr, mph.alpha)
                self.update_cloud_dry(k, EnvVar, mph.T, mph.H, mph.qt, mph.ql, mph.qv)

                if in_Env:
                    self.Hvar_rain_dt[k]   = 0.
                    self.QTvar_rain_dt[k]  = 0.
                    self.HQTcov_rain_dt[k] = 0.

        return

    cdef void sommeria_deardorff(self, EnvironmentVariables EnvVar):
        # this function follows the derivation in
        # Sommeria and Deardorff 1977: Sub grid scale condensation in models of non-precipitating clouds.
        # J. Atmos. Sci., 34, 344-355.
        cdef:
            Py_ssize_t gw = self.Gr.gw
            double Lv, Tl, q_sl, beta1, lambda1, alpha1, sigma1, Q1, R, C0, C1, C2, C2_THL, qv

        if EnvVar.H.name == 'thetal':
            with nogil:
                for k in xrange(gw, self.Gr.nzg-gw):
                    Lv = latent_heat(EnvVar.T.values[k])
                    cp = cpd
                    # paper notation used below
                    Tl = EnvVar.H.values[k]*exner_c(self.Ref.p0_half[k])
                    q_sl = qv_star_t(self.Ref.p0[k], Tl) # using the qv_star_c function instead of the approximation in eq. (4) in SD
                    beta1 = 0.622*Lv**2/(Rd*cp*Tl**2) # eq. (8) in SD
                    #q_s = q_sl*(1+beta1*EnvVar.QT.values[k])/(1+beta1*q_sl) # eq. (7) in SD
                    lambda1 = 1/(1+beta1*q_sl) # text under eq. (20) in SD
                    # check the pressure units - mb vs pa
                    alpha1 = (self.Ref.p0[k]/100000.0)**0.286*0.622*Lv*q_sl/Rd/Tl**2 # eq. (14) and eq. (6) in SD
                    # see if there is another way to calculate dq/dT from scmapy
                    sigma1 = EnvVar.QTvar.values[k]-2*alpha1*EnvVar.HQTcov.values[k]+alpha1**2*EnvVar.Hvar.values[k] # eq. (18) in SD , with r from (11)
                    Q1 = (EnvVar.QT.values[k]-q_sl)/sigma1 # eq. (17) in SD
                    R = 0.5*(1+erf(Q1/sqrt(2.0))) # approximation in eq. (16) in SD
                    #R1 = 0.5*(1+Q1/1.6) # approximation in eq. (22) in SD
                    C0 = 1.0+0.61*q_sl-alpha1*lambda1*EnvVar.THL.values[k]*(Lv/cp/Tl*(1.0+0.61*q_sl)-1.61) # eq. (37) in SD
                    C1 = (1.0-R)*(1+0.61*q_sl)+R*C0 # eq. (42a) in SD
                    C2 = (1.0-R)*0.61+R*(C0*Lv/cp/Tl-1.0) # eq. (42b) in SD
                    C2_THL = C2*EnvVar.THL.values[k] # defacto the coefficient in eq(41) is C2*THL
                    # the THVvar is given as a function of THVTHLcov and THVQTcov from eq. (41) in SD.
                    # these covariances with THL are obtained by substituting w for THL or QT in eq. (41),
                    # i.e. applying eq. (41) twice. The resulting expression yields: C1^2*THL_var+2*C1*C2*THL_var*QT_var+C2^2**QT_var
                    EnvVar.THVvar.values[k] = C1**2*EnvVar.Hvar.values[k] + 2*C1*C2_THL*EnvVar.HQTcov.values[k]+ C2_THL**2*EnvVar.QTvar.values[k]
                    # equation (19) exact form for QL
                    EnvVar.QL.values[k] = 1.0/(1.0+beta1*q_sl)*(R*(EnvVar.QT.values[k]-q_sl)+sigma1/sqrt(6.14)*exp(-((EnvVar.QT.values[k]-q_sl)*(EnvVar.QT.values[k]-q_sl)/(2.0*sigma1*sigma1))))
                    EnvVar.T.values[k] = Tl + Lv/cp*EnvVar.QL.values[k] # should this be the differnece in ql - would it work for evaporation as well ?
                    EnvVar.CF.values[k] = R
                    qv = EnvVar.QT.values[k] - EnvVar.QL.values[k]
                    alpha = alpha_c(self.Ref.p0_half[k], EnvVar.T.values[k], EnvVar.QT.values[k], qv)
                    EnvVar.B.values[k] = buoyancy_c(self.Ref.alpha0_half[k], alpha)
                    EnvVar.THL.values[k] = t_to_thetali_c(self.Ref.p0_half[k], EnvVar.T.values[k], EnvVar.QT.values[k],
                                                          EnvVar.QL.values[k], 0.0)

                    self.qt_dry[k] = EnvVar.QT.values[k]
                    self.th_dry[k] = EnvVar.T.values[k]/exner_c(self.Ref.p0_half[k])
                    self.t_cloudy[k] = EnvVar.T.values[k]
                    self.qv_cloudy[k] = EnvVar.QT.values[k] - EnvVar.QL.values[k]
                    self.qt_cloudy[k] = EnvVar.QT.values[k]
                    self.th_cloudy[k] = EnvVar.T.values[k]/exner_c(self.Ref.p0_half[k])

                    #using the approximation in eq. (25) in SD, noting that in the paper there is a typo in the first
                    # condition and 1.6 there should be -1.6
                    # if Q1<-1.6:
                    #     EnvVar.QL.values[k] = 0.0*lambda1*sigma1
                    # elif Q1>-1.6 and Q1<1.6:
                    #     EnvVar.QL.values[k] = ((Q1+1.6)**2/6.4)*lambda1*sigma1
                    # elif Q1>1.6:
                    #     EnvVar.QL.values[k] = Q1*lambda1*sigma1
                    #

        elif EnvVar.H.name == 's':
            sys.exit('EDMF_Environment: Sommeria Deardorff is not defined for using entropy as thermodyanmic variable')
        return

    cpdef satadjust(self, EnvironmentVariables EnvVar, bint in_Env):#, TimeStepping TS):

        if EnvVar.EnvThermo_scheme == 'sa_mean':
            self.eos_update_SA_mean(EnvVar, in_Env)
        elif EnvVar.EnvThermo_scheme == 'sa_quadrature':
            self.eos_update_SA_sgs(EnvVar, in_Env)#, TS)
        elif EnvVar.EnvThermo_scheme == 'sommeria_deardorff':
            self.sommeria_deardorff(EnvVar)
        else:
            sys.exit('EDMF_Environment: Unrecognized EnvThermo_scheme. Possible options: sa_mean, sa_quadrature, sommeria_deardorff')

        return
