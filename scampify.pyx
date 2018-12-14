import numpy as np
from netCDF4 import Dataset

from Grid cimport Grid
from ReferenceState cimport ReferenceState
cimport TimeStepping

cimport EDMF_Updrafts
cimport EDMF_Environment
cimport EDMF_Rain
from Variables cimport GridMeanVariables

from thermodynamic_functions cimport  *
from microphysics_functions cimport  *

import matplotlib.pyplot as plt
import pylab
from matplotlib import ticker

import cython

# helper functions for reading LES netcdf files
def read_profiles(ncdf_dataset, var_list):
    p_dict = {}

    for var in var_list:
        p_dict[var] = np.array(LES_data["profiles/" + var][:, :])

    p_dict['z_half'] = np.array(LES_data["profiles/z_half"][:])
    p_dict['z']      = np.array(LES_data["profiles/z"][:])
    p_dict['time']   = np.array(LES_data["profiles/t"][:])

    return p_dict

def read_reference(ncdf_dataset, ref_list):
    r_dict = {}

    for var in ref_list:
        r_dict[var] = np.array(LES_data["reference/" + var][:])

    return r_dict

def helper_var_reader(name1, name2, name3):
    var1 = np.array(LES_data["profiles/" + name1][:, :])
    var2 = np.array(LES_data["profiles/" + name2][:, :])
    var3 = np.array(LES_data["profiles/" + name3][:, :])

    return var1 - var2 * var3

def read_varcovar(ncdf_dataset):

    v_dict = {}

    v_dict['Hvar']   = helper_var_reader("env_thetali2", "env_thetali", "env_thetali")
    v_dict['QTvar']  = helper_var_reader("env_qt2", "env_qt", "env_qt")
    v_dict['HQTcov'] = helper_var_reader("env_qt_thetali", "env_thetali", "env_qt")

    return v_dict

# user-defined parameters for SCM (TODO - read them from file)
namelist = {}
namelist['turbulence'] = {}
namelist['turbulence']['scheme'] = 'EDMF_PrognosticTKE'
namelist['turbulence']['EDMF_PrognosticTKE'] = {}
namelist['turbulence']['EDMF_PrognosticTKE']['calc_scalar_var']        = True
namelist['turbulence']['EDMF_PrognosticTKE']['updraft_number']         = 1     #TODO - 5
namelist['turbulence']['EDMF_PrognosticTKE']['entrainment']            = "suselj"
namelist['turbulence']['EDMF_PrognosticTKE']['mixing_length']          = "sbl"
namelist['turbulence']['EDMF_PrognosticTKE']["use_sommeria_deardorff"] = False
namelist['turbulence']['EDMF_PrognosticTKE']['use_local_micro']        = True

namelist['thermodynamics'] = {}
namelist['thermodynamics']['thermal_variable'] = 'thetal'
namelist['thermodynamics']['saturation'] = {}
namelist['thermodynamics']['saturation'] = 'sa_quadrature'

namelist['microphysics'] = {}
namelist['microphysics']['rain_model']          = True
namelist['microphysics']['rain_const_area']     = True
namelist['microphysics']['max_supersaturation'] = 0.0001 #0.1 # 1e-4

namelist['grid'] = {}
namelist['grid']['dims'] = 1
namelist['grid']['dz'] = 40. # 40  in SCM 40 in LES
namelist['grid']['gw'] = 0   #  2  in SCM  7 in LES
namelist['grid']['nz'] = 150 # 100 in SCM 150 in LES

namelist['time_stepping'] = {}
namelist['time_stepping']['t_max'] = 86400.
namelist['time_stepping']['dt'] = 1.  # TODO - interpolate between data points

paramlist = {}
paramlist['turbulence'] = {}
paramlist['turbulence']['EDMF_PrognosticTKE'] = {}
paramlist['turbulence']['EDMF_PrognosticTKE']['surface_area']       = 0.1  # 0.2  # 0.1
paramlist['turbulence']['EDMF_PrognosticTKE']['max_area_factor']    = 10.  # 1.0  # 9.9
paramlist['turbulence']['EDMF_PrognosticTKE']['tke_diss_coeff']     = 1.0  # 0.22 # 2.
paramlist['turbulence']['EDMF_PrognosticTKE']['tke_ed_coeff']       = 0.1  # 0.17 # 0.1
paramlist['turbulence']['EDMF_PrognosticTKE']["detrainment_factor"] = 0.75  # 1.
paramlist['turbulence']['EDMF_PrognosticTKE']["entrainment_factor"] = 0.75  # 1.

time = namelist['time_stepping']['t_max']

LES_data = Dataset('../SCAMPY_tests/plots/RICO_comparison/Stats_Rico_pycles_Pr_Ar.nc')

var_list = ['ql_mean',         'env_ql',         'updraft_ql',\
            'qt_mean',         'env_qt',         'updraft_qt',\
            'thetali_mean',    'env_thetali',    'updraft_thetali',\
            'temperature_mean','env_temperature','updraft_temperature',\
            'qrain_mean',\
                               'env_fraction',   'updraft_fraction',\
            'w_mean']

ref_list = ['p0_full', 'p0', 'alpha0_full', 'alpha0', 'rho0_full', 'rho0']

p_dict = read_profiles(LES_data, var_list)
v_dict = read_varcovar(LES_data)
r_dict = read_reference(LES_data, ref_list)

# TODO - interpolate in space and time between scampy and pycles?
#print p_dict['time'][0], p_dict['time'][1], p_dict['time'][-1]
#print p_dict['time']

cdef class Scampify1d:

    def __init__(self):
        self.Gr  = Grid(namelist)
        self.Ref = ReferenceState(self.Gr)
        self.TS  = TimeStepping.TimeStepping(namelist)
        self.GMV = GridMeanVariables(namelist, self.Gr, self.Ref)

        self.Ref.p0          = r_dict['p0_full']
        self.Ref.p0_half     = r_dict['p0']
        self.Ref.alpha0      = r_dict['alpha0_full']
        self.Ref.alpha0_half = r_dict['alpha0']
        self.Ref.rho0        = r_dict['rho0_full']
        self.Ref.rho0_half   = r_dict['rho0']

        # init GMV
        self.GMV.W.values   = p_dict["w_mean"][0,:]
        self.GMV.QT.values  = p_dict["qt_mean"][0,:]
        self.GMV.THL.values = p_dict["thetali_mean"][0,:]
        self.GMV.T.values   = p_dict["temperature_mean"][0,:]

        # init Rain
        self.rain_var = EDMF_Rain.RainVariables(namelist, self.Gr)
        self.rain_phs = EDMF_Rain.RainPhysics(self.Gr, self.Ref)

        # init Environment
        self.env_var = EDMF_Environment.EnvironmentVariables(namelist, self.Gr)
        self.env_var.QT.values     = p_dict['env_qt'][0,:]
        self.env_var.H.values      = p_dict['env_thetali'][0,:]
        self.env_var.EnvArea.values= p_dict['env_fraction'][0,:]
        self.env_var.Hvar.values   = v_dict['Hvar'][0,:]
        self.env_var.QTvar.values  = v_dict['QTvar'][0,:]
        self.env_var.HQTcov.values = v_dict['HQTcov'][0,:]
        self.env_thr  = EDMF_Environment.EnvironmentThermodynamics(namelist, self.Gr, self.Ref, self.env_var, self.rain_var)

        # init Updrafts
        self.upd_var = EDMF_Updrafts.UpdraftVariables(1, namelist, paramlist, self.Gr)
        # cant assign nupy array to memory view slice
        # https://mail.python.org/pipermail/cython-devel/2014-December/004288.html
        cdef double[:] tmp_qt       = p_dict['updraft_qt'][0,:]
        cdef double[:] tmp_thetali  = p_dict['updraft_thetali'][0,:]
        cdef double[:] tmp_area     = p_dict['updraft_fraction'][0,:]
        self.upd_var.QT.values[0,:]   = tmp_qt
        self.upd_var.H.values[0,:]    = tmp_thetali
        self.upd_var.Area.values[0,:] = tmp_area
        self.upd_thr  = EDMF_Updrafts.UpdraftThermodynamics(1, self.Gr, self.Ref, self.upd_var, self.rain_var)

    def do_environment(self):

        if namelist['thermodynamics']['saturation'] == 'sa_quadrature':
            self.env_thr.eos_update_SA_sgs(self.env_var, True)
        else:
            self.env_thr.eos_update_SA_mean(self.env_var, True)

    cpdef do_updrafts(self):
        cdef:
            Py_ssize_t k
            eos_struct sa

        for k in xrange(0, self.Gr.nzg):
            if self.upd_var.Area.values[0,k] != 0:
                sa = eos(
                    self.upd_thr.t_to_prog_fp,
                    self.upd_thr.prog_to_t_fp,
                    self.Ref.p0_half[k],
                    self.upd_var.QT.values[0,k],
                    self.upd_var.H.values[0,k]
                )
                self.upd_var.QL.values[0,k] = sa.ql
                self.upd_var.T.values[0,k] = sa.T

                # TODO
                #self.upd_thr.compute_update_combined_local_thetal(
                #    self.Ref.p0_half[k],
                #    self.upd_var.T.values[0,k],
                #    &self.upd_var.QT.values[0,k],
                #    &self.upd_var.QL.values[0,k],
                #    &self.upd_var.QR.values[0,k],
                #    &self.upd_var.H.values[0,k],
                #    0,
                #    k
                #)

    #TODO - do rain

    def run(self):

        self.TS.t = 0.

        while self.TS.t <= self.TS.t_max:
            self.read_in_LES_data()
            self.do_updrafts()
            #self.do_environment()
            #self.do_rain()
            #self.update_mean()
            #self.output()
            self.TS.t += self.TS.dt

            print self.TS.t
            print np.min(self.upd_var.QL.values[0,:]), np.max(self.upd_var.QL.values[0,:])
            print np.min(self.env_var.QL.values[:]), np.max(self.env_var.QL.values[:])
            print "-----------------------------------------------------"

        return

    #def plot_all(self):

    #    env_ql = np.array(self.env_var.QL.values[:])
    #    upd_ql = np.array(self.upd_var.QL.values[0, :])
    #    area   = np.array(self.upd_var.Area.values[0, :])
    #    ql_mean = (area * upd_ql + (1.-area) * env_ql) * 1e3

    #    fig, ax_arr = plt.subplots(3, 3)
    #    fig.set_figheight(10)
    #    fig.set_figwidth(12)

    #    ax_arr[0, 0].plot(p_dict['env_thetali'][time,:], p_dict['z_half'], c='darkorange', label='env', lw=3)
    #    ax_arr[0, 0].plot(p_dict['updraft_thetali'][time,:], p_dict['z_half'], c='blue', label='upd', lw=3)
    #    ax_arr[0, 0].set(xlabel="thetali [K]", ylabel='height [m]')
    #    ax_arr[0, 0].set_xlim(285, 310)
    #    ax_arr[0, 0].legend(frameon=False)
    #    ax_arr[0, 0].grid(True)

    #    ax_arr[0, 1].plot(p_dict['env_qt'][time,:]*1e3 , p_dict['z_half'], c='darkorange', label='env', lw=3)
    #    ax_arr[0, 1].plot(p_dict['updraft_qt'][time,:]*1e3 , p_dict['z_half'], c='blue', label='upd', lw=3)
    #    ax_arr[0, 1].set(xlabel="qt [g/kg]")
    #    ax_arr[0, 1].tick_params(labelleft='off')
    #    ax_arr[0, 1].legend(frameon=False)
    #    ax_arr[0, 1].grid(True)

    #    ax_arr[0, 2].plot(p_dict['updraft_fraction'][time,:]*100, p_dict['z_half'], c='blue', label='updraft frac.', lw=3)
    #    ax_arr[0, 2].set(xlabel="upd fraction [%]")
    #    ax_arr[0, 2].tick_params(labelleft='off')
    #    ax_arr[0, 2].legend(frameon=False)
    #    ax_arr[0, 2].grid(True)

    #    ax_arr[1, 0].plot(v_dict['Hvar'][time,:], p_dict['z_half'], c='darkorange', lw=3)
    #    ax_arr[1, 0].set(xlabel="Hvar", ylabel='height [m]')
    #    ax_arr[1, 0].set_xlim(-0.2, 4.)
    #    ax_arr[1, 0].legend(frameon=False)
    #    ax_arr[1, 0].grid(True)

    #    ax_arr[1, 1].plot(v_dict['QTvar'][time,:], p_dict['z_half'], c='darkorange', lw=3)
    #    ax_arr[1, 1].set(xlabel="QTvar")
    #    ax_arr[1, 1].tick_params(labelleft='off')
    #    ax_arr[1, 1].set_xlim(-1e-7, 2.4e-6)
    #    ax_arr[1, 1].xaxis.set_major_locator(ticker.MultipleLocator(1.2e-6))
    #    ax_arr[1, 1].legend(frameon=False)
    #    ax_arr[1, 1].grid(True)

    #    ax_arr[1, 2].plot(v_dict['HQTcov'][time,:], p_dict['z_half'], c='darkorange', lw=3)
    #    ax_arr[1, 2].set(xlabel="HQTcov")
    #    ax_arr[1, 2].set_xlim(-3e-3, 1e-4)
    #    ax_arr[1, 2].xaxis.set_major_locator(ticker.MultipleLocator(1e-3))
    #    ax_arr[1, 2].tick_params(labelleft='off')
    #    ax_arr[1, 2].legend(frameon=False)
    #    ax_arr[1, 2].grid(True)

    #    ax_arr[2, 0].plot(np.array(self.upd_var.QL.values[0,:])*1e3, p_dict['z_half'], label='mock SCM', c='limegreen', lw=3)
    #    ax_arr[2, 0].plot(p_dict['updraft_ql'][time,:]*1e3, p_dict['z_half'], label='LES', c='blue', lw=3)
    #    ax_arr[2, 0].set(xlabel="upd ql [g/kg]", ylabel='height [m]')
    #    ax_arr[2, 0].set_xlim(-1e-1, .6)
    #    ax_arr[2, 0].legend(frameon=False)
    #    ax_arr[2, 0].grid(True)

    #    ax_arr[2, 1].plot(np.array(self.env_var.QL.values[:])*1e3, p_dict['z_half'], label='mock SCM', c='limegreen', lw=3)
    #    ax_arr[2, 1].plot(p_dict['env_ql'][time,:]*1e3, p_dict['z_half'], label='LES', c='darkorange', lw=3)
    #    ax_arr[2, 1].set(xlabel="env ql [g/kg]")
    #    ax_arr[2, 1].set_xlim(-1e-1, .6)
    #    ax_arr[2, 1].tick_params(labelleft='off')
    #    ax_arr[2, 1].legend(frameon=False)
    #    ax_arr[2, 1].grid(True)

    #    ax_arr[2, 2].plot(ql_mean, p_dict['z_half'], label='mock SCM', c='limegreen', lw=3)
    #    ax_arr[2, 2].plot(p_dict['ql_mean'][time,:]*1e3, p_dict['z_half'], label='LES', c='black', lw=3)
    #    ax_arr[2, 2].set_xlim(-1e-1, .6)
    #    ax_arr[2, 2].set(xlabel="mean ql [g/kg]")
    #    ax_arr[2, 2].tick_params(labelleft='off')
    #    ax_arr[2, 2].legend(frameon=False)
    #    ax_arr[2, 2].grid(True)

    #    plt.savefig("tests/scampify/output/Rico_scoff.pdf")

