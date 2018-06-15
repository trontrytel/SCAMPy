import numpy as np
from netCDF4 import Dataset

from Grid cimport Grid
from ReferenceState cimport ReferenceState
cimport EDMF_Updrafts
cimport EDMF_Environment

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
    p_dict['time']   = np.array(LES_data["profiles/t"][:]) / 60. / 60.

    return p_dict

def read_reference(ncdf_dataset):
    r_dict = {}
    if casename == 'Bomex':
        var_list = ['p0',      'p0_half', 'alpha0',      'alpha0_half', 'rho0',      'rho0_half']
    elif casename == 'DYCOMS_RF01' or casename == 'Bomex_100_20':
        var_list = ['p0_full', 'p0',      'alpha0_full', 'alpha0',      'rho0_full', 'rho0']
    else:
        print "wrong casename", casename
        assert(False)

    for var in var_list:
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

# user-defined parameters for SCM
# TODO - read them from netcdf LES files
namelist = {}
namelist['turbulence'] = {}
namelist['turbulence']['scheme'] = 'EDMF_PrognosticTKE'
namelist['turbulence']['EDMF_PrognosticTKE'] = {}
namelist['turbulence']['EDMF_PrognosticTKE']['use_scalar_var'] = {}
namelist['turbulence']['EDMF_PrognosticTKE']['use_scalar_var'] = True
namelist['thermodynamics'] = {}
namelist['thermodynamics']['thermal_variable'] = 'thetal'
namelist['thermodynamics']['saturation'] = {}
namelist['thermodynamics']['saturation'] = 'sa_quadrature'
namelist['condensation'] = {}
namelist['condensation']['quadrature_order'] = 100

paramlist = {}
paramlist['turbulence'] = {}
paramlist['turbulence']['EDMF_PrognosticTKE'] = {}
paramlist['turbulence']['EDMF_PrognosticTKE']['surface_area'] = 0.1
paramlist['turbulence']['updraft_microphysics'] = {}
paramlist['turbulence']['updraft_microphysics']['max_supersaturation'] = 44.#0.01

namelist_gr = {}
namelist_gr['grid'] = {}
casename = 'Bomex'
#casename = 'Bomex_100_20'
#casename = 'DYCOMS_RF01'
time = 200

if casename == 'Bomex' or casename == 'Bomex_100_20':
    namelist_gr['grid']['dz'] = 20. # 40 in SCM 20 in LES
    namelist_gr['grid']['gw'] = 0   #  2 in SCM  7 in LES
    namelist_gr['grid']['nz'] = 150 # 75 in SCM 75 in LES
elif casename == 'DYCOMS_RF01':
    namelist_gr['grid']['dz'] = 5.
    namelist_gr['grid']['gw'] = 0
    namelist_gr['grid']['nz'] = 300
else:
    print "not a valid casename"
    assert(False)

LES_data = Dataset('tests/scampify/LES_data/Stats.'+casename+'.nc')
var_list = ['ql_mean', 'env_ql', 'updraft_ql', 'env_qt', 'env_thetali', 'updraft_qt',\
            'updraft_thetali', 'updraft_fraction', 'env_qt_thetali']
p_dict = read_profiles(LES_data, var_list)
v_dict = read_varcovar(LES_data)
r_dict = read_reference(LES_data)

cdef class Scampify1d:

    def __init__(self):
        self.Gr = Grid(namelist_gr)
        self.ref_state = ReferenceState(self.Gr)

    def initialize(self):
        if casename == 'Bomex':
            self.ref_state.p0          = r_dict['p0']
            self.ref_state.p0_half     = r_dict['p0_half']
            self.ref_state.alpha0      = r_dict['alpha0']
            self.ref_state.alpha0_half = r_dict['alpha0_half']
            self.ref_state.rho0        = r_dict['rho0']
            self.ref_state.rho0_half   = r_dict['rho0_half']
        elif casename == 'DYCOMS_RF01' or casename == 'Bomex_100_20':
            self.ref_state.p0          = r_dict['p0_full']
            self.ref_state.p0_half     = r_dict['p0']
            self.ref_state.alpha0      = r_dict['alpha0_full']
            self.ref_state.alpha0_half = r_dict['alpha0']
            self.ref_state.rho0        = r_dict['rho0_full']
            self.ref_state.rho0_half   = r_dict['rho0']
        else:
            print "wrong casename"
            assert(False)


        self.env_var = EDMF_Environment.EnvironmentVariables(namelist, self.Gr)
        self.env_var.QT.values = p_dict['env_qt'][time,:]
        self.env_var.H.values  = p_dict['env_thetali'][time,:]
        self.env_var.Hvar.values   = v_dict['Hvar'][time,:]
        self.env_var.QTvar.values  = v_dict['QTvar'][time,:]
        self.env_var.HQTcov.values = v_dict['HQTcov'][time,:]

        self.upd_var = EDMF_Updrafts.UpdraftVariables(1, namelist, paramlist, self.Gr)
        # cant assign nupy array to memory view slice
        # https://mail.python.org/pipermail/cython-devel/2014-December/004288.html
        cdef double[:] tmp_qt       = p_dict['updraft_qt'][time,:]
        cdef double[:] tmp_thetali  = p_dict['updraft_thetali'][time,:]
        cdef double[:] tmp_area     = p_dict['updraft_fraction'][time,:]
        self.upd_var.QT.values[0,:]   = tmp_qt
        self.upd_var.H.values[0,:]    = tmp_thetali
        self.upd_var.Area.values[0,:] = tmp_area
        self.upd_var.QR.values[0,:]   = 0.0   #TODO

        self.env_thr = EDMF_Environment.EnvironmentThermodynamics(namelist, paramlist, self.Gr, self.ref_state, self.env_var)
        self.upd_thr = EDMF_Updrafts.UpdraftThermodynamics(1, self.Gr, self.ref_state, self.upd_var)
        self.upd_mcr = EDMF_Updrafts.UpdraftMicrophysics(paramlist, 1, self.Gr, self.ref_state)

    def do_environment(self):

        #self.env_thr.eos_update_SA_mean(self.env_var, True)
        self.env_thr.eos_update_SA_sgs(self.env_var, True)

    cpdef do_updrafts(self):
        cdef:
            Py_ssize_t k
            eos_struct sa

        for k in xrange(0, self.Gr.nzg):
            if self.upd_var.Area.values[0,k] != 0:
                sa = eos(
                    self.upd_thr.t_to_prog_fp,
                    self.upd_thr.prog_to_t_fp,
                    self.ref_state.p0_half[k],
                    self.upd_var.QT.values[0,k],
                    self.upd_var.H.values[0,k]
                )
                self.upd_var.QL.values[0,k] = sa.ql
                self.upd_var.T.values[0,k] = sa.T

                self.upd_mcr.compute_update_combined_local_thetal(
                    self.ref_state.p0_half[k],
                    self.upd_var.T.values[0,k],
                    &self.upd_var.QT.values[0,k],
                    &self.upd_var.QL.values[0,k],
                    &self.upd_var.QR.values[0,k],
                    &self.upd_var.H.values[0,k],
                    0,
                    k
                )

    def plot_all(self):

        env_ql = np.array(self.env_var.QL.values[:])
        upd_ql = np.array(self.upd_var.QL.values[0, :])
        area   = np.array(self.upd_var.Area.values[0, :])
        ql_mean = (area * upd_ql + (1.-area) * env_ql) * 1e3

        fig, ax_arr = plt.subplots(3, 3)
        fig.set_figheight(10)
        fig.set_figwidth(12)

        ax_arr[0, 0].plot(p_dict['env_thetali'][time,:], p_dict['z_half'], c='darkorange', label='env', lw=3)
        ax_arr[0, 0].plot(p_dict['updraft_thetali'][time,:], p_dict['z_half'], c='blue', label='upd', lw=3)
        ax_arr[0, 0].set(xlabel="thetali [K]", ylabel='height [m]')
        ax_arr[0, 0].set_xlim(285, 310)
        ax_arr[0, 0].legend(frameon=False)
        ax_arr[0, 0].grid(True)

        ax_arr[0, 1].plot(p_dict['env_qt'][time,:]*1e3 , p_dict['z_half'], c='darkorange', label='env', lw=3)
        ax_arr[0, 1].plot(p_dict['updraft_qt'][time,:]*1e3 , p_dict['z_half'], c='blue', label='upd', lw=3)
        ax_arr[0, 1].set(xlabel="qt [g/kg]")
        ax_arr[0, 1].tick_params(labelleft='off')
        ax_arr[0, 1].legend(frameon=False)
        ax_arr[0, 1].grid(True)

        ax_arr[0, 2].plot(p_dict['updraft_fraction'][time,:]*100, p_dict['z_half'], c='blue', label='updraft frac.', lw=3)
        ax_arr[0, 2].set(xlabel="upd fraction [%]")
        ax_arr[0, 2].tick_params(labelleft='off')
        ax_arr[0, 2].legend(frameon=False)
        ax_arr[0, 2].grid(True)

        ax_arr[1, 0].plot(v_dict['Hvar'][time,:], p_dict['z_half'], c='darkorange', lw=3)
        ax_arr[1, 0].set(xlabel="Hvar", ylabel='height [m]')
        ax_arr[1, 0].set_xlim(-0.01, 0.22)
        ax_arr[1, 0].legend(frameon=False)
        ax_arr[1, 0].grid(True)

        ax_arr[1, 1].plot(v_dict['QTvar'][time,:], p_dict['z_half'], c='darkorange', lw=3)
        ax_arr[1, 1].set(xlabel="QTvar")
        ax_arr[1, 1].tick_params(labelleft='off')
        ax_arr[1, 1].set_xlim(-1e-8, 4*1e-7)
        ax_arr[1, 1].xaxis.set_major_locator(ticker.MultipleLocator(1e-7))
        #ax_arr[1, 1].xaxis.set_major_locator(ticker.MultipleLocator(1e-6))
        ax_arr[1, 1].legend(frameon=False)
        ax_arr[1, 1].grid(True)

        ax_arr[1, 2].plot(v_dict['HQTcov'][time,:], p_dict['z_half'], c='darkorange', lw=3)
        ax_arr[1, 2].set(xlabel="HQTcov")
        ax_arr[1, 2].set_xlim(-3e-4, 1e-5)
        #ax_arr[1, 2].xaxis.set_major_locator(ticker.MultipleLocator(1e-3))
        ax_arr[1, 2].xaxis.set_major_locator(ticker.MultipleLocator(1e-4))
        ax_arr[1, 2].tick_params(labelleft='off')
        ax_arr[1, 2].legend(frameon=False)
        ax_arr[1, 2].grid(True)

        ax_arr[2, 0].plot(np.array(self.upd_var.QL.values[0,:])*1e3, p_dict['z_half'], label='mock SCM', c='limegreen', lw=3)
        ax_arr[2, 0].plot(p_dict['updraft_ql'][time,:]*1e3, p_dict['z_half'], label='LES', c='blue', lw=3)
        ax_arr[2, 0].set(xlabel="upd ql [g/kg]", ylabel='height [m]')
        ax_arr[2, 0].set_xlim(-1e-1, 1.8)
        ax_arr[2, 0].legend(frameon=False)
        ax_arr[2, 0].grid(True)

        ax_arr[2, 1].plot(np.array(self.env_var.QL.values[:])*1e3, p_dict['z_half'], label='mock SCM', c='limegreen', lw=3)
        ax_arr[2, 1].plot(p_dict['env_ql'][time,:]*1e3, p_dict['z_half'], label='LES', c='darkorange', lw=3)
        ax_arr[2, 1].set(xlabel="env ql [g/kg]")
        ax_arr[2, 1].set_xlim(-1e-4, 8e-3)
        ax_arr[2, 1].tick_params(labelleft='off')
        ax_arr[2, 1].legend(frameon=False)
        ax_arr[2, 1].grid(True)

        ax_arr[2, 2].plot(ql_mean, p_dict['z_half'], label='mock SCM', c='limegreen', lw=3)
        ax_arr[2, 2].plot(p_dict['ql_mean'][time,:]*1e3, p_dict['z_half'], label='LES', c='black', lw=3)
        ax_arr[2, 2].set_xlim(-1e-3, 0.015)
        ax_arr[2, 2].set(xlabel="mean ql [g/kg]")
        ax_arr[2, 2].tick_params(labelleft='off')
        ax_arr[2, 2].legend(frameon=False)
        ax_arr[2, 2].grid(True)

        plt.savefig("tests/scampify/output/Bomex_mock_SCM_sgs_quad_100.pdf")

