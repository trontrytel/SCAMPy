import numpy as np
from netCDF4 import Dataset

from Grid cimport Grid
from ReferenceState cimport ReferenceState
cimport TimeStepping

from Cases import CasesFactory
from Cases import CasesBase

cimport EDMF_Updrafts
cimport EDMF_Environment
cimport EDMF_Rain
from Variables cimport GridMeanVariables
from NetCDFIO cimport NetCDFIO_Stats

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

namelist["meta"] = {}
namelist["meta"]["casename"] = "Rico"
namelist["meta"]["simname"] = "Rico"
namelist["meta"]["uuid"] = "6666"

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
namelist['grid']['gw'] = 2   #  2  in SCM  7 in LES
namelist['grid']['nz'] = 150 # 100 in SCM 150 in LES

namelist['time_stepping'] = {}
namelist['time_stepping']['t_max'] = 86400.
namelist['time_stepping']['dt'] = 60.  # TODO - interpolate between data points

namelist["output"] = {}
namelist["output"]["output_root"] = "./"

namelist["stats_io"] = {}
namelist["stats_io"]["frequency"] = 1.
namelist["stats_io"]["stats_dir"] = "stats"

paramlist = {}

paramlist["meta"] = {}
paramlist["meta"]["casename"] = "Rico"

paramlist['turbulence'] = {}
paramlist['turbulence']['EDMF_PrognosticTKE'] = {}
paramlist['turbulence']['EDMF_PrognosticTKE']['surface_area']       = 0.1  # 0.2  # 0.1
paramlist['turbulence']['EDMF_PrognosticTKE']['max_area_factor']    = 10.  # 1.0  # 9.9
paramlist['turbulence']['EDMF_PrognosticTKE']['tke_diss_coeff']     = 1.0  # 0.22 # 2.
paramlist['turbulence']['EDMF_PrognosticTKE']['tke_ed_coeff']       = 0.1  # 0.17 # 0.1
paramlist['turbulence']['EDMF_PrognosticTKE']["detrainment_factor"] = 0.75  # 1.
paramlist['turbulence']['EDMF_PrognosticTKE']["entrainment_factor"] = 0.75  # 1.
paramlist['turbulence']['EDMF_PrognosticTKE']["pressure_buoy_coeff"] = 0.3333333333333333
paramlist['turbulence']['EDMF_PrognosticTKE']["pressure_drag_coeff"] = 0.375
paramlist['turbulence']['EDMF_PrognosticTKE']["pressure_plume_spacing"] = 500.0

paramlist['turbulence']["Ri_bulk_crit"] = 0.0
paramlist['turbulence']["prandtl_number"] = 1.0

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
        self.it  = 0

        self.Gr    = Grid(namelist)
        self.Ref   = ReferenceState(self.Gr)
        self.TS    = TimeStepping.TimeStepping(namelist)
        self.GMV   = GridMeanVariables(namelist, self.Gr, self.Ref)
        self.Stats = NetCDFIO_Stats(namelist, paramlist, self.Gr)
        self.Case  = CasesFactory(namelist, paramlist)

        #TODO - check if needed
        self.Ref.p0           = np.zeros(self.Gr.nzg)
        self.Ref.p0_half      = np.zeros(self.Gr.nzg)
        self.Ref.alpha0       = np.zeros(self.Gr.nzg)
        self.Ref.alpha0_half  = np.zeros(self.Gr.nzg)
        self.Ref.rho0         = np.zeros(self.Gr.nzg)
        self.Ref.rho0_half    = np.zeros(self.Gr.nzg)

        for idx in range(self.Gr.gw, self.Gr.nzg - self.Gr.gw):
            self.Ref.p0[idx]          = r_dict['p0_full'][idx - self.Gr.gw]
            self.Ref.p0_half[idx]     = r_dict['p0'][idx - self.Gr.gw]
            self.Ref.alpha0[idx]      = r_dict['alpha0_full'][idx - self.Gr.gw]
            self.Ref.alpha0_half[idx] = r_dict['alpha0'][idx - self.Gr.gw]
            self.Ref.rho0[idx]        = r_dict['rho0_full'][idx - self.Gr.gw]
            self.Ref.rho0_half[idx]   = r_dict['rho0'][idx - self.Gr.gw]

        # init GMV
        for idx in range(self.Gr.gw, self.Gr.nzg - self.Gr.gw):
            self.GMV.W.values[idx]    = p_dict["w_mean"][0,idx - self.Gr.gw]
            self.GMV.QT.values[idx]   = p_dict["qt_mean"][0,idx - self.Gr.gw]
            self.GMV.THL.values[idx]  = p_dict["thetali_mean"][0,idx - self.Gr.gw]
            self.GMV.T.values[idx]    = p_dict["temperature_mean"][0,idx - self.Gr.gw]

        # init Rain
        self.rain_var  = EDMF_Rain.RainVariables(namelist, self.Gr)
        self.rain_phs  = EDMF_Rain.RainPhysics(self.Gr, self.Ref)

        # init Environment
        self.env_var = EDMF_Environment.EnvironmentVariables(namelist, self.Gr)

        for idx in range(self.Gr.gw, self.Gr.nzg - self.Gr.gw):
            self.env_var.QT.values[idx]       = p_dict['env_qt'][0,idx - self.Gr.gw]
            self.env_var.H.values[idx]        = p_dict['env_thetali'][0,idx - self.Gr.gw]
            self.env_var.EnvArea.values[idx]  = p_dict['env_fraction'][0,idx - self.Gr.gw]
            self.env_var.Hvar.values[idx]     = v_dict['Hvar'][0,idx - self.Gr.gw]
            self.env_var.QTvar.values[idx]    = v_dict['QTvar'][0,idx - self.Gr.gw]
            self.env_var.HQTcov.values[idx]   = v_dict['HQTcov'][0,idx - self.Gr.gw]
        self.env_thr  = EDMF_Environment.EnvironmentThermodynamics(namelist, self.Gr, self.Ref, self.env_var, self.rain_var)

        # init Updrafts
        self.upd_var = EDMF_Updrafts.UpdraftVariables(1, namelist, paramlist, self.Gr)
        # cant assign nupy array to memory view slice
        # https://mail.python.org/pipermail/cython-devel/2014-December/004288.html
        cdef double[:] tmp_qt       = p_dict['updraft_qt'][0,:]
        cdef double[:] tmp_thetali  = p_dict['updraft_thetali'][0,:]
        cdef double[:] tmp_area     = p_dict['updraft_fraction'][0,:]
        self.upd_var.QT.values[0, self.Gr.gw : -self.Gr.gw]   = tmp_qt
        self.upd_var.H.values[0, self.Gr.gw : -self.Gr.gw]    = tmp_thetali
        self.upd_var.Area.values[0, self.Gr.gw : -self.Gr.gw] = tmp_area
        self.upd_thr  = EDMF_Updrafts.UpdraftThermodynamics(1, self.Gr, self.Ref, self.upd_var, self.rain_var)

        # init output
        self.Case.initialize_reference(self.Gr, self.Ref, self.Stats)
        self.Case.initialize_profiles(self.Gr, self.GMV, self.Ref)
        self.GMV.initialize_io(self.Stats)
        self.upd_var.initialize_io(self.Stats)
        self.env_var.initialize_io(self.Stats)
        self.rain_var.initialize_io(self.Stats)

    def read_in_LES_data(self, int it):

        # read updrafts
        # cant assign nupy array to memory view slice
        # https://mail.python.org/pipermail/cython-devel/2014-December/004288.html
        cdef double[:] tmp_qt       = p_dict['updraft_qt'][it,:]
        cdef double[:] tmp_thetali  = p_dict['updraft_thetali'][it,:]
        cdef double[:] tmp_area     = p_dict['updraft_fraction'][it,:]
        self.upd_var.QT.values[0,   self.Gr.gw : -self.Gr.gw] = tmp_qt
        self.upd_var.H.values[0,    self.Gr.gw : -self.Gr.gw] = tmp_thetali
        self.upd_var.Area.values[0, self.Gr.gw : -self.Gr.gw] = tmp_area

        for idx in range(self.Gr.gw, self.Gr.nzg - self.Gr.gw):
            # read environment
            self.env_var.QT.values[idx]      = p_dict['env_qt'][it,idx-self.Gr.gw]
            self.env_var.H.values[idx]       = p_dict['env_thetali'][it,idx-self.Gr.gw]
            self.env_var.EnvArea.values[idx] = p_dict['env_fraction'][it,idx-self.Gr.gw]
            self.env_var.Hvar.values[idx]    = v_dict['Hvar'][it,idx-self.Gr.gw]
            self.env_var.QTvar.values[idx]   = v_dict['QTvar'][it,idx-self.Gr.gw]
            self.env_var.HQTcov.values[idx]  = v_dict['HQTcov'][it,idx-self.Gr.gw]

            # read GMV
            self.GMV.W.values[idx]    = p_dict["w_mean"][it,idx-self.Gr.gw]
            self.GMV.QT.values[idx]   = p_dict["qt_mean"][it,idx-self.Gr.gw]
            self.GMV.THL.values[idx]  = p_dict["thetali_mean"][it,idx-self.Gr.gw]
            self.GMV.T.values[idx]    = p_dict["temperature_mean"][it,idx-self.Gr.gw]

    def do_environment(self):

        if namelist['thermodynamics']['saturation'] == 'sa_quadrature':
            self.env_thr.eos_update_SA_sgs(self.env_var, self.rain_var)
        else:
            self.env_thr.eos_update_SA_mean(self.env_var, self.rain_var)

    cpdef do_updrafts(self):
        cdef:
            Py_ssize_t k
            eos_struct sa

        for k in xrange(self.Gr.gw, self.Gr.nzg - self.Gr.gw):
            if self.upd_var.Area.values[0,k] != 0:

                # saturation adjustment
                sa = eos(
                    self.upd_thr.t_to_prog_fp,
                    self.upd_thr.prog_to_t_fp,
                    self.Ref.p0_half[k],
                    self.upd_var.QT.values[0,k],
                    self.upd_var.H.values[0,k]
                )

                # autoconversion
                mph = microphysics(
                    sa.T,
                    sa.ql,
                    self.Ref.p0_half[k],
                    self.upd_var.QT.values[0,k],
                    self.upd_var.Area.values[0,k],
                    self.rain_var.max_supersaturation,
                    True
                )

                # update updraft variables
                self.upd_var.QL.values[0,k] = mph.ql
                self.upd_var.T.values[0,k]  = mph.T
                self.upd_var.QT.values[0,k] = mph.qt
                self.upd_var.H.values[0,k]  = mph.thl

                #TODO - do rain
                #if self.Rain.rain_model:
                #    #TODO - after adding rain evaporation think on the correct order of updates
                #    #       and on the correct source terms for QT and H
                #    self.UpdThermo.update_UpdRain(
                #        &self.UpdVar.Area.new[i,k],
                #        &self.Rain.Upd_QR.new[k],
                #        &self.Rain.Upd_RainArea.new[k], mph.qr,
                #        self.Rain.rain_area_value, i, k
                #    )

    def run(self):

        self.TS.t = 0.
        self.it   = 0

        while self.TS.t <= self.TS.t_max:
            self.read_in_LES_data(self.it)
            self.do_updrafts()
            self.do_environment()
            #self.do_rain()
            self.Stats.open_files()
            self.Stats.write_simulation_time(self.TS.t)
            #self.Stats.write_profile('updraft_ql', self.upd_var.QL.values[0, self.Gr.gw : self.Gr.nzg-self.Gr.gw])
            self.GMV.io(self.Stats)
            self.upd_var.io(self.Stats, self.Ref)
            self.env_var.io(self.Stats, self.Ref)
            #self.rain_var.io(self.Stats, self.Ref)
            self.Stats.close_files()
            self.TS.t += self.TS.dt
            self.it += 1

            #print self.TS.t
            #print self.it
            #print np.min(self.upd_var.QL.values[0,:]), np.max(self.upd_var.QL.values[0,:])
            #print np.min(self.env_var.QL.values[:]), np.max(self.env_var.QL.values[:])
            #print "-----------------------------------------------------"

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

