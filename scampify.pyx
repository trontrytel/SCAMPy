import cython
import numpy as np
import math as mt
from netCDF4 import Dataset

from Cases import CasesFactory
from Cases import CasesBase
from Grid  cimport Grid
from ReferenceState cimport ReferenceState
from Variables cimport GridMeanVariables
from NetCDFIO cimport NetCDFIO_Stats

cimport TimeStepping
cimport EDMF_Updrafts
cimport EDMF_Environment
cimport EDMF_Rain

from thermodynamic_functions cimport  *
from microphysics_functions cimport  *

cdef class Scampify1d:

    def __init__(self, namelist, paramlist, scampifylist):

        self.it    = 0 # TODO - get rid of me
        self.norm  = scampifylist["les_stats_freq"] / namelist["time_stepping"]["dt"]

        # create instances of grid, reference state, ...
        self.Gr    = Grid(namelist)
        self.Ref   = ReferenceState(self.Gr)
        self.GMV   = GridMeanVariables(namelist, self.Gr, self.Ref)
        self.Case  = CasesFactory(namelist, paramlist)

        self.TS    = TimeStepping.TimeStepping(namelist)
        self.Stats = NetCDFIO_Stats(namelist, paramlist, self.Gr)

        return

    # helper functions for reading LES netcdf files
    def read_profiles(self, LES_data, var_list):
        p_dict = {}

        for var in var_list:
            p_dict[var] = np.array(LES_data["profiles/" + var][:, :])

        p_dict['z_half'] = np.array(LES_data["profiles/z_half"][:])
        p_dict['z']      = np.array(LES_data["profiles/z"][:])
        p_dict['time']   = np.array(LES_data["profiles/t"][:])

        return p_dict

    def read_reference(self, LES_data, ref_list):
        r_dict = {}

        for var in ref_list:
            r_dict[var] = np.array(LES_data["reference/" + var][:])

        return r_dict

    def helper_var_reader(self, name1, name2, name3, LES_data):
        var1 = np.array(LES_data["profiles/" + name1][:, :])
        var2 = np.array(LES_data["profiles/" + name2][:, :])
        var3 = np.array(LES_data["profiles/" + name3][:, :])

        return var1 - var2 * var3

    def read_varcovar(self, LES_data):

        v_dict = {}

        v_dict['Hvar']   = self.helper_var_reader("env_thetali2", "env_thetali", "env_thetali", LES_data)
        v_dict['QTvar']  = self.helper_var_reader("env_qt2", "env_qt", "env_qt", LES_data)
        v_dict['HQTcov'] = self.helper_var_reader("env_qt_thetali", "env_thetali", "env_qt", LES_data)

        return v_dict

    def initialize(self, namelist, paramlist, scampifylist):

        nc_file  = scampifylist["nc_file"]
        LES_data = Dataset(nc_file)
        if namelist['meta']['casename'] in ['Bomex']:
           half_flag = True
        else:
           half_flag = False

        var_list = [\
            'ql_mean',         'env_ql',         'updraft_ql',\
            'qt_mean',         'env_qt',         'updraft_qt',\
            'thetali_mean',    'env_thetali',    'updraft_thetali',\
            'temperature_mean',\
                               'env_fraction',   'updraft_fraction',\
            'w_mean']
            #'qrain_mean',\

        #if half_flag:
        #    ref_list = ['p0'     , 'p0_half', 'alpha0'     , 'alpha0_half', 'rho0'     , 'rho0_half']
        #else:
        ref_list = ['p0_full', 'p0'     , 'alpha0_full', 'alpha0'     , 'rho0_full', 'rho0']

        self.r_dict = self.read_reference(LES_data, ref_list)
        self.p_dict = self.read_profiles( LES_data, var_list)
        self.v_dict = self.read_varcovar( LES_data)

        # read in reference profiles
        for idx in range(self.Gr.gw, self.Gr.nzg - self.Gr.gw):
            #if half_flag:
            #    self.Ref.p0[idx]          = self.r_dict['p0'         ][idx - self.Gr.gw]
            #    self.Ref.p0_half[idx]     = self.r_dict['p0_half'    ][idx - self.Gr.gw]
            #    self.Ref.alpha0[idx]      = self.r_dict['alpha0'     ][idx - self.Gr.gw]
            #    self.Ref.alpha0_half[idx] = self.r_dict['alpha0_half'][idx - self.Gr.gw]
            #    self.Ref.rho0[idx]        = self.r_dict['rho0'       ][idx - self.Gr.gw]
            #    self.Ref.rho0_half[idx]   = self.r_dict['rho0_half'  ][idx - self.Gr.gw]
            #else:
            self.Ref.p0[idx]          =self.r_dict['p0_full'    ][idx - self.Gr.gw]
            self.Ref.p0_half[idx]     =self.r_dict['p0'         ][idx - self.Gr.gw]
            self.Ref.alpha0[idx]      =self.r_dict['alpha0_full'][idx - self.Gr.gw]
            self.Ref.alpha0_half[idx] =self.r_dict['alpha0'     ][idx - self.Gr.gw]
            self.Ref.rho0[idx]        =self.r_dict['rho0_full'  ][idx - self.Gr.gw]
            self.Ref.rho0_half[idx]   =self.r_dict['rho0'       ][idx - self.Gr.gw]

        # read in initial GMV
        for idx in range(self.Gr.gw, self.Gr.nzg - self.Gr.gw):
            self.GMV.W.values[idx]    = self.p_dict["w_mean"          ][0,idx - self.Gr.gw]
            self.GMV.QT.values[idx]   = self.p_dict["qt_mean"         ][0,idx - self.Gr.gw]
            self.GMV.THL.values[idx]  = self.p_dict["thetali_mean"    ][0,idx - self.Gr.gw]
            self.GMV.T.values[idx]    = self.p_dict["temperature_mean"][0,idx - self.Gr.gw]

        # init Environment, Updraft and Rain
        self.env_var   = EDMF_Environment.EnvironmentVariables(namelist, self.Gr)
        self.upd_var   = EDMF_Updrafts.UpdraftVariables(1, namelist, paramlist, self.Gr)
        self.rain_var  = EDMF_Rain.RainVariables(namelist, self.Gr)

        # and their functions
        self.env_thr   = EDMF_Environment.EnvironmentThermodynamics(namelist, self.Gr, self.Ref, self.env_var, self.rain_var)
        self.upd_thr   = EDMF_Updrafts.UpdraftThermodynamics(1, self.Gr, self.Ref, self.upd_var, self.rain_var)
        self.rain_phs  = EDMF_Rain.RainPhysics(self.Gr, self.Ref)

        # read in initial updraft and environment properties
        for idx in range(self.Gr.gw, self.Gr.nzg - self.Gr.gw):
            self.env_var.QT.values[idx]       = self.p_dict['env_qt'          ][0,idx - self.Gr.gw]
            self.env_var.H.values[idx]        = self.p_dict['env_thetali'     ][0,idx - self.Gr.gw]
            self.env_var.EnvArea.values[idx]  = self.p_dict['env_fraction'    ][0,idx - self.Gr.gw]
            self.env_var.Hvar.values[idx]     = self.v_dict['Hvar'            ][0,idx - self.Gr.gw]
            self.env_var.QTvar.values[idx]    = self.v_dict['QTvar'           ][0,idx - self.Gr.gw]
            self.env_var.HQTcov.values[idx]   = self.v_dict['HQTcov'          ][0,idx - self.Gr.gw]
            self.upd_var.QT.values[0,   idx]  = self.p_dict['updraft_qt'      ][0,idx - self.Gr.gw]
            self.upd_var.H.values[0,    idx]  = self.p_dict['updraft_thetali' ][0,idx - self.Gr.gw]
            self.upd_var.Area.values[0, idx]  = self.p_dict['updraft_fraction'][0,idx - self.Gr.gw]

        # initialize output
        self.Case.initialize_reference(self.Gr, self.Ref, self.Stats)
        self.Case.initialize_profiles(self.Gr, self.GMV, self.Ref)
        self.GMV.initialize_io(self.Stats)
        self.upd_var.initialize_io(self.Stats)
        self.env_var.initialize_io(self.Stats)
        self.rain_var.initialize_io(self.Stats)

    def read_in_LES_data(self, int it):

        for idx in range(self.Gr.gw, self.Gr.nzg - self.Gr.gw):
            # read updrafts
            self.upd_var.QT.values[0,   idx] = self.p_dict['updraft_qt'      ][it, idx-self.Gr.gw]
            self.upd_var.H.values[0,    idx] = self.p_dict['updraft_thetali' ][it, idx-self.Gr.gw]
            self.upd_var.Area.values[0, idx] = self.p_dict['updraft_fraction'][it, idx-self.Gr.gw]

            # read environment
            self.env_var.QT.values[idx]      = self.p_dict['env_qt'      ][it, idx-self.Gr.gw]
            self.env_var.H.values[idx]       = self.p_dict['env_thetali' ][it, idx-self.Gr.gw]
            self.env_var.EnvArea.values[idx] = self.p_dict['env_fraction'][it, idx-self.Gr.gw]
            self.env_var.Hvar.values[idx]    = self.v_dict['Hvar'        ][it, idx-self.Gr.gw]
            self.env_var.QTvar.values[idx]   = self.v_dict['QTvar'       ][it, idx-self.Gr.gw]
            self.env_var.HQTcov.values[idx]  = self.v_dict['HQTcov'      ][it, idx-self.Gr.gw]

            #print idx, self.env_var.H.values[idx], self.env_var.QT.values[idx], self.env_var.EnvArea.values[idx]

            # read GMV
            self.GMV.W.values[idx]    = self.p_dict["w_mean"          ][it, idx-self.Gr.gw]
            self.GMV.QT.values[idx]   = self.p_dict["qt_mean"         ][it, idx-self.Gr.gw]
            self.GMV.THL.values[idx]  = self.p_dict["thetali_mean"    ][it, idx-self.Gr.gw]
            self.GMV.T.values[idx]    = self.p_dict["temperature_mean"][it, idx-self.Gr.gw]

    def do_environment(self):

        if self.env_var.EnvThermo_scheme == 'sa_quadrature':
            self.env_thr.eos_update_SA_sgs(self.env_var, self.rain_var)
        else:
            self.env_thr.eos_update_SA_mean(self.env_var, self.rain_var)

    cpdef do_updrafts(self):
        cdef:
            Py_ssize_t k
            eos_struct sa

        for k in xrange(self.Gr.gw, self.Gr.nzg - self.Gr.gw):

            if self.upd_var.Area.values[0,k] > 1e-20:

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
                    False
                )

                # update updraft variables
                self.upd_var.QL.values[0,k] = mph.ql
                self.upd_var.T.values[0, k] = mph.T
                self.upd_var.QT.values[0,k] = mph.qt
                self.upd_var.H.values[0, k] = mph.thl

            else:
                self.upd_var.Area.values[0,k] = 0.
                self.upd_var.QL.values[0,k] = 0.
                self.upd_var.T.values[0, k] = self.GMV.T.values[k]
                self.upd_var.QT.values[0,k] = self.GMV.QT.values[k]
                self.upd_var.H.values[0, k] = self.GMV.H.values[k]

                #TODO - update rain variables

    def run(self):

        self.TS.t = 0.

        print "scampify!"

        while self.TS.t <= self.TS.t_max:
            print self.it, self.norm, mt.floor(self.it / self.norm)
            self.read_in_LES_data(mt.floor(self.it / self.norm)) #TODO - change to interpolating
            self.upd_var.set_means(self.GMV)
            self.do_updrafts()
            self.do_environment()
            #self.do_rain()
            self.upd_var.set_means(self.GMV)

            # TODO - calculate the GMV QL and QR
            for idx in range(self.Gr.gw, self.Gr.nzg - self.Gr.gw):

                #self.upd_var.QL.bulkvalues[idx] = self.upd_var.QL.values[0, idx]
                self.upd_var.Area.bulkvalues[idx] = self.upd_var.Area.values[0, idx]

                self.GMV.QT.values[idx] = self.upd_var.QT.bulkvalues[idx] * self.upd_var.Area.bulkvalues[idx] +\
                                          self.env_var.QT.values[idx]     * self.env_var.EnvArea.values[idx]

                self.GMV.QL.values[idx] = self.upd_var.QL.bulkvalues[idx] * self.upd_var.Area.bulkvalues[idx] +\
                                          self.env_var.QL.values[idx]     * self.env_var.EnvArea.values[idx]

                self.GMV.QR.values[idx] = 0.#self.upd_var.QR.bulkvalues[idx] * self.upd_var.Area.bulkvalues[idx] +\
                #self.env_var.QR.values[idx]     * self.env_var.EnvArea.values[idx]


            self.Stats.open_files()
            self.Stats.write_simulation_time(self.TS.t)
            self.GMV.io(self.Stats)
            self.upd_var.io(self.Stats, self.Ref)
            self.env_var.io(self.Stats, self.Ref)
            self.rain_var.io(self.Stats, self.Ref)
            self.Stats.close_files()

            self.TS.t += self.TS.dt
            self.it += 1

        print "done"

        return
