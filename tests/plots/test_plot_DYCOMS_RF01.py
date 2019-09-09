import sys
sys.path.insert(0, "./")
sys.path.insert(0, "../")

import os
import subprocess
import json
import warnings

from netCDF4 import Dataset

import pytest
import numpy as np

import main as scampy
import common as cmn
import plot_scripts as pls

@pytest.fixture(scope="module")
def sim_data(request):

    # generate namelists and paramlists
    setup = cmn.simulation_setup('DYCOMS_RF01')

    setup['namelist']['grid']['dims'] = 1
    setup['namelist']['grid']['dz'] = 5.
    setup['namelist']['grid']['gw'] = 2
    setup['namelist']['grid']['nz'] = 300

    setup['namelist']['time_stepping']['t_max'] = 14400
    setup['namelist']['time_stepping']['dt'] = 2

    setup['namelist']["stats_io"]["frequency"] = 60.

    # additional parameters for offline runs
    scampifylist = {}
    scampifylist["offline"] = False
    scampifylist["nc_file"] = '../SCAMPY_tests/plots/DYCOMS_comparison/DYCOMS_SA_LES_new_tracers/Stats.DYCOMS_RF01.Restart_3.nc'
    scampifylist["les_stats_freq"] = 10.

    if scampifylist["offline"]:
        # run scampy offline
        print "offline run"
        scampy.main_scampify(setup["namelist"], setup["paramlist"], scampifylist)
    else:
        # run scampy online
        print "online run"
        scampy.main1d(setup["namelist"], setup["paramlist"])

    # simulation results
    sim_data = Dataset(setup["outfile"], 'r')

    # remove netcdf files after tests
    request.addfinalizer(cmn.removing_files)

    return sim_data

def test_plot_DYCOMS_RF01(sim_data):
    """
    plot DYCOMS_RF01 quicklook profiles
    """
    data_to_plot = cmn.read_data_avg(sim_data, n_steps=100)

    pls.plot_mean(data_to_plot,   "DYCOMS_RF01_quicklook.pdf")
    pls.plot_drafts(data_to_plot, "DYCOMS_RF01_quicklook_drafts.pdf")

def test_plot_var_covar_DYCOMS_RF01(sim_data):
    """
    plot DYCOMS_RF01 quicklook profiles
    """
    data_to_plot = cmn.read_data_avg(sim_data, n_steps=100, var_covar=True)

    pls.plot_var_covar_mean(data_to_plot,       "DYCOMS_RF01_var_covar_mean.pdf")
    pls.plot_var_covar_components(data_to_plot, "DYCOMS_RF01_var_covar_components.pdf")

def test_plot_timeseries_DYCOMS(sim_data):
    """
    plot timeseries
    """
    data_to_plot = cmn.read_data_srs(sim_data)

    pls.plot_timeseries(data_to_plot, "DYCOMS")

def test_plot_timeseries_1D_DYCOMS_RF01(sim_data):
    """
    plot DYCOMS_RF01 1D timeseries
    """
    data_to_plot = cmn.read_data_timeseries(sim_data)

    pls.plot_timeseries_1D(data_to_plot, "DYCOMS_RF01_timeseries_1D.pdf")

def test_DYCOMS_RF01_radiation(sim_data):
    """
    plots DYCOMS_RF01
    """
    import matplotlib as mpl
    mpl.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
    import matplotlib.pyplot as plt

    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    plt_data = cmn.read_data_avg(sim_data,     n_steps=100, var_covar=False)
    rad_data = cmn.read_rad_data_avg(sim_data, n_steps=100)

    plots = []
    # loop over simulation and reference data for t=0 and t=-1
    x_lab  = ['longwave radiative flux [W/m2]', 'dTdt [K/day]',       'QT [g/kg]',         'QL [g/kg]']
    legend = ["lower right",                    "lower left",         "lower left",        "lower right"]
    line   = ['--',                             '--',                 '-',                 '-']
    plot_y = [rad_data["rad_flux"],             rad_data["rad_dTdt"], plt_data["qt_mean"], plt_data["ql_mean"]]
    plot_x = [rad_data["z"],                    plt_data["z_half"],   plt_data["z_half"],  plt_data["z_half"]]
    color  = ["palegreen",                      "forestgreen"]
    label  = ["ini",                            "end"        ]

    for plot_it in range(4):
        plots.append(plt.subplot(2,2,plot_it+1))
                              #(rows, columns, number)
        for it in range(2):
            plots[plot_it].plot(plot_y[plot_it][it], plot_x[plot_it], '.-', color=color[it], label=label[it])
        plots[plot_it].legend(loc=legend[plot_it])
        plots[plot_it].set_xlabel(x_lab[plot_it])
        plots[plot_it].set_ylabel('z [m]')
    plots[2].set_xlim([1, 10])
    plots[3].set_xlim([-0.1, 0.5])

    plt.savefig("plots/output/DYCOMS_RF01_radiation.pdf")
    plt.clf()

