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
import common_scampify as cmn
import plot_scripts_scampify as pls

@pytest.fixture(scope="module")
def sim_data(request):

    # generate namelists and paramlists
    cmn.removing_files
    setup = cmn.simulation_setup('DYCOMS_RF01')

    setup['namelist']['grid']['gw'] = 3
    setup['namelist']['grid']['dz'] = 5.
    setup['namelist']['grid']['nz'] = 300

    setup['namelist']['time_stepping']['t_max'] = 14400
    setup["namelist"]['time_stepping']['dt'] = 4.

    setup['namelist']['microphysics']['rain_model'] = 'clima_1m'
    setup['namelist']['thermodynamics']['sgs'] = 'quadrature'
    #setup["namelist"]["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]="moisture_deficit"

    # additional parameters for offline runs
    scampifylist = {}
    scampifylist["offline"] = False
    scampifylist["les_stats_freq"] = 60.
    scampifylist["les_outfile"] = setup["les_outfile"]

    subprocess.call("python setup.py build_ext --inplace", shell=True, cwd='../')
    if scampifylist["offline"]:
        # run scampy offline
        print "offline run"
        scampy.main_scampify(setup["namelist"], setup["paramlist"], scampifylist)
    else:
        # run scampy online
        print "online run"
        scampy.main1d(setup["namelist"], setup["paramlist"])

    # simulation results
    scm_data = Dataset(setup["scm_outfile"], 'r')
    les_data = Dataset(setup["les_outfile"], 'r')

    sim_data = {}
    sim_data["scm_data"] = scm_data
    sim_data["les_data"] = les_data
    sim_data["scampifylist"] = scampifylist

    # remove netcdf files after tests
    #request.addfinalizer(cmn.removing_files)

    return sim_data

def test_plot_timeseries_DYCOMS_RF01(sim_data):
    """
    plot DYCOMS_RF01 timeseries
    """
    localpath = os.getcwd()
    try:
        os.mkdir(localpath + "/scampify_plots/output/DYCOMS_RF01_drizzle/")
    except:
        print('DYCOMS_RF01_drizzle folder exists')

    scm_data_to_plot = cmn.read_scm_data_srs(sim_data["scm_data"])
    les_data_to_plot = cmn.read_les_data_srs(sim_data["les_data"])

    scm_data_to_plot_timesrs = cmn.read_scm_data_timeseries(sim_data["scm_data"])
    les_data_to_plot_timesrs = cmn.read_les_data_timeseries(sim_data["les_data"])

    cb_min = [285,   285,   288.75,   0,      0,  9,   0,   0,    0,     0,     0,      0, -0.16,    0,   0]
    cb_max = [307.5, 307.5, 289.5, 12.5, 12.5, 11, 0.8, 0.8, 0.64, 0.002, 0.002, 0.0012,     0, 0.24, 1.5]
    t0 = 3
    t1 = 4
    folder = "scampify_plots/output/DYCOMS_RF01_drizzle/"
    case = "drizzle_DYCOMS_RF01_"

    pls.plot_humidities(scm_data_to_plot, les_data_to_plot, t0, t1, case+"humidities.pdf", folder=folder)
    pls.plot_cloud_rain_components(scm_data_to_plot, les_data_to_plot, t0, t1, case + "cloud_rain_comp.pdf", folder=folder)
    pls.plot_updraft_properties(scm_data_to_plot, les_data_to_plot, t0, t1, case+"updraft_properties.pdf", folder=folder)
    pls.plot_timeseries_1D(scm_data_to_plot_timesrs, les_data_to_plot_timesrs, folder=folder)
    pls.plot_timeseries(scm_data_to_plot, les_data_to_plot, cb_min, cb_max, folder=folder)