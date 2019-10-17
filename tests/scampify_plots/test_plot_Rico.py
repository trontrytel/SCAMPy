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
    setup = cmn.simulation_setup('Rico')

    setup["namelist"]['microphysics']['rain_model'] = 'clima_1m'
    setup['namelist']['thermodynamics']['sgs'] = 'quadrature'
    #setup["namelist"]["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]="moisture_deficit"

    setup['namelist']['grid']['gw'] = 3
    setup['namelist']['grid']['nz'] = 150
    setup['namelist']['grid']['dz'] = 40

    setup["namelist"]['time_stepping']['dt'] = 10.
    setup["namelist"]['time_stepping']["t_max"] = 86400.0

    #setup["paramlist"]['turbulence']['EDMF_PrognosticTKE']['entrainment_factor'] = 0.075   # 0.15
    #setup["paramlist"]['turbulence']['EDMF_PrognosticTKE']['detrainment_factor'] = 1.75      # 2
    #setup["paramlist"]['turbulence']['EDMF_PrognosticTKE']['entrainment_erf_const'] = 2   # 2

    # additional parameters for offline runs
    scampifylist = {}
    scampifylist["offline"] = False
    scampifylist["les_stats_freq"] = 100.
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

    # remove netcdf file after tests
    #request.addfinalizer(cmn.removing_files)

    return sim_data

def test_plot_timeseries_Rico(sim_data):
    """
    plot Rico timeseries
    """
    localpath = os.getcwd()
    try:
        os.mkdir(localpath + "/scampify_plots/output/Rico/")
    except:
        print('Rico folder exists')

    scm_data_to_plot = cmn.read_scm_data_srs(sim_data["scm_data"])
    les_data_to_plot = cmn.read_les_data_srs(sim_data["les_data"])

    scm_data_to_plot_timesrs = cmn.read_scm_data_timeseries(sim_data["scm_data"])
    les_data_to_plot_timesrs = cmn.read_les_data_timeseries(sim_data["les_data"])

    cb_min = [296, 296, 297,  0,    0,   7.5,    0,    0,   0,    0,     0,   0,  -0.1,    0, 0]
    cb_max = [332, 332, 305, 17.5, 17.5, 18,  0.05, 0.02, 2.8, 0.007, 0.004, 0.56,   0, 0.24, 5]
    t0 = 22
    t1 = 24
    folder = "scampify_plots/output/Rico/"
    case = "Rico_"

    pls.plot_humidities(scm_data_to_plot, les_data_to_plot, t0, t1, case+"humidities.pdf", folder=folder)
    pls.plot_cloud_rain_components(scm_data_to_plot, les_data_to_plot, t0, t1, case + "cloud_rain_comp.pdf", folder=folder)
    pls.plot_updraft_properties(scm_data_to_plot, les_data_to_plot, t0, t1, case+"updraft_properties.pdf", folder=folder)
    pls.plot_timeseries_1D(scm_data_to_plot_timesrs, les_data_to_plot_timesrs, folder=folder)
    pls.plot_timeseries(scm_data_to_plot, les_data_to_plot, cb_min, cb_max, folder=folder)
