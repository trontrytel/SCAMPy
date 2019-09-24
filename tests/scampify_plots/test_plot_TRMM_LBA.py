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
    setup = cmn.simulation_setup('TRMM_LBA')
    setup['namelist']['thermodynamics']['sgs'] = 'quadrature'
    setup['namelist']['microphysics']['rain_model'] = True
    setup["namelist"]["turbulence"]["EDMF_PrognosticTKE"]["entrainment"]="moisture_deficit"

    setup["paramlist"]['turbulence']['EDMF_PrognosticTKE']['entrainment_factor'] = 0.2       # 0.15
    setup["paramlist"]['turbulence']['EDMF_PrognosticTKE']['detrainment_factor'] = 3         # 2
    setup["paramlist"]['turbulence']['EDMF_PrognosticTKE']['entrainment_erf_const'] = 2      # 2

    # additional parameters for offline runs
    scampifylist = {}
    scampifylist["offline"] = False
    scampifylist["les_stats_freq"] = 10.

    subprocess.call("python setup.py build_ext --inplace", shell=True, cwd='../')
    if scampifylist["offline"]:
        # run scampy offline
        print "offline run"
        scampy.main_scampify(setup["namelist"], setup["paramlist"], scampifylist)
    else:
        # run scampy online
        print "online run"
        #scampy.main1d(setup["namelist"], setup["paramlist"])

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

def test_plot_timeseries_TRMM_LBA(sim_data):
    """
    plot TRMM_LBA timeseries
    """
    # make directory
    localpath = os.getcwd()
    try:
        os.mkdir(localpath + "/scampify_plots/output/TRMM_LBA/")
    except:
        print('TRMM_LBA folder exists')

    scm_data_to_plot = cmn.read_scm_data_srs(sim_data["scm_data"])
    les_data_to_plot = cmn.read_les_data_srs(sim_data["les_data"])

    scm_data_to_plot_timesrs = cmn.read_scm_data_timeseries(sim_data["scm_data"])
    les_data_to_plot_timesrs = cmn.read_les_data_timeseries(sim_data["les_data"])

    pls.plot_humidities(
        scm_data_to_plot,
        les_data_to_plot,
        5,
        6,
        "TRMM_LBA_humidities.pdf",
        folder="scampify_plots/output/TRMM_LBA/"
    )

    pls.plot_updraft_properties(
        scm_data_to_plot,
        les_data_to_plot,
        5,
        6,
        "TRMM_LBA_updraft_properties.pdf",
        folder="scampify_plots/output/TRMM_LBA/"
    )

    pls.plot_timeseries_1D(
        scm_data_to_plot_timesrs,
        les_data_to_plot_timesrs,
        folder="scampify_plots/output/TRMM_LBA/"
    )

    cb_min = [280, 280, 294,  0,  0,  0,    0,     0,   0,    0,     0,   0, -0.35,    0,    0]
    cb_max = [370, 370, 348, 20, 20, 20, 0.09, 0.028, 2.8, 0.09, 0.064, 2.4,     0, 0.28, 10.5]

    pls.plot_timeseries(
        scm_data_to_plot,
        les_data_to_plot,
        cb_min,
        cb_max,
        folder="scampify_plots/output/TRMM_LBA/"
    )
