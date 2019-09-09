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
    cmn.removing_files
    setup = cmn.simulation_setup('DYCOMS_RF01')

    setup['namelist']['grid']['dims'] = 1
    setup['namelist']['grid']['dz'] = 5.
    setup['namelist']['grid']['gw'] = 2
    setup['namelist']['grid']['nz'] = 300

    setup['namelist']['time_stepping']['t_max'] = 14400
    setup['namelist']['time_stepping']['dt'] = 2

    setup['namelist']["stats_io"]["frequency"] = 60.

    setup['namelist']['thermodynamics']['sgs'] = 'quadrature'
    setup['namelist']['microphysics']['rain_model']          = True
    setup['namelist']['microphysics']['max_supersaturation'] = 0.05

    # additional parameters for offline runs
    scampifylist = {}
    scampifylist["offline"] = True
    scampifylist["nc_file"] = '../../SCAMPY_tests/plots/DYCOMS_comparison/DYCOMS_SA_LES_new_tracers/Stats.DYCOMS_RF01.Restart_3.nc'
    scampifylist["les_stats_freq"] = 10.

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
    sim_data = Dataset(setup["outfile"], 'r')

    # remove netcdf files after tests
    request.addfinalizer(cmn.removing_files)

    return sim_data

def test_plot_timeseries_DYCOMS_RF01(sim_data):
    """
    plot DYCOMS_RF01 timeseries
    """
    # make directory
    localpath = os.getcwd()
    try:
        os.mkdir(localpath + "/plots/output/DYCOMS_RF01/")
    except:
        print('DYCOMS_RF01 folder exists')
    try:
        os.mkdir(localpath + "/plots/output/DYCOMS_RF01/all_variables/")
    except:
        print('DYCOMS_RF01/all_variables folder exists')
    les_data = Dataset(localpath + '/les_data/DYCOMS_RF01.nc', 'r')
    data_to_plot = cmn.read_data_srs(sim_data)
    les_data_to_plot = cmn.read_les_data_srs(les_data)

    pls.plot_humidities(data_to_plot, les_data_to_plot,3,4,         "DYCOMS_RF01_humidities.pdf",         folder="plots/output/DYCOMS_RF01/")
