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
    setup = cmn.simulation_setup('GATE_III')

    # run scampy
    subprocess.call("python setup.py build_ext --inplace", shell=True, cwd='../')
    scampy.main1d(setup["namelist"], setup["paramlist"])

    # simulation results
    sim_data = Dataset(setup["outfile"], 'r')

    # remove netcdf file after tests
    request.addfinalizer(cmn.removing_files)

    return sim_data

# @pytest.mark.skip(reason="deep convection not working with current defaults")
def test_plot_timeseries_GATE_III(sim_data):
    """
    plot GATE_III timeseries
    """
    # make directory
    localpath = os.getcwd()
    try:
        os.mkdir(localpath + "/plots/output/GATE_III/")
    except:
        print('GATE_III folder exists')
    les_data = Dataset(localpath + '/les_data/GATE_III.nc', 'r')
    data_to_plot = cmn.read_data_srs(sim_data)
    les_data_to_plot = cmn.read_les_data_srs(les_data)

    pls.plot_timeseries(data_to_plot, les_data_to_plot,          folder="plots/output/GATE_III/")
    pls.plot_mean(data_to_plot, les_data_to_plot,5,6,            folder="plots/output/GATE_III/")
    pls.plot_closures(data_to_plot, les_data_to_plot,5,6,        "GATE_III_closures.pdf", folder="plots/output/GATE_III/")
    pls.plot_var_covar_mean(data_to_plot, les_data_to_plot, 5,6, "GATE_III_var_covar_mean.pdf", folder="plots/output/GATE_III/")
    pls.plot_var_covar_components(data_to_plot,5,6,              "GATE_III_var_covar_components.pdf", folder="plots/output/GATE_III/")
    pls.plot_tke_components(data_to_plot, les_data_to_plot, 5,6, "GATE_III_tke_components.pdf", folder="plots/output/GATE_III/")
    pls.plot_tke_breakdown(data_to_plot, les_data_to_plot, 5,6,  "GATE_III_tke_breakdown.pdf", folder="plots/output/GATE_III/")

# @pytest.mark.skip(reason="deep convection not working with current defaults")
def test_plot_timeseries_1D_GATE_III(sim_data):
    """
    plot GATE_III 1D timeseries
    """
    # localpath = os.getcwd()
    # try:
    #     os.mkdir(localpath + "/plots/output/GATE_III/")
    # except:
    #     print('GATE_III folder exists')
    les_data = Dataset(localpath + '/les_data/GATE_III.nc', 'r')
    data_to_plot = cmn.read_data_timeseries(sim_data)
    les_data_to_plot = cmn.read_les_data_timeseries(les_data)

    pls.plot_timeseries_1D(data_to_plot,  les_data_to_plot, folder="plots/output/GATE_III/")