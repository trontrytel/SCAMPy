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
    setup = cmn.simulation_setup('TRMM_LBA')

    # run scampy
    subprocess.call("python setup.py build_ext --inplace", shell=True, cwd='../')
    scampy.main1d(setup["namelist"], setup["paramlist"])

    # simulation results
    sim_data = Dataset(setup["outfile"], 'r')

    # remove netcdf file after tests
    request.addfinalizer(cmn.removing_files)

    return sim_data

def test_plot_timeseries_TRMM_LBA(sim_data):
    """
    plot TRMM_LBA timeseries
    """
    # make directory
    localpath = os.getcwd()
    try:
        os.mkdir(localpath + "/plots/output/TRMM_LBA/")
    except:
        print('TRMM_LBA folder exists')
    try:
        os.mkdir(localpath + "/plots/output/TRMM_LBA/all_variables/")
    except:
        print('TRMM_LBA/all_variables folder exists')

    if (os.path.exists(localpath + "/les_data/TRMM_LBA.nc")):
        les_data = Dataset(localpath + "/les_data/TRMM_LBA.nc", 'r')
    else:
        url_ = "https://www.dropbox.com/s/snhxbzxt4btgiis/TRMM_LBA.nc?dl=0"
        os.system("wget -O "+localpath+"/les_data/TRMM_LBA.nc "+url_)
        les_data = Dataset(localpath + "/les_data/TRMM_LBA.nc", 'r')

    data_to_plot = cmn.read_data_srs(sim_data)
    les_data_to_plot = cmn.read_les_data_srs(les_data)

    pls.plot_closures(data_to_plot, les_data_to_plot,5,6,           "TRMM_LBA_closures.pdf",           folder="plots/output/TRMM_LBA/")
    pls.plot_humidities(data_to_plot, les_data_to_plot,5,6,         "TRMM_LBA_humidities.pdf",         folder="plots/output/TRMM_LBA/")
    pls.plot_updraft_properties(data_to_plot, les_data_to_plot,5,6, "TRMM_LBA_updraft_properties.pdf", folder="plots/output/TRMM_LBA/")
    pls.plot_tke_components(data_to_plot, les_data_to_plot, 5,6,    "TRMM_LBA_tke_components.pdf",     folder="plots/output/TRMM_LBA/")

    pls.plot_timeseries(data_to_plot, les_data_to_plot,          folder="plots/output/TRMM_LBA/all_variables/")
    pls.plot_mean(data_to_plot, les_data_to_plot,5,6,            folder="plots/output/TRMM_LBA/all_variables/")
    pls.plot_var_covar_mean(data_to_plot, les_data_to_plot, 5,6, "TRMM_LBA_var_covar_mean.pdf", folder="plots/output/TRMM_LBA/all_variables/")
    pls.plot_var_covar_components(data_to_plot,5,6,              "TRMM_LBA_var_covar_components.pdf", folder="plots/output/TRMM_LBA/all_variables/")
    pls.plot_tke_breakdown(data_to_plot, les_data_to_plot, 5,6,  "TRMM_LBA_tke_breakdown.pdf", folder="plots/output/TRMM_LBA/all_variables/")

def test_plot_timeseries_1D_TRMM_LBA(sim_data):
    """
    plot TRMM_LBA 1D timeseries
    """
    localpath = os.getcwd()
    try:
        os.mkdir(localpath + "/plots/output/TRMM_LBA/")
        print()
    except:
        print('TRMM_LBA folder exists')
    try:
        os.mkdir(localpath + "/plots/output/TRMM_LBA/all_variables/")
    except:
        print('TRMM_LBA/all_variables folder exists')

    if (os.path.exists(localpath + "/les_data/TRMM_LBA.nc")):
        les_data = Dataset(localpath + "/les_data/TRMM_LBA.nc", 'r')
    else:
        url_ = "https://www.dropbox.com/s/snhxbzxt4btgiis/TRMM_LBA.nc?dl=0"
        os.system("wget -O "+localpath+"/les_data/TRMM_LBA.nc "+url_)
        les_data = Dataset(localpath + "/les_data/TRMM_LBA.nc", 'r')

    data_to_plot = cmn.read_data_timeseries(sim_data)
    les_data_to_plot = cmn.read_les_data_timeseries(les_data)
    data_to_plot_ = cmn.read_data_srs(sim_data)
    les_data_to_plot_ = cmn.read_les_data_srs(les_data)

    pls.plot_main_timeseries(data_to_plot, les_data_to_plot, data_to_plot_, les_data_to_plot_, "TRMM_LBA_main_timeseries.pdf", folder="plots/output/TRMM_LBA/")
    pls.plot_timeseries_1D(data_to_plot,  les_data_to_plot,  folder="plots/output/TRMM_LBA/all_variables/")