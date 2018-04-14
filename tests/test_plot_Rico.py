import sys
sys.path.insert(0, "./")

import os
import subprocess
import json
import warnings

from netCDF4 import Dataset

import pytest
import pprint as pp
import numpy as np

import main as scampy
import plot_scripts as pls

@pytest.fixture(scope="module")
def sim_data(request):

    # generate namelists and paramlists
    setup = pls.simulation_setup('Rico')
    # chenge the defaults  
    #setup["namelist"]['stats_io']['frequency'] = setup["namelist"]['time_stepping']['t_max']
    setup['namelist']['turbulence']['EDMF_PrognosticTKE']['use_similarity_diffusivity'] = False
    setup["namelist"]['turbulence']['EDMF_PrognosticTKE']['use_local_micro'] = True
    setup['namelist']['turbulence']['EDMF_PrognosticTKE']['use_scalar_var'] = True
    setup['namelist']['turbulence']['sgs'] = {}
    setup['namelist']['turbulence']['sgs']['use_prescribed_scalar_var'] = True
    setup['namelist']['turbulence']['sgs']['prescribed_QTvar'] = 0.5 * 1e-7
    setup['namelist']['turbulence']['sgs']['prescribed_Hvar'] = 0.01
    setup['namelist']['turbulence']['sgs']['prescribed_HQTcov'] = -1e-3
    setup['paramlist']['turbulence']['updraft_microphysics']['max_supersaturation'] = 0.01 #0.1

    setup['namelist']['thermodynamics']['saturation'] = 'sa_quadrature'        

    #print " "
    #print "namelist"
    #print pp.pprint(setup["namelist"])
    #print " "
    #print "paramlist"
    #print pp.pprint(setup["paramlist"])

    # run scampy
    scampy.main1d(setup["namelist"], setup["paramlist"])
    
    # simulation results 
    sim_data = Dataset(setup["outfile"], 'r')

    # remove netcdf file after tests
    request.addfinalizer(pls.removing_files)

    return sim_data

def test_plot_Rico(sim_data):
    """
    plot Rico profiles
    """
    data_to_plot = pls.read_data_avg(sim_data, 100)

    pls.plot_mean(data_to_plot,   "Rico_quicklook.pdf")
    pls.plot_drafts(data_to_plot, "Rico_quicklook_drafts.pdf")

def test_plot_timeseries_Rico(sim_data):
    """
    plot timeseries
    """
    data_to_plot = pls.read_data_srs(sim_data)

    pls.plot_timeseries(data_to_plot, "Rico")

