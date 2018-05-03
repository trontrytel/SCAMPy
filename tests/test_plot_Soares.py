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
    setup = pls.simulation_setup('Soares')
    # chenge the defaults  
    #setup["namelist"]['stats_io']['frequency'] = setup["namelist"]['time_stepping']['t_max']
    setup['namelist']['turbulence']['EDMF_PrognosticTKE']['use_scalar_var'] = True
    setup['namelist']['turbulence']['sgs'] = {}
    setup['namelist']['turbulence']['sgs']['use_prescribed_scalar_var'] = False
    setup['namelist']['turbulence']['sgs']['prescribed_QTvar'] = 0.5 * 1e-7
    setup['namelist']['turbulence']['sgs']['prescribed_Hvar'] = 0.01
    setup['namelist']['turbulence']['sgs']['prescribed_HQTcov'] = -1e-3

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

def test_plot_Soares(sim_data):
    """
    plot Soares profiles
    """
    data_to_plot = pls.read_data_avg(sim_data, 100)

    pls.plot_mean(data_to_plot,   "Soares_quicklook.pdf")
    pls.plot_drafts(data_to_plot, "Soares_quicklook_drafts.pdf")

def test_plot_timeseries_Soares(sim_data):
    """
    plot Soares timeseries
    """
    data_to_plot = pls.read_data_srs(sim_data)

    pls.plot_timeseries(data_to_plot, "Soares")

def test_plot_var_covar_Soares(sim_data):                                 
    """                                                                        
    plot Soares var covar profiles                                        
    """                                                                        
    data_to_plot = pls.read_data_avg(sim_data, 100)                            

    pls.plot_var_covar_mean(data_to_plot,       "Soares_var_covar_mean.pdf")                                                                                      
    pls.plot_var_covar_components(data_to_plot, "Soares_var_covar_comp.pdf")        
