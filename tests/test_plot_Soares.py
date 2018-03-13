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
    setup["namelist"]['stats_io']['frequency'] = setup["namelist"]['time_stepping']['t_max']

    print " "
    print "namelist"
    print pp.pprint(setup["namelist"])
    print " "
    print "paramlist"
    print pp.pprint(setup["paramlist"])

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
    ref_data = Dataset("tests/reference/Soares/stats/Stats.Soares.nc", 'r')

    data_to_plot = pls.read_data(sim_data, ref_data)

    pls.plot_mean(data_to_plot,   "Soares_quicklook.pdf")
    pls.plot_drafts(data_to_plot, "Soares_quicklook_drafts.pdf")

