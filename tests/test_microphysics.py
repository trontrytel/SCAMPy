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

# list of possible test cases
case_list = ['Soares', 'DYCOMS_RF01', 'Bomex', 'life_cycle_Tan2018', 'Rico', 'GATE_III', 'ARM_SGP', 'TRMM_LBA']

@pytest.fixture(scope="module")
def data(request):

    # dictionary where simulation results will be stored
    data = {}

    # loop over all test cases
    for case in case_list:

        # generate namelist and paramlist
        setup = pls.simulation_setup(case)
        # change the default parameters
        setup["namelist"]['stats_io']['frequency'] = setup["namelist"]['time_stepping']['t_max']

        #print " "
        #print "namelist"
        #print pp.pprint(setup["namelist"])
        #print " "
        #print "paramlist"
        #print pp.pprint(setup["paramlist"])

        # run scampy
        scampy.main1d(setup["namelist"], setup["paramlist"])
        
        # simulation results 
        data[case] = Dataset(setup["outfile"], 'r')

    request.addfinalizer(pls.removing_files)

    return data

@pytest.mark.parametrize("case", case_list)
def test_mean_qt(data, case, eps = 1e-2):
    """
    Check if the mean qt is equal to updraft_area * updraft_qt + (1 - updraft_area) * env_q
    """
    #print data.groups
    #z_dim = data["profiles"].dimensions['z'].size
    #t_dim = data["profiles"].dimensions['t'].size
    #rho_ref  = np.array(data["reference/rho0"][:])

    # read in the data
    qt_mean  = np.array(data[str(case)]["profiles/qt_mean"][:,:])
    udr_area = np.array(data[str(case)]["profiles/updraft_area"][:,:])
    udr_qt   = np.array(data[str(case)]["profiles/updraft_qt"][:,:])
    env_qt   = np.array(data[str(case)]["profiles/env_qt"][:,:])

    # calculate the mean qt basing un updraft and environmet means
    tmp_mean = udr_area * udr_qt + (1 - udr_area) * env_qt

    #print "case = ", case
    #for it in xrange(len(tmp_mean[-1,:])):
    #    print tmp_mean[-1, it], " vs ", qt_mean[-1, it]

    for idx in [1, -1]:
        assert np.allclose(qt_mean[idx,:], tmp_mean[idx,:], rtol = eps),\
            "mean != updraft_area * updraft_mean + (1-updraft_area) * env_mean"
