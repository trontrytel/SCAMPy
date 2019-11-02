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
import pprint as pp

import main as scampy
import common_scampify as cmn
import plot_scripts_scampify as pls

import matplotlib as mpl
mpl.use('Agg')  # To allow plotting when display is off
import matplotlib.pyplot as plt
from matplotlib import ticker

params = {}
cases = ["DYCOMS_RF01", "Rico", "TRMM_LBA"]

# parameters specific to each case
params["DYCOMS_RF01"] = \
    {"gw": 3, "dz": 5., "nz": 300, "dt": 4., "t_max": 14400, "les_stats_freq": 60.,\
     "cb_min": [285, 285, 288.75, 0, 0, 9, 0, 0, 0, 0, 0, 0, -0.16, 0, 0],\
     "cb_max": [307.5, 307.5, 289.5, 12.5, 12.5, 11, 0.8, 0.8, 0.64, 0.002, 0.002, 0.0012, 0, 0.24, 1.5],\
     "t0": 3, "t1": 4, "case_name": "drizzle_DYCOMS_RF01",\
     "ql_max": 0.4, "qr_max": 1.1e-4, "ql_min": -1e-2, "qr_min": -1e-5, "z_min": -0.02, "z_max": 1.2,\
     "les": "scampify_plots/data_les/DYCOMS_RF01_drizzle/CLIMA_1M_micro/Stats.DYCOMS_RF01.nc"}
params["Rico"] = \
    {"gw": 3, "dz": 40., "nz": 150, "dt": 10., "t_max": 86400, "les_stats_freq": 100.,\
     "cb_min": [296, 296, 297, 0, 0, 7.5, 0, 0, 0, 0, 0, 0, -0.1, 0, 0],\
     "cb_max": [332, 332, 305, 17.5, 17.5, 18,  0.05, 0.02, 2.8, 0.007, 0.004, 0.56, 0, 0.24, 5] ,\
     "t0": 22, "t1": 24, "case_name": "Rico",\
     "ql_max": 0.015, "qr_max": 3e-3, "ql_min": -1e-3, "qr_min": -1e-4, "z_min": -0.1, "z_max": 3.1,\
     "les": "scampify_plots/data_les/Rico/CLIMA_1M_micro/Stats.Rico.nc"}
params["TRMM_LBA"] = \
    {"gw": 3, "dz": 100., "nz": 220, "dt": 10., "t_max": 21600, "les_stats_freq": 100.,\
     "cb_min": [280, 280, 294, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.35, 0, 0],\
     "cb_max": [370, 370, 348, 20, 20, 20, 0.09, 0.028, 2.8, 0.09, 0.064, 0.06, 0, 0.28, 10.5],\
     "t0": 5, "t1": 6, "case_name": "TRMM_LBA",\
     "ql_max": 0.07, "qr_max": 0.07, "ql_min": -1e-3, "qr_min": -1e-3, "z_min": -0.5, "z_max": 15,\
     "les": "scampify_plots/data_les/TRMM_LBA/CLIMA_1M_micro/Stats.TRMM_LBA.nc"}

@pytest.mark.parametrize("case", cases)
@pytest.mark.parametrize("sgs", ['mean', 'quadrature'])
@pytest.mark.parametrize("mode", [True, False])
def test_plot_offline_individual(case, sgs, mode):

    # generate namelists and paramlists
    subprocess.call("rm *.in", shell=True)
    setup = cmn.simulation_setup(case)

    setup['namelist']['microphysics']['rain_model'] = 'clima_1m'
    setup['namelist']['thermodynamics']['sgs'] = sgs

    setup['namelist']['grid']['gw'] = params[case]["gw"]
    setup['namelist']['grid']['dz'] = params[case]["dz"]
    setup['namelist']['grid']['nz'] = params[case]["nz"]

    setup["namelist"]['time_stepping']['dt'] =params[case]["dt"]
    setup['namelist']['time_stepping']['t_max'] = params[case]["t_max"]

    # additional parameters for offline runs
    scampifylist = {}
    scampifylist["offline"] = mode
    scampifylist["les_stats_freq"] = params[case]["les_stats_freq"]
    scampifylist["les_outfile"] = setup["les_outfile"]

    # construct the filepath
    if setup['namelist']['thermodynamics']['sgs'] == 'mean':
        quad = "of"
    if setup['namelist']['thermodynamics']['sgs'] == 'quadrature':
        quad = "on"
    if scampifylist["offline"] == False:
        model = "online"
    if scampifylist["offline"] == True:
        model = "offline"
    fpath = "scampify_plots/data_scm_"+model+"/quad_"+quad

    # remove old simulation files
    subprocess.call("rm -r " + fpath + "/Tests.*", shell=True)

    # set folder for the new simulation files
    setup["namelist"]["output"]["output_root"]=fpath + "/Tests."

    #subprocess.call("python setup.py build_ext --inplace", shell=True, cwd='../')
    if scampifylist["offline"]:
        scampy.main_scampify(setup["namelist"], setup["paramlist"], scampifylist)
    else:
        scampy.main1d(setup["namelist"], setup["paramlist"])

    setup["scm_outfile"] = fpath + "/" + setup["scm_outfile"][1:]
    params[case]["s"+model[0:2]+"_q"+quad] = setup["scm_outfile"]

    # simulation results
    scm_data = Dataset(setup["scm_outfile"], 'r')
    les_data = Dataset(setup["les_outfile"], 'r')
    sim_data = {}
    sim_data["scm_data"] = scm_data
    sim_data["les_data"] = les_data
    sim_data["scampifylist"] = scampifylist
    scm_data_to_plot = cmn.read_scm_data_srs(sim_data["scm_data"])
    les_data_to_plot = cmn.read_les_data_srs(sim_data["les_data"])
    scm_data_to_plot_timesrs = cmn.read_scm_data_timeseries(sim_data["scm_data"])
    les_data_to_plot_timesrs = cmn.read_les_data_timeseries(sim_data["les_data"])

    # plot setup
    cb_min = params[case]["cb_min"]
    cb_max = params[case]["cb_max"]
    t0 = params[case]["t0"]
    t1 = params[case]["t1"]

    try:
        os.mkdir(fpath + "/plots/")
    except:
        pass
    folder = fpath + "/plots/" + params[case]["case_name"] + "/"
    try:
        os.mkdir(folder)
    except:
        pass

    # plot results
    pls.plot_humidities(scm_data_to_plot, les_data_to_plot, params[case], "humidities.pdf", folder=folder)
    pls.plot_cloud_rain_components(scm_data_to_plot, les_data_to_plot, params[case], "cloud_rain_comp.pdf", folder=folder)
    pls.plot_updraft_properties(scm_data_to_plot, les_data_to_plot, params[case], "updraft_properties.pdf", folder=folder)
    pls.plot_timeseries_1D(scm_data_to_plot_timesrs, les_data_to_plot_timesrs, folder=folder)
    pls.plot_timeseries(scm_data_to_plot, les_data_to_plot, params[case], folder=folder)

def test_plot_offline_all():

    # dictionary with simulation results
    rd = {}
    for case in cases:
        rd[case] = {}

        # les specific data
        rd[case]["les"] = {}
        rd[case]["les"]["data"] = cmn.read_les_data_srs(Dataset(params[case]["les"], 'r'))
        rd[case]["les"]["t0"] = int(np.where(rd[case]["les"]["data"]["t"] > params[case]["t0"] * 3600.)[0][0])
        rd[case]["les"]["t1"] = int(np.where(params[case]["t1"] * 3600.0 <= rd[case]["les"]["data"]["t"])[0][0])
        rd[case]["les"]["upd_area"] = rd[case]["les"]["data"]["updraft_fraction"]
        rd[case]["les"]["env_area"] = 1. - rd[case]["les"]["upd_area"]
        rd[case]["les"]["qr_upd"] = np.multiply(rd[case]["les"]["upd_area"], rd[case]["les"]["data"]["updraft_qr"])
        rd[case]["les"]["qr_env"] = np.multiply(rd[case]["les"]["env_area"], rd[case]["les"]["data"]["env_qr"])

        # scm specific data
        for model in ["son_qon", "son_qof", "sof_qon", "sof_qof"]:
            rd[case][model] = {}
            rd[case][model]["data"] = cmn.read_scm_data_srs(Dataset(params[case][model], 'r'))
            rd[case][model]["t0"] = int(np.where(rd[case][model]["data"]["t"] > params[case]["t0"] * 3600.)[0][0])
            rd[case][model]["t1"] = int(np.where(params[case]["t1"] * 3600.0 <= rd[case][model]["data"]["t"])[0][0])
            rd[case][model]["upd_area"] = rd[case][model]["data"]["updraft_area"]
            rd[case][model]["env_area"] = 1. - rd[case][model]["upd_area"]
            rd[case][model]["qr_upd"] = rd[case][model]["data"]["updraft_qr"]
            rd[case][model]["qr_env"] = rd[case][model]["data"]["env_qr"]

        # all data
        for model in ["les", "son_qon", "son_qof", "sof_qon", "sof_qof"]:
            rd[case][model]["ql_mean"] = rd[case][model]["data"]["ql_mean"]
            rd[case][model]["qr_mean"] = rd[case][model]["data"]["qr_mean"]
            rd[case][model]["ql_upd"]  = np.multiply(rd[case][model]["upd_area"], rd[case][model]["data"]["updraft_ql"])
            rd[case][model]["ql_env"]  = np.multiply(rd[case][model]["env_area"], rd[case][model]["data"]["env_ql"])
            rd[case][model]["z_half"]  = rd[case][model]["data"]["z_half"] / 1000.

        # plot results
        pls.plot_les_online_offline(case, rd[case], params[case])
