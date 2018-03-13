import sys
sys.path.insert(0, "./")
sys.path.insert(0, "./tests")

import os
import subprocess
import json
import warnings

from netCDF4 import Dataset
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

def simulation_setup(case):
    """
    generate namelist and paramlist files for scampy
    choose the name of the output folder
    """
    # Filter annoying Cython warnings that serve no good purpose.
    # see https://stackoverflow.com/questions/40845304/runtimewarning-numpy-dtype-size-changed-may-indicate-binary-incompatibility
    warnings.filterwarnings("ignore", message="numpy.dtype size changed")
    warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

    # simulation related parameters
    os.system("python generate_namelist.py " + case)
    file_case = open(case + '.in').read()
    # turbulence related parameters
    os.system("python generate_paramlist.py " +  case)
    file_params = open('paramlist_' + case + '.in').read()

    namelist  = json.loads(file_case)
    paramlist = json.loads(file_params)

    namelist['output']['output_root'] = "./Tests."
    namelist['meta']['uuid'] = case
    # TODO - copied from NetCDFIO
    # ugly way to know the name of the folder where the data is saved
    uuid = str(namelist['meta']['uuid'])
    outpath = str(
        os.path.join(
            namelist['output']['output_root'] +
            'Output.' +
            namelist['meta']['simname'] +
            '.' +
            uuid[len(uuid )-5:len(uuid)]
        )
    )
    outfile = outpath + "/stats/Stats." + case + ".nc"

    res = {"namelist"  : namelist,
           "paramlist" : paramlist,
           "outfile"   : outfile}
    return res


def removing_files():
    """
    remove the folder with netcdf files from tests 
    """
    #TODO - think of something better
    cmd = "rm -r Tests.Output.*"
    subprocess.call(cmd , shell=True)


def read_data(sim_data, ref_data):
    """
    Read in the data from netcdf files into a dictionary that can be used for quicklook plots
    
    dict[var] = [reference var[t=0], reference var[t=end], current var[t=0], current var[t=end]
    """
    variables = ["temperature_mean", "thetal_mean", "qt_mean", "ql_mean", "buoyancy_mean", "u_mean", "v_mean", "tke_mean",\
                 "updraft_buoyancy", "updraft_area", "env_qt", "updraft_qt", "env_ql", "updraft_ql", "updraft_w", "env_w"]
    time = [0,        -1,       0,        -1]
    sim  = [ref_data, ref_data, sim_data, sim_data]

    # read the data (both simulation and reference)
    data_to_plot = {"z_half" : np.array(sim_data["profiles/z_half"][:])}
    for var in variables:
        data_to_plot[var] = []
        for it in xrange(4):
            if ("buoyancy" in var):
                data_to_plot[var].append(np.array(sim[it]["profiles/" + var][time[it], :]) * 10000) #cm2/s3
            elif ("qt" in var or "ql" in var):
                data_to_plot[var].append(np.array(sim[it]["profiles/" + var][time[it], :]) * 1000)  #g/kg
            elif ("p0" in var):
                data_to_plot[var].append(np.array(sim[it]["reference/" + var][time[it], :]) * 1000)  #g/kg
            else:
                data_to_plot[var].append(np.array(sim[it]["profiles/" + var][time[it], :]))

    return data_to_plot


def read_rad_data(sim_data, ref_data):
    """
    Read in the radiation forcing data from netcdf files into a dictionary that can be used for quicklook plots
    
    dict[var] = [reference var[t=0], reference var[t=end], current var[t=0], current var[t=end]
    """
    variables = ["rad_flux", "rad_dTdt"]
    time = [0,        -1,       0,        -1]
    sim  = [ref_data, ref_data, sim_data, sim_data]

    # read the data (both simulation and reference)
    rad_data = {"z" : np.array(sim_data["profiles/z"][:])}
    for var in variables:
        rad_data[var] = []
        for it in xrange(4):
            if ("rad_dTdt" in var):
                rad_data[var].append(np.array(sim[it]["profiles/" + var][time[it], :]) * 60 * 60 * 24) # K/day
            else:
                rad_data[var].append(np.array(sim[it]["profiles/" + var][time[it], :]))

    return rad_data


def plot_mean(data, title, folder="tests/output/"):
    """
    Plots mean profiles from Scampy (current test run and reference simulation - EDMF_BulkSteady)
    """
    # customize defaults
    mpl.rc('lines', linewidth=3, markersize=8)

    plt.figure(1, figsize=(18,14))
    mpl.rc('lines', linewidth=4, markersize=10)
    mpl.rcParams.update({'font.size': 18})
    plots = []
    # iteration over plots
    x_lab  = ['T [K]',                   'THL [K]',           'buoyancy [cm2/s3]',   'QT [g/kg]',     'QL [g/kg]',     'TKE']
    plot_x = [ data["temperature_mean"], data["thetal_mean"], data["buoyancy_mean"], data["qt_mean"], data["ql_mean"], data["tke_mean"]]
    # iteration over currnt vs reference simulation, t=0, t=-1
    color  = ["palegreen", "forestgreen", "gold",    "orangered"]
    label  = ["ref ini",   "ref end",     "sim ini", "sim end"]

    for plot_it in range(6):
        plots.append(plt.subplot(2,3,plot_it+1))
                               #(rows, columns, number)
        plots[plot_it].set_xlabel(x_lab[plot_it])
        plots[plot_it].set_ylabel('z [m]')
        plots[plot_it].set_ylim([0, data["z_half"][-1] + (data["z_half"][1] - data["z_half"][0]) * 0.5])
        plots[plot_it].grid(True)
        #for it in xrange(2,4,1):  # only current version
        for it in xrange(4):
            plots[plot_it].plot(plot_x[plot_it][it], data["z_half"], '.-', color=color[it], label=label[it])

    plots[0].legend(loc='upper right')
    #plots[2].set_xlim([1, 10])
    #plots[3].set_xlim([-0.1, 0.5])
    #plots[4].set_xlim([-50, 350])
    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()


def plot_drafts(data, title, folder="tests/output/"):
    """
    Plots updraft and environment profiles from Scampy (current and reference simulation - EDMF_Bulksteady)
    """
    # customize defaults
    mpl.rc('lines', linewidth=3, markersize=8)

    plt.figure(1, figsize=(18,14))
    mpl.rc('lines', linewidth=4, markersize=10)
    mpl.rcParams.update({'font.size': 18})
    plots = []
    # iteration over plots
    x_lab    = ["QT [g/kg]",        "QL [g/kg]",         "w [m/s]",         "updraft buoyancy [cm2/s3]", "updraft area [%]"]
    plot_upd = [data["updraft_qt"], data["updraft_ql"],  data["updraft_w"], data["updraft_buoyancy"],    data["updraft_area"]]
    plot_env = [data["env_qt"],     data["env_ql"],      data["env_w"]]
    plot_mean= [data["qt_mean"],    data["ql_mean"]]
    # iteration over current vs reference simulation
    color_mean= ["plum",        "purple"]
    color_env = ["lightsalmon", "red"]
    color_upd = ["deepskyblue", "blue"]
    label_mean= ["mean ref",    "mean sim"]
    label_env = ["env ref",     "env sim"]
    label_upd = ["upd ref",     "upd sim"]

    for plot_it in xrange(5):
        plots.append(plt.subplot(3,2,plot_it+1))
                               #(rows, columns, number)
        plots[plot_it].set_xlabel(x_lab[plot_it])
        plots[plot_it].set_ylabel('z [m]')
        plots[plot_it].set_ylim([0, data["z_half"][-1] + (data["z_half"][1] - data["z_half"][0]) * 0.5])
        plots[plot_it].grid(True)
        #for it in xrange(1,2,1): # only current version
        for it in xrange(2):
            #plot updrafts
            if (plot_it < 4):
                plots[plot_it].plot(plot_upd[plot_it][2 * it + 1], data["z_half"], ".-", color=color_upd[it], label=label_upd[it])
            if (plot_it == 4):
                plots[plot_it].plot(plot_upd[plot_it][2 * it + 1] * 100, data["z_half"], ".-", color=color_upd[it], label=label_upd[it])
            if (plot_it in [0, 1, 2]):
                # plot environment
                plots[plot_it].plot(plot_env[plot_it ][2 * it + 1], data["z_half"], ".-", color=color_env[it], label=label_env[it])
            if (plot_it in [0, 1]):
                # plot mean
                plots[plot_it].plot(plot_mean[plot_it][2 * it + 1], data["z_half"], ".-", color=color_mean[it], label=label_mean[it])


    plots[1].legend(loc='lower right')
    #plots[0].set_xlim([1, 10])
    #plots[5].set_xlim([-50, 350])
    #plots[5].set_xlim([-0.1, 0.5])
    plt.savefig(folder + title)
    plt.clf()
