import sys
sys.path.insert(0, "./")
sys.path.insert(0, "./tests")

import os
import subprocess
import json
import warnings

from netCDF4 import Dataset
import numpy as np
import pprint as pp

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


def read_data_avg(sim_data, n_steps):
    """
    Read in the data from netcdf file into a dictionary that can be used for quicklook plots
    """
    variables = ["temperature_mean", "thetal_mean", "qt_mean", "ql_mean", "buoyancy_mean", "u_mean", "v_mean", "tke_mean",\
                 "updraft_buoyancy", "updraft_area", "env_qt", "updraft_qt", "env_ql", "updraft_ql", "updraft_w", "env_w"]

    # read the data
    data_to_plot = {"z_half" : np.array(sim_data["profiles/z_half"][:])}

    time = [0, -1]
    for var in variables:
        data_to_plot[var] = []
        for it in xrange(2):
            if ("buoyancy" in var):
                data_to_plot[var].append(np.array(sim_data["profiles/" + var][time[it], :]) * 10000) #cm2/s3
            elif ("qt" in var or "ql" in var):
                data_to_plot[var].append(np.array(sim_data["profiles/" + var][time[it], :]) * 1000)  #g/kg
            elif ("p0" in var):
                data_to_plot[var].append(np.array(sim_data["reference/" + var][time[it], :]) * 100)  #hPa
            else:
                data_to_plot[var].append(np.array(sim_data["profiles/" + var][time[it], :]))

    # add averaging over last n_steps timesteps
    if(n_steps > 0):
        for var in variables:
            for time_it in xrange(-2, -1*n_steps-1, -1):
                if ("buoyancy" in var):
                    data_to_plot[var][1] += np.array(sim_data["profiles/" + var][time_it, :]) * 10000  #cm2/s3
                elif ("qt" in var or "ql" in var):
                    data_to_plot[var][1] += np.array(sim_data["profiles/" + var][time_it, :]) * 1000   #g/kg
                elif ("p0" in var):
                    data_to_plot[var][1] += np.array(sim_data["reference/" + var][time_it, :]) * 100   #hPa
                else:
                    data_to_plot[var][1] += np.array(sim_data["profiles/" + var][time_it, :])

            data_to_plot[var][1] /= n_steps

    return data_to_plot


def read_rad_data_avg(sim_data, n_steps):
    """
    Read in the radiation forcing data from netcdf files into a dictionary that can be used for quicklook plots
    """
    variables = ["rad_flux", "rad_dTdt"]

    time = [0, -1]
    rad_data = {"z" : np.array(sim_data["profiles/z"][:])}
    for var in variables:
        rad_data[var] = []
        for it in xrange(2):
            if ("rad_dTdt" in var):
                rad_data[var].append(np.array(sim_data["profiles/" + var][time[it], :]) * 60 * 60 * 24) # K/day
            else:
                rad_data[var].append(np.array(sim_data["profiles/" + var][time[it], :]))

    # add averaging over last n_steps timesteps
    if(n_steps > 0):
        for var in variables:
            for time_it in xrange(-2, -1*n_steps-1, -1):
                if ("rad_dTdt" in var):
                    rad_data[var][1] += np.array(sim_data["profiles/" + var][time_it, :] * 60 * 60 * 24) # K/day
                else:
                    rad_data[var][1] += np.array(sim_data["profiles/" + var][time_it, :])

            rad_data[var][1] /= n_steps

    return rad_data

def read_data_srs(sim_data):
    """
    Read in the data from netcdf file into a dictionary that can be used for quicklook timeseries plots
    """
    variables = ["temperature_mean", "thetal_mean", "qt_mean", "ql_mean", "buoyancy_mean", "u_mean", "v_mean", "tke_mean",\
                 "updraft_buoyancy", "updraft_area", "env_qt", "updraft_qt", "env_ql", "updraft_ql", "updraft_w", "env_w"]

    # read the data
    data_to_plot = {"z_half" : np.array(sim_data["profiles/z_half"][:]), "t" : np.array(sim_data["profiles/t"][:])}

    for var in variables:
        data_to_plot[var] = []
        if ("buoyancy" in var):
            data_to_plot[var] = np.array(sim_data["profiles/"  + var][:, :]) * 10000 #cm2/s3
        elif ("qt" in var or "ql" in var):
            data_to_plot[var] = np.array(sim_data["profiles/"  + var][:, :]) * 1000  #g/kg
        elif ("p0" in var):
            data_to_plot[var] = np.array(sim_data["reference/" + var][:, :]) * 100   #hPa
        else:
            data_to_plot[var] = np.array(sim_data["profiles/"  + var][:, :])

    return data_to_plot

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
    color  = ["palegreen", "forestgreen"]
    label  = ["ini", "end"]

    for plot_it in range(6):
        plots.append(plt.subplot(2,3,plot_it+1))
                               #(rows, columns, number)
        plots[plot_it].set_xlabel(x_lab[plot_it])
        plots[plot_it].set_ylabel('z [m]')
        plots[plot_it].set_ylim([0, data["z_half"][-1] + (data["z_half"][1] - data["z_half"][0]) * 0.5])
        plots[plot_it].grid(True)
        for it in xrange(2):
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
    color_mean= "purple"
    color_env = "red"
    color_upd = "blue"
    label_mean= "mean"
    label_env = "env"
    label_upd = "upd"

    for plot_it in xrange(5):
        plots.append(plt.subplot(3,2,plot_it+1))
                               #(rows, columns, number)
        plots[plot_it].set_xlabel(x_lab[plot_it])
        plots[plot_it].set_ylabel('z [m]')
        plots[plot_it].set_ylim([0, data["z_half"][-1] + (data["z_half"][1] - data["z_half"][0]) * 0.5])
        plots[plot_it].grid(True)
        #plot updrafts
        if (plot_it < 4):
            plots[plot_it].plot(plot_upd[plot_it][1], data["z_half"], ".-", color=color_upd, label=label_upd)
        if (plot_it == 4):
            plots[plot_it].plot(plot_upd[plot_it][1] * 100, data["z_half"], ".-", color=color_upd, label=label_upd)
        if (plot_it in [0, 1, 2]):
            # plot environment
            plots[plot_it].plot(plot_env[plot_it ][1], data["z_half"], ".-", color=color_env, label=label_env)
        if (plot_it in [0, 1]):
            # plot mean
            plots[plot_it].plot(plot_mean[plot_it][1], data["z_half"], ".-", color=color_mean, label=label_mean)


    plots[1].legend(loc='upper right')
    #plots[0].set_xlim([1, 10])
    #plots[5].set_xlim([-50, 350])
    #plots[5].set_xlim([-0.1, 0.5])
    plt.savefig(folder + title)
    plt.clf()

def plot_timeseries(plt_data, case):

    output_folder="tests/output/"
    plt.figure(1, figsize=(20,14))
    mpl.rcParams.update({'font.size': 18})

    z_half  = plt_data["z_half"]
    time    = plt_data["t"] / 60. / 60.

    mean_ql  = np.transpose(plt_data["ql_mean"])
    mean_qt  = np.transpose(plt_data["qt_mean"])
    mean_qv  = mean_qt - mean_ql
    mean_tke = np.transpose(plt_data["tke_mean"])

    updr_buo  = np.transpose(plt_data["updraft_buoyancy"])
    updr_qt   = np.transpose(plt_data["updraft_qt"])
    updr_ql   = np.transpose(plt_data["updraft_ql"])
    updr_qv   = updr_qt - updr_ql
    updr_w    = np.transpose(plt_data["updraft_w"])
    updr_area = np.transpose(plt_data["updraft_area"])

    env_qt = np.transpose(plt_data["env_qt"])
    env_ql = np.transpose(plt_data["env_ql"])
    env_qv = env_qt - env_ql
    env_w  = np.transpose(plt_data["env_w"])

    x_lab = ["mean tke", "mean qt", "mean_qv", "mean_ql"]
    fig = plt.figure(1)
    ax = []
    for plot_it in range(4):
        ax.append(fig.add_subplot(2,2,plot_it+1))
                               #(rows, columns, number)
        ax[plot_it].set_xlabel('t [hrs]')
        ax[plot_it].set_ylabel('z [m]')

    plot0 = ax[0].pcolormesh(time, z_half, mean_tke)
    fig.colorbar(plot0, ax=ax[0], label='mean tke [m2/s2]')
    plot1 = ax[1].pcolormesh(time, z_half, mean_qt)
    fig.colorbar(plot1, ax=ax[1], label='mean qt [g/kg]')
    plot2 = ax[2].pcolormesh(time, z_half, mean_qv)
    fig.colorbar(plot2, ax=ax[2], label='mean qv [g/kg]')
    plot3 = ax[3].pcolormesh(time, z_half, mean_ql)
    fig.colorbar(plot3, ax=ax[3], label='mean ql [g/kg]')
    plt.savefig(output_folder + case + "_timeseries_mean.png")
    plt.clf()

    x_lab = ["env w", "env qt", "env qv", "env ql"]
    fig = plt.figure(1)
    ax = []
    for plot_it in range(4):
        ax.append(fig.add_subplot(2,2,plot_it+1))
                               #(rows, columns, number)
        ax[plot_it].set_xlabel('t [hrs]')
        ax[plot_it].set_ylabel('z [m]')

    plot0 = ax[0].pcolormesh(time, z_half, env_w)
    fig.colorbar(plot0, ax=ax[0], label='env w [m/2]')
    plot1 = ax[1].pcolormesh(time, z_half, env_qt)
    fig.colorbar(plot1, ax=ax[1], label='env qt [g/kg]')
    plot2 = ax[2].pcolormesh(time, z_half, env_qv)
    fig.colorbar(plot2, ax=ax[2], label='env qv [g/kg]')
    plot3 = ax[3].pcolormesh(time, z_half, env_ql)
    fig.colorbar(plot3, ax=ax[3], label='env ql [g/kg]')
    plt.savefig(output_folder + case + "_timeseries_env.png")
    plt.clf()

    x_lab = ["updr buo", "updr area", "updr w", "updr ql"]
    fig = plt.figure(1)
    ax = []
    for plot_it in range(4):
        ax.append(fig.add_subplot(2,2,plot_it+1))
                               #(rows, columns, number)
        ax[plot_it].set_xlabel('t [hrs]')
        ax[plot_it].set_ylabel('z [m]')

    plot0 = ax[0].pcolormesh(time, z_half, updr_buo)
    fig.colorbar(plot0, ax=ax[0], label='updr buo [cm2/s3]')
    plot1 = ax[1].pcolormesh(time, z_half, updr_area)
    fig.colorbar(plot1, ax=ax[1], label='updr area ')
    plot2 = ax[2].pcolormesh(time, z_half, updr_w)
    fig.colorbar(plot2, ax=ax[2], label='updr w [m/s]')
    plot3 = ax[3].pcolormesh(time, z_half, updr_ql)
    fig.colorbar(plot3, ax=ax[3], label='updr ql [g/kg]')
    plt.savefig(output_folder + case + "_timeseries_updr.png")
    plt.clf()

