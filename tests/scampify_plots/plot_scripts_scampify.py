import sys
sys.path.insert(0, "./")
sys.path.insert(0, "../")

import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # To allow plotting when display is off
import matplotlib.pyplot as plt
from matplotlib import ticker

def plot_les_online_offline(case, data, params):
    """
    Plots updraft and environment water realted profiles from
    les, offline scm and online scm

    Input:
    case   - les case name
    data   - dictionary with simulation data
    params - dictionary with simulation parameters
    """
    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(16)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)
    cd1 = {"mean" : "black", "env" : "red", "upd" : "blue"}
    cd2 = {"mean" : "grey",  "env" : "darksalmon", "upd" : "skyblue"}

    plt.subplot(2, 3, 1)
    for sdm in ["mean", "env", "upd"]:
        plt.plot(np.nanmean(data["les"]["ql_"+sdm][:, data["les"]["t0"] : data["les"]["t1"]], axis=1), data["les"]["z_half"], '-', c=cd1[sdm], label='les_'+sdm)
    plt.xlabel("LES ql [g/kg]")
    plt.ylabel("z [km]")
    plt.xlim([params["ql_min"], params["ql_max"]])
    plt.ylim([params["z_min"], params["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 2)
    for sdm in ["mean", "env", "upd"]:
        plt.plot(np.nanmean(data["sof_qon"]["ql_"+sdm][:, data["sof_qon"]["t0"] : data["sof_qon"]["t1"]], axis=1), data["sof_qon"]["z_half"], '-', c=cd2[sdm], label='scm_sgs_'+sdm)
        plt.plot(np.nanmean(data["sof_qof"]["ql_"+sdm][:, data["sof_qof"]["t0"] : data["sof_qof"]["t1"]], axis=1), data["sof_qof"]["z_half"], '-', c=cd1[sdm], label='scm_mean_'+sdm)
    plt.xlabel("SCM-off ql [g/kg]")
    plt.xlim([params["ql_min"], params["ql_max"]])
    plt.ylim([params["z_min"], params["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 3)
    for sdm in ["mean", "env", "upd"]:
        plt.plot(np.nanmean(data["son_qon"]["ql_"+sdm][:, data["son_qon"]["t0"] : data["son_qon"]["t1"]], axis=1), data["son_qon"]["z_half"], '-', c=cd2[sdm], label='scm_sgs_'+sdm)
        plt.plot(np.nanmean(data["son_qof"]["ql_"+sdm][:, data["son_qof"]["t0"] : data["son_qof"]["t1"]], axis=1), data["son_qof"]["z_half"], '-', c=cd1[sdm], label='scm_mean_'+sdm)
    plt.xlabel("SCM-on ql [g/kg]")
    plt.xlim([params["ql_min"], params["ql_max"]])
    plt.ylim([params["z_min"], params["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 4)
    for sdm in ["mean"]:
        plt.plot(np.nanmean(data["les"]["qr_"+sdm][:, data["les"]["t0"] : data["les"]["t1"]], axis=1), data["les"]["z_half"], '-', c=cd1[sdm], label='les_'+sdm)
    plt.xlabel("LES qr [g/kg]")
    plt.ylabel("z [km]")
    plt.xlim([params["qr_min"], params["qr_max"]])
    plt.ylim([params["z_min"], params["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 5)
    for sdm in ["mean"]:
        plt.plot(np.nanmean(data["les"]["qr_"+sdm][:, data["les"]["t0"] : data["les"]["t1"]], axis=1), data["les"]["z_half"], '-', c="silver", label='les_'+sdm)
    for sdm in ["mean", "env", "upd"]:
        plt.plot(np.nanmean(data["sof_qon"]["qr_"+sdm][:, data["sof_qon"]["t0"] : data["sof_qon"]["t1"]], axis=1), data["sof_qon"]["z_half"], '-', c=cd2[sdm], label='scm_sgs_'+sdm)
        plt.plot(np.nanmean(data["sof_qof"]["qr_"+sdm][:, data["sof_qof"]["t0"] : data["sof_qof"]["t1"]], axis=1), data["sof_qof"]["z_half"], '-', c=cd1[sdm], label='scm_mean_'+sdm)
    plt.xlabel("SCM-off qr [g/kg]")
    plt.xlim([params["qr_min"], params["qr_max"]])
    plt.ylim([params["z_min"], params["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 6)
    for sdm in ["mean"]:
        plt.plot(np.nanmean(data["les"]["qr_"+sdm][:, data["les"]["t0"] : data["les"]["t1"]], axis=1), data["les"]["z_half"], '-', c="silver", label='les')
    for sdm in ["mean"]:
        plt.plot(np.nanmean(data["son_qon"]["qr_"+sdm][:, data["son_qon"]["t0"] : data["son_qon"]["t1"]], axis=1), data["son_qon"]["z_half"], '-', c=cd2[sdm], label='scm_sgs')
        plt.plot(np.nanmean(data["son_qof"]["qr_"+sdm][:, data["son_qof"]["t0"] : data["son_qof"]["t1"]], axis=1), data["son_qof"]["z_half"], '-', c=cd1[sdm], label='scm_mean')
    for sdm in ["env", "upd"]:
        plt.plot(np.nanmean(data["son_qon"]["qr_"+sdm][:, data["son_qon"]["t0"] : data["son_qon"]["t1"]], axis=1), data["son_qon"]["z_half"], '-', c=cd2[sdm], label='scm_sgs_'+sdm)
        plt.plot(np.nanmean(data["son_qof"]["qr_"+sdm][:, data["son_qof"]["t0"] : data["son_qof"]["t1"]], axis=1), data["son_qof"]["z_half"], '-', c=cd1[sdm], label='scm_mean_'+sdm)
    plt.xlabel("SCM-on qr [g/kg]")
    plt.xlim([params["qr_min"], params["qr_max"]])
    plt.ylim([params["z_min"], params["z_max"]])
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0.5)
    plt.grid(True)

    plt.tight_layout(rect=[0,0,0.8,1])
    plt.savefig("scampify_plots/subdomains_comp_" + case + ".pdf")
    plt.clf()


def plot_humidities(scm_data, les_data, params, title, folder="scampify_plots/output/"):
    """
    Plots updraft and environment water realted profiles from scm and les

    Input:
    scm_data - scm stats file
    les_data - les stats file
    params   - dictionary with simulation parameters
    title    - name for the created plot
    folder   - folder where to save the created plot
    """
    tmin = params["t0"]
    tmax = params["t1"]

    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin*3600.)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0 <= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax*3600.0 <= les_data["t"]))[0][0])

    scm_data["qv_mean"] = scm_data["qt_mean"]    - scm_data["ql_mean"]
    scm_data["upd_qv"]  = scm_data["updraft_qt"] - scm_data["updraft_ql"]
    scm_data["env_qv"]  = scm_data["env_qt"]     - scm_data["env_ql"]

    les_data["qv_mean"] = les_data["qt_mean"]    - les_data["ql_mean"]
    les_data["upd_qv"]  = les_data["updraft_qt"] - les_data["updraft_ql"]
    les_data["env_qv"]  = les_data["env_qt"]     - les_data["env_ql"]

    les_data["env_qr"] = np.multiply(les_data["env_qr"], (1. - les_data["updraft_fraction"]))
    les_data["updraft_qr"] = np.multiply(les_data["updraft_qr"], les_data["updraft_fraction"])

    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    data_arr = ["qv_mean", "upd_qv", "env_qv",\
                "ql_mean", "updraft_ql", "env_ql",\
                "qr_mean", "updraft_qr", "env_qr"]

    label_arr = ["mean_qv [g/kg]", "updraft qv [g/kg]", "env qv [g/kg]",\
                 "mean ql [g/kg]", "updraft ql [g/kg]", "env ql [g/kg]",\
                 "mean qr [g/kg]", "updraft qr [g/kg]", "env qr [g/kg]"]

    for it in range(9):
        plt.subplot(3, 3, it+1)
        plt.plot(np.nanmean(les_data[data_arr[it]][:, t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='gray',      label='les', lw = 2)
        plt.plot(np.nanmean(scm_data[data_arr[it]][:, t0_scm : t1_scm], axis=1), scm_data["z_half"] / 1000., "-", c="royalblue", label='scm', lw = 4)
        plt.xlabel(label_arr[it])
        plt.ylim([params["z_min"], params["z_max"]])
        plt.grid(True)
        if it in [0,3,6]:
            plt.ylabel("z [km]")
        if it == 8:
            plt.legend()

    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()

def plot_cloud_rain_components(scm_data, les_data, params, title, folder="scampify_plots/output/"):
    """
    Plots updraft and environment water realted profiles from scm and les

    Input:
    scm_data - scm stats file
    les_data - les stats file
    params   - dictionary with simulation parameters
    title    - name for the created plot
    folder   - folder where to save the created plot
    """
    tmin = params["t0"]
    tmax = params["t1"]

    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin*3600.)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0 <= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax*3600.0 <= les_data["t"]))[0][0])

    scm_upd_area = scm_data["updraft_area"]
    scm_env_area = 1. - scm_upd_area
    les_upd_area = les_data["updraft_fraction"]
    les_env_area = 1. - les_upd_area

    scm_upd_ql = np.multiply(scm_upd_area, scm_data["updraft_ql"])
    scm_upd_qr = scm_data["updraft_qr"]
    scm_env_ql = np.multiply(scm_env_area, scm_data["env_ql"])
    scm_env_qr = scm_data["env_qr"]

    les_upd_ql = np.multiply(les_upd_area, les_data["updraft_ql"])
    les_upd_qr = np.multiply(les_upd_area, les_data["updraft_qr"])
    les_env_ql = np.multiply(les_env_area, les_data["env_ql"])
    les_env_qr = np.multiply(les_env_area, les_data["env_qr"])

    les_mean_ql = les_data["ql_mean"]
    les_mean_qr = les_data["qr_mean"]
    scm_mean_ql = scm_data["ql_mean"]
    scm_mean_qr = scm_data["qr_mean"]

    fig = plt.figure(1)
    fig.set_figheight(10)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    data_arr = ["qv_mean", "upd_qv", "env_qv",\
                "ql_mean", "updraft_ql", "env_ql",\
                "qr_mean", "updraft_qr", "env_qr"]

    label_arr = ["mean_qv [g/kg]", "updraft qv [g/kg]", "env qv [g/kg]",\
                 "mean ql [g/kg]", "updraft ql [g/kg]", "env ql [g/kg]",\
                 "mean qr [g/kg]", "updraft qr [g/kg]", "env qr [g/kg]"]

    plt.subplot(1, 2, 1)
    plt.plot(np.nanmean(les_mean_ql[:, t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='grey',       label='les_mean', lw = 2)
    plt.plot(np.nanmean(les_env_ql[:,  t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='darksalmon', label='les_env', lw = 2)
    plt.plot(np.nanmean(les_upd_ql[:,  t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='skyblue',    label='les_upd', lw = 2)

    plt.plot(np.nanmean(scm_mean_ql[:, t0_scm : t1_scm], axis=1), scm_data["z_half"] / 1000., '-', c='black', label='scm_mean', lw = 2)
    plt.plot(np.nanmean(scm_env_ql[:,  t0_scm : t1_scm], axis=1), scm_data["z_half"] / 1000., '-', c='red',   label='scm_env', lw = 2)
    plt.plot(np.nanmean(scm_upd_ql[:,  t0_scm : t1_scm], axis=1), scm_data["z_half"] / 1000., '-', c='blue',  label='scm_upd', lw = 2)

    plt.locator_params(axis='x', nbins=3)
    plt.xlabel("ql [g/kg]")
    plt.ylabel("z [km]")
    plt.xlim([params["ql_min"], params["ql_max"]])
    plt.ylim([params["z_min"], params["z_max"]])
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(np.nanmean(les_mean_qr[:, t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='grey',       label='les_mean', lw = 2)
    plt.plot(np.nanmean(les_env_qr[:,  t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='darksalmon', label='les_env', lw = 2)
    plt.plot(np.nanmean(les_upd_qr[:,  t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='skyblue',    label='les_upd', lw = 2)

    plt.plot(np.nanmean(scm_mean_qr[:, t0_scm : t1_scm], axis=1), scm_data["z_half"] / 1000., '-', c='black', label='scm_mean', lw = 2)
    plt.plot(np.nanmean(scm_env_qr[:,  t0_scm : t1_scm], axis=1), scm_data["z_half"] / 1000., '-', c='red',   label='scm_env', lw = 2)
    plt.plot(np.nanmean(scm_upd_qr[:,  t0_scm : t1_scm], axis=1), scm_data["z_half"] / 1000., '-', c='blue',  label='scm_upd', lw = 2)

    plt.locator_params(axis='x', nbins=3)
    plt.xlabel("qr [g/kg]")
    plt.ylabel("z [km]")
    plt.xlim([params["qr_min"], params["qr_max"]])
    plt.ylim([params["z_min"], params["z_max"]])
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()


def plot_updraft_properties(scm_data, les_data, tmin, tmax, title, folder="scampify_plots/output/"):
    """
    Plots updraft and environment profiles from scm and les

    Input:
    data   - scm stats file
    les    - les stats file
    tmin   - lower bound for time mean
    tmax   - upper bound for time mean
    title  - name for the created plot
    folder - folder where to save the created plot
    """
    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin * 3600.)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin * 3600.)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600. <= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax*3600. <= les_data["t"]))[0][0])

    les_data["massflux"]  = np.multiply(les_data["updraft_fraction"], les_data["updraft_w"])
    scm_data["massflux"]  = np.multiply(scm_data["updraft_area"],     scm_data["updraft_w"])

    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    les_data_arr = ["updraft_fraction", "updraft_w",       "massflux",\
                    "thetali_mean",     "updraft_thetali", "env_thetali"]

    scm_data_arr = ["updraft_area",     "updraft_w",       "massflux",\
                    "thetal_mean",      "updraft_thetal",  "env_thetal"]

    label_arr = ["updraft fraction", "updraft w [m/s]",    "massflux [kg/m^2/s]",\
                 "thetal mean [K]",  "updraft thetal [K]", "env thetal [K]"]

    for it in range(6):
      plt.subplot(2,3,it+1)
      plt.plot(np.nanmean(les_data[les_data_arr[it]][:, t0_les : t1_les], axis=1), les_data["z_half"]/1000., '-', c='gray',      label='les', lw = 4)
      plt.plot(np.nanmean(scm_data[scm_data_arr[it]][:, t0_scm : t1_scm], axis=1), scm_data["z_half"]/1000., "-", c="royalblue", label='scm', lw = 2)
      plt.xlabel(label_arr[it])
      plt.grid(True)
      if it in [0,3]:
          plt.ylabel("z [km]")

    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()

def plot_timeseries_1D(scm_data, les_data, folder="scampify_plots/output/"):
    """
    Plots timeseries from Scampy

    Input:
    scm_data - scm stats file
    les_data - les stats file
    tmin     - lower bound for time mean
    tmax     - upper bound for time mean
    folder   - folder where to save the created plot
    """
    # cloud timeseries
    plot_scm_y = [scm_data["lwp_mean"],\
                  scm_data["cloud_cover_mean"],\
                  scm_data["rwp_mean"],\
                  scm_data["cloud_top_mean"]/1e3, scm_data["cloud_base_mean"]/1e3]
    plot_les_y = [les_data["lwp_mean"],\
                  les_data["cloud_cover_mean"],\
                  les_data["rwp_mean"],\
                  les_data["cloud_top_mean"]/1e3, les_data["cloud_base_mean"]/1e3]
    y_lab      = ['lwp', 'cloud_cover', 'rwp', 'CB, CT [km]']

    fig = plt.figure(1)
    for plot_it in range(4):
      plt.subplot(2,2,plot_it+1)
      plt.plot(les_data["t"][1:]/3600., plot_les_y[plot_it][1:], '-', color="gray", label="LES", lw=3)
      plt.plot(scm_data["t"][1:]/3600., plot_scm_y[plot_it][1:], '-', color="b",    label="SCM", lw=3)
      if plot_it == 3:
        plt.plot(les_data["t"][1:]/3600., plot_les_y[4][1:], '-', color="gray", lw=3)
        plt.plot(scm_data["t"][1:]/3600., plot_scm_y[4][1:], '-', color="b", lw=3)
      plt.legend()
      plt.grid(True)
      plt.xlim([0, scm_data["t"][-1]/3600.])
      plt.xlabel('time [h]')
      plt.ylabel(y_lab[plot_it])
    plt.tight_layout()
    plt.savefig(folder + "timeseries_cloud_properties.pdf")
    plt.clf()

def plot_timeseries(scm_data, les_data, cb_min, cb_max, folder="scampify_plots/output/"):
    """
    Plots the time series of Scampy simulations

    Input:
    data   - scm stats file
    les    - les stats file
    tmin   - lower bound for time mean
    tmax   - upper bound for time mean
    folder - folder where to save the created plot
    """

    scm_z_half = scm_data["z_half"] / 1000.
    scm_time   = scm_data["t"] / 3600.

    les_z_half = les_data["z_half"] / 1000.
    les_time   = les_data["t"] / 3600.

    les_data["env_qr"] = np.multiply(les_data["env_qr"], (1. - les_data["updraft_fraction"]))
    les_data["updraft_qr"] = np.multiply(les_data["updraft_qr"], les_data["updraft_fraction"])

    les_vars  = ["thetali_mean", "env_thetali",      "updraft_thetali",\
                 "qt_mean",      "env_qt",           "updraft_qt",\
                 "ql_mean",      "env_ql",           "updraft_ql",\
                 "qr_mean",      "env_qr",           "updraft_qr",\
                 "env_w",        "updraft_fraction", "updraft_w"]

    scm_vars  = ["thetal_mean",  "env_thetal",       "updraft_thetal",\
                 "qt_mean",      "env_qt",           "updraft_qt",\
                 "ql_mean",      "env_ql",           "updraft_ql",\
                 "qr_mean",      "env_qr",           "updraft_qr",\
                 "env_w",        "updraft_area",     "updraft_w"]

    labels    = ["mean thl [K]",   "env thl [K]",   "updr thl [K]",\
                 "mean qt [g/kg]", "env qt [g/kg]", "updr qt [g/kg]",\
                 "mean ql [g/kg]", "env ql [g/kg]", "updr ql [g/kg]",\
                 "mean qr [g/kg]", "env qr [g/kg]", "updr qr [g/kg",\
                 "env w [m/s]",    "updr area [%]", "updr w [m/s]"]

    fig_name =  ["contour_thl_mean", "contour_env_thl",  "contour_upd_thl",\
                 "contour_qt_mean",  "contour_env_qt",   "contour_upd_qt",\
                 "contour_ql_mean",  "contour_env_ql",   "contour_upd_ql",\
                 "contour_qr_mean",  "contour_env_qr",   "contour_upd_qr",\
                 "contour_env_w",    "contour_upd_area", "contour_upd_w"]

    # initial condition for env_thetal starts with zeros
    scm_data["env_thetal"][:, 0] = scm_data["env_thetal"][:, 1]

    for plot_it in range(len(labels)):

        fig = plt.figure(fig_name[plot_it])
        fig.set_figheight(12)
        fig.set_figwidth(14)
        mpl.rcParams.update({'font.size': 18})
        mpl.rc('lines', linewidth=4, markersize=10)

        scm_field = scm_data[scm_vars[plot_it]]
        les_field = les_data[les_vars[plot_it]]

        if ("updraft" in scm_vars[plot_it]):
            a_scm = scm_data["updraft_area"]
            scm_field[np.where(a_scm==0.0)] = np.nan
            scm_field[np.where(np.isnan(a_scm))] = np.nan

        if ("updraft" in les_vars[plot_it]):
            a_les = les_data["updraft_fraction"]
            les_field[np.where(a_les==0.0)] = np.nan
            les_field[np.where(np.isnan(a_les))] = np.nan

        levels = np.linspace(cb_min[plot_it], cb_max[plot_it], 11)
        cmap = "RdBu_r"

        plt.subplot(211)
        cntrf = plt.contourf(les_time, les_z_half, les_field, cmap=cmap, levels=levels, vmin=cb_min[plot_it], vmax=cb_max[plot_it])
        cbar = plt.colorbar(cntrf)
        cbar.set_label(labels[plot_it])
        plt.ylabel('height [km]')
        plt.ylim([0, np.max(scm_data["z_half"]/1000.)])
        plt.grid(True)

        plt.subplot(212)
        cntrf = plt.contourf(scm_time, scm_z_half, scm_field, cmap=cmap, levels=levels, vmin=cb_min[plot_it], vmax=cb_max[plot_it])
        cbar = plt.colorbar(cntrf)
        cbar.set_label(labels[plot_it])
        plt.xlabel('time [h]')
        plt.ylabel('height [km]')
        plt.ylim([0, np.max(scm_data["z_half"]/1000.)])
        plt.grid(True)

        plt.tight_layout()
        plt.savefig(folder + fig_name[plot_it]+".pdf")
        plt.clf()
        plt.close()
