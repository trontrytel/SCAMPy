import sys
sys.path.insert(0, "./")
sys.path.insert(0, "../")

import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # To allow plotting when display is off
import matplotlib.pyplot as plt
from matplotlib import ticker

def plot_humidities(scm_data, les_data, tmin, tmax, title, folder="scampify_plots/output/"):
    """
    Plots updraft and environment water realted profiles from scm and les

    Input:
    scm_data - scm stats file
    les_data - les stats file
    tmin     - lower bound for time mean
    tmax     - upper bound for time mean
    title    - name for the created plot
    folder   - folder where to save the created plot
    """
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
        plt.grid(True)
        if it in [1,4,7]:
            plt.ylabel("z [km]")

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
      if it in [1,4]:
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
    plot_scm_y = [scm_data["cloud_cover_mean"], scm_data["lwp_mean"], scm_data["rwp_mean"],\
                  scm_data["lhf"], scm_data["shf"], scm_data["cloud_top_mean"], scm_data["cloud_base_mean"]]

    plot_les_y = [les_data["cloud_cover_mean"], les_data["lwp_mean"], les_data["rwp_mean"],\
                  les_data["lhf"], les_data["shf"], les_data["cloud_top_mean"], les_data["cloud_base_mean"]]

    y_lab      = ['cloud cover', 'lwp', 'rwp', 'lhf', 'shf', 'CB, CT [km]']

    fig_name   = ['cloud_cover', 'liquid_water_path', 'rain_water_path',\
                  'latent_heat_flux',  'sensible_heat_flux', 'cloud_base_top']

    # iteration over plots
    plots = []
    for plot_it in range(6):

        fig = plt.figure(1)
        plt.xlabel('time [h]')
        plt.ylabel(y_lab[plot_it])
        plt.xlim([0, scm_data["t"][-1]/3600.])
        plt.grid(True)

        if plot_it < 5:
            plt.plot(les_data["t"][1:]/3600., plot_les_y[plot_it][1:], '-', color="gray", lw = 4, label="LES")
            plt.plot(scm_data["t"][1:]/3600., plot_scm_y[plot_it][1:], '-', color="b",    lw = 2, label="SCM")
            plt.legend()

        else:
            plt.plot(les_data["t"][1:]/3600., les_data["cloud_base_mean"][1:]/1000., '-', color="gray",      label="CB_les", lw = 4)
            plt.plot(les_data["t"][1:]/3600., les_data["cloud_top_mean"][1:]/1000.,  '-', color="gray",      label="CT_les", lw = 4)
            plt.plot(scm_data["t"][1:]/3600., scm_data["cloud_base_mean"][1:]/1000., '-', color="crimson",   label="CB_scm", lw = 2)
            plt.plot(scm_data["t"][1:]/3600., scm_data["cloud_top_mean"][1:]/1000.,  '-', color="royalblue", label="CT_scm", lw = 2)

        plt.tight_layout()
        plt.savefig(folder + fig_name[plot_it]+".pdf")
        plt.clf()

def plot_timeseries(scm_data, les_data, folder="scampify_plots/output/"):
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

        plt.subplot(211)
        plt.contourf(les_time, les_z_half, les_field, cmap='RdBu_r')
        plt.colorbar()
        plt.ylabel('height [km]')
        #plt.ylim([0, np.max(scm_data["z_half"])])

        plt.subplot(212)
        plt.contourf(scm_time, scm_z_half, scm_field, cmap='RdBu_r')
        plt.colorbar()
        plt.xlabel('time [h]')
        plt.ylabel('height [km]')
        #plt.ylim([0, np.max(scm_data["z_half"])])

        plt.tight_layout()
        plt.savefig(folder + fig_name[plot_it]+".pdf")
        plt.clf()
        plt.close()
