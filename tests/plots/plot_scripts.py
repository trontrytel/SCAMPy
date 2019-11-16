import sys
sys.path.insert(0, "./")
sys.path.insert(0, "../")

import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # To allow plotting when display is off
import matplotlib.pyplot as plt
from matplotlib import ticker

def plot_mean_prof(scm_data, les_data, tmin, tmax, folder="plots/output/"):
    """
    Plots mean profiles from Scampy
    Input:
    scm_data - scm stats file
    les_data - les stats file
    tmin     - lower bound for time mean
    tmax     - upper bound for time mean
    folder   - folder where to save the created plot
    """
    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.0)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0<= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax<= les_data["t"]))[0][0])

    qv_mean_scm = np.array(scm_data["qt_mean"]) - np.array(scm_data["ql_mean"])
    qv_mean_les = np.array(les_data["qt_mean"]) - np.array(les_data["ql_mean"])

    x_labels  =  [r'$q_{t} mean [\mathrm{g/kg}]$',
                  r'$q_{l} mean [\mathrm{g/kg}]$',
                  r'$q_{i} mean [\mathrm{g/kg}]$',
                  r'$q_{r} mean [\mathrm{g/kg}]$',
                  r'$q_{v} mean [\mathrm{g/kg}]$',
                  r'$q_{s} mean [\mathrm{g/kg}]$',
                  r'$\theta_{l} [\mathrm{K}]$',
                  r'$TKE [\mathrm{m^2/s^2}]$',
                  'u [m/s]',
                  'v [m/s]',
                  r'$\bar{w}_{upd} [\mathrm{m/s}]$',
                  r'$\bar{b}_{upd} [\mathrm{m/s^2}]$',
                  "updraft area [%]",
                  r'$\bar{q}_{l,upd} [\mathrm{g/kg}]$',
                  r'$\bar{q}_{i,upd} [\mathrm{g/kg}]$',
                  r'$\bar{q}_{r,upd} [\mathrm{g/kg}]$',
                  r'$\bar{q}_{s,upd} [\mathrm{g/kg}]$',
                  r'$\bar{q}_{l,env} [\mathrm{g/kg}]$',
                  r'$\bar{q}_{i,env} [\mathrm{g/kg}]$',
                  r'$\bar{q}_{r,env} [\mathrm{g/kg}]$',
                  r'$\bar{q}_{s,env} [\mathrm{g/kg}]$']

    fig_name  =  ["mean_qt", "mean_ql", "mean_qi", "mean_qr", "mean_qv", "mean_qs",\
                  "mean_thetal", "mean_TKE", "mean_u", "mean_v",\
                  "updraft_w", "updraft_buoyancy", "updraft_area",\
                  "updraft_ql", "updraft_qi", "updraft_qr", "updraft_qs",\
                  "env_ql", "env_qi", "env_qr", "env_qs"]

    plot_x_scm = [scm_data["qt_mean"], scm_data["ql_mean"], scm_data["qi_mean"],\
                  scm_data["qr_mean"], qv_mean_scm, scm_data["qs_mean"],\
                  scm_data["thetal_mean"], scm_data["tke_mean"],\
                  scm_data["u_mean"], scm_data["v_mean"], scm_data["updraft_w"],\
                  scm_data["updraft_buoyancy"], scm_data["updraft_area"],\
                  scm_data["updraft_ql"], scm_data["updraft_qi"],\
                  scm_data["updraft_qr"], scm_data["updraft_qs"],\
                  scm_data["env_ql"], scm_data["env_qi"],\
                  scm_data["env_qr"], scm_data["env_qs"]]

    plot_x_les = [les_data["qt_mean"], les_data["ql_mean"], les_data["qi_mean"],\
                  les_data["qr_mean"], qv_mean_les, les_data["qs_mean"],\
                  les_data["thetali_mean"], les_data["tke_mean"],\
                  les_data["u_translational_mean"], les_data["v_translational_mean"],\
                  les_data["updraft_w"], les_data["updraft_buoyancy"],\
                  les_data["updraft_fraction"],\
                  les_data["updraft_ql"], les_data["updraft_qi"],\
                  les_data["updraft_qr"], les_data["updraft_qs"],\
                  les_data["env_ql"], les_data["env_qi"],\
                  les_data["env_qr"], les_data["env_qs"]]

    plots = []
    for plot_it in range(len(x_labels)):
        fig = plt.figure(fig_name[plot_it])
        plt.xlabel(x_labels[plot_it])
        plt.ylabel('height [km]')
        plt.ylim([0,\
                  scm_data["z_half"][-1]/1000.0 +\
                  (scm_data["z_half"][1]/1000.0 - scm_data["z_half"][0]/1000.0) * 0.5\
                 ])
        plt.grid(True)
        plt.plot(np.nanmean(plot_x_les[plot_it][:, t0_les:t1_les],axis=1), les_data["z_half"],     '-', color='k', label='les', linewidth = 2)
        plt.plot(np.nanmean(plot_x_scm[plot_it][:, t0_scm:t1_scm],axis=1), scm_data["z_half"]/1e3, '-', color = '#157CC7', label='scm', linewidth = 2)

        plt.legend()
        plt.tight_layout()
        plt.savefig(folder + "profile_" + fig_name[plot_it]+".pdf")
        plt.clf()

def plot_closures(scm_data, les_data, tmin, tmax, title, folder="plots/output/"):
    """
    Plots updraft and environment profiles from Scampy
    Input:
    scm_data - scm stats file
    les_data - les stats file
    tmin     - lower bound for time mean
    tmax     - upper bound for time mean
    title    - name for the created plot
    folder   - folder where to save the created plot
    """
    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.0)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0<= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax<= les_data["t"]))[0][0])

    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    scm_vars = [np.nanmean(scm_data["eddy_diffusivity"][:, t0_scm : t1_scm], axis=1),\
                np.nanmean(scm_data["mixing_length"][:,t0_scm : t1_scm] / 1e3, axis=1),\
                np.nanmean(scm_data["nh_pressure"][:,  t0_scm : t1_scm] /\
                           scm_data["updraft_area"][:, t0_scm : t1_scm], axis=1\
                          ) / scm_data["rho_half"][:],\
                np.nanmean(scm_data["turbulent_entrainment"][:, t0_scm : t1_scm], axis=1),\
                np.nanmean(scm_data["updraft_buoyancy"][:, t0_scm : t1_scm], axis=1),\
                np.nanmean(scm_data["entrainment_sc"][:, t0_scm : t1_scm], axis=1)]

    x_lab = ["eddy_diffusivity", "mixing_length [km]", "non hydro pressure [Pa]",\
             "turbulent_entrainment", "buoyancy [m/s^2]", "entr and detr [1/m]"]

    for it in range(6):
        plt.subplot(2,3,it+1)
        if it < 4:
            plt.plot(scm_vars[it], scm_data["z_half"]/1e3, "-", c="royalblue", lw=3)

        if it == 2:
            plt.plot(np.nanmean(-les_data["updraft_ddz_p_alpha"][:, t0_les : t1_les], axis=1),\
                     les_data["z_half"], '-', color='gray', label='les', lw=3)
        if it == 4:
            plt.plot(scm_vars[it], scm_data["z_half"]/1e3, "-", c="royalblue", lw=3, label="b_upd")
            plt.plot(np.nanmean(scm_data["b_mix"][:, t0_scm : t1_scm],axis=1),\
                     scm_data["z_half"]/1e3, "-", color="darkorange", label="b_mix", lw=3)
            plt.legend()
        if it == 5:

            xmax = 0.015#np.min([np.max(scm_data["detrainment_sc"]), 0.05])
            if xmax == 0.0:
                xmax = np.max(scm_data["detrainment_sc"])

            plt.plot(scm_vars[it], scm_data["z_half"]/1e3, "-", c="royalblue", lw=3, label="entr")
            plt.plot(np.nanmean(scm_data["detrainment_sc"][:, t0_scm : t1_scm], axis=1),\
                     scm_data["z_half"]/1e3, "-", color="darkorange", label="detr", lw=3)
            plt.xlim([-0.0001,xmax])
            plt.legend()

        plt.xlabel(x_lab[it])
        plt.ylabel("z [km]")
        plt.grid(True)

    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()

def plot_tke_comp(scm_data, les_data, tmin, tmax, title, folder="plots/output/"):
    """
    Plots updraft and environment profiles from Scampy
    Input:
    scm_data - scm stats file
    les_data - les stats file
    tmin     - lower bound for time mean
    tmax     - upper bound for time mean
    title    - name for the created plot
    folder   - folder where to save the created plot
    """
    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.0)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0<= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax<= les_data["t"]))[0][0])

    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    x_lab =  ["tke_advection", "tke_buoy", "tke_dissipation", "tke_pressure",\
              "tke_transport","tke_shear"]

    plot_vars =  [scm_data["tke_advection"], scm_data["tke_buoy"],\
                  scm_data["tke_dissipation"], scm_data["tke_pressure"],\
                  scm_data["tke_transport"], scm_data["tke_shear"]]

    plot_x_les = [les_data["tke_prod_A"], les_data["tke_prod_B"],\
                  les_data["tke_prod_D"], les_data["tke_prod_P"],\
                  les_data["tke_prod_T"], les_data["tke_prod_S"]]

    xmax = 5*np.max(np.nanmean(scm_data["tke_entr_gain"][3:, t0_scm:t1_scm], axis=1))

    plots = []
    for plot_it in range(6):
        plots.append(plt.subplot(2,3,plot_it+1))
                               #(rows, columns, number)
        plots[plot_it].set_ylabel('z [km]')
        plots[plot_it].grid(True)
        if plot_it<6:
            # plots[plot_it].plot(np.nanmean(plot_x_les[plot_it][:, t0_les:t1_les],axis=1),\
            #                     les_data["z_half"], '-', color='gray', label='les', lw=3)
            plots[plot_it].plot(np.nanmean(plot_vars[plot_it][:, t0_scm:t1_scm],axis=1),\
                                scm_data["z_half"]/1e3, "-", color="royalblue", label='les', lw=3)
            plots[plot_it].set_xlabel(x_lab[plot_it])
            plots[plot_it].set_ylim([0, np.max(scm_data["z_half"]/1000.0)])
        else:
            plots[plot_it].plot(np.nanmean(scm_data["tke_entr_gain"][:, t0_scm:t1_scm],axis=1),\
                                scm_data["z_half"]/1e3, "-", color="royalblue", label="tke entr", lw=3)
            plots[plot_it].plot(np.nanmean(scm_data["tke_detr_loss"][:, t0_scm:t1_scm],axis=1),\
                                scm_data["z_half"]/1e3, "-", color="darkorange", label="tke detr", lw=3)
            plots[plot_it].set_xlabel('tke entr detr [1/m]')
            plots[plot_it].set_xlim([-1e-4, xmax])
            plots[plot_it].set_ylim([0, np.max(scm_data["z_half"]/1000.0)])
            plots[plot_it].legend()

    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()

def plot_spec_hum(scm_data, les_data, tmin, tmax, title, folder="plots/output/"):
    """
    Plots updraft and environment profiles from Scampy
    Input:
    scm_data - scm stats file
    les_data - les stats file
    tmin     - lower bound for time mean
    tmax     - upper bound for time mean
    title    - name for the created plot
    folder   - folder where to save the created plot
    """
    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.0)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0<= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax<= les_data["t"]))[0][0])

    scm_data["qv_mean"] = scm_data["qt_mean"]    - scm_data["ql_mean"]
    scm_data["upd_qv"]  = scm_data["updraft_qt"] - scm_data["updraft_ql"]
    scm_data["env_qv"]  = scm_data["env_qt"]     - scm_data["env_ql"]
    les_data["qv_mean"] = les_data["qt_mean"]    - les_data["ql_mean"]
    les_data["upd_qv"]  = les_data["updraft_qt"] - les_data["updraft_ql"]
    les_data["env_qv"]  = les_data["env_qt"]     - les_data["env_ql"]

    var = ["qv_mean", "upd_qv", "env_qv",\
           "ql_mean", "updraft_ql", "env_ql",\
           "qr_mean", "updraft_qr", "env_qr"]

    lab = ["mean qv [g/kg]", "updraft qv [g/kg]", "env qv [g/kg]",\
           "mean ql [g/kg]", "updraft ql [g/kg]", "env ql [g/kg]",\
           "mean qr [g/kg]", "updraft qr [g/kg]", "env qr [g/kg]"]

    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    for it in range(9):
        plt.subplot(3,3,it+1)
        plt.grid(True)
        plt.xlabel(lab[it])
        plt.plot(np.nanmean(les_data[var[it]][:, t0_les:t1_les],axis=1),\
                 les_data["z_half"], '-', color='gray', label='les', lw=3)
        plt.plot(np.nanmean(scm_data[var[it]][:, t0_scm:t1_scm],axis=1),\
                 scm_data["z_half"]/1e3, "-", color="royalblue", label='scm', lw=3)
        if it in [0,3,6]:
            plt.ylabel("z [km]")
        if it == 3:
           plt.plot(np.nanmean(les_data["qi_mean"][:, t0_les:t1_les],axis=1),\
                    les_data["z_half"], '-', color='silver', label='les ice', lw=3)
           plt.plot(np.nanmean(scm_data["qi_mean"][:, t0_scm:t1_scm],axis=1),\
                    scm_data["z_half"]/1e3, "-", color="cyan", label='scm ice', lw=3)
        if it == 4:
           plt.plot(np.nanmean(les_data["updraft_qi"][:, t0_les:t1_les],axis=1),\
                    les_data["z_half"], '-', color='silver', label='les ice', lw=3)
           plt.plot(np.nanmean(scm_data["updraft_qi"][:, t0_scm:t1_scm],axis=1),\
                    scm_data["z_half"]/1e3, "-", color="cyan", label='scm ice', lw=3)
        if it == 5:
           plt.plot(np.nanmean(les_data["env_qi"][:, t0_les:t1_les],axis=1),\
                    les_data["z_half"], '-', color='silver', label='les ice', lw=3)
           plt.plot(np.nanmean(scm_data["env_qi"][:, t0_scm:t1_scm],axis=1),\
                    scm_data["z_half"]/1e3, "-", color="cyan", label='scm ice', lw=3)
        if it == 6:
           plt.plot(np.nanmean(les_data["qs_mean"][:, t0_les:t1_les],axis=1),\
                    les_data["z_half"], '-', color='silver', label='les snow', lw=3)
           plt.plot(np.nanmean(scm_data["qs_mean"][:, t0_scm:t1_scm],axis=1),\
                    scm_data["z_half"]/1e3, "-", color="cyan", label='scm snow', lw=3)
        if it == 7:
           plt.plot(np.nanmean(les_data["updraft_qs"][:, t0_les:t1_les],axis=1),\
                    les_data["z_half"], '-', color='silver', label='les snow', lw=3)
           plt.plot(np.nanmean(scm_data["updraft_qs"][:, t0_scm:t1_scm],axis=1),\
                    scm_data["z_half"]/1e3, "-", color="cyan", label='scm snow', lw=3)
        if it == 8:
           plt.plot(np.nanmean(les_data["env_qs"][:, t0_les:t1_les],axis=1),\
                    les_data["z_half"], '-', color='silver', label='les snow', lw=3)
           plt.plot(np.nanmean(scm_data["env_qs"][:, t0_scm:t1_scm],axis=1),\
                    scm_data["z_half"]/1e3, "-", color="cyan", label='scm snow', lw=3)

    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()

def plot_upd_prop(scm_data, les_data, tmin, tmax, title, folder="plots/output/"):
    """
    Plots updraft and environment profiles from Scampy
    Input:
    scm_data - scm stats file
    les_data - les stats file
    tmin     - lower bound for time mean
    tmax     - upper bound for time mean
    title    - name for the created plot
    folder   - folder where to save the created plot
    """
    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.0)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0<= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax<= les_data["t"]))[0][0])

    les_data["massflux"]  = np.multiply(les_data["updraft_fraction"], les_data["updraft_w"])

    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    scm_var = ["updraft_area","updraft_w","massflux",\
               "u_mean", "thetal_mean","qt_mean"]

    les_var = ["updraft_fraction", "updraft_w", "massflux",\
               "u_translational_mean", "thetali_mean","qt_mean"]

    lab = ["updraft fraction", "updraft w [m/s]", "massflux [kg/m^2/s]",\
           "horizontal velocities [m/s]", "thetal mean [K]", "qt mean [g/kg]"]

    for it in range(6):
        plt.subplot(2,3,it+1)
        plt.grid(True)
        plt.plot(np.nanmean(les_data[les_var[it]][:, t0_les:t1_les], axis=1),\
                 les_data["z_half"], '-', color='gray', label='les', lw=3)
        plt.plot(np.nanmean(scm_data[scm_var[it]][:, t0_scm:t1_scm], axis=1),\
                 scm_data["z_half"]/1e3, "-", color="royalblue", label='scm', lw=3)
        plt.xlabel(lab[it])
        if it in [0,3]:
            plt.ylabel("z [km]")
        if it == 3:
            plt.plot(np.nanmean(les_data["v_translational_mean"][:,t0_les:t1_les],axis=1),
                     les_data["z_half"], '--', color='gray', label='v-les', lw=3)
            plt.plot(np.nanmean(scm_data["v_mean"][:, t0_scm:t1_scm], axis=1),\
                     scm_data["z_half"]/1e3, "-", color="darkorange", label='v-scm', lw=3)
            plt.legend()

    plt.savefig(folder + title)
    plt.clf()

def plot_tke_break(scm_data, les_data, tmin, tmax, title, folder="plots/output/"):
    """
    Plots updraft and environment profiles from Scampy
    Input:
    scm_data - scm stats file
    les_data - les stats file
    tmin     - lower bound for time mean
    tmax     - upper bound for time mean
    title    - name for the created plot
    folder   - folder where to save the created plot
    """
    # customize defaults
    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.0)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0<= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax<= les_data["t"]))[0][0])

    fig = plt.figure(1)
    fig.set_figheight(8)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    col = ["royalblue", "darkorange", "k", "darkgreen", "red", "purple"]

    scm_var = ["tke_advection","tke_buoy","tke_dissipation","tke_pressure",\
               "tke_transport","tke_shear"]

    les_var = ["tke_prod_A", "tke_prod_B", "tke_prod_D", "tke_prod_P",\
               "tke_prod_T", "tke_prod_S"]

    plt.subplot(121)
    for it in range(6):
        plt.plot(np.nanmean(scm_data[scm_var[it]][:, t0_scm:t1_scm], axis=1),\
                 scm_data["z_half"]/1e3, "-", color=col[it],  label=scm_var[it],\
                 lw=3)
    plt.ylim([0, np.max(scm_data["z_half"]/1e3)])
    plt.xlabel('tke componenets scm')
    plt.ylabel('height [km]')
    plt.legend()

    plt.subplot(122)
    for it in range(6):
        plt.plot(np.nanmean(les_data[les_var[it]][:, t0_les:t1_les], axis=1),\
                 les_data["z_half"], "-", color=col[it],  label=les_var[it],\
                 lw=3)
    plt.ylim([0, np.max(les_data["z_half"])])
    plt.xlabel('tke componenets les')
    plt.legend()

    plt.savefig(folder + title)
    plt.clf()

def plot_cvar_mean(scm_data, les_data, tmin, tmax, title, folder="plots/output/"):
    """
    Plots variance and covariance profiles from Scampy
    Input:
    scm_data - scm stats file
    les_data - les stats file
    tmin     - lower bound for time mean
    tmax     - upper bound for time mean
    title    - name for the created plot
    folder   - folder where to save the created plot
    """
    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.0)[0][0])
    t0_les = int(np.where(np.array(les_data["t"]) > tmin)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0<= scm_data["t"]))[0][0])
    t1_les = int(np.where(np.array(tmax<= les_data["t"]))[0][0])

    fig = plt.figure(1)
    fig.set_figheight(8)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    x_lab         = ["Hvar",       "QTvar",       "HQTcov"]
    plot_var_mean = ["Hvar_mean",  "QTvar_mean",  "HQTcov_mean"]
    plot_var_env  = ["env_Hvar",   "env_QTvar",   "env_HQTcov"]

    plots = []
    for plot_it in range(3):
        plots.append(plt.subplot(1,3,plot_it+1))
                               #(rows, columns, number)
        plots[plot_it].set_xlabel(x_lab[plot_it])
        plots[plot_it].set_ylabel('height [km]')
        plots[plot_it].set_ylim([0,\
                                 scm_data["z_half"][-1]/1000.0 +\
                                 (scm_data["z_half"][1]/1000.0 - scm_data["z_half"][0]/1000.0) * 0.5\
                                ])
        plots[plot_it].grid(True)
        plots[plot_it].xaxis.set_major_locator(ticker.MaxNLocator(2))

        plots[plot_it].plot(np.nanmean(les_data[plot_var_env[plot_it]][:, t0_les:t1_les], axis=1),\
                            les_data["z_half"], "-", label= 'les', c="gray", lw=4)
        plots[plot_it].plot(np.nanmean(scm_data[plot_var_mean[plot_it]][:,t0_scm:t1_scm], axis=1),\
                            scm_data["z_half"]/1e3, "-", label=plot_var_mean[plot_it], c="crimson", lw=3)
        plots[plot_it].plot(np.nanmean(scm_data[plot_var_env[plot_it]][:, t0_scm:t1_scm], axis=1),\
                            scm_data["z_half"]/1e3, "-", label=plot_var_env[plot_it],  c="forestgreen", lw=3)

    plots[0].legend()
    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()

def plot_cvar_comp(scm_data, tmin, tmax, title, folder="plots/output/"):
    """
    Plots variance and covariance components profiles from Scampy
    Input:
    scm_data   - scm stats file
    tmin   - lower bound for time mean
    tmax   - upper bound for time mean
    title  - name for the created plot
    folder - folder where to save the created plot
    """
    t0_scm = int(np.where(np.array(scm_data["t"]) > tmin*3600.0)[0][0])
    t1_scm = int(np.where(np.array(tmax*3600.0<= scm_data["t"]))[0][0])

    fig = plt.figure(1)
    fig.set_figheight(8)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', linewidth=4, markersize=10)

    plot_Hvar_c   = ["Hvar_dissipation",   "Hvar_entr_gain",   "Hvar_detr_loss",   "Hvar_shear",   "Hvar_rain"]
    plot_QTvar_c  = ["QTvar_dissipation",  "QTvar_entr_gain",  "QTvar_detr_loss",  "QTvar_shear",  "QTvar_rain"]
    plot_HQTcov_c = ["HQTcov_dissipation", "HQTcov_entr_gain", "HQTcov_detr_loss", "HQTcov_shear", "HQTcov_rain"]
    color_c       = ['darkgreen',          'purple',           'purple',           'darkorange',    'royalblue']

    x_lab         = ["Hvar",      "QTvar",      "HQTcov"]
    plot_var_data = [plot_Hvar_c, plot_QTvar_c, plot_HQTcov_c]

    plots = []
    for plot_it in range(3):
        plots.append(plt.subplot(1,3,plot_it+1))
                               #(rows, columns, number)
        plots[plot_it].set_xlabel(x_lab[plot_it])
        plots[plot_it].set_ylabel('height [km]')
        plots[plot_it].set_ylim([0,\
                                 scm_data["z_half"][-1]/1e3 +\
                                 (scm_data["z_half"][1]/1e3 - scm_data["z_half"][0]/1e3) * 0.5\
                                ])
        plots[plot_it].grid(True)
        plots[plot_it].xaxis.set_major_locator(ticker.MaxNLocator(2))

        for var in range(5):
            plots[plot_it].plot(np.nanmean(scm_data[plot_var_data[plot_it][var]][:, t0_scm:t1_scm], axis=1),\
                                scm_data["z_half"]/1e3, "-", label=plot_Hvar_c[var], c=color_c[var])

    plots[0].legend()
    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()

def plot_main(scm_srs, les_srs, scm_data, les_data, title,\
              cb_min, cb_max, folder="plots/output/"):
    """
    Plots the time series of Scampy simulations
    Input:
    scm_srs  - scm timeseries file
    les_srs  - les timeseries file
    scm_data - scm stats file
    les_data - les stats file
    title    - figure title
    cb_min   - lower colorbar limit for ql_mean and upd_w contour plots
    cb_min   - upper colorbar limit for ql_mean and upd_w contour plots
    folder   - folder where to save the created figure
    """

    scm_z_half = scm_data["z_half"]/1000.0
    scm_time   = scm_data["t"] /3600.0
    les_z_half = les_data["z_half"]
    les_time   = les_data["t"]

    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 16})
    mpl.rc('lines', linewidth=4, markersize=10)

    cmap = "RdBu_r"

    les_var = ["ql_mean", "updraft_w"]
    les_tit = ["LES ql mean [g/kg]", "LES upd w [m/s]"]
    for it in range(2):
        plt.subplot(3,2,it+1)
        levels = np.linspace(cb_min[it], cb_max[it], 11)
        cntrf = plt.contourf(les_time, les_z_half, les_data[les_var[it]],\
                             cmap=cmap, levels=levels, vmin=cb_min[it], vmax=cb_max[it])
        cbar = plt.colorbar(cntrf)
        plt.ylim([0, np.max(scm_z_half)])
        plt.ylabel('height [km]')
        plt.title(les_tit[it])

    scm_var = ["ql_mean", "updraft_w"]
    scm_tit = ["SCM ql mean [g/kg]", "SCM upd w [m/s]"]
    for it in range(2):
        plt.subplot(3,2,it+3)
        levels = np.linspace(cb_min[it], cb_max[it], 11)
        cntrf = plt.contourf(scm_time, scm_z_half, scm_data[scm_var[it]],\
                             cmap=cmap, levels=levels, vmin=cb_min[it], vmax=cb_max[it])
        cbar = plt.colorbar(cntrf)
        plt.ylim([0, np.max(scm_z_half)])
        plt.xlabel('time [h]')
        plt.ylabel('height [km]')
        plt.title(scm_tit[it])

    var = ["lwp_mean", "rwp_mean"]
    lab = ["lwp", "rwp"]
    for it in range(2):
        plt.subplot(3,2,it+5)
        plt.plot(les_srs["t"][1:],        les_srs[var[it]][1:], '-', c="gray", lw=3)
        plt.plot(scm_srs["t"][1:]/3600.0, scm_srs[var[it]][1:], '-', c="royalblue", lw=3)
        plt.xlim([0, scm_srs["t"][-1]/3600.0])
        plt.xlabel('time [h]')
        plt.ylabel(lab[it])
        plt.grid(True)

    plt.tight_layout()
    plt.savefig(folder + title)
    plt.clf()
    plt.close()

def plot_1D(scm_data, les_data, case, folder="plots/output/"):
    """
    Plots timeseries from Scampy
    Input:
    scm_data - scm stats file
    les_data - les stats file
    case     - case name
    folder   - folder where to save the created plot
    """
    fig = plt.figure(1)
    fig.set_figheight(12)
    fig.set_figwidth(14)
    mpl.rcParams.update({'font.size': 18})
    mpl.rc('lines', lw=3, markersize=10)

    # surface fluxes
    plot_scm_y = [scm_data["lhf"], scm_data["shf"]]
    plot_les_y = [les_data["lhf"], les_data["shf"]]
    y_lab = ["LHF", "SHF"]

    fig = plt.figure(1)
    for plot_it in range(2):
        plt.subplot(2,1,plot_it+1)
        plt.plot(les_data["t"][1:], plot_les_y[plot_it][1:], '-', color="gray", lw=3, label="LES")
        plt.plot(scm_data["t"][1:]/3600., plot_scm_y[plot_it][1:], '-', color="b", lw=3, label="SCM")
        plt.ylabel(y_lab[plot_it])
        plt.xlim([0, scm_data["t"][-1]/3600.])
        plt.grid(True)
    plt.xlabel('time [h]')
    plt.tight_layout()
    plt.savefig(folder + case + "surface_heat_fluxes.pdf")
    plt.clf()

    # cloud timeseries
    plot_scm_y = [scm_data["lwp_mean"],\
                  scm_data["cloud_cover_mean"],\
                  scm_data["rwp_mean"],\
                  scm_data["cloud_top_mean"]/1e3,
                  scm_data["cloud_base_mean"]/1e3]
    plot_les_y = [les_data["lwp_mean"],\
                  les_data["cloud_cover_mean"],\
                  les_data["rwp_mean"],\
                  les_data["cloud_top_mean"]/1e3,
                  les_data["cloud_base_mean"]/1e3]

    y_lab      = ['lwp, iwp', 'cloud_cover', 'rwp, swp', 'CB, CT [km]']

    fig = plt.figure(1)
    for plot_it in range(4):
        plt.subplot(2,2,plot_it+1)
        plt.plot(les_data["t"][1:], plot_les_y[plot_it][1:], '-', color="gray", label="LES", lw=3)
        plt.plot(scm_data["t"][1:]/3600., plot_scm_y[plot_it][1:], '-', color="b", label="SCM", lw=3)
        if plot_it == 0:
            plt.plot(les_data["t"][1:], les_data["iwp_mean"][1:], '-', color="silver", lw=3)
            plt.plot(scm_data["t"][1:]/3600., scm_data["iwp_mean"][1:], '-', color="cyan", lw=3)
        if plot_it == 2:
            plt.plot(les_data["t"][1:], les_data["swp_mean"][1:], '-', color="silver", lw=3)
            plt.plot(scm_data["t"][1:]/3600., scm_data["swp_mean"][1:], '-', color="cyan", lw=3)
        if plot_it == 3:
            plt.plot(les_data["t"][1:], plot_les_y[4][1:], '-', color="gray", lw=3)
            plt.plot(scm_data["t"][1:]/3600., plot_scm_y[4][1:], '-', color="b", lw=3)
        plt.legend()
        plt.grid(True)
        plt.xlim([0, scm_data["t"][-1]/3600.])
        plt.xlabel('time [h]')
        plt.ylabel(y_lab[plot_it])
    plt.tight_layout()
    plt.savefig(folder + case + "timeseries_cloud_properties.pdf")
    plt.clf()

    # separation radius
    fig = plt.figure(1)
    plt.plot(scm_data["t"][1:]/3600., scm_data["rd"][1:], '-', color="b", lw=3, label="SCM")
    plt.xlim([0, scm_data["t"][-1]/3600.])
    plt.xlabel('time [h]')
    plt.ylabel("plume separation radius [m]")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(folder + case + "plume_separation_radius.pdf")
    plt.clf()

def plot_contour_t(scm_data, les_data, fixed_cbar, cb_min_t, cb_max_t, folder="plots/output/"):
    """
    Plots the time series of Scampy simulations
    Input:
    scm_data - scm stats file
    les_data - les stats file
    fixed_cbar - bool flag, True if you want to plot with specified colorbar range
    cb_min_t - min values for colorbar
    cb_max_t - max_values for colorbar
    folder - folder where to save the created plot
    """

    scm_z_half = scm_data["z_half"]/1000.0
    scm_time   = scm_data["t"] /3600.0
    les_z_half = les_data["z_half"]
    les_time   = les_data["t"]

    les_vars  = ["thetali_mean", "env_thetali", "updraft_thetali",\
                 "ql_mean", "env_ql", "updraft_ql",\
                 "qi_mean", "env_qi", "updraft_qi",\
                 "qr_mean", "env_qr", "updraft_qr",\
                 "qt_mean", "env_qt", "updraft_qt",\
                 "qs_mean", "env_qs", "updraft_qs",\
                 "env_w", "updraft_w", "u_translational_mean", "v_translational_mean",\
                 "updraft_fraction", "updraft_buoyancy", "tke_mean",\
                 "massflux_h", "diffusive_flux_h", "total_flux_h",\
                 "massflux_qt", "diffusive_flux_qt", "total_flux_qt"]

    scm_vars  = ["thetal_mean", "env_thetal", "updraft_thetal",\
                 "ql_mean", "env_ql", "updraft_ql",\
                 "qi_mean", "env_qi", "updraft_qi",\
                 "qr_mean", "env_qr", "updraft_qr",\
                 "qt_mean", "env_qt", "updraft_qt",\
                 "qs_mean", "env_qs", "updraft_qs",\
                 "env_w", "updraft_w", "u_mean", "v_mean",\
                 "updraft_area", "updraft_buoyancy", "tke_mean",\
                 "massflux_h", "diffusive_flux_h", "total_flux_h",\
                 "massflux_qt", "diffusive_flux_qt", "total_flux_qt"]

    labels    = ["mean thl [K]", "env thl [K]", "updr thl [K]",\
                 "mean ql [g/kg]", "env ql [g/kg]", "updr ql [g/kg]",\
                 "mean qi [g/kg]", "env qi [g/kg]", "updr qi [g/kg]",\
                 "mean qr [g/kg]", "env qr [g/kg]", "updr qr [g/kg",\
                 "mean qt [g/kg]", "env qt [g/kg]", "updr qt [g/kg]",\
                 "mean qs [g/kg]", "env qs [g/kg]", "updr qs [g/kg]",\
                 "env w [m/s]", "updr w [m/s]", "u [m/s]", "v [m/s]",\
                 "updr area [%]", "updr buoyancy [m/s^2]", "mean TKE [m2/s2]",\
                 "massflux_h [kg*K/ms^2]", "diffusive_flux_h [kg*K/ms^2]","total_flux_h [kg*K/ms^2]",\
                 "massflux_qt [g*/ms^2]", "diffusive_flux_qt [g*/ms^2]", "total_flux_qt [g*/ms^2]"]

    fig_name =  ["mean_thl", "env_thl", "upd_thl",\
                 "mean_ql", "env_ql", "upd_ql",\
                 "mean_qi", "env_qi", "upd_qi",\
                 "mean_qr", "env_qr", "upd_qr",\
                 "mean_qt", "env_qt", "upd_qt",\
                 "mean_qs", "env_qs", "upd_qs",\
                 "env_w", "upd_w", "mean_u", "mean_v",\
                 "upd_area", "upd_buoyancy", "mean_TKE",\
                 "edmf_massflux_h", "edmf_diffusive_flux_h", "edmf_total_flux_h",\
                 "edmf_massflux_qt", "edmf_diffusive_flux_qt", "edmf_total_flux_qt"]

    for plot_it in range(len(labels)):
        fig = plt.figure(fig_name[plot_it])
        fig.set_figheight(12)
        fig.set_figwidth(14)
        mpl.rcParams.update({'font.size': 18})
        mpl.rc('lines', linewidth=4, markersize=10)
        cmap = "RdBu_r"

        # the initial condition for env thetal in scampy starts with zeros
        if scm_vars[plot_it]=="env_thetal":
            scm_data[scm_vars[plot_it]][:,0] = scm_data[scm_vars[plot_it]][:,1]

        scm_field = scm_data[scm_vars[plot_it]]
        les_field = les_data[les_vars[plot_it]]
        a_scm = scm_data['updraft_area']
        a_les = les_data['updraft_fraction']

        if ("updraft" in scm_vars[plot_it]):
            scm_field[np.where(a_scm==0.0)] = np.nan
            scm_field[np.where(np.isnan(a_scm))] = np.nan

        if ("updraft" in les_vars[plot_it]):
            les_field[np.where(a_les==0.0)] = np.nan
            les_field[np.where(np.isnan(a_les))] = np.nan

        plt.subplot(211)
        if fixed_cbar:
            levels = np.linspace(cb_min_t[plot_it], cb_max_t[plot_it], 11)
            cntrf = plt.contourf(les_time, les_z_half, les_field, cmap=cmap,\
                                 levels=levels, vmin=cb_min_t[plot_it], vmax=cb_max_t[plot_it])
            cbar = plt.colorbar(cntrf)
            cbar.set_label(labels[plot_it])
        else:
            plt.contourf(les_time, les_z_half, les_field, cmap=cmap)
            plt.colorbar()
        plt.ylim([0,np.max(scm_data["z_half"]/1000.0)])
        plt.ylabel('height [km]')
        plt.grid(True)

        plt.subplot(212)
        if fixed_cbar:
            levels = np.linspace(cb_min_t[plot_it], cb_max_t[plot_it], 11)
            cntrf = plt.contourf(scm_time, scm_z_half, scm_field, cmap=cmap,\
                                 levels=levels, vmin=cb_min_t[plot_it], vmax=cb_max_t[plot_it])
            cbar = plt.colorbar(cntrf)
            cbar.set_label(labels[plot_it])
        else:
            plt.contourf(scm_time, scm_z_half, scm_field, cmap=cmap)
            plt.colorbar()
        plt.ylim([0,np.max(scm_data["z_half"]/1000.0)])
        plt.xlabel('time [h]')
        plt.ylabel('height [km]')
        plt.grid(True)

        plt.tight_layout()
        plt.savefig(folder + "contour_" + fig_name[plot_it]+".pdf")
        plt.clf()
        plt.close()
