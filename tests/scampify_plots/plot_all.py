import sys
sys.path.insert(0, "./")
sys.path.insert(0, "../")

import os
import subprocess
import json
import warnings
from netCDF4 import Dataset
import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # To allow plotting when display is off
import matplotlib.pyplot as plt
from matplotlib import ticker

import common_scampify as cmn
import plot_scripts_scampify as pls

# nd - input dictionary
nd = {}
cases = ["Dycoms", "Rico", "TRMM"]
for case in cases:
  nd[case] = {}

# son - single column model online run
# sof - single column model offline run
# qon - with quadratures
# qof - without quadratures
localpath = os.getcwd()
nd["Dycoms"]["les"]     = localpath + "/data_les/DYCOMS_RF01_drizzle/CLIMA_1M_micro/Stats.DYCOMS_RF01.nc"
nd["Rico"]["les"]       = localpath + "/data_les/Rico/CLIMA_1M_micro/Stats.Rico.nc"
nd["TRMM"]["les"]       = localpath + "/data_les/TRMM_LBA/CLIMA_1M_micro/Stats.TRMM_LBA.nc"
nd["Dycoms"]["son_qon"] = localpath + "/data_scm_online/quad_on/Tests.Output.DYCOMS_RF01._RF01/stats/Stats.DYCOMS_RF01.nc"
nd["Rico"]["son_qon"]   = localpath + "/data_scm_online/quad_on/Tests.Output.Rico.o/stats/Stats.Rico.nc"
nd["TRMM"]["son_qon"]   = localpath + "/data_scm_online/quad_on/Tests.Output.TRMM_LBA.M_LBA/stats/Stats.TRMM_LBA.nc"
nd["Dycoms"]["sof_qon"] = localpath + "/data_scm_offline/quad_on/Tests.Output.DYCOMS_RF01._RF01/stats/Stats.DYCOMS_RF01.nc"
nd["Rico"]["sof_qon"]   = localpath + "/data_scm_offline/quad_on/Tests.Output.Rico.o/stats/Stats.Rico.nc"
nd["TRMM"]["sof_qon"]   = localpath + "/data_scm_offline/quad_on/Tests.Output.TRMM_LBA.M_LBA/stats/Stats.TRMM_LBA.nc"
nd["Dycoms"]["son_qof"] = localpath + "/data_scm_online/quad_of/Tests.Output.DYCOMS_RF01._RF01/stats/Stats.DYCOMS_RF01.nc"
nd["Rico"]["son_qof"]   = localpath + "/data_scm_online/quad_of/Tests.Output.Rico.o/stats/Stats.Rico.nc"
nd["TRMM"]["son_qof"]   = localpath + "/data_scm_online/quad_of/Tests.Output.TRMM_LBA.M_LBA/stats/Stats.TRMM_LBA.nc"
nd["Dycoms"]["sof_qof"] = localpath + "/data_scm_offline/quad_of/Tests.Output.DYCOMS_RF01._RF01/stats/Stats.DYCOMS_RF01.nc"
nd["Rico"]["sof_qof"]   = localpath + "/data_scm_offline/quad_of/Tests.Output.Rico.o/stats/Stats.Rico.nc"
nd["TRMM"]["sof_qof"]   = localpath + "/data_scm_offline/quad_of/Tests.Output.TRMM_LBA.M_LBA/stats/Stats.TRMM_LBA.nc"

# t0 - t1 - time range for averaging
nd["Dycoms"]["t0"] = 3
nd["Dycoms"]["t1"] = 4
nd["Rico"]["t0"] = 22
nd["Rico"]["t1"] = 24
nd["TRMM"]["t0"] = 5
nd["TRMM"]["t1"] = 6

# xlim
nd["Dycoms"]["ql_max"] = 0.4
nd["Dycoms"]["qr_max"] = 1.1e-4
nd["Dycoms"]["ql_min"] = -1e-2
nd["Dycoms"]["qr_min"] = -1e-5
nd["Rico"]["ql_max"] = 0.015
nd["Rico"]["qr_max"] = 3e-3
nd["Rico"]["ql_min"] = -1e-3
nd["Rico"]["qr_min"] = -1e-4
nd["TRMM"]["ql_max"] = 0.07
nd["TRMM"]["qr_max"] = 0.07
nd["TRMM"]["ql_min"] = -1e-3
nd["TRMM"]["qr_min"] = -1e-3
# ylim
nd["Dycoms"]["z_min"] = -0.02
nd["Dycoms"]["z_max"] = 1.2
nd["Rico"]["z_min"] = -0.1
nd["Rico"]["z_max"] = 3.1
nd["TRMM"]["z_min"] = -0.5
nd["TRMM"]["z_max"] = 15

# rd - data dictionary
rd = {}

for case in cases:

    rd[case] = {}

    # les specific data
    rd[case]["les"] = {}
    rd[case]["les"]["data"] = cmn.read_les_data_srs(Dataset(nd[case]["les"], 'r'))
    rd[case]["les"]["t0"] = int(np.where(rd[case]["les"]["data"]["t"] > nd[case]["t0"] * 3600.)[0][0])
    rd[case]["les"]["t1"] = int(np.where(nd[case]["t1"] * 3600.0 <= rd[case]["les"]["data"]["t"])[0][0])
    rd[case]["les"]["upd_area"] = rd[case]["les"]["data"]["updraft_fraction"]
    rd[case]["les"]["env_area"] = 1. - rd[case]["les"]["upd_area"]
    rd[case]["les"]["qr_upd"] = np.multiply(rd[case]["les"]["upd_area"], rd[case]["les"]["data"]["updraft_qr"])
    rd[case]["les"]["qr_env"] = np.multiply(rd[case]["les"]["env_area"], rd[case]["les"]["data"]["env_qr"])

    # scm specific data
    for model in ["son_qon", "son_qof", "sof_qon", "sof_qof"]:
        rd[case][model] = {}
        rd[case][model]["data"] = cmn.read_scm_data_srs(Dataset(nd[case][model], 'r'))
        rd[case][model]["t0"] = int(np.where(rd[case][model]["data"]["t"] > nd[case]["t0"] * 3600.)[0][0])
        rd[case][model]["t1"] = int(np.where(nd[case]["t1"] * 3600.0 <= rd[case][model]["data"]["t"])[0][0])
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

fig = plt.figure(1)
fig.set_figheight(12)
fig.set_figwidth(14)
mpl.rcParams.update({'font.size': 18})
mpl.rc('lines', linewidth=4, markersize=10)
cd1 = {"mean" : "black", "env" : "red", "upd" : "blue"}
cd2 = {"mean" : "grey",  "env" : "darksalmon", "upd" : "skyblue"}

for case in cases:

    plt.subplot(2, 3, 1)
    for sdm in ["mean", "env", "upd"]:
        plt.plot(np.nanmean(rd[case]["les"]["ql_"+sdm][:, rd[case]["les"]["t0"] : rd[case]["les"]["t1"]], axis=1), rd[case]["les"]["z_half"], '-', c=cd1[sdm], label='les_'+sdm)
    plt.xlabel("LES ql [g/kg]")
    plt.ylabel("z [km]")
    plt.xlim([nd[case]["ql_min"], nd[case]["ql_max"]])
    plt.ylim([nd[case]["z_min"], nd[case]["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 2)
    for sdm in ["mean", "env", "upd"]:
        plt.plot(np.nanmean(rd[case]["sof_qon"]["ql_"+sdm][:, rd[case]["sof_qon"]["t0"] : rd[case]["sof_qon"]["t1"]], axis=1), rd[case]["sof_qon"]["z_half"], '-', c=cd2[sdm], label='scm_sgs_'+sdm)
        plt.plot(np.nanmean(rd[case]["sof_qof"]["ql_"+sdm][:, rd[case]["sof_qof"]["t0"] : rd[case]["sof_qof"]["t1"]], axis=1), rd[case]["sof_qof"]["z_half"], '-', c=cd1[sdm], label='scm_mean_'+sdm)
    plt.xlabel("SCM-off ql [g/kg]")
    plt.xlim([nd[case]["ql_min"], nd[case]["ql_max"]])
    plt.ylim([nd[case]["z_min"], nd[case]["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 3)
    for sdm in ["mean", "env", "upd"]:
        plt.plot(np.nanmean(rd[case]["son_qon"]["ql_"+sdm][:, rd[case]["son_qon"]["t0"] : rd[case]["son_qon"]["t1"]], axis=1), rd[case]["son_qon"]["z_half"], '-', c=cd2[sdm], label='scm_sgs_'+sdm)
        plt.plot(np.nanmean(rd[case]["son_qof"]["ql_"+sdm][:, rd[case]["son_qof"]["t0"] : rd[case]["son_qof"]["t1"]], axis=1), rd[case]["son_qof"]["z_half"], '-', c=cd1[sdm], label='scm_mean_'+sdm)
    plt.xlabel("SCM-on ql [g/kg]")
    plt.xlim([nd[case]["ql_min"], nd[case]["ql_max"]])
    plt.ylim([nd[case]["z_min"], nd[case]["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 4)
    for sdm in ["mean"]:
        plt.plot(np.nanmean(rd[case]["les"]["qr_"+sdm][:, rd[case]["les"]["t0"] : rd[case]["les"]["t1"]], axis=1), rd[case]["les"]["z_half"], '-', c=cd1[sdm], label='les_'+sdm)
    plt.xlabel("LES qr [g/kg]")
    plt.ylabel("z [km]")
    plt.xlim([nd[case]["qr_min"], nd[case]["qr_max"]])
    plt.ylim([nd[case]["z_min"], nd[case]["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 5)
    for sdm in ["mean"]:
        plt.plot(np.nanmean(rd[case]["les"]["qr_"+sdm][:, rd[case]["les"]["t0"] : rd[case]["les"]["t1"]], axis=1), rd[case]["les"]["z_half"], '-', c="silver", label='les_'+sdm)
    for sdm in ["mean", "env", "upd"]:
        plt.plot(np.nanmean(rd[case]["sof_qon"]["qr_"+sdm][:, rd[case]["sof_qon"]["t0"] : rd[case]["sof_qon"]["t1"]], axis=1), rd[case]["sof_qon"]["z_half"], '-', c=cd2[sdm], label='scm_sgs_'+sdm)
        plt.plot(np.nanmean(rd[case]["sof_qof"]["qr_"+sdm][:, rd[case]["sof_qof"]["t0"] : rd[case]["sof_qof"]["t1"]], axis=1), rd[case]["sof_qof"]["z_half"], '-', c=cd1[sdm], label='scm_mean_'+sdm)
    plt.xlabel("SCM-off qr [g/kg]")
    plt.xlim([nd[case]["qr_min"], nd[case]["qr_max"]])
    plt.ylim([nd[case]["z_min"], nd[case]["z_max"]])
    plt.grid(True)

    plt.subplot(2, 3, 6)
    for sdm in ["mean"]:
        plt.plot(np.nanmean(rd[case]["les"]["qr_"+sdm][:, rd[case]["les"]["t0"] : rd[case]["les"]["t1"]], axis=1), rd[case]["les"]["z_half"], '-', c="silver", label='les_'+sdm)
    for sdm in ["mean", "env", "upd"]:
        plt.plot(np.nanmean(rd[case]["son_qon"]["qr_"+sdm][:, rd[case]["son_qon"]["t0"] : rd[case]["son_qon"]["t1"]], axis=1), rd[case]["son_qon"]["z_half"], '-', c=cd2[sdm], label='scm_sgs_'+sdm)
        plt.plot(np.nanmean(rd[case]["son_qof"]["qr_"+sdm][:, rd[case]["son_qof"]["t0"] : rd[case]["son_qof"]["t1"]], axis=1), rd[case]["son_qof"]["z_half"], '-', c=cd1[sdm], label='scm_mean_'+sdm)
    plt.xlabel("SCM-on qr [g/kg]")
    plt.xlim([nd[case]["qr_min"], nd[case]["qr_max"]])
    plt.ylim([nd[case]["z_min"], nd[case]["z_max"]])
    plt.grid(True)

    plt.tight_layout()
    plt.savefig("subdomains_comp_" + case + ".pdf")
    plt.clf()
