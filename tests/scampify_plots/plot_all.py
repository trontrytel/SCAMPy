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

localpath = os.getcwd()

cases = ["Dycoms", "Rico", "TRMM"]

nd = {}
nd["Dycoms"] = {}
nd["Rico"] = {}
nd["TRMM"] = {}

nd["Dycoms"]["les"] = localpath + "/data_les/DYCOMS_RF01_drizzle/CLIMA_1M_micro/Stats.DYCOMS_RF01.nc"
nd["Rico"]["les"]   = localpath + "/data_les/Rico/CLIMA_1M_micro/Stats.Rico.nc"
nd["TRMM"]["les"]   = localpath + "/data_les/TRMM_LBA/CLIMA_1M_micro/Stats.TRMM_LBA.nc"

nd["Dycoms"]["son"] = localpath + "/data_scm_online/Tests.Output.DYCOMS_RF01._RF01/stats/Stats.DYCOMS_RF01.nc"
nd["Rico"]["son"]   = localpath + "/data_scm_online/Tests.Output.Rico.o/stats/Stats.Rico.nc"
nd["TRMM"]["son"]   = localpath + "/data_scm_online/Tests.Output.TRMM_LBA.M_LBA/stats/Stats.TRMM_LBA.nc"

nd["Dycoms"]["sof"] = localpath + "/data_scm_offline/Tests.Output.DYCOMS_RF01._RF01/stats/Stats.DYCOMS_RF01.nc"
nd["Rico"]["sof"]   = localpath + "/data_scm_offline/Tests.Output.Rico.o/stats/Stats.Rico.nc"
nd["TRMM"]["sof"]   = localpath + "/data_scm_offline/Tests.Output.TRMM_LBA.M_LBA/stats/Stats.TRMM_LBA.nc"

nd["Dycoms"]["t0"] = 3
nd["Dycoms"]["t1"] = 4
nd["Rico"]["t0"] = 22
nd["Rico"]["t1"] = 24
nd["TRMM"]["t0"] = 5
nd["TRMM"]["t1"] = 6

nd["Dycoms"]["ql_max"] = 0.4
nd["Dycoms"]["qr_max"] = 1e-3
nd["Dycoms"]["ql_min"] = -1e-2
nd["Dycoms"]["qr_min"] = -1e-5

nd["Rico"]["ql_max"] = 0.015
nd["Rico"]["qr_max"] = 0.005
nd["Rico"]["ql_min"] = -1e-3
nd["Rico"]["qr_min"] = -1e-3

nd["TRMM"]["ql_max"] = 0.07
nd["TRMM"]["qr_max"] = 0.07
nd["TRMM"]["ql_min"] = -1e-3
nd["TRMM"]["qr_min"] = -1e-3

case = "Dycoms"

les_data = cmn.read_les_data_srs(Dataset(nd[case]["les"], 'r'))
son_data = cmn.read_scm_data_srs(Dataset(nd[case]["son"], 'r'))
sof_data = cmn.read_scm_data_srs(Dataset(nd[case]["sof"], 'r'))

t0_les = int(np.where(np.array(les_data["t"]) > nd[case]["t0"] * 3600.)[0][0])
t0_son = int(np.where(np.array(son_data["t"]) > nd[case]["t0"] * 3600.)[0][0])
t0_sof = int(np.where(np.array(sof_data["t"]) > nd[case]["t0"] * 3600.)[0][0])

t1_les = int(np.where(np.array(nd[case]["t1"] * 3600.0 <= les_data["t"]))[0][0])
t1_son = int(np.where(np.array(nd[case]["t1"] * 3600.0 <= son_data["t"]))[0][0])
t1_sof = int(np.where(np.array(nd[case]["t1"] * 3600.0 <= sof_data["t"]))[0][0])

les_upd_area = les_data["updraft_fraction"]
les_env_area = 1. - les_upd_area
son_upd_area = son_data["updraft_area"]
son_env_area = 1. - son_upd_area
sof_upd_area = sof_data["updraft_area"]
sof_env_area = 1. - sof_upd_area

son_upd_ql = np.multiply(son_upd_area, son_data["updraft_ql"])
son_upd_qr = son_data["updraft_qr"]
son_env_ql = np.multiply(son_env_area, son_data["env_ql"])
son_env_qr = son_data["env_qr"]

sof_upd_ql = np.multiply(sof_upd_area, sof_data["updraft_ql"])
sof_upd_qr = sof_data["updraft_qr"]
sof_env_ql = np.multiply(sof_env_area, sof_data["env_ql"])
sof_env_qr = sof_data["env_qr"]

les_upd_ql = np.multiply(les_upd_area, les_data["updraft_ql"])
les_upd_qr = np.multiply(les_upd_area, les_data["updraft_qr"])
les_env_ql = np.multiply(les_env_area, les_data["env_ql"])
les_env_qr = np.multiply(les_env_area, les_data["env_qr"])

les_mean_ql = les_data["ql_mean"]
les_mean_qr = les_data["qr_mean"]
son_mean_ql = son_data["ql_mean"]
son_mean_qr = son_data["qr_mean"]
sof_mean_ql = sof_data["ql_mean"]
sof_mean_qr = sof_data["qr_mean"]


fig = plt.figure(1)
fig.set_figheight(12)
fig.set_figwidth(14)
mpl.rcParams.update({'font.size': 18})
mpl.rc('lines', linewidth=4, markersize=10)

plt.subplot(2, 3, 1)
plt.plot(np.nanmean(les_mean_ql[:, t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='black', label='les_mean')
plt.plot(np.nanmean(les_env_ql[:,  t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='red',   label='les_env')
plt.plot(np.nanmean(les_upd_ql[:,  t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='blue',  label='les_upd')
plt.xlabel("LES ql [g/lg]")
plt.ylabel("z [km]")
plt.xlim([nd[case]["ql_min"], nd[case]["ql_max"]])
plt.grid(True)

plt.subplot(2, 3, 4)
plt.plot(np.nanmean(les_mean_qr[:, t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='black', label='les_mean')
plt.plot(np.nanmean(les_env_qr[:,  t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='red',   label='les_env')
plt.plot(np.nanmean(les_upd_qr[:,  t0_les : t1_les], axis=1), les_data["z_half"] / 1000., '-', c='blue',  label='les_upd')
plt.xlabel("LES qr [g/lg]")
plt.ylabel("z [km]")
plt.xlim([-1e-3, nd[case]["qr_max"]])
plt.xlim([nd[case]["qr_min"], nd[case]["qr_max"]])
plt.grid(True)

plt.subplot(2, 3, 2)
plt.plot(np.nanmean(son_mean_ql[:, t0_son : t1_son], axis=1), son_data["z_half"] / 1000., '-', c='black', label='scm_on_mean')
plt.plot(np.nanmean(son_env_ql[:,  t0_son : t1_son], axis=1), son_data["z_half"] / 1000., '-', c='red',   label='scm_on_env')
plt.plot(np.nanmean(son_upd_ql[:,  t0_son : t1_son], axis=1), son_data["z_half"] / 1000., '-', c='blue',  label='scm_on_upd')
plt.xlabel("SCM-on ql [g/lg]")
plt.ylabel("z [km]")
plt.xlim([nd[case]["ql_min"], nd[case]["ql_max"]])
plt.grid(True)

plt.subplot(2, 3, 5)
plt.plot(np.nanmean(son_mean_qr[:, t0_son : t1_son], axis=1), son_data["z_half"] / 1000., '-', c='black', label='scm_on_mean')
plt.plot(np.nanmean(son_env_qr[:,  t0_son : t1_son], axis=1), son_data["z_half"] / 1000., '-', c='red',   label='scm_on_env')
plt.plot(np.nanmean(son_upd_qr[:,  t0_son : t1_son], axis=1), son_data["z_half"] / 1000., '-', c='blue',  label='scm_on_upd')
plt.xlabel("SCM-on qr [g/lg]")
plt.ylabel("z [km]")
plt.xlim([nd[case]["qr_min"], nd[case]["qr_max"]])
plt.grid(True)

plt.subplot(2, 3, 3)
plt.plot(np.nanmean(sof_mean_ql[:, t0_sof : t1_sof], axis=1), sof_data["z_half"] / 1000., '-', c='black', label='scm_off_mean')
plt.plot(np.nanmean(sof_env_ql[:,  t0_sof : t1_sof], axis=1), sof_data["z_half"] / 1000., '-', c='red',   label='scm_off_env')
plt.plot(np.nanmean(sof_upd_ql[:,  t0_sof : t1_sof], axis=1), sof_data["z_half"] / 1000., '-', c='blue',  label='scm_off_upd')
plt.xlabel("SCM-off ql [g/lg]")
plt.ylabel("z [km]")
plt.xlim([nd[case]["ql_min"], nd[case]["ql_max"]])
plt.grid(True)

plt.subplot(2, 3, 6)
plt.plot(np.nanmean(sof_mean_qr[:, t0_sof : t1_sof], axis=1), sof_data["z_half"] / 1000., '-', c='black', label='scm_off_mean')
plt.plot(np.nanmean(sof_env_qr[:,  t0_sof : t1_sof], axis=1), sof_data["z_half"] / 1000., '-', c='red',   label='scm_off_env')
plt.plot(np.nanmean(sof_upd_qr[:,  t0_sof : t1_sof], axis=1), sof_data["z_half"] / 1000., '-', c='blue',  label='scm_off_upd')
plt.xlabel("SCM-off qr [g/lg]")
plt.xlim([nd[case]["qr_min"], nd[case]["qr_max"]])
plt.ylabel("z [km]")
plt.grid(True)

plt.tight_layout()
plt.savefig("Dycoms_les_on_off.pdf")
plt.clf()

