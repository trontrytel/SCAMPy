import os
import subprocess
import json
import warnings
import pprint as pp
from netCDF4 import Dataset
import numpy as np

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
    os.system("python ../generate_namelist.py " + case)
    file_case = open(case + '.in').read()
    namelist  = json.loads(file_case)
    namelist['output']['output_root'] = "./Tests."
    namelist['meta']['uuid'] = case
    os.system("python ../generate_paramlist.py " +  case)
    file_params = open('paramlist_' + case + '.in').read()
    paramlist = json.loads(file_params)

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
    scm_outfile = outpath + "/stats/Stats." + case + ".nc"

    # reference LES simulations
    les_outpath = "scampify_plots/data_les/"
    if case == "TRMM_LBA":
        les_outfile = les_outpath + "TRMM_LBA/CLIMA_1M_micro/Stats.TRMM_LBA.nc"
    elif case == "Rico":
        les_outfile = les_outpath + "Rico/CLIMA_1M_micro/Stats.Rico.nc"
    elif case == "DYCOMS_RF01":
        les_outfile = les_outpath + "DYCOMS_RF01_drizzle/CLIMA_1M_micro/Stats.DYCOMS_RF01.nc"
    else:
       print("No LES reference data")

    res = {"namelist"      : namelist,
           "paramlist"     : paramlist,
           "scm_outfile"   : scm_outfile,
           "les_outfile"   : les_outfile}

    return res

def read_scm_data_srs(scm_data):
    """
    Read in the data from netcdf file into a dictionary that can be used for timeseries of profiles plots

    Input:
    sim_data  - netcdf Dataset with simulation results
    """
    variables = ["temperature_mean", "thetal_mean", "qt_mean", "ql_mean", "qr_mean",\
                 "updraft_area", "updraft_thetal", "updraft_w", "updraft_qt", "updraft_ql", "updraft_qr",\
                 "env_thetal", "env_w", "env_qt", "env_ql", "env_qr",\
                 #"Hvar_mean", "QTvar_mean", "HQTcov_mean", "env_Hvar", "env_QTvar", "env_HQTcov"\
                ]

    data = {"z_half" : np.array(scm_data["profiles/z_half"][:]),\
            "t"      : np.array(scm_data["profiles/t"][:])}

    for var in variables:
        data[var] = []
        if ("qt" in var or "ql" in var or "qr" in var):
            try:
                data[var] = np.transpose(np.array(scm_data["profiles/"  + var][:, :])) * 1000 #g/kg
            except:
                data[var] = np.transpose(np.array(scm_data["profiles/w_mean" ][:, :])) * 0  # have zeros when no ql or qr
        else:
            data[var] = np.transpose(np.array(scm_data["profiles/"  + var][:, :]))

    return data


def read_les_data_srs(les_data):
    """
    Read in the data from netcdf file into a dictionary that can be used for timeseries of profiles plots

    Input:
    les_data - netcdf Dataset with specific fileds taken from LES stats file
    """
    variables = ["temperature_mean", "thetali_mean", "qt_mean", "ql_mean", "qr_mean",\
                 "updraft_fraction", "updraft_thetali", "updraft_w", "updraft_qt", "updraft_ql", "updraft_qr",\
                 "env_thetali", "env_w", "env_qt", "env_ql", "env_qr",\
                 #"Hvar_mean", "QTvar_mean", "HQTcov_mean", "env_Hvar", "env_QTvar", "env_HQTcov"\
                ]

    data = {"z_half" : np.array(les_data["profiles/z_half"][:]),\
            "t"      : np.array(les_data["profiles/t"][:])}

    for var in variables:
        data[var] = []
        if ("qt" in var or "ql" in var or "qr" in var):
            try:
                data[var] = np.transpose(np.array(les_data["profiles/"  + var][:, :])) * 1000  #g/kg
            except:
                data[var] = np.transpose(np.array(les_data["profiles/w_mean" ][:, :])) * 0  # have zeros when no ql or qr
        else:
            data[var] = np.transpose(np.array(les_data["profiles/"  + var][:, :]))

    return data

def read_scm_data_timeseries(scm_data):
    """
    Read in the 1D data from netcdf file into a dictionary that can be used for timeseries plots

    Input:
    sim_data - netcdf Dataset with simulation results
    """
    variables = ["cloud_cover_mean", "cloud_base_mean", "cloud_top_mean",\
                 "lwp_mean", "rwp_mean"]

    data = {"z_half" : np.array(scm_data["profiles/z_half"][:]),\
            "t"      : np.array(scm_data["profiles/t"][:])}

    maxz = np.max(data['z_half'])
    for var in variables:
        data[var] = []
        data[var] = np.array(scm_data["timeseries/" + var][:])

    CT = np.array(scm_data["timeseries/cloud_top_mean"][:])
    CT[np.where(CT<=0.0)] = np.nan
    data["cloud_top_mean"] = CT

    CB = np.array(scm_data["timeseries/cloud_base_mean"][:])
    CB[np.where(CB>=maxz)] = np.nan
    data["cloud_base_mean"] = CB

    return data

def read_les_data_timeseries(les_data):
    """
    Read in the 1D data from netcdf file into a dictionary that can be used for timeseries plots

    Input:
    les_data - netcdf Dataset with specific fileds taken from LES stats file
    """
    variables = ["cloud_fraction", "cloud_base", "cloud_top",\
                 "lwp_mean", "rwp_mean"]

    data = {"z_half_les" : np.array(les_data["profiles/z_half"][:]),\
            "t"          : np.array(les_data["profiles/t"][:])}

    maxz = np.max(data['z_half_les'])

    data["shf"] = np.array(les_data["timeseries/shf_surface_mean"][:])
    data["lhf"] = np.array(les_data["timeseries/lhf_surface_mean"][:])
    data["lwp_mean"] = np.array(les_data["timeseries/lwp"][:])
    data["rwp_mean"] = np.zeros_like(data["lwp_mean"]) # TODO

    CF = np.array(les_data["timeseries/cloud_fraction"][:])
    CF[np.where(CF<=0.0)] = np.nan
    data["cloud_cover_mean"] = CF

    CT = np.array(les_data["timeseries/cloud_top"][:])
    CT[np.where(CT<=0.0)] = np.nan
    data["cloud_top_mean"] = CT

    CB = np.array(les_data["timeseries/cloud_base"][:])
    CB[np.where(CB>maxz)] = np.nan
    data["cloud_base_mean"] = CB

    return data