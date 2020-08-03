#!/bin/env python

"""
Run script for the KNMI radar rainfall nowcasts.

More info: https://pysteps.github.io/
"""
import os
import datetime
import netCDF4
import numpy as np
import pprint
import sys
import time

start = time.time()

import pysteps as stp

os.chdir('c:\\Users\\imhof002\\pysteps')

data_source = "knmi"

# Verification settings - Some of these settings are not necessary for the run, they are used for verification steps.
verification = {
    "experiment_name"   : "pysteps_radar_RAC_1",
    "overwrite"         : True,            # to recompute nowcasts
    "v_thresholds"      : [0.1, 1.0, 5.0],       # [mm/h]                 
    "v_leadtimes"       : [10, 20, 30, 45, 60],     # [min]
    "v_accu"            : None,             # [min]
    "seed"              : 42,               # for reproducibility
    "doplot"            : True,            # save figures
    "dosaveresults"     : True              # save verification scores to csv
}

# Forecast settings
forecast = {
    "n_lead_times"      : 36,       # timesteps per nowcast
    "r_threshold"       : 0.1,      # rain/no rain threshold [mm/h]
    "unit"              : "mm/h",   # mm/h or dBZ
    "transformation"    : "dB",     # None or dB 
    "adjust_domain"     : None      # None or square
}

# The experiment set-up, this includes tuneable parameters
experiment = {
    ## the events           event start     event end       update cycle  data source
    "data"              : [("201106100815", "201106111100", 15,           "knmi"), # 2 - 20110611
                           ("201108060815", "201108071100", 15,           "knmi"), # 3 - 20110807
                           ("201108070815", "201108081100", 15,           "knmi")], # 4 - 20110808
#                           ("201108080815", "201108091100", 15,           "knmi"), # 5 - 20110809
#                           ("201108130815", "201108141100", 15,           "knmi"), # 6 - 20110814
#                           ("201108180815", "201108181100", 15,           "knmi"), # 7 - 20110819
#                           ("201108220815", "201108231100", 15,           "knmi"), # 8 - 20110823
#                           ("201108250815", "201108261100", 15,           "knmi"), # 9 - 20110826
#                           ("201108270815", "201108281100", 15,           "knmi"), # 10 - 20110828
#                           ("201109060815", "201109071100", 15,           "knmi"), # 11 - 20110907
#                           ("201109070815", "201109081100", 15,           "knmi")], # 12 - 20110908
                           # ("201109101830", "201109110245", 15,           "knmi"), # 1 - 20110911
    
    ## the methods
    "oflow_method"      : ["lucaskanade"],      # lucaskanade, darts
    "adv_method"        : ["semilagrangian"],   # semilagrangian, eulerian
    "nwc_method"        : ["steps"],
    "noise_method"      : ["nonparametric"],    # parametric, nonparametric, ssft
    "decomp_method"     : ["fft"],
    
    ## the parameters
    "n_ens_members"     : [20],
    "ar_order"          : [2],
    "n_cascade_levels"  : [8],
    "noise_adjustment"  : ['auto'],
    "conditional"       : [False],
    "precip_mask"       : [True],
    "mask_method"       : ["incremental"],      # obs, incremental, sprog
    "prob_matching"     : ["cdf"],
}

###########################################################
# Get the correct paths
###########################################################

root_path = stp.rcparams.data_sources[data_source]["root_path"]
path_fmt = stp.rcparams.data_sources[data_source]["path_fmt"]
fn_pattern = stp.rcparams.data_sources[data_source]["fn_pattern"]
fn_ext = stp.rcparams.data_sources[data_source]["fn_ext"]
importer_name = stp.rcparams.data_sources[data_source]["importer"]
importer_kwargs = stp.rcparams.data_sources[data_source]["importer_kwargs"]
timestep = stp.rcparams.data_sources[data_source]["timestep"]

# Conditional parameters
## parameters that can be directly related to other parameters
def cond_pars(pars):
    for key in list(pars):
        if key == "oflow_method":
            if pars[key].lower() == "darts":  pars["n_prvs_times"] = 9
            else:                             pars["n_prvs_times"] = 3
        elif key.lower() == "n_cascade_levels":
            if pars[key] == 1 : pars["bandpass_filter"] = "uniform"
            else:               pars["bandpass_filter"] = "gaussian"
        elif key.lower() == "nwc_method":
            if pars[key] == "extrapolation" : pars["n_ens_members"] = 1
    return pars
    
# Prepare the list of all parameter sets of the verification
parsets = [[]]
for _, items in experiment.items():
    parsets = [parset+[item] for parset in parsets for item in items]

# Now loop all parameter sets
for n, parset in enumerate(parsets):
    
    # Build parameter set
    
    p = {}
    for m, key in enumerate(experiment.keys()):
        p[key] = parset[m]
    ## apply conditional parameters
    p = cond_pars(p)
    ## include all remaining parameters
    p.update(verification)
    p.update(forecast)
    
    print("************************")
    print("* Parameter set %02d/%02d: *" % (n+1, len(parsets)))
    print("************************")
    
    pprint.pprint(p)
    
    # If necessary, build path to results
    path_to_experiment = os.path.join(stp.rcparams.outputs["path_outputs"], p["experiment_name"])
    # subdir with event date
    path_to_nwc = os.path.join(path_to_experiment, '-'.join([p["data"][0], p["data"][3]]))
    for key, item in p.items():
		# include only variables that change
        if len(experiment.get(key,[None])) > 1 and key.lower() is not "data":
            path_to_nwc = os.path.join(path_to_nwc, '-'.join([key, str(item)]))
    try:
        os.makedirs(path_to_nwc)
    except FileExistsError:
        pass
        
    # **************************************************************************
    # NOWCASTING
    # ************************************************************************** 
    
    # Loop forecasts within given event using the prescribed update cycle interval
  
    if p["v_accu"] is None:
        p["v_accu"] = timestep
    
    # Loop forecasts for given event
    startdate   = datetime.datetime.strptime(p["data"][0], "%Y%m%d%H%M")
    enddate     = datetime.datetime.strptime(p["data"][1], "%Y%m%d%H%M")
    countnwc = 0
    while startdate + datetime.timedelta(minutes = p["n_lead_times"]*timestep) <= enddate:
        try:
            
            # filename of the nowcast netcdf
            outfn = os.path.join(path_to_nwc, "%s_nowcast.netcdf" % startdate.strftime("%Y%m%d%H%M"))
        
            ## check if results already exists
            run_exist = False
            if os.path.isfile(outfn):
                fid = netCDF4.Dataset(outfn, 'r')
                if fid.dimensions["time"].size == p["n_lead_times"]:
                    run_exist = True
                    if p["overwrite"]:
                        os.remove(outfn)
                        run_exist = False    
                else:
                    os.remove(outfn)
                    
            if run_exist:
                print("Nowcast %s_nowcast already exists in %s" % (startdate.strftime("%Y%m%d%H%M"),path_to_nwc))
    
            else:
                countnwc += 1
                print("Computing the nowcast (%02d) ..." % countnwc)
                
                print("Starttime: %s" % startdate.strftime("%Y%m%d%H%M"))
                
                ## redirect stdout to log file
                logfn =  os.path.join(path_to_nwc, "%s_log.txt" % startdate.strftime("%Y%m%d%H%M")) 
                print("Log: %s" % logfn)
                orig_stdout = sys.stdout
                f = open(logfn, 'w')
                sys.stdout = f
                
                print("*******************")
                print("* %s *****" % startdate.strftime("%Y%m%d%H%M"))
                print("* Parameter set : *")
                pprint.pprint(p)
                print("*******************")
                
                print("--- Start of the run : %s ---" % (datetime.datetime.now()))
                
                ## time
                t0 = time.time()
            
                # Read inputs
                print("Read the data...")
                
                ## find radar field filenames
                input_files = stp.io.find_by_date(startdate, root_path, path_fmt, fn_pattern,
                                                  fn_ext, timestep, p["n_prvs_times"])
                
        
                ## read radar field files
                importer    = stp.io.get_method(importer_name, "importer")
                R, _, metadata = stp.io.read_timeseries(input_files, importer, **importer_kwargs)
                metadata0 = metadata.copy()
                metadata0["shape"] = R.shape[1:]
                
                # Prepare input files
                print("Prepare the data...")
                
                ## if requested, make sure we work with a square domain
                reshaper = stp.utils.get_method(p["adjust_domain"])
                R, metadata = reshaper(R, metadata)
        
                ## if necessary, convert to rain rates [mm/h]    
                converter = stp.utils.get_method("mm/h")
                R, metadata = converter(R, metadata)
                
                ## threshold the data
                R[R < p["r_threshold"]] = 0.0
                metadata["threshold"] = p["r_threshold"]
                
                ## convert the data
                converter = stp.utils.get_method(p["unit"])
                R, metadata = converter(R, metadata)
                    
                ## transform the data
                transformer = stp.utils.get_method(p["transformation"])
                R, metadata = transformer(R, metadata)
                
                ## set NaN equal to zero
                R[~np.isfinite(R)] = metadata["zerovalue"]
                
                # Compute motion field
                oflow_method = stp.motion.get_method(p["oflow_method"])
                UV = oflow_method(R)
                
                ####
                # Perform the nowcast       
                ####
        
                ## define the callback function to export the nowcast to netcdf
                converter   = stp.utils.get_method("mm/h")
                def export(X):
                    ## convert to mm/h
                    X,_ = converter(X, metadata)
                    # readjust to initial domain shape
                    X,_ = reshaper(X, metadata, inverse=True)
                    # export to netcdf
                    stp.io.export_forecast_dataset(X, exporter)
                
                ## initialize netcdf file
                incremental = "timestep" if p["nwc_method"].lower() == "steps" else None
                exporter = stp.io.initialize_forecast_exporter_netcdf(outpath = path_to_nwc, outfnprefix = startdate.strftime("%Y%m%d%H%M"), startdate = startdate,
                                  timestep = timestep, n_timesteps = p["n_lead_times"], shape = metadata0["shape"], 
                                  n_ens_members = p["n_ens_members"], metadata = metadata0, incremental=incremental)
                
                ## start the nowcast
                nwc_method = stp.nowcasts.get_method(p["nwc_method"])
                R_fct = nwc_method(R, UV, p["n_lead_times"], p["n_ens_members"],
                                p["n_cascade_levels"], kmperpixel=metadata["xpixelsize"]/1000, 
                                timestep=timestep, R_thr=metadata["threshold"], 
                                extrap_method=p["adv_method"], 
                                decomp_method=p["decomp_method"], 
                                bandpass_filter_method=p["bandpass_filter"], 
                                noise_method=p["noise_method"], 
                                noise_stddev_adj=p["noise_adjustment"],
                                ar_order=p["ar_order"],conditional=p["conditional"], 
                                probmatching_method=p["prob_matching"], 
                                mask_method=p["mask_method"], 
                                callback=export, 
                                return_output=False, seed=p["seed"],
                                measure_time = True)
                
                ## save results
                stp.io.close_forecast_files(exporter)
                R_fct = None
                
                # save log
                print("--- End of the run : %s ---" % (datetime.datetime.now()))
                print("--- Total time : %s seconds ---" % (time.time() - t0))
                sys.stdout = orig_stdout
                f.close()
                
            # next forecast
            startdate += datetime.timedelta(minutes = p["data"][2])

        except ValueError:
            print('No nowcast for this date')
            
            # next forecast
            startdate += datetime.timedelta(minutes = p["data"][2])
    
        
end = time.time()
print('Total process took', (end - start)/3600.0, 'hours')
