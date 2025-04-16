import xarray as xr # version used here: 2022.11.0 (python 3.8.15)
import numpy as np  # version used here: 1.23.4 (python 3.8.15)
import pandas as pd # Version used here: 1.5.2 (python 3.8.15)
from scipy import interpolate # Version used here: 1.10.0 (python 3.8.15)

def crh(datain_path, dataout_path, test, version):
    in_file   = 'output_1Dset_' + test + '_1O_' + version + '.nc' # Just to retrieve number of columns and vertical levels.
    ecrad_out = xr.open_dataset(datain_path + in_file)
    param = ecrad_out.sizes['column']     # Test parameter
    profi = ecrad_out.sizes['half_level'] # Vertical Profile
    optsh = 4

    ## Heat capacity [J kg^-1 K^-1]:
    cp = 1.08*10**(3)
    ## Volumetric heat capacity as ICON evaluates on model levels, not pl:
    # cv = 0.718*10**3
    ## Gravity [m s^-2]:
    g = 9.8

    ## Cloud radiative heating rates matrices [K d^-1]:
    crh_sw  = np.zeros((optsh, profi, param))
    crh_lw  = np.zeros((optsh, profi, param))
    crh_net = np.zeros((optsh, profi, param))

    for opt in range(1, optsh + 1): # per optical scheme
        file      = 'output_1Dset_' + test + '_' + str(opt) + 'O_' + version + '.nc'
        ecrad_out = xr.open_dataset(datain_path + file)
        for i in range(param):
            pres = ecrad_out.pressure_hl.values[i] # Pressure [Pa]

            # Net Fluxes [W m-2]
            lw_net_flux_cloudy = ecrad_out.flux_dn_lw.values[i] - ecrad_out.flux_up_lw.values[i]
            lw_net_flux_clear  = ecrad_out.flux_dn_lw_clear.values[i] - ecrad_out.flux_up_lw_clear.values[i]
            sw_net_flux_cloudy = ecrad_out.flux_dn_sw.values[i] - ecrad_out.flux_up_sw.values[i]
            sw_net_flux_clear  = ecrad_out.flux_dn_sw_clear.values[i] - ecrad_out.flux_up_sw_clear.values[i]

            # Cloud Radiative Effect [W m-2]
            lw_cre = lw_net_flux_cloudy - lw_net_flux_clear
            sw_cre = sw_net_flux_cloudy - sw_net_flux_clear

            # By using the np.gradient function:
            crh_sw[opt-1, :, i] = -(g/cp)*(np.gradient(sw_cre, pres))*86400 # K day^(-1)
            crh_lw[opt-1, :, i] = -(g/cp)*(np.gradient(lw_cre, pres))*86400 # K day^(-1)

        crh_net[opt-1, :, :] = crh_sw[opt-1, :, :] + crh_lw[opt-1, :, :]
        
    # Creation of out .nc file. Include data coordinates in the future.
    ds = xr.Dataset(
#         coords = dict(
#             scheme   = (['scheme'], np.array('Fu','Yi13','Baran14','Baran16'),
#             level    = (['level'], altitude from climatology),
#             variable = (['variable'], test variable array),
#         ),
        data_vars = dict(
            crh_sw  = (['scheme', 'level', 'variable'], crh_sw),
            crh_lw  = (['scheme', 'level', 'variable'], crh_lw),
            crh_net = (['scheme', 'level', 'variable'], crh_net),
        ),
        attrs = dict(
            title = 'CRH 1D set from ' + test,
            description = 'SW, LW and Net CRH matrices for ' + test, 
        ),
    )

    out_file = 'CRH_1Dset_' + test + '_' + version + '.nc'
    ds.to_netcdf(dataout_path + out_file)
    
def crh3levels(datain_path, dataout_path, profile_file, test, version):
    file      = 'output_1Dset_' + test + '_236_1O_' + version + '.nc' # Just to retrieve number of columns and vertical levels.
    ecrad_out = xr.open_dataset(datain_path + file)
    param = ecrad_out.sizes['column']     # Test parameter
    profi = ecrad_out.sizes['half_level'] # Vertical Profile
    
    tropical_profile = pd.read_csv(profile_file, sep='\s+ ', engine = 'python')
    temp_int         = interpolate.interp1d(tropical_profile['pressure (hPa)'].iloc[::-1]*100,
                                            tropical_profile['temperature (K)'].iloc[::-1])
    temperature_hl   = temp_int(ecrad_out.pressure_hl.values)
    temp             = pd.Series(temperature_hl[0])
    tropopause_i     = temp.argmin()
    
    if test == 'test3': optsh = 4
    else: optsh = 2

    ## Heat capacity [J kg^-1 K^-1]:
    cp = 1.08*10**(3)
    ## Volumetric heat capacity as ICON evaluates on model levels, not pl:
    # cv = 0.718*10**3
    ## Gravity [m s^-2]:
    g = 9.8
    
    crh_sw  = np.zeros((optsh, profi, param))
    crh_lw  = np.zeros((optsh, profi, param))
    crh_net = np.zeros((optsh, profi, param))
    
    for pert in ['236','218','201']:

        ## Cloud radiative heating rates matrices [K d^-1]:
        globals()['crh_sw_' + pert]  = np.zeros((optsh, profi, param))
        globals()['crh_lw_' + pert]  = np.zeros((optsh, profi, param))
        globals()['crh_net_' + pert] = np.zeros((optsh, profi, param))

        for opt in range(1, optsh + 1): # per optical scheme
            file      = 'output_1Dset_' + test + '_' + pert + '_' + str(opt) + 'O_' + version + '.nc'
            ecrad_out = xr.open_dataset(datain_path + file)
            for i in range(param):
                pres = ecrad_out.pressure_hl.values[i] # Pressure [Pa]

                # Net Fluxes [W m-2]
                lw_net_flux_cloudy = ecrad_out.flux_dn_lw.values[i] - ecrad_out.flux_up_lw.values[i]
                lw_net_flux_clear  = ecrad_out.flux_dn_lw_clear.values[i] - ecrad_out.flux_up_lw_clear.values[i]
                sw_net_flux_cloudy = ecrad_out.flux_dn_sw.values[i] - ecrad_out.flux_up_sw.values[i]
                sw_net_flux_clear  = ecrad_out.flux_dn_sw_clear.values[i] - ecrad_out.flux_up_sw_clear.values[i]

                # Cloud Radiative Effect [W m-2]
                lw_cre = lw_net_flux_cloudy - lw_net_flux_clear
                sw_cre = sw_net_flux_cloudy - sw_net_flux_clear

                # By using the np.gradient function:
                eval('crh_sw_' + pert)[opt-1, :, i] = -(g/cp)*(np.gradient(sw_cre, pres))*86400 # K day^(-1)
                eval('crh_lw_' + pert)[opt-1, :, i] = -(g/cp)*(np.gradient(lw_cre, pres))*86400 # K day^(-1)

            eval('crh_net_' + pert)[opt-1, :, :] = eval('crh_sw_' + pert)[opt-1, :, :] + eval('crh_lw_' + pert)[opt-1, :, :]
            
        pert_i = (temp.iloc[tropopause_i:] - int(pert)).abs().argmin() # iloc position
        i      = tropopause_i + pert_i
#         if pert == 218: bottom_step = 5
#         else: bottom_step = 4
#         if pert == 236: top_step = 3
#         else: top_step = 4
        crh_sw[:, i-4:i+5, :]  = eval('crh_sw_'+str(pert))[:, i-4:i+5, :]
        crh_lw[:, i-4:i+5, :]  = eval('crh_lw_'+str(pert))[:, i-4:i+5, :]
        crh_net[:, i-4:i+5, :] = eval('crh_net_'+str(pert))[:, i-4:i+5, :]
   
    # Creation of out .nc file. Include data coordinates in the future.
    ds = xr.Dataset(
#         coords = dict(
#             scheme   = (['scheme'], np.array('Fu','Yi13','Baran14','Baran16'),
#             level    = (['level'], altitude from climatology),
#             variable = (['variable'], test variable array),
#         ),
        data_vars = dict(
            crh_sw  = (['scheme', 'level', 'variable'], crh_sw),
            crh_lw  = (['scheme', 'level', 'variable'], crh_lw),
            crh_net = (['scheme', 'level', 'variable'], crh_net),
        ),
        attrs = dict(
            title = 'CRH 1D set from ' + test,
            description = 'SW, LW and Net CRH matrices for ' + test + '. These profiles include three independent cloud layers at 201, 218 and 236 K. Radiative calculations where carried out independently.', 
        ),
    )

    out_file = 'CRH_1Dset_' + test + '_' + version + '.nc'
    ds.to_netcdf(dataout_path + out_file)
    
def crh_diff(in_data, test):
    for comp in ['sw','lw','net']:
        globals()['crhd_' + comp] = np.zeros((3, in_data.sizes['level'], in_data.sizes['variable']))
        if test == 'test4': scheme_range = [1]
        else: scheme_range = [1, 2, 3]
        for scheme in scheme_range: # Yi13 - Fu, Baran2016 - Fu, Baran2014 - Fu
            eval('crhd_' + comp)[scheme - 1] = in_data['crh_' + comp].values[scheme, :, :] - in_data['crh_' + comp].values[0, :, :]
    crh_diff = xr.Dataset(
        data_vars = dict(
            crhd_sw  = (['inter_scheme', 'level', 'variable'], crhd_sw), 
            crhd_lw  = (['inter_scheme', 'level', 'variable'], crhd_lw),
            crhd_net = (['inter_scheme', 'level', 'variable'], crhd_net),
        ),
        attrs = dict(
            title = 'Interscheme CRH difference - 1D set from ' + test,
            description = 'SW, LW and Net Interscheme CRH difference matrices for ' + test, 
        ),
    )
    return crh_diff

def crh_rdiff(crh_diff, crh_data, test):
    for comp in ['sw', 'lw', 'net']:
        globals()['crhrd_' + comp] = np.zeros((3, crh_diff.sizes['level'], crh_diff.sizes['variable']))
        ref = np.where(crh_data['crh_' + comp].values[0, :, :] == 0, 10**5,
                       crh_data['crh_' + comp].values[0, :, :])
        for inter_scheme in [1, 2, 3]: # Yi13 - Fu, Baran2016 - Fu, Baran2014 - Fu
            eval('crhrd_' + comp)[inter_scheme - 1] = crh_diff['crhd_' + comp].values[inter_scheme - 1, :, :]/ref
    crh_rdiff = xr.Dataset(
        data_vars = dict(
            crhrd_sw  = (['inter_scheme', 'level', 'variable'], crhrd_sw), 
            crhrd_lw  = (['inter_scheme', 'level', 'variable'], crhrd_lw),
            crhrd_net = (['inter_scheme', 'level', 'variable'], crhrd_net),
        ),
        attrs = dict(
            title = 'Interscheme CRH RELATIVE difference - 1D set from ' + test,
            description = 'SW, LW and Net Interscheme CRH RELATIVE difference matrices for ' + test, 
        ),
    )
    return crh_rdiff
    
def colorbar_range(in_data, test, diff):
    MaxMin_list = []
    for var in list(in_data.keys()):
        MaxMin_list.append(abs(in_data[var].values.max()))
        MaxMin_list.append(abs(in_data[var].values.min())) 
    MaxMax = max(MaxMin_list)
    pos_tick_list = []
    neg_tick_list = []
    if MaxMax > 0 and MaxMax < 1:
        MaxMax = MaxMax - (MaxMax % .05) + .05
        if MaxMax % 0.1 == 0: extra = 0
        else: extra = .05
        linthresh = .01
        linscale  = .1
        list01 = np.round(np.linspace(0.02, 0.09, 8), 2)
        list1  = np.round(np.linspace(0.2, MaxMax + extra, int((MaxMax + extra)*10 - 1)), 1)
        pos_tick_list = list(list01) + list(list1)
        neg_tick_list = list(-list1[::-1]) + list(-list01[::-1])
        maj_tickbar = [-0.1, 0, 0.1]
        min_tickbar = neg_tick_list + pos_tick_list
        vmin = -(MaxMax + extra)
        vmax = MaxMax + extra
    elif MaxMax < 10:
        MaxMax = MaxMax - (MaxMax % .5) + .5
        if MaxMax % 1 == 0: extra = 0
        else: extra = .5
        linthresh = .1
        linscale  = .1
        list01 = np.round(np.linspace(0.02, 0.09, 8), 2)
        list1  = np.round(np.linspace(0.2, 0.9, 9), 1)
        list2  = np.round(np.linspace(2, MaxMax + extra, int((MaxMax + extra) - 1)), 0)
        pos_tick_list = list(list1) + list(list2)
        neg_tick_list = list(-list2[::-1]) + list(-list1[::-1])
#         pos_tick_list = list(list01) + list(list1) + list(list2)
#         neg_tick_list = list(-list2[::-1]) + list(-list1[::-1]) + list(-list01[::-1])
        maj_tickbar = [-1, 0, 1]
        min_tickbar = neg_tick_list + pos_tick_list
        vmin = -(MaxMax + extra)
        vmax = MaxMax + extra
    elif MaxMax < 100:
        MaxMax = MaxMax - (MaxMax % 5) + 5
        #if MaxMax > 50: extra = 5
        if MaxMax % 10 == 0: extra = 0
        else: extra = 5
        linthresh = 1
        linscale  = .1
        list01 = np.round(np.linspace(0.02, 0.09, 8), 2)
        list1  = np.round(np.linspace(0.2, 0.9, 8), 2)
        list2  = np.round(np.linspace(2, 9, 8), 2)
        list3  = np.round(np.linspace(20, MaxMax + extra, int((MaxMax + extra)*0.1 - 1)), 2)
        if MaxMax > 20:
            if test == 'test3' and diff:
                pos_tick_list = list(list01) + list(list1) + list(list2) + list(list3)
                neg_tick_list = list(-list3[::-1]) + list(-list2[::-1]) + list(-list1[::-1]) + list(-list01[::-1])
                maj_tickbar = [-10, -1, -0.1, 0, 0.1, 1, 10]
            else:
                pos_tick_list = list(list2) + list(list3)
                neg_tick_list = list(-list3[::-1]) + list(-list2[::-1])
                maj_tickbar = [-10, -1, 0, 1, 10]
        else:
            pos_tick_list = list(list1) + list(list2) + list(list3)
            neg_tick_list = list(-list3[::-1]) + list(-list2[::-1]) + list(-list1[::-1])
            maj_tickbar = [-10, -1, 0, 1, 10]
        min_tickbar = neg_tick_list + pos_tick_list
        vmin = -(MaxMax + extra)
        vmax = MaxMax + extra
    else:
        MaxMax = MaxMax - (MaxMax % 10) + 10
        if MaxMax % 10 == 0: extra = 0
        else: extra = 10
        linthresh = .1
        linscale  = .1
        list1 = np.round(np.linspace(0.2, 0.9, 8), 1)
        list2 = np.round(np.linspace(2, 9, 8), 2)
        list3 = np.round(np.linspace(20, 90, 8), 2)
        list4 = np.round(np.linspace(200, MaxMax + extra, int((MaxMax + extra)*0.1 - 1)), 2)
        pos_tick_list = list(list1) + list(list2) + list(list3)# + list(list4)
        #neg_tick_list = list(-list4[::-1]) + list(-list3[::-1]) + list(-list2[::-1]) + list(-list1[::-1])
        neg_tick_list = list(-list3[::-1]) + list(-list2[::-1]) + list(-list1[::-1])
        maj_tickbar = [-100, -10, -1, 0, 1, 10, 100]
        min_tickbar = neg_tick_list + pos_tick_list
        vmin = -(MaxMax + extra)
        vmax = MaxMax + extra
    return vmin, vmax, linthresh, linscale, maj_tickbar, min_tickbar

def cloud_range(datain_path, test, version, in_data, diff):
    file     = 'input_1Dset_' + test + '_' + version + '.nc'
    input_qi = xr.open_dataset(datain_path + file).q_ice.values.T
    # Adding an extra fake number at cloud bottom to avoid half level pressure differences:
    for col in range(len(input_qi[0,:])):
        for row in range(len(input_qi[:,0])):
            if input_qi[::-1, col][row] != 0: 
                input_qi[::-1, col][row-1] = 1 
                break
                
    var_name = list(in_data.keys())[0][:-2]  
    if diff: opt_range = 3
    else: opt_range = 4
    for comp in ['sw', 'lw', 'net']:
        globals()[comp + '_ranges'] = np.zeros((opt_range, 2)) # min and max
    for opt in range(opt_range):
        for comp in ['sw','lw','net']:
            eval(comp + '_ranges')[opt, 0] = in_data[var_name + comp].values[opt, :-1, :].min(where = (input_qi > 0), initial = 100)
            eval(comp + '_ranges')[opt, 1] = in_data[var_name + comp].values[opt, :-1, :].max(where = (input_qi > 0), initial = -100)
    
    return sw_ranges, lw_ranges, net_ranges

def cloud_range_3levels(datain_path, test, version, in_data, diff):
    for pert in ['236','218','201']:
        file = 'input_1Dset_' + test + '_' + pert + '_' + version + '.nc'
        globals()['input_qi_' + pert] = xr.open_dataset(datain_path + file).q_ice.values.T
        # Adding an extra fake number at cloud bottom to avoid half level pressure differences:
        for col in range(len(eval('input_qi_' + pert)[0, :])):
            for row in range(len(eval('input_qi_' + pert)[:, 0])):
                if eval('input_qi_' + pert)[::-1, col][row] != 0: 
                    eval('input_qi_' + pert)[::-1, col][row - 2: row] = 1 
                    break
    input_qi = input_qi_236 + input_qi_218 + input_qi_201
    size = in_data.variable.size
    if test == 'test3' and size < 2000: input_qi = input_qi[:, :size]
    if test == 'test4':
        re = np.arange(0.1, 100.1, 0.1)
        lim5  = np.where(re == min(re, key = lambda x: abs(x - 10)))[0][0]
        lim60 = np.where(re == min(re, key = lambda x: abs(x - 60)))[0][0]
        input_qi = input_qi[:, lim5:lim60+1]

    var_name  = list(in_data.keys())[0][:-2]
    if diff: opt_range = in_data.inter_scheme.size
    else: opt_range = in_data.scheme.size
#     if diff: opt_range = 3
#     else: opt_range = 4
#     if test == 'test4': opt_range = 2
    for comp in ['sw', 'lw', 'net']:
        globals()[comp + '_ranges'] = np.zeros((opt_range, 2)) # min and max
    for opt in range(opt_range):
        for comp in ['sw','lw','net']:
            eval(comp + '_ranges')[opt, 0] = in_data[var_name + comp].values[opt, :-1, :].min(where = (input_qi > 0), initial = 100)
            eval(comp + '_ranges')[opt, 1] = in_data[var_name + comp].values[opt, :-1, :].max(where = (input_qi > 0), initial = -100)
    
    return sw_ranges, lw_ranges, net_ranges