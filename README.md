Data and R code for analyzing the stream water isotopes of the Pituffik Peninsula, Greenland

This repository contains the data and code needed to perform analyses and modelling on the variability of stream water isotopic composition (d18O, d2H, dxs) across the Pituffik Peninsula, Greenland. Water samples were collected from June through August 2018 and in July 2019.
The file pituffik_stream_iso_analyses_script.r contains the R code to perform all analyses and plotting. It calls upon datafiles included here. Packages needed to run to code are all listed within the provided code.
The primary datafile is pituffik_stream_iso_2018_2019.csv. This data contains every water sample collected and isotopically analyzed. Data columns for this file are: 
sample_id: The unique sample identification for each water sample.
date: The date that the water was sampled on.
latitude: The latitude of the water sample site, in decimal degrees.
longitude: The longitude of the water sample site, in decimal degrees.
elevation: The elevation of the water sample site, in meters above sea level, extracted from ArcticDEM based on geographic coordinate.
d18O: The oxygen stable isotopic ratio of the water sample in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰).
d2H: The hydrogen stable isotopic ratio of the water sample in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰).
dxs: The deuterium excess of the water sample, calculated as dxs = 8*d2H - d18O. Multiply by 1000 to express in permil (‰).
site_name: The common name for the sampling site.
type: The type category of water source sampled.
type_number: The type category of water source sampled, expressed as numeric categories where 1 = lake, 2 = pool, 3 = stream, 4 = surface flow, 5 = snow or ice, 6 = precipitation rain event, and 7 = precipitation snow event. 
catchment_name: The name of the catchment basin that the sample is located within.
alt_catchment_name: The alternate or secondary name of the watershed basin that the sample is located within. 
stream_group: The name of the larger stream group that the sample belongs to.
minor_stream_group: The name of the stream subgroup that the sample is part of.
qc_flag: A marker indicating samples that were flagged (1) during quality checking for total loss, vial cracking, and/or evaporative water loss.

The datafile pituffik_stream_sources_iso_2018_2019.csv contains the potential source water samples that were collected and isotopically analyzed. Data columns for this file are: 
sample_id: The unique sample identification for each water sample.
date: The date that the water was sampled on.
latitude: The latitude of the water sample site, in decimal degrees.
longitude: The longitude of the water sample site, in decimal degrees.
elevation: The elevation of the water sample site, in meters above sea level, extracted from ArcticDEM based on geographic coordinate.
d18O: The oxygen stable isotopic ratio of the water sample in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰).
d2H: The hydrogen stable isotopic ratio of the water sample in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰).
dxs: The deuterium excess of the water sample, calculated as dxs = 8*d2H - d18O. Multiply by 1000 to express in permil (‰).
site_name: The common name for the sampling site.
catchment_name: The name of the catchment basin that the sample is located within.
alt_catchment_name: The alternate or secondary name of the watershed basin that the sample is located within. 
sourcetype: The major category of stream water source.
sourcetype_detailed: Further detail on the type of stream water source that was sampled.
qc_flag: A marker indicating samples that were flagged (1) during quality checking for total loss, vial cracking, and/or evaporative water loss.

The datafile pituffik_sample_sources_iterations.csv contains the output of the simmr mixing model. This output contains contributions of stream source waters to individual stream water samples. This file contains the output of 100 iterations of the model. Data columns for this file are: 
sample_id: The unique sample identification for each water sample.
stream_group: The name of the larger stream group that the sample belongs to.
minor_stream_group: The name of the stream subgroup that the sample is part of.
date: Date the sample was taken.
glacial_frac: The mean of the model output contribution distribution attributed to the glacial source.
glacial_frac_sd: The standard deviation of the model output contribution distribution attributed to the glacial source.
snowpack_frac: The mean of the model output contribution distribution attributed to the seasonal snowpack source.
snowpack_frac_sd: The standard deviation of the model output contribution distribution attributed to the seasonal snowpack source.
prcp_act_frac: The mean of the model output contribution distribution attributed to the rain+active layer source.
prcp_act_frac_sd: The standard deviation of the model output contribution distribution attributed to the rain+active layer source.
lentic_frac: The mean of the model output contribution distribution attributed to the glacial source.
lentic_frac_sd: The standard deviation of the model output contribution distribution attributed to the glacial source.

The datafile pituffik_catchment_sources_iterations.csv contains the output of the simmr mixing model. This output contains contributions of stream source waters to stream water samples aggregated at the catchment level. This file contains the output of 100 iterations of the model. Data columns for this file are: 
stream_group: The name of the catchment that the samples were aggregated into.
glacial_frac: The mean of the model output contribution distribution attributed to the glacial source.
glacial_frac_sd: The standard deviation of the model output contribution distribution attributed to the glacial source.
snowpack_frac: The mean of the model output contribution distribution attributed to the seasonal snowpack source.
snowpack_frac_sd: The standard deviation of the model output contribution distribution attributed to the seasonal snowpack source.
prcp_act_frac: The mean of the model output contribution distribution attributed to the rain+active layer source.
prcp_act_frac_sd: The standard deviation of the model output contribution distribution attributed to the rain+active layer source.
lentic_frac: The mean of the model output contribution distribution attributed to the glacial source.
lentic_frac_sd: The standard deviation of the model output contribution distribution attributed to the glacial source.

The file pfk_iso_wx_day.csv contains daily weather and water vapor isotope data taken at Pituffik Space Base from August 2017 through May 2020. Data previously published in Akers et al., 2020 (https://doi.org/10.5194/acp-20-13929-2020), and more information on data sourcing and collection can be found there.
Data columns for this file are: 
daybreak: The day that the rowdata covers. 
wind_az: Mean azimuth of wind, in degrees. 
wind_sp: Mean wind speed, in m/s. 
H2O: Specific humidity of ambient air, measured on a Picarro L2130-i, in ppmv. d18O: Mean oxygen stable isotopic ratio of ambient air in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰). 
d2H: Mean hydrogen stable isotopic ratio of ambient air in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰). dxs: Mean deuterium excess of ambient air, calculated as dxs = 8*d2H - d18O. Multiply by 1000 to express in permil (‰). 
tavg: Mean of 10 min ambient air temperature readings at SMT site, in °C. 
tdew: Mean of 10 min ambient air dew point temperature readings at SMT site, in °C. 
pres: Mean of 10 min station air pressure readings at SMT site, in mb.
rh_ice: Mean of 10 min relative humidity with respect to ice readings at SMT site, in percent. 
nao: Daily NAO index. 
ao: Daily AO index. 
aao: Daily AAO index. 
pna: Daily PNA index. 
seaice: Baffin Bay sea ice extent, in km². Eness: Katabatic deviation of wind, defined as |100° - wind azimuth|, in degrees.
mo: The month of the rowdata.
sn: The season of the rowdata. 
tmax_af: Maximum air temperature measured by USAF at THU airport, in °C. 
tmin_af: Minimum air temperature measured by USAF at THU airport, in °C. 
tavg_af: Mean air temperature measured by USAF at THU airport, in °C. 
tdew_af: Mean dew point temperature measured by USAF at THU airport, in °C. 
wind_pk_az_af: Azimuth of maximum wind speed measured by USAF at THU airport, in degrees. 
wind_pk_sp_af: Maximum wind speed measured by USAF at THU airport, in m/s. 
prcp_snow: Precipitation amount of snow measured by USAF at THU airport, in mm. 
prcp_H2O: Liquid precipitation amount measured by USAF at THU airport, in mm. 
snow_depth: Depth of snow on ground measured by USAF at THU airport, in cm. 
prcp_occur: Precipitation recorded during day by USAF at THU airport, including trace amounts.

The file pituffik_pet_1981_2023.csv contains daily potential evapotranspiration (PET) rates extracted for the Pituffik region from 1981 to 2023 from the global database provided by Singer et al., 2021 (https://doi.org/10.1038/s41597-021-01003-9). 
Data columns for this file are: 
pet: Daily potential evapotranspiration in mm / day. 
yr: Year that PET data falls within. 
doy: Day of year that PET falls on. 
date: Date of daily PET value. 
mo: Month that PET data falls within. 
sn: Season that PET data falls within.

The file THU_L_melt_data.csv contains melt proxy data from the THU_L PROMICE weather station. Data from the Programme for Monitoring of the Greenland Ice Sheet (PROMICE) are provided and supported by the Geological Survey of Denmark and Greenland (GEUS) at http://www.promice.dk. and provided by How et al., 2022 (https://doi.org/10.22008/FK2/IW73UU).
Data columns for this file are: 
date: Date of daily recorded value. 
snow_height: Height of the snowpack above the ice surface in m.
z_ice_surf: Surface elevation of the ice sheet at the station, given as m where 0 was the elevation at original installation.
Z_ice_surf_yr_anom: Surface elevation of the ice sheet at the station, given as m, but the surface elevation is set at 0 on 01 Jan to show the annual surface loss due to ablation.
yr: Year that PET data falls within. 
doy: Day of year that PET falls on.




