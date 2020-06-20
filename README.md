# enso_diabatic_fluxes
Analysis Scripts and Data used for the publication:

Huguenin, M. F., Holmes, R. M., & England, M. H. (2020). Diabatic Contributions to Warm Water Volume Variability over ENSO Events

# Packages and Functions

I use the following packages which are publicly available online:

- the [CSIRO netCDF/OPeNDAP](http://www.marine.csiro.au/sw/matlab-netcdf.html) interface to import .nc data with the getnc() function

- the [cbrewer](https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) package for custom colours. Copyright (c) 2015, Charles Robert. All rights reserved.

- the [m_map](https://www.eoas.ubc.ca/~rich/map.html) package to visualize spatial maps. Pawlowicz, R., 2019. "M_Map: A mapping package for MATLAB", version 1.4k, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.

- the [TEOS-10](http://www.teos-10.org/publications.htm) package which contains oceanographic functions for Matlab, e.g. the conversion from potential to conservative temperature. Calculations using conservative temperature are virtually identical to those using potential temperature. IOC,  SCOR  and  IAPSO,  2010: The  international  thermodynamic  equation  of  seawater –  2010:  Calculation  and use of thermodynamic properties.  Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

# List of Figures and Analysis Scripts

[EOF_based_regression_analysis_of_wind_stress_anomalies.m](EOF_based_regression_analysis_of_wind_stress_anomalies.m): describes the regression analysis of wind stress anomalies in the equatorial Pacific as in [McGregor et al., 2014](https://doi.org/10.1002/2014JC010286).

[nanmonmean.m](nanmonmean.m): function to create the monthly mean of equatorial transect data

[testmap.m](testmap.m): function to quickly create spatial maps of data with testmap(lon, lat, data) where all input has dimensions [1440 1080]. I use this function often in the analysis scripts for quick checks

[boom.m](boom.m): function that closes all active figures

----

__Fig. 1__: Visualize data from the regression analysis ([EOF_based_regression_analysis_of_wind_stress_anomalies.m](EOF_based_regression_analysis_of_wind_stress_anomalies.m)) stored in the data/workspace_regression_patterns_PC1_equals_N34_rev2.mat workspace with Fig1_wind_stress_regression_enso.m

__Fig. 2__: Import data from data/workspace_regression_patterns_PC1_equal_nino34_rev2.mat and workspace_EXP1_and_EXP2_polynomial_PC_composites_symmetric.mat for visualization in [Fig2_idealized_symmetric_timeseries.m](Fig2_idealized_symmetric_timeseries.m)

__Fig. 3__: Using [Fig3_timeseries_nino34_and_wwv_and_taux_and_nsfc.m](Fig3_timeseries_nino34_and_wwv_and_taux_and_nsfc.m) to plot climate indices, new calculations for the upper ocean 0-2000 m heat calculations are imported from 

__Fig. 4 & 7__: Creating equatorial transects of vertical mixing/surface forcing with [Fig4_and_Fig7_equatorial_transects_climatology](Fig4_and_Fig7_equatorial_transects_climatology). The first part extracts the data first and saves it in a .mat file while the second part describes the plotting routine. The .mat workspaces are located in data/equatorial_transect_workspace/.. The file nanmonmean.m is used to calculate monthly mean values for the transects.

__Fig. 5__: [Fig5_timeseries_wwv_transport_clim.m](Fig5_timeseries_wwv_transport_clim.m) used for plotting the climatological warm water volume (WWV) budget

----

__Fig. 6__: Calculating and plotting the WWV budget for the idealized El Niño and La Niña events with [Fig6_and_Fig8_time_series_wwv_transport_short.m](Fig6_and_Fig8_time_series_wwv_transport_short.m)

__Fig. 8__: [Fig9_discharge_and_recharge_schematics.pptx](Fig8_discharge_and_recharge_schematics.pptx) shows the schematics with data extracted from ..... drawn in powerpoint

__Fig. 9__: Using [Fig9_time_series_nino34_and_wwv.m](Fig9_time_series_nino34_and_wwv.m) to plot time series of observed and simulated N34 values and WWV anomalies.

----

__Fig. 10__: Plotting equatorial isotherm distribution with [Fig10_20_degrees_isotherm_depth.m](Fig10_20_degrees_isotherm_depth.m)

__Fig. 11__: Calculation of the WWV budget over the 1979-2016 period with [Fig11_and_Fig13_time_series_wwv_transport_and_composites.m](Fig11_and_Fig13_time_series_wwv_transport_and_composites.m)

__Fig. 12__: Creating equatorial transects as in Fig. 4 and Fig. 7 with [Fig12_equatorial_transects_during_LN1988.m](Fig12_equatorial_transects_during_LN1988.m) The first part of the script saves the equatorial data in data/ACCESS-OM2_output*varsat_Eq.mat which are then loaded in and plotted

__Fig. 13__: Create composite time series of the three strong El Niño and La Niña events using the WWV budget time series from Fig. 12 and the script [Fig11_and_Fig13_time_series_wwv_transport.m](Fig11_and_Fig13_time_series_wwv_transport.m)

__Fig. 14__: First I visualize as a bar plot the values saved in data/WMT_time_series_1979-2016.mat with python ([Fig14_visualizing_script_barplot.py](Fig14_visualizing_script_barplot.py) and then add the individual bar plots together in [Fig14_barplot_wwv_contributing_terms.pptx](Fig15_barplot_wwv_contributing_terms.pptx). Fig14_visualizing_script_barplot.py acts as an example, the scripts for the five other events are, except for the hard-coded data, identical. I choose this approach with a separate script for each event as it was the fastest that worked the first time.

__Fig. A1__: Visualize regression patterns from the plotting_all_PC1_and_PC2_patterns_rev1.mat workspace with [FigA1_anomaly_fields_ENSO_regression](FigA1_anomaly_fields_ENSO_regression)


[OLD_Fig14_lag_regression_of_wwv_terms_and_N34.m](OLD_Fig14_lag_regression_of_wwv_terms_and_N34.m) is the script creating a lag regression of the WWV budget terms. This figure is not in the manuscript anymore.

# Data folder
- __equatorial_transect_workspace__:
	+ data from the 1979-2016 ACCESS-OM2 simulation
	+ equatorial_slices_pnEXP1 # data for the idealized El Niño
	+ equatorial_slices_pnEXP2 # data for the idealized La Niña

- __WMT_time_series_1979-2016.mat__: the WWV budget terms as time series for the 1979-2016 period
- __plotting_all_PC1_and_PC2_patterns_rev1.mat__: data for all spatial regression patterns
- __workspace_EXP1_and_EXP2_polynomial_PC_composites_symmetric.mat__: saved variables for Fig. 3
- __workspace_regression_patterns_PC1_equal_nino34_rev2.mat__: data for the spatial patterns and time series for Fig. 1
- __OHC_clim_2000.mat__ global uper 0-2000 m climatological ocean heat content time series used for Fig. 3.
- __OHC_pnEXP1_and_pnEXP2_2000.mat__ global uper 0-2000 m ocean heat content time series from the idealised simulations used for Fig. 3.


# Model output data

The model output is quite big(114 GB for the idealised simulations and 1.2 TB for the 1979-2016 hindcast simulation). To get access to the model output, please send me an email to m.huguenin-virchaux@unsw.edu.au.

- __idealized El Niño__: H:/Maurice_ENSO_Data/pn_EXP1_composite_nino_windstress/ on my personal hard drive
- __idealized La Niña__: H:/Maurice_ENSO_Data/pn_EXP2_composite_nina_windstress/ on my personal hard drive
- __1979-2016 simulation__: mv7494@gadi.nci.org.au:/g/data/e14/mv7494/access-om2/archive/025deg_jra55_iaf/ 

