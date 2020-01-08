# enso_diabatic_fluxes
Scripts used for my Msc thesis/PhD project: __Diabatic Contributions to Warm Water Volume Variability over ENSO Events__

Maurice F. Huguenin(1,2), Ryan M. Holmes (2,3) and Matthew H. England (2)

1 Institute for Atmospheric and Climate Science, ETH Zurich, 8092 Zurich, Switzerland
 
2 Climate Change Research Centre and ARC Centre of Excellence for Climate Extremes, University of New South Wales, New South Wales 2052, Australia 

3 School of Mathematics and Statistics, University of New South Wales, New South Wales 2052, Australia 

----

# Packages and Functions

I use the following packages which are publicly available online:

- the [CSIRO netCDF/OPeNDAP](http://www.marine.csiro.au/sw/matlab-netcdf.html) interface to import .nc data with the getnc() function

- the [cbrewer](https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) package for custom colours. Copyright (c) 2015, Charles Robert. All rights reserved.

- the [m_map](https://www.eoas.ubc.ca/~rich/map.html) package to visualize spatial maps. Pawlowicz, R., 2019. "M_Map: A mapping package for MATLAB", version 1.4k, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.

- the [TEOS-10](http://www.teos-10.org/publications.htm) package which contains oceanographic functions for Matlab, e.g. the conversion from potential to conservative temperature. Calculations using conservative temperature are virtually identical to those using potential temperature. IOC,  SCOR  and  IAPSO,  2010: The  international  thermodynamic  equation  of  seawater –  2010:  Calculation  and use of thermodynamic properties.  Intergovernmental Oceanographic Commission, Manuals and Guides No. 56, UNESCO (English), 196 pp.

# List of Figures and their Analysis Scripts

__EOF_based_regression_analysis_of_wind_stress_anomalies.m__: describes the regression analysis of wind stress anomalies in the equatorial Pacific as in [McGregor et al., 2014](https://doi.org/10.1002/2014JC010286).

__nanmonmean.m__: function to create the monthly mean of equatorial transect data

__testmap.m__: function to quickly create spatial maps of data with testmap(lon, lat, data) where all input has dimensions [1440 1080]. I use this function often in the analysis scripts for quick checks

__boom.m__: function that closes all active figures

----

__Fig. 1__: Visualize data from the regression analysis (EOF_based_regression_analysis_of_wind_stress_anomalies.m) stored in the data/workspace_regression_patterns_PC1_equals_N34_rev2.mat workspace with Fig1_wind_stress_regression_enso.m

__Fig. 2__: Import data from workspace_regression_patterns_PC1_equal_nino34_rev2.mat and workspace_EXP1_and_EXP2_polynomial_PC_composites_symmetric.mat for visualization in Fig2_idealized_symmetric_timeseries.m

__Fig. 3__: Visualize regression patterns from the plotting_all_PC1_and_PC2_patterns_rev1.mat workspace as in Fig. 1. __CANNOT FIND THE SCRIPT ANYMORE.__

__Fig. 4__: Using Fig4_timeseries_nino34_and_wwv_and_taux_and_nsfc.m to plot climate indices

__Fig. 5__: Fig5_timeseries_wwv_transport_clim.m used for plotting the climatological warm water volume (WWV) budget

----

__Fig. 6__: Creating equatorial transects of vertical mixing/surface forcing with Fig6_and_Fig8_equatorial_transects_of_water_mass_transformation_velocities.m. The first part extracts the data first and saves it in a .mat file while the second part describes the plotting routine. The .mat workspaces are located in data/equatorial_transect_workspace/.. The file nanmonmean.m is used to calculate monthly mean values for the transects.

__Fig. 7__: Calculating and plotting the WWV budget for the idealized El Niño and La Niña events with Fig7_time_series_wwv_transport_short.m

__Fig. 8__: As for Fig. 6 but using the model output of the idealized El Niño and La Niña simulations

__Fig. 9__: Fig9_discharge_and_recharge_schematics.pptx shows the schematics with data extracted from ..... drawn in powerpoint

__Fig. 10__: Plotting equatorial isotherm distribution with Fig10_20_degrees_isotherm_depth.m

----

__Fig. 11__: Using Fig11_time_series_nino34_and_wwv_and_sst_and_ohc.m to plot time series of observed and simulated N34 values and WWV anomalies.

__Fig. 12__: Calculation of the WWV budget over the 1979-2016 period with Fig12_time_series_wwv_transport.m

__Fig. 13__: Creating equatorial transects as in Fig. 6 and Fig. 8 with Fig13_equatorial_transects_during_LN1988.m The first part of the script saves the equatorial data in data/ACCESS-OM2_output*varsat_Eq.mat which are then loaded in and plotted

__Fig. 14__: Create lag regression with Fig14_lag_regression_of_wwv_terms_and_N34.m while loading in the saved time series from the data/WMT_time_series_1979-2016.mat workspace

__Fig. 15__: First I visualize as a bar plot the values saved in data/WMT_time_series_1979-2016.mat with python (Fig15_visualizing_script_barplot.py) and then add the individual bar plots together in Fig15_barplot_wwv_contributing_terms.pptx. Fig15_visualizing_script_barplot.py acts as an example, the scripts for the five other events are, except for the hard-coded data, identical. I choose this approach with a separate script for each event as it was the fastest that worked the first time.


# Data folder
```
data
?   ??? equatorial_transect_workspace          # data for 
?   ?   ??? ACCESS-OM2_output*_varsat_Eq.mat
?   ?   ??? equatorial_slices_pnEXP1 # data for the ideal. El Nino
?   ?   ??? equatorial_slices_pnEXP2 # data for the ideal. La Nina
Figs. 6, 8 and 13
?   ??? WMT_time_series_1979-2016.mat # data for time series of Fig. 12 and tabular data for Fig. 15
?   ??? plotting_all_PC1_and_PC2_patterns_rev1.mat # data for all spatial regression patterns (Figs. 1 and 3)
?   ??? workspace_EXP1_and_EXP2_polynomial_PC_composites_symmetric.mat # variables for Fig. 3
?   ??? workspace_regression_patterns_PC1_equal_nino34_rev2.mat # data for spatial patterns and time series for Fig. 1
```

# Model output data

- __idealized El Niño__: G:/Maurice_ENSO_Data/pn_EXP1_composite_nino_windstress/ on cube
- __idealized La Niña__: G:/Maurice_ENSO_Data/pn_EXP2_composite_nina_windstress/ on cube
- __1979-2016 simulation__: /g/data/e14/mv7494/access-om2/archive/025deg_jra55_iaf/ on gadi

