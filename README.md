# enso_diabatic_fluxes
Scripts used for my Msc thesis/PhD project: __Diabatic Contributions to Warm Water Volume Variability over ENSO Events__

Maurice F. Huguenin(1,2), Ryan M. Holmes (2,3) and Matthew H. England (2)

1 Institute for Atmospheric and Climate Science, ETH Zurich, 8092 Zurich, Switzerland
 
2 Climate Change Research Centre and ARC Centre of Excellence for Climate Extremes, University of New South Wales, New South Wales 2052, Australia 

3 School of Mathematics and Statistics, University of New South Wales, New South Wales 2052, Australia 

----

To import and visualize the model output in Matlab, I use the following packages available online:

- the [CSIRO netCDF/OPeNDAP](http://www.marine.csiro.au/sw/matlab-netcdf.html) interface to import .nc data with the getnc() function

- the [cbrewer](https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab) package for custom colours. Copyright (c) 2015, Charles Robert. All rights reserved.

- the [m_map](https://www.eoas.ubc.ca/~rich/map.html) package to visualize spatial maps. Pawlowicz, R., 2019. "M_Map: A mapping package for MATLAB", version 1.4k, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.




# List of Figures and their Analysis Scripts

__Fig. 1__: Fig1_wind_stress_regression_enso.m, import and visualize data from /data/workspace_regression_patterns_PC1_equals_N34_rev2.mat workspace

__Fig. 2__: Fig2_idealized_symmetric_timeseries.m

__Fig. 3__: visualize regression patterns from the plotting_all_PC1_and_PC2_patterns_rev1.mat workspace

__Fig. 4__: 

__Fig. 5__: 

----

__Fig. 6__: 

__Fig. 7__: 

__Fig. 8__: 

__Fig. 9__: Fig9_discharge_and_recharge_schematics.pptx shows the schematics with data extracted from ..... drawn in powerpoint

__Fig. 10__: 

----

__Fig. 11__: 

__Fig. 12__: 

__Fig. 13__: 

__Fig. 14__: Create lag regression with Fig14_lag_regression_of_wwv_terms_and_N34.m while loading in the saved time series from the data/WMT_time_series_1979-2016.mat workspace

__Fig. 15__: First I visualize as a bar plot the values saved in data/WMT_time_series_1979-2016.mat with python (Fig15_visualizing_script_barplot.py) and then add the individual bar plots together in Fig15_barplot_wwv_contributing_terms.pptx. Fig15_visualizing_script_barplot.py acts as an example, the scripts for the five other events are, except for the hard-coded data, identical. I choose this approach with a separate script for each event as it was the fastest that worked the first time.


# Data folder

- __plotting_all_PC1_and_PC2_patterns_rev1.mat__ contains data for all spatial regression patterns
- __WMT_time_series_1979-2016.mat__ contains the time series of Fig. 12 as well as the table for the contributions of all terms over the discharge and recharge periods
- __workspace_regression_patterns_PC1_equal_nino34_rev2.mat__ contains the spatial patterns and time series for Fig. 1
- __workspace_EXP1_and_EXP2_polynomial_PC_composites_symmetric.mat__ contains the saved variables for Fig. 3
