%% Plotting Anomaly Fields from Regression Analysis

RdBu_short = cbrewer('div', 'RdBu', 20, 'PCHIP');
RdYlBu = cbrewer('div', 'RdYlBu', 60, 'PCHIP');

% load('workspace_regression_patterns_PC1_equal_nino34_rev4.mat'); 
load('plotting_all_PC1_and_PC2_patterns_rev1.mat')


% [works] subplots of air temperature
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%'Specific humidity [kg kg$^{-1}$]',...
%    'Downward short-wave [W m$^{-2}$]','Sea level pressure [Pa]'];
figure('units', 'pixels', 'position', [0 0 1920 1080]);

subplot(2,2,1); % AIR TEMPERATURE
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colormap(flipud(RdBu_short));
m_proj('Equidistant Cylindrical','lat', [-90 90], 'lon',[23 383]);
h=m_pcolor(lon, lat, t2m_PC1); 
set(h,'linestyle','none'); hold on;
h=m_pcolor(lon+360, lat, t2m_PC1); 
set(h,'linestyle','none'); hold on;

m_coast('color',[.1 .1 .1]); % black coastline
m_grid('box', 'on', 'xtick', 7, 'ytick', 7, 'tickdir', 'in', ...
    'yaxislocation', 'left', 'Fontname', 'Times New Roman', ...
    'Fontsize', 18, 'color', RdYlBu(60,:), 'linewidth', .5);

col_limit = [-1.5 1.5];
 set(gca, 'clim', [col_limit(1) col_limit(2)]); hold on; % set colour limit

h3 = colorbar('color', [.1922 .2118 .5843], 'box', 'on', ...
    'tickdirection', 'in', 'Fontname', 'Times New Roman', ...
    'Fontsize', 18);
h = ylabel(h3, '[°C]', 'color', RdYlBu(60,:), ...
    'Fontname', 'Times New Roman', 'Fontsize', 18);
set(h3, 'YTick', linspace(col_limit(1), col_limit(2), 5)); % set limit of colourbar

text([-3.75 -3.75],[2.25 2.25],'a) Air temperature [°C]', 'fontname', 'Times New Roman', ...
       'fontsize', 20, 'color', 'k');

subplot(2,2,2); % SPECIFIC HUMIDITY
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colormap(flipud(RdBu_short));
m_proj('Equidistant Cylindrical','lat', [-90 90], 'lon',[23 383]);
h=m_pcolor(lon, lat, qv_PC1); 
set(h,'linestyle','none'); hold on;
h=m_pcolor(lon+360, lat, qv_PC1); 
set(h,'linestyle','none'); hold on;

m_coast('color',[.1 .1 .1]); % black coastline
m_grid('box', 'on', 'xtick', 7, 'ytick', 7, 'tickdir', 'in', ...
    'yaxislocation', 'left', 'Fontname', 'Times New Roman', ...
    'Fontsize', 18, 'color', RdYlBu(60,:), 'linewidth', .5);

col_limit = [-1e-3,1e-3];
 set(gca, 'clim', [col_limit(1) col_limit(2)]); hold on; % set colour limit

h3 = colorbar('color', [.1922 .2118 .5843], 'box', 'on', ...
    'tickdirection', 'in', 'Fontname', 'Times New Roman', ...
    'Fontsize', 18);
h = ylabel(h3, '[kg kg^{-1}]', 'color', RdYlBu(60,:), ...
    'Fontname', 'Times New Roman', 'Fontsize', 18);
set(h3, 'YTick', linspace(col_limit(1), col_limit(2), 5)); % set limit of colourbar

text([-3.75 -3.75],[2.25 2.25],'b) Specific humidity [kg kg^{-1}]', 'fontname', 'Times New Roman', ...
       'fontsize', 20, 'color', 'k');


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Subplots (c, d)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

subplot(2,2,3); % Downward short-wave
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colormap(flipud(RdBu_short));
m_proj('Equidistant Cylindrical','lat', [-90 90], 'lon',[23 383]);
h=m_pcolor(lon, lat, ssrd_PC1); 
set(h,'linestyle','none'); hold on;
h=m_pcolor(lon+360, lat, ssrd_PC1); 
set(h,'linestyle','none'); hold on;

m_coast('color',[.1 .1 .1]); % black coastline
m_grid('box', 'on', 'xtick', 7, 'ytick', 7, 'tickdir', 'in', ...
    'yaxislocation', 'left', 'Fontname', 'Times New Roman', ...
    'Fontsize', 18, 'color', RdYlBu(60,:), 'linewidth', .5);

col_limit = [-25 25];
 set(gca, 'clim', [col_limit(1) col_limit(2)]); hold on; % set colour limit

h3 = colorbar('color', [.1922 .2118 .5843], 'box', 'on', ...
    'tickdirection', 'in', 'Fontname', 'Times New Roman', ...
    'Fontsize', 18);
h = ylabel(h3, '[W m^{-2}]', 'color', RdYlBu(60,:), ...
    'Fontname', 'Times New Roman', 'Fontsize', 18);
set(h3, 'YTick', linspace(col_limit(1), col_limit(2), 5)); % set limit of colourbar

text([-3.75 -3.75],[2.25 2.25],'c) Downward short-wave [W m^{-2}]', 'fontname', 'Times New Roman', ...
       'fontsize', 20, 'color', 'k');

subplot(2,2,4); % Downward short-wave
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colormap(flipud(RdBu_short));
m_proj('Equidistant Cylindrical','lat', [-90 90], 'lon',[23 383]);
h=m_pcolor(lon, lat, msl_PC1/10); 
set(h,'linestyle','none'); hold on;
h=m_pcolor(lon+360, lat, msl_PC1/10); 
set(h,'linestyle','none'); hold on;

[cs, h] = m_contour(lon,lat,msl_PC1/10, [-40:10:40]);
set(h, 'color', [.3 .3 .3]);
clabel(cs, h);


m_coast('color',[.1 .1 .1]); % black coastline
m_grid('box', 'on', 'xtick', 7, 'ytick', 7, 'tickdir', 'in', ...
    'yaxislocation', 'left', 'Fontname', 'Times New Roman', ...
    'Fontsize', 18, 'color', RdYlBu(60,:), 'linewidth', .5);

col_limit = [-40 40];
 set(gca, 'clim', [col_limit(1) col_limit(2)]); hold on; % set colour limit

h3 = colorbar('color', [.1922 .2118 .5843], 'box', 'on', ...
    'tickdirection', 'in', 'Fontname', 'Times New Roman', ...
    'Fontsize', 18);
h = ylabel(h3, '[hPa]', 'color', RdYlBu(60,:), ...
    'Fontname', 'Times New Roman', 'Fontsize', 18);
set(h3, 'YTick', linspace(col_limit(1), col_limit(2), 5)); % set limit of colourbar

text([-3.75 -3.75],[2.25 2.25],'d) Sea level pressure [hPa]', 'fontname', 'Times New Roman', ...
       'fontsize', 20, 'color', 'k');


% finish fancy plot and print
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
p1 = 'C:\Users\Maurice Huguenin\Desktop\';
print('-dpng','-r300', [p1 'regression_maps']);



    
    











