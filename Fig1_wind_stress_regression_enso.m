%% Era-Interim plotting first two EOF of wind stress patterns

% Dimensions:
%            LAT  = 480             % 0.75x0.75 degree grid
%            LON  = 241
%            TIME = 13880  (UNLIMITED)
%            MSL  = mean sea level 456 [Pa]
%            bnds = 2

% empty world map
empty = nan(480,241);
years_2 = linspace(1979,2016,456);
RdBu_short = cbrewer('div', 'RdBu', 21, 'PCHIP');
RdBu_short = cbrewer('div', 'RdBu', 21, 'PCHIP');
RdYlBu = cbrewer('div', 'RdYlBu', 60, 'PCHIP');
antarctica = nature;

load('workspace_regression_patterns_PC1_equal_nino34_rev4.mat'); 


%% [works, 13.11.17, 19:42 AEDT] subplots of first two EOFs
  
figure(3);
colormap(flipud(RdBu_short));
m_proj('miller','lat', [-30 30], 'lon',[100 300]);
h=m_pcolor(lon, lat, EOFs(:,:,1)*1.0e2); 
set(h,'linestyle','none'); hold on;
h=m_pcolor(lon+360, lat, EOFs(:,:,1)*1.0e2); 
set(h,'linestyle','none'); hold on;


m_coast('color',[.1 .1 .1]);
m_grid('box', 'on', 'xtick', 5, 'ytick', 5, 'tickdir', 'in', ...
    'yaxislocation', 'left', 'Fontname', 'Times New Roman', ...
    'Fontsize', 10, 'color', RdYlBu(60,:), 'linewidth', .5);

% % drawing quiver arrows (on top of Nino3.4 area)
h5 = m_vec(.15, lon, lat, arrows(:,:,1), arrows(:,:,2), [0 0 0], ...
    'headlength',3.5, 'shaftwidth', .9, 'Edgecolor', [1 1 1], ...
    'linewidth', .5); % plot white area around arrows
h6 = m_vec(.15, lon, lat, arrows(:,:,1), arrows(:,:,2), [0 0 0], ...
    'headlength',3.5, 'shaftwidth', .9, 'Edgecolor', [1 1 1]);

set(gca, 'clim', [-3 3]); hold on; % set colour limit

h3 = colorbar('color', [.1922 .2118 .5843], 'box', 'on', ...
    'tickdirection', 'in', 'Fontname', 'Times New Roman', ...
    'Fontsize', 10);
h = ylabel(h3, '[10^{-2} N m^{-2}]', 'color', RdYlBu(60,:), ...
    'Fontname', 'Times New Roman', 'Fontsize', 10);
set(h3, 'YTick', linspace(-3, 3, 5)); % set limit of colourbar

% draw rectangular area in Niño3.4 region: [-5 5], [170W, 120W]
m_line([190, 190],[-5, 5],'color', [0 0 0], 'linewidth', .5); % m_line([lat],[lon])
m_line([240, 240],[-5, 5], 'color', [0 0 0], 'linewidth', .5); % m_line([lat],[lon])
m_line([190, 240],[-5, -5], 'color', [0 0 0], 'linewidth', .5); % m_line([lat],[lon])
m_line([190, 240],[5, 5], 'color', [0 0 0], 'linewidth', .5); % m_line([lat],[lon])
m_text(225.5, -2.3, 'Niño3.4','color', [0 0 0], 'fontsize', 10, ...
    'Fontname', 'Times New Roman', 'Fontsize',6); hold on;


% finish fancy plot and print
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = 'E:\2017 Herbstsemester\Master Seminar 2\';
print('-dpng','-r1000', [directory 'EOF1_wind_stress_anomalies_v11']);


figure(4);
colormap(flipud(RdBu_short));
m_proj('miller','lat', [-30 30], 'lon',[100 300]);
h=m_pcolor(lon, lat, EOFs(:,:,3)*1.0e2); set(h,'linestyle','none'); hold on;
h=m_pcolor(lon+360, lat, EOFs(:,:,3)*1.0e2); set(h,'linestyle','none'); hold on;
set(h,'linestyle','none');
hold on;
m_coast('color',[.1 .1 .1]);

m_grid('box', 'on', 'xtick', 5, 'ytick', 5, 'tickdir', 'in', 'yaxislocation', 'left', ...
    'Fontname', 'Times New Roman', 'Fontsize', 10, 'color', RdYlBu(60,:));

% drawing quiver arrows (on top of Nino3.4 area)
h5 = m_vec(.15, lon, lat, arrows(:,:,3), arrows(:,:,4), [0 0 0], ...
    'headlength',3.5, 'shaftwidth', .9, 'Edgecolor', [1 1 1], 'linewidth', .5);
h6 = m_vec(.15, lon, lat, arrows(:,:,3), arrows(:,:,4), [0 0 0], ...
    'headlength',3.5, 'shaftwidth', .9, 'Edgecolor', [1 1 1]);

set(gca, 'clim', [-3 3]); hold on;

h3 = colorbar('color', [.1922 .2118 .5843], 'box', 'on', 'tickdirection', 'in', ....
        'Fontname', 'Times New Roman', 'Fontsize', 10);
h = ylabel(h3, '[10^{-2} N m^{-2}]', 'color', RdYlBu(60,:), ...
    'Fontname', 'Times New Roman', 'Fontsize',10);
set(h3, 'YTick', linspace(-3, 3, 5));

% set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
% directory = 'E:\2017 Herbstsemester\Master Seminar 2\';
% print('-dpng','-r1000', [directory 'EOF2_wind_stress_anomalies_v11']);


%% [works] subplots of first two PCs
% load('workspace_plotting_testmap_rev5.mat'); 
Reds = cbrewer('seq', 'Reds', 16, 'PCHIP');
Blues = cbrewer('seq', 'Blues', 16, 'PCHIP');

figure(2);
% str1 = {'r = 0.748'};
% text(2007, 1.25, str1, 'Fontname', 'Times New Roman', 'Fontsize', 10, 'color', RdYlBu(1,:));

% strong el nino events
xbars = [years_2(37) years_2(60)];    
h5 = patch([xbars(1) xbars(1), xbars(2) xbars(2)], [-5 5 5 -5], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [years_2(216) years_2(239)];   
h5 = patch([xbars(1) xbars(1), xbars(2) xbars(2)], [-5 5 5 -5], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [years_2(432) years_2(455)];  
h5 = patch([xbars(1) xbars(1), xbars(2) xbars(2)], [-5 5 5 -5], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;

% strong la nina events
xbars = [years_2(110) years_2(130)];    
h5 = patch([xbars(1) xbars(1), xbars(2) xbars(2)], [-5 5 5 -5], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [years_2(335) years_2(354)];   
h5 = patch([xbars(1) xbars(1), xbars(2) xbars(2)], [-5 5 5 -5], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [years_2(368) years_2(391)];   
h5 = patch([xbars(1) xbars(1), xbars(2) xbars(2)], [-5 5 5 -5], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;

h3 = plot(linspace(1979,2016,456), PC(1,:), 'color', [0 0 0], 'linewidth', .7); hold on;
h4 = plot(linspace(1979,2016,456), PC(2,:), 'color', RdYlBu(1,:), 'linewidth', .7); hold on;


% create fancy plot
hXLabel = xlabel('Year');
hYLabel = ylabel('Amplitude');
set([hXLabel], 'FontName'   , 'Times New Roman', 'Fontsize', 9);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 9);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'          , ...
  'YMinorTick'      , 'off'          , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -2.:1:2.5       , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'LineWidth'       , .5);
  % 'GridLineStyle'   , '--'          , ...


ylim([-2.5 2.5]);
xlim([1979 2017]);
pbaspect([7 1 1]);                        % aspect ratios: x, y, z
legend([h3 h4], 'Niño3.4', 'PC2', ...
    'location', 'southoutside', 'orientation', 'horizontal');
clear hTitle hXLabel hYLabel h1 h2;

% set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
% directory = 'E:\2017 Herbstsemester\Master Seminar 2\';
% print('-dpng','-r1000', [directory 'PC_wind_stress_scores_v14']);


%% [works] Oceanic Nino Index Era-Interim
load('workspace_plotting_testmap_rev5.mat'); 
load('plotting_NOAA_ERSST_nino34_index.mat', 'nino');
Reds = cbrewer('seq', 'Reds', 16, 'PCHIP');
Blues = cbrewer('seq', 'Blues', 16, 'PCHIP');

% calculating 3-month running mean for ONI
nino(2,:) = movmean(nino, 3);               % nino(1,:) = monthly mean values
                                            % nino(2,:) = 3-month running mean
figure(5);
% function anomaly which colours areas above and below zero
[h1, h2] = anomaly(years_2, nino(2,:), 'edgecolor', 'none', 'topcolor', RdYlBu(1,:), ...
    'bottomcolor', RdYlBu(60,:)); hold on; 

% plot the +/-0.5 degree threshold lines
hline1 = line([1979 2016], [.5 .5]); % draw threshold line for nino
set(hline1, 'color', RdYlBu(1,:), 'linestyle', '--','linewidth', .5);
hline2 = line([1979 2016], [-.5 -.5]); % draw threshold line for nino
set(hline2, 'color', RdYlBu(60,:), 'linestyle', '--','linewidth', .5);

% str1 = {'El Niño'};
% text(1999, 2.5, str1, 'Fontname', 'Times New Roman', 'Fontsize', 10, 'color', RdYlBu(5,:));
% % drawing arrow indicating nino event with the drawArrow function
% drawArrow([1998.8; 2.5],[1998; 1.7], [0 0 0]);
% str2 = {'La Niña'};
% text(2011.5, -1.5, str2, 'Fontname', 'Times New Roman', 'Fontsize', 10, 'color', RdYlBu(60,:));
% drawArrow([2011.3; -1.5],[2010.5; -1.1], [0 0 0]);


% create fancy plot
hXLabel = xlabel('Year');
hYLabel = ylabel('ONI [\circC]');
set([hXLabel], 'FontName'   , 'Times New Roman', 'Fontsize', 9);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 9);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'          , ...
  'YMinorTick'      , 'off'          , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -3:1:4       , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'LineWidth'       , .5);
  % 'GridLineStyle'   , '--'          , ...


ylim([-2 2]);
xlim([1979 2016]);
pbaspect([7 1 1]);                        % aspect ratios: x, y, z
% legend([h3 h4], 'PC1', 'PC2', ...
%     'location', 'southwest', 'orientation', 'horizontal');
clear hTitle hXLabel hYLabel h1 h2;

set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = 'E:\2017 Frühlingssemester\Master Seminar 1\ProposalTemplate\figures\';
print('-dpng','-r1000', [directory 'ONI_index_from_NOAA']);


% correlation coefficients
r_pc_nino34 = corrcoef(scores_tau_x(1,:), nino(1,:));
    % r_pc_nino34 = 0.748 (as in McGregor et al., 2014
r_pc_oni = corrcoef(scores_tau_x(1,:), nino(2,:));
    % r_pc_oni    = 0.770
    
    











