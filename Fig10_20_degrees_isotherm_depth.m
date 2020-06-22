%% Purpose: Plotting the 20 degree Celsius isotherm depth during
%           (1) OBS, (2) climatological, (3) El Nino, (4) La Nina 
%           (5) realistic El Nino, (6) realistic La Nina
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %                         
%                        Maurice Huguenin-Virchaux                        %                         
%                     m.huguenin-virchaux@unsw.edu.au                     %                         
%                          25.09.2019,11:32 AEDT                          %                         
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %  
first = 'output020'; % first output in time period, output020 = year 1979
first_num = 20;
ts_len = 57; % last output entry, i.e. ts_len = 25 -> output025


%% [2.3 sec.] load in workspace
tic;
% load workspace
f1 = 'pEXP9601_timeseries_wwv_transport_20S_20N';
f2 = 'pEXP9601_timeseries_wwv_transport_8S_8N';
f3 = 'pEXP9601_timeseries_wwv_transport_5S_5N';
base = '/srv/ccrc/data67/z5180028/MSC_thesis_access_output/'; % new data from ACCESS stored in here
basem = '/srv/ccrc/data15/z5180028/MSC_thesis_mom_output/pEXP9601_restart000_windstress/'; % new data from ACCESS stored in here
base2 = '/srv/ccrc/data15/z5180028/MSC_thesis_mom_output/EXP0_control_run/' % old data from MOM
load('/srv/ccrc/data15/z5180028/workspace_EXP1_analysis_composite_ninos_rev3.mat');
antarctica = nature;
RdBu_short  = cbrewer('div', 'RdBu', 21, 'PCHIP');
RdYlGn  = cbrewer('div', 'RdYlGn', 60, 'PCHIP');
Reds  = cbrewer('seq', 'Reds', 16, 'PCHIP');
Blues  = cbrewer('seq', 'Blues', 16, 'PCHIP');

lon = lon(91:811,478:518);
lat = lat(91:811,478:518);

clear ans antarctica basin_mask bathymetry filename2017 wwv_mask wwv_mask_2S

clear lev_bnds % delete vertical depth levels from MOM025
% load in vertical depth levels from ACCESS
lev_bnds = getnc([base first '/ocean/ocean.nc'],'st_ocean');   
         % [month, depth, lat, lon], [end], [stride]

toc; 


%% [38 seconds] calculate 20 degree isotherm for ACCESS for December 1997
tic;    
% load in potential temperature
thetao = squeeze(nanmean(permute(getnc([base first '/ocean/ocean.nc'], ...
           'temp', [1,1,498,91], [12,-1,499,811], [1,1,1,1]), [4 3 2 1]),2));   
         % [month, depth, lat, lon], [end], [stride]
         % take mean over the two equatorial grid cells -0.125 and 0.125
size(thetao) 
         
for i = (first_num+1):ts_len % load in subsequent years and concatenate
    if i < 10
        a = ['0' num2str(i)];
    else
        a = num2str(i);
    end
    thetao = cat(3,thetao, squeeze(nanmean(permute(getnc([base 'output0' a '/ocean/ocean.nc'], ...
           'temp', [1,1,498,91], [12,-1,499,811], [1,1,1,1]), [4 3 2 1]),2)));     
end

fprintf('~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ \n');
fprintf(['thetao:' first 'to output0' num2str(i)])
size(thetao) % ok, it gets concatenated correctly    

% calculate the mean over SON 1997 and SON 1998 for both the El Nino and
% La Nina as well as calculating the climatology over the full period
isotherm_iEN = squeeze(nanmean(thetao(:,:,225:227),3)) - 273.15;
isotherm_iLN = squeeze(nanmean(thetao(:,:,237:239),3)) - 273.15;
isotherm_iclim = squeeze(nanmean(thetao,3)) - 273.15;


[Xt,Zt] = ndgrid(lon(:,1),lev_bnds); % prepare grid


% climatology from MOM025
    [c4,h4] = contour(Xt,-Zt,isotherm_iclim,[1:35],'color',[.6 .6 .6],'linewidth',.5); hold on;
    [c1,h1] = contour(Xt,-Zt,isotherm_iEN,[20 20],'color','r','linewidth',2.5); hold on;
    [c2,h2] = contour(Xt,-Zt,isotherm_iLN,[20 20],'color','b','linewidth',2.5); hold on;
    [c3,h3] = contour(Xt,-Zt,isotherm_iclim,[20 20],'color','k','linewidth',2.5); hold on;
% clabel(c4,h4)
ylim([-300 0]);
xlim([-220 -80]);

toc;

% clear workspace
clearvars -except isotherm_iEN isotherm_iLN isotherm_iclim RdYlBu ...
    RdBu lon lat lev_bnds Xt Zt basem;

% ok great, now I have the climatological, El Nino and La Nina thermocline
% depth-longitude transects in the ACCESS model

% what I need now: OBS, i.e. WOA data


%% [9.5 seconds] calculate 20 degree isotherm for MOM025 for December 1997
tic;    
% load in potential temperature
thetao = squeeze(nanmean(permute(getnc([basem 'output000/ocean.nc'], ...
           'temp', [1,1,498,91], [12,-1,499,811], [1,1,1,1]), [4 3 2 1]),2));   
         % [month, depth, lat, lon], [end], [stride]
         % take mean over the two equatorial grid cells -0.125 and 0.125
size(thetao) 
         
for i = 1:5 % load in subsequent years and concatenate
    if i < 10
        a = ['0' num2str(i)];
    else
        a = num2str(i);
    end
    thetao = cat(3,thetao, squeeze(nanmean(permute(getnc([basem 'output0' a '/ocean.nc'], ...
           'temp', [1,1,498,91], [12,-1,499,811], [1,1,1,1]), [4 3 2 1]),2)));     
end

fprintf('~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ \n');
fprintf(['thetao: output000 to output00' num2str(i)])
size(thetao) % ok, it gets concatenated correctly    

% calculate the mean over SON 1997 and SON 1998 for both the El Nino and
% La Nina as well as calculating the climatology over the full period
isotherm_iENm = squeeze(nanmean(thetao(:,:,21:23),3));
isotherm_iLNm = squeeze(nanmean(thetao(:,:,33:35),3));
isotherm_iclimm = squeeze(nanmean(thetao,3));

load('/srv/ccrc/data15/z5180028/workspace_EXP1_analysis_composite_ninos_rev3.mat', 'lev_bnds');

[Xm,Zm] = ndgrid(lon(:,1),lev_bnds); % prepare grid


% climatology from MOM025
    [c4,h4] = contour(Xm,-Zm,isotherm_iclimm,[1:35],'color',[.6 .6 .6],'linewidth',.5); hold on;
    [c1,h1] = contour(Xm,-Zm,isotherm_iENm,[20 20],'color','r','linewidth',2.5); hold on;
    [c2,h2] = contour(Xm,-Zm,isotherm_iLNm,[20 20],'color','b','linewidth',2.5); hold on;
    [c3,h3] = contour(Xm,-Zm,isotherm_iclimm,[20 20],'color','k','linewidth',2.5); hold on;
% clabel(c4,h4)
ylim([-300 0]);
xlim([-220 -80]);

toc;

% clear workspace
clearvars -except isotherm_iENm isotherm_iLNm isotherm_iclimm RdYlBu ...
    RdBu lon lat lev_bnds Xm Zm basem Xt Zt isotherm_iclim isotherm_iEN isotherm_iLN;

% ok great, now I have the climatological, El Nino and La Nina thermocline
% depth-longitude transects in the ACCESS model

% what I need now: OBS, i.e. WOA data


%% [0.19 seconds] Import WOA equatorial 20 degree isotherm data

tic;
base2 = '/srv/ccrc/data15/z5180028/';
f1 = 'WOA13/woa13_decav_t00_04v2_wwv_region.nc';

lon = getnc([base2 f1], 'lon');
lat = getnc([base2 f1], 'lat');
depth = getnc([base2 f1], 'depth');

[Xw,Zw] = ndgrid(lon,depth);


[lat, lon] = meshgrid(lat, lon);
WOA = permute(getnc([base2 f1], 't_an'), [3 2 1]); % [lon, lat, depth]
% WOA = [lon, lat, depth]

% testmap(lon, lat, WOA(:,:,1))

WOA = squeeze(nanmean(WOA(:,20:21,:),2)); % take equatorial transect as the -0.125 degree grid cell
[c1,h1] = contour(Xw,-Zw,WOA,[20 20], 'color', 'k','linewidth',4.5); hold on;
ylim([-300 0]);
toc;


%% [20 seconds] Import ORA-S5 equatorial 20 degree isotherm data
tic;
% create colours
RdBu = flipud(cbrewer('div', 'RdBu', 20, 'PCHIP'));
RdYlBu = cbrewer('div', 'RdBu', 60, 'PCHIP');

base2 = '/srv/ccrc/data15/z5180028/';
f1 = 'ORA-S5/votemper_ORAs5_merged_r1x1_wwv_region_1996-2001.nc';
f2 = 'ORA-S5/votemper_ORAs5_merged_r1x1_wwv_region_1979-2017.nc';

lon = getnc([base2 f1], 'lon');
lat = getnc([base2 f1], 'lat');
depth = getnc([base2 f1], 'deptht');
[Xo,Zo] = ndgrid(lon,depth);
[lat, lon] = meshgrid(lat, lon);


% [lat, lon] = meshgrid(lat, lon);
ORA = permute(getnc([base2 f1], 'votemper'), [4 3 2 1]); % [lon, lat, depth]
ORA_clim = permute(getnc([base2 f2], 'votemper'), [4 3 2 1]); % [lon, lat, depth]
size(ORA) % [lon, lat, depth, time]




% ... and over the EN 1997 period of SON
ORA_EN97 = squeeze(nanmean(squeeze(nanmean(ORA(:,5,:,24:26),2)),3));
% ... and over the LN 1998 period of SON
ORA_LN98 = squeeze(nanmean(squeeze(nanmean(ORA(:,5,:,33),2)),3));

% testmap(lon, lat, WOA(:,:,1))
[c2,h1] = contour(Xw,-Zw,WOA,[2:30], ...
        'color', [.83 .83 .83],'linewidth',1); hold on;
[c2,h1] = contour(Xw,-Zw,WOA,[20 20], ...
        'color', [0 0 0],'linewidth',3); hold on;
[c1,h2] = contour(Xo,-Zo,ORA_EN97,[20 20], ...
        'color', 'blue','linewidth',3, 'linestyle', '--'); hold on;
[c2,h3] = contour(Xo,-Zo,ORA_LN98,[20 20], ...
        'color', 'red','linewidth',3, 'linestyle', '--'); hold on;
ylim([-300 0]);
% xlim([-220 -80]);
h66 = legend([h1 h2 h3], 'WOA clim.', ...
    'ORAS5 EN 97', ...
    'ORAS5 LN 98', ...
    'location', 'southeast', 'orientation', 'vertical');
set(h66, 'interpreter', 'latex', 'fontsize', 12, 'textcolor', RdYlBu(60,:));
toc;


%% [7.6 seconds] plotting routine
tic;
figure;
warning off
set(gcf,'Position',[1          36        1000         350]);

i=4
% El Nino 1997 
    [c4,h4] = contour(Xt,-Zt,isotherm_iEN,[20 20],'color',RdYlBu(8,:),'linewidth',2.5); hold on;
% La Nina 1998
    [c5,h5] = contour(Xt,-Zt,isotherm_iLN,[20 20],'color',RdYlBu(52,:),'linewidth',2.5); hold on;
% El Nino 1997 from ORAS5
    [c6,h6] = contour(Xo-350,-Zo,ORA_EN97,[20 20], ...
    'color', RdYlBu(8,:),'linewidth',1.5, 'linestyle', '--'); hold on;
% La Nina 1998 from ORAS5
    [c7,h7] = contour(Xo-350,-Zo,ORA_LN98,[20 20], ...
    'color', RdYlBu(52,:),'linewidth',1.5, 'linestyle', '--'); hold on;
    [c2,h2] = contour(Xw-356,-Zw,WOA,[20 20], ...
        'color', 'k','linewidth', 1.5, 'linestyle', '--'); hold on;
    [c1,h1] = contour(Xt,-Zt,isotherm_iclim,[20 20], ...
        'color', 'k','linewidth', 2.5); hold on;
    
ylim([-200 0]);
xlim([-220 -80]);

% set(gca,'XTickLabel',{'','$160^{\circ}$E','','$160^{\circ}$W','', ...
%                          '$120^{\circ}$W', '', '$80^{\circ}$W'}, ...
%                          'ticklabelinterpreter', 'latex', 'FontSize', 16);
set(gca,'XTickLabel',{'','$160^{\circ}$E',...                        
'','$160^{\circ}$W','','$120^{\circ}$W','','$80^{\circ}$W'}, ...
'ticklabelinterpreter', 'latex', 'FontSize', 16);

    set(gca, ...
      'Box'             , 'off'         , ...
      'TickDir'         , 'out'         , ...
      'TickLength'      , [.01 .01]     , ...
      'XMinorTick'      , 'off'         , ...
      'YMinorTick'      , 'off'         , ...
      'YGrid'           , 'on'          , ...
      'YTick'           , -5000:50:0   , ...     % define grid, every 1000 m a grid line
      'XColor'          , [.1922 .2118 .5843]  , ...
      'YColor'          , [.1922 .2118 .5843]  , ...
      'ticklabelinterpreter', 'latex'   , ...
      'LineWidth'       , 1.25);
      % 'GridLineStyle'   , '--'          , ...

set(gca,'layer','top'); % bring grid to the top layer
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 14);


% add custom labels
% - El Nino
text(-160, -180, 'La Ni\~na', 'interpreter', 'latex', 'Fontsize', 16, ...
 'color', RdYlBu(52,:));
% - La Nina
text(-140, -180, 'El Ni\~no', 'interpreter', 'latex', 'Fontsize', 16, ...
 'color', RdYlBu(8,:));
% - CLIM
text(-120, -180, 'Climatology', 'interpreter', 'latex', 'Fontsize', 16, ...
 'color', 'k');



h77 = legend([h2 h1 h6 h4 h7 h5], ...
    'WOA13 v2', ...
    'Hindcast simulation', ...
    'ORA-S5 1997', 'Hindcast simulation 1997', ...
    'ORA-S5 1998', 'Hindcast simulation 1998', ...
    'location', 'northwest', 'orientation', 'vertical');
set(h77, 'interpreter', 'latex', 'fontsize', 12, 'textcolor', RdYlBu(60,:), ...
    'EdgeColor', [.83 .83 .83]);


xlabel('Longitude','interpreter', 'latex', 'FontSize', 16);
ylabel('Depth [m]','interpreter', 'latex', 'FontSize', 16);

% save plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = '/home/z5180028/MSC_thesis/access_figures/';
print('-dpng','-r500', [directory 'thermocline_depth_transects_including_access_2']);

toc;

