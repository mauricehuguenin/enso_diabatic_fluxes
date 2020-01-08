%% Purpose: Plotting time series of ocean diagnostics for my 1979-2016 
%           ACCESS-OM2-025 JRA55-do iaf run
%           a) GMSSTa and OHCa
%           b) Niño3.4 and WWVa
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %                         
%                        Maurice Huguenin-Virchaux                        %                         
%                     m.huguenin-virchaux@unsw.edu.au                     %                         
%                          18.09.2019, 17:27 AEST                         %                         
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 

first = 'output020'; % first output in time period, output020 = year 1979
first_num = 20;
ts_len = 57; % last output entry, i.e. ts_len = 25 -> output025


%% [1.4 seconds] load in workspace
tic;
% load workspace
f1 = 'pEXP9601_timeseries_wwv_transport_20S_20N';
f2 = 'pEXP9601_timeseries_wwv_transport_8S_8N';
f3 = 'pEXP9601_timeseries_wwv_transport_5S_5N';
base = '/srv/ccrc/data67/z5180028/MSC_thesis_access_output/'; % new data from ACCESS stored in here
base2 = '/srv/ccrc/data15/z5180028/MSC_thesis_mom_output/EXP0_control_run/' % old data from MOM
load('/srv/ccrc/data15/z5180028/workspace_EXP1_analysis_composite_ninos_rev3.mat');
load('/home/z5180028/MSC_thesis/access_matlab_scripts/workspace_regression_patterns_PC1_equal_nino34_rev2.mat', 'PC1');



antarctica = nature;
RdBu_short  = cbrewer('div', 'RdBu', 21, 'PCHIP');
RdYlGn  = cbrewer('div', 'RdYlGn', 60, 'PCHIP');
Reds  = cbrewer('seq', 'Reds', 16, 'PCHIP');
Blues  = cbrewer('seq', 'Blues', 16, 'PCHIP');


wwv_mask = wwv_mask(91:811,478:518);
wwv_mask_2S = wwv_mask_2S(91:811,478:518);
wwv_mask_2S_borneo = wwv_mask_2S_borneo(91:811,478:518);
lon = lon(91:811,478:518);
lat = lat(91:811,478:518);
% testmap(lon, lat, wwv_mask_2S_borneo);

clear ans antarctica basin_mask bathymetry filename2017 wwv_mask wwv_mask_2S

toc;         


%% [7.5 seconds] load in N34 value 
tic;
% load in potential temperature
thetao = squeeze(permute(getnc([base first '/ocean/ocean_snap.nc'], ...
            'temp', [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
         % [month, depth, lat, lon], [end], [stride]
         % size(thetao) % = [721, 41, 50, 12)
         
for i = (first_num+1):ts_len % load in subsequent years and concatenate
    if i < 10
        a = ['0' num2str(i)];
    else
        a = num2str(i);
    end
    thetao = cat(3,thetao, squeeze(permute(getnc([base 'output0' a '/ocean/ocean_snap.nc'], ...
           'temp', [1,1,478,440], [12,1,518,640], [1,1,1,1]), [3 2 1])));     
end

fprintf('~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ \n');
fprintf(['thetao:' first 'to output0' num2str(i)])
size(thetao) % ok, it gets concatenated correctly    


nino = squeeze(nanmean(squeeze(nanmean(thetao, 1)),1));
% plot(nino); hold on;

% create anomaly for simulated nino
clear ninoc;
for t = 1:12 % loop over months
    ninoc(t)         =  nanmean(nino(t:12:end)); % calculate mean over all Jan/Feb/etc.
end

% repeat climatology n-times to make it as long as the original
ninoc         =  repmat(ninoc,1,(ts_len-first_num+1));

ninoa = nino - ninoc; % create anomaly but only select Jan 1979-Dec 2016 period

plot(ninoa);



toc;


%% [0.00/44/899 seconds load in WWV anomaly
tic;
% load in previously calculated WWV time series from the script
% where I calculated the warm water volume budget


directory = '/home/z5180028/MSC_thesis/access_figures/';
load([directory 'WMT_time_series_1979-2016.mat'], 'wwv');
% remove first data entry as that one was from the 1st snap file
% so that I have once again 1x456 data entries. The extention with the one
% data point at the beginning was necessary for the derivative in the other
% script
wwv = wwv(2:457);
% plot(wwv)


% create anomalies
clear wwvc
for t = 1:12 % loop over months
    wwvc(t)         =  nanmean(wwv(t:12:end)); % calculate mean over all Jan/Feb/etc.
end

% repeat climatology n-times to make it as long as the wwv data
wwvc         =  repmat(wwvc,1,(ts_len-first_num+1));


% calculate anomalies
wwva = wwv - wwvc;


plot(ninoa*1e14, 'color', [0 0 0]); hold on; plot(wwva, 'color', [203,24,30]/255);
toc;


%% [0.28/2.45 seconds] load in OHC anomaly
tic;
% load in potential temperature
ocean_heat = getnc([base first '/ocean/ocean_scalar.nc'],...
                 'total_ocean_heat'); 
         % [month, depth, lat, lon], [end], [stride]
         % size(thetao) % = [721, 41, 50, 12)
         
for i = (first_num+1):ts_len % load in subsequent years and concatenate
    if i < 10
        a = ['0' num2str(i)];
    else
        a = num2str(i);
    end
    ocean_heat = cat(1,ocean_heat, getnc([base 'output0' a '/ocean/ocean_scalar.nc'], ...
           'total_ocean_heat'));     
end

fprintf('~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ \n');
fprintf(['ocean_heat:' first 'to output0' num2str(i)])
size(ocean_heat) % ok, it gets concatenated correctly    


% unit is in [1e25 J] as shown in metadata of ocean_scalar.nc
ocean_heat = ocean_heat * 1e25;

% create anomaly for simulated ocean heat content
clear ocean_heatc;
for t = 1:12 % loop over months
    ocean_heatc(t)         =  nanmean(ocean_heat(t:12:end)); % calculate mean over all Jan/Feb/etc.
end

% repeat climatology n-times to make it as long as the original
ocean_heatc         =  permute(repmat(ocean_heatc,1,(ts_len-first_num+1)), [2 1]);

ocean_heata = ocean_heat - ocean_heatc; % create anomaly but only select Jan 1979-Dec 2016 period

plot(ocean_heata);




toc;


%% [0.29/0.79 seconds] load in GMSST anomaly
tic;
% load in potential temperature
temp = getnc([base first '/ocean/ocean_scalar.nc'],...
                 'temp_surface_ave'); 
         
for i = (first_num+1):ts_len % load in subsequent years and concatenate
    if i < 10
        a = ['0' num2str(i)];
    else
        a = num2str(i);
    end
    temp = cat(1,temp, getnc([base 'output0' a '/ocean/ocean_scalar.nc'], ...
           'temp_surface_ave'));     
end

fprintf('~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ \n');
fprintf(['temp:' first 'to output0' num2str(i)])
size(temp) % ok, it gets concatenated correctly    



% create anomaly for simulated ocean heat content
clear tempc;
for t = 1:12 % loop over months
    tempc(t)         =  nanmean(temp(t:12:end)); % calculate mean over all Jan/Feb/etc.
end

% repeat climatology n-times to make it as long as the original
tempc         =  permute(repmat(tempc,1,(ts_len-first_num+1)), [2 1]);

tempa = temp - tempc; % create anomaly but only select Jan 1979-Dec 2016 period

plot(tempa);

toc;


%% create running mean before plotting

% creating 5-month running mean of all time series as in Meinen and McPhaden 2000
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
ninoa = filter(b,1,ninoa);
wwva =filter(b,1,wwva);
ocean_heata = filter(b,1,ocean_heata);
tempa = filter(b,1,tempa);


subplot(2,2,1);plot(ninoa); hold on; title('Niño3.4') 
subplot(2,2,2);plot(wwva); hold on; title('WWV') 
subplot(2,2,3);plot(ocean_heata); hold on; title('OHC') 
subplot(2,2,4);plot(tempa); hold on; title('GMSST')  

toc;


%% preparing NOAA_PMEL_WWVa.mat dataset (observations) & clean up workspace
srv2 = '/home/z5180028/MSC_thesis/access_matlab_scripts/'
load([srv2 'workspace_NOAA_PMEL_WWVa.mat']);

wwv_obs(1:12) = nan;
wwv_obs(13:(ts_len-first_num+1)*12) = NOAA_PMEL_WWVa(1:(ts_len-first_num+1)*12-12,2);


plot(wwva); hold on; plot(wwv_obs);

corrcoef(wwva, wwv_obs)

clearvars -except ninoa wwva ocean_heata tempa RdYlBu Reds Blues PC1 ...
    first_num ts_len wwv_obs directory


%% plotting routine with two subplots
boom;
directory = '/home/z5180028/MSC_thesis/access_figures/';
% load([directory 'Diagnostics_time_series_n34_and_wwv_and_ohc_and_gmsst.mat']);

%           a) GMSSTa and OHCa
%           b) N34 and WWV as well as observed N34 and observed WWV

ticks = linspace(1,length(wwva),length(wwva));
figure('units', 'pixels', 'position', [0 0 2300 800]);

subplot(2,1,1)
% -------------------------------------------------------------------------
[hAx,h1,h2] = plotyy(ticks,tempa,ticks,ocean_heata);
set(hAx(1),'YLim',[-.5 .5],'YTick',[-.5:.25:.5]);
set(hAx(2),'YLim',[-10e22 10e22],'YTick',[-10e22:5e22:10e22]);

set(h1, 'color', RdYlBu(60,:), 'linewidth', 2.5); hold on;
set(h2, 'color', RdYlBu(15,:), 'linewidth', 2.5); hold on;
set(gca, 'xtick',[]);

ylabel(hAx(1),['[$^{\circ}$C]'],'interpreter','latex', 'Fontsize', 21,...
    'color', RdYlBu(60,:), 'linewidth', 2.5); hold on;

ylabel(hAx(2), '[J]', 'interpreter', 'latex','FontSize',21, ...
    'color', RdYlBu(13,:), 'linewidth', 2.5); hold on;

% set([hYLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(60,:));

    set(gca, ...
      'Box'             , 'off'         , ...
      'TickDir'         , 'out'         , ...
      'TickLength'      , [.01 .01]     , ...
      'XMinorTick'      , 'off'         , ...
      'YMinorTick'      , 'off'         , ...
      'YGrid'           , 'on'          , ...
      'XColor'          , RdYlBu(60,:)  , ...
      'YColor'          , RdYlBu(60,:)  , ...
      'YTick'           , -5:.25:5  , ...     % define grid, every 1000 m a grid line
      'ticklabelinterpreter', 'latex'   , ...
      'LineWidth'       , 1.25);
      % 'GridLineStyle'   , '--'          , ...
set(gca,'fontsize',21)
secondAxes = findobj('Type','axes','Color','none')
set(secondAxes,'fontsize',21, 'ticklabelinterpreter', 'latex');

text(-30, 0.63, 'a)', 'interpreter', 'latex', 'Fontsize', 20, ...
     'color', [0 0 0]);

set(hAx,'XLim',[1 length(wwva)],'XTick',[0:24:length(wwva)])
set(hAx,'xticklabel',({'','','','',...
    '','','','', '','','','','',...
    '','','','','',''}))

set(hAx,'Ycolor', RdYlBu(60,:));

% legend
h99 = legend([h1 h2], 'GMSST anomaly', 'OHC anomaly', ...
    'location', 'northwest', 'orientation', 'horizontal');
set(h99, 'interpreter', 'latex', 'fontsize', 14, 'textcolor', RdYlBu(60,:), ...
    'EdgeColor', [.83 .83 .83]);
 
subplot(2,1,2)
% -------------------------------------------------------------------------
clear hAx
[hAx,h3,h4] = plotyy(ticks,ninoa,ticks,wwva);
set(hAx(1),'YLim',[-5 5] ,'YTick',[-.5:.25:.5]);
set(hAx(2),'YLim',[-4e14 4e14],'YTick',[-4e14:2e14:4e14]);

%hold both axes
hold(hAx(1)); hold(hAx(2));
%plot additional data to the axes
h5 = plot(hAx(1),ticks,PC1(1:length(wwva)),'color', RdYlBu(50,:), 'linewidth', 1.5, 'linestyle', '--');
h6 = plot(hAx(2),ticks,wwv_obs,'color', RdYlBu(10,:), 'linewidth', 1.5, 'linestyle', '--');



set(h3, 'color', RdYlBu(50,:), 'linewidth', 2.5); hold on;
set(h4, 'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;



ylabel(hAx(1),['[$^{\circ}$C]'],'interpreter','latex', 'Fontsize', 21,...
    'color', RdYlBu(50,:), 'linewidth', 2.5); hold on;

ylabel(hAx(2), '[m$^3$]', 'interpreter', 'latex','FontSize',21, ...
    'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;

hXLabel = xlabel('Year');
set([hXLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(60,:));

    set(gca, ...
      'Box'             , 'off'         , ...
      'TickDir'         , 'out'         , ...
      'TickLength'      , [.01 .01]     , ...
      'XMinorTick'      , 'off'         , ...
      'YMinorTick'      , 'off'         , ...
      'YGrid'           , 'on'          , ...
      'XColor'          , RdYlBu(60,:)  , ...
      'YColor'          , RdYlBu(60,:)  , ...
      'ticklabelinterpreter', 'latex'   , ...
      'YTick'           , -10:2.5:10  , ...     % define grid, every 1000 m a grid line
      'LineWidth'       , 1.25);
      % 'GridLineStyle'   , '--'          , ...
set(gca,'fontsize',21)
secondAxes = findobj('Type','axes','Color','none')
set(secondAxes,'fontsize',21, 'ticklabelinterpreter', 'latex');

text(-30, 6, 'b)', 'interpreter', 'latex', 'Fontsize', 20, ...
     'color', [0 0 0]);

set(hAx,'XLim',[1 length(wwva)],'XTick',[0:24:length(wwva)])

% evaluate correlation:
a = corrcoef(ninoa, PC1(1:length(wwva)));
a(2) % correlation nino vs. nino_obs
b = corrcoef(wwva(13:end), wwv_obs(13:end));
b(2) % correlation nino vs. nino_obs

set(hAx,'Ycolor', RdYlBu(60,:));
set(hAx,'xticklabel',({'','1981','','1985',...
    '','1989','','1993', '','1997','','2001','',...
    '2005','','2009','','2013','','2017'}))


text(16, -1.5, ['r = ' num2str(round(a(2),2))], 'interpreter', 'latex', 'Fontsize', 21, ...
     'color', RdYlBu(50,:));
text(250, -3.5, ['r = ' num2str(round(b(2),2))], 'interpreter', 'latex', 'Fontsize', 21, ...
     'color', RdYlBu(10,:));

% legend
h99 = legend([h3 h5 h4 h6], 'N34 anomaly   ', 'N34 anomaly NOAA ERSST v4   ', 'WWV anomaly   ', 'WWV anomaly NOAA PMEL   ', ...
    'location', 'northwest', 'orientation', 'horizontal');
set(h99, 'interpreter', 'latex', 'fontsize', 14, 'textcolor', RdYlBu(60,:), ...
    'EdgeColor', [.83 .83 .83]);

% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = '/home/z5180028/MSC_thesis/access_figures/';
print('-dpng','-r300', [directory 'Diagnostics_timeseries_1979-2016']);

% save([directory 'Diagnostics_time_series_n34_and_wwv_and_ohc_and_gmsst.mat']);

%% plotting routine only the second subplot
boom;
directory = '/home/z5180028/MSC_thesis/access_figures/';
% load([directory 'Diagnostics_time_series_n34_and_wwv_and_ohc_and_gmsst.mat']);

%           a) GMSSTa and OHCa
%           b) N34 and WWV as well as observed N34 and observed WWV

ticks = linspace(1,length(wwva),length(wwva));
figure('units', 'pixels', 'position', [0 0 2300 600]);
 
% -------------------------------------------------------------------------
clear hAx
[hAx,h3,h4] = plotyy(ticks,ninoa,ticks,wwva);
set(hAx(1),'YLim',[-5 5] ,'YTick',[-.5:.25:.5]);
set(hAx(2),'YLim',[-4e14 4e14],'YTick',[-4e14:2e14:4e14]);

%hold both axes
hold(hAx(1)); hold(hAx(2));
% plot additional data to the axes, this time the observed values
% the observed N34 is the same as the PC1 time series from the workspace
% that I used in the regression analysis, hence the name 'PC1'
h5 = plot(hAx(1),ticks,PC1(1:length(wwva)),'color', RdYlBu(50,:), 'linewidth', 1.5, 'linestyle', '--');
h6 = plot(hAx(2),ticks,wwv_obs,'color', RdYlBu(10,:), 'linewidth', 1.5, 'linestyle', '--');


% set colours for time series
set(h3, 'color', RdYlBu(50,:), 'linewidth', 2.5); hold on;
set(h4, 'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;



ylabel(hAx(1),['[$^{\circ}$C]'],'interpreter','latex', 'Fontsize', 21,...
    'color', RdYlBu(50,:), 'linewidth', 2.5); hold on;

ylabel(hAx(2), '[m$^3$]', 'interpreter', 'latex','FontSize',21, ...
    'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;

hXLabel = xlabel('Year');
set([hXLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(60,:));

    set(gca, ...
      'Box'             , 'off'         , ...
      'TickDir'         , 'out'         , ...
      'TickLength'      , [.01 .01]     , ...
      'XMinorTick'      , 'off'         , ...
      'YMinorTick'      , 'off'         , ...
      'YGrid'           , 'on'          , ...
      'XColor'          , RdYlBu(60,:)  , ...
      'YColor'          , RdYlBu(60,:)  , ...
      'ticklabelinterpreter', 'latex'   , ...
      'YTick'           , -10:2.5:10  , ...     % define grid, every 1000 m a grid line
      'LineWidth'       , 1.25);
      % 'GridLineStyle'   , '--'          , ...
set(gca,'fontsize',21)
secondAxes = findobj('Type','axes','Color','none')
set(secondAxes,'fontsize',21, 'ticklabelinterpreter', 'latex');

text(-30, 6, 'b)', 'interpreter', 'latex', 'Fontsize', 20, ...
     'color', [0 0 0]);

set(hAx,'XLim',[1 length(wwva)],'XTick',[0:24:length(wwva)])

% evaluate correlation:
a = corrcoef(ninoa, PC1(1:length(wwva)));
a(2) % correlation nino vs. nino_obs
b = corrcoef(wwva(13:end), wwv_obs(13:end));
b(2) % correlation nino vs. nino_obs

set(hAx,'Ycolor', RdYlBu(60,:));
set(hAx,'xticklabel',({'','1981','','1985',...
    '','1989','','1993', '','1997','','2001','',...
    '2005','','2009','','2013','','2017'}))


text(16, -1.5, ['r = ' num2str(round(a(2),2))], 'interpreter', 'latex', 'Fontsize', 21, ...
     'color', RdYlBu(50,:));
text(250, -3.5, ['r = ' num2str(round(b(2),2))], 'interpreter', 'latex', 'Fontsize', 21, ...
     'color', RdYlBu(10,:));
text(300, -3.5, ['Year'], 'interpreter', 'latex', 'Fontsize', 21, ...
     'color', RdYlBu(60,:));


% legend
h99 = legend([h3 h5 h4 h6], 'N34 anomaly   ', 'N34 anomaly NOAA ERSST v4   ', 'WWV anomaly   ', 'WWV anomaly NOAA PMEL   ', ...
    'location', 'northwest', 'orientation', 'horizontal');
set(h99, 'interpreter', 'latex', 'fontsize', 14, 'textcolor', RdYlBu(60,:), ...
    'EdgeColor', [.83 .83 .83]);

% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = '/home/z5180028/MSC_thesis/access_figures/';
print('-dpng','-r300', [directory 'Diagnostics_timeseries_1979-2016_second_subplot']);

% save([directory 'Diagnostics_time_series_n34_and_wwv_and_ohc_and_gmsst.mat']);
