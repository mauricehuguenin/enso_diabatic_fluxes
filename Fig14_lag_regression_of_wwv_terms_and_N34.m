%% Purpose: Create Warm Water Volume (WWV) budget for the period 1960-1985
%           to check if we have all the needed diagnostics for the full
%           ACCESS-OM2-025 JRA55-do iaf run 1979-2018
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %                         
%                        Maurice Huguenin-Virchaux                        %                         
%                     m.huguenin-virchaux@unsw.edu.au                     %                         
%                          29.07.2019, 08:48 AEST                         %                         
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
directory = '/home/z5180028/MSC_thesis/access_figures/';
load([directory 'WMT_time_series_1979-2016.mat']);
clearvars -except ...
    dV horizontal ITF mass forcing mixing neutral skew numerical ...
    dVc horizontalc ITFc massc forcingc mixingc neutralc skewc...
    dVa horizontala ITFa massa forcinga mixinga neutrala skewa ...
    dVaf horizontalaf ITFaf massaf forcingaf mixingaf neutralaf skewaf ...
    antarctica wwv EN LN ...
    RdYlBu RdYlGn Reds Blues;

% load in observed N34 from NOAA OI SST v4, i.e. my first principle
% component, PC1
load('/home/z5180028/MSC_thesis/access_matlab_scripts/workspace_regression_patterns_PC1_equal_nino34_rev2.mat', 'PC1');
nino = PC1; clear PC1;


%% calculate lag regression on N34 ~~~~~~~~
t = 1:456;
tL = length(t);
ts1 = nino;
% ts2 = (horizontal - mean(horizontal)) / std(horizontal);
% ts3 = (ITF - mean(ITF)) / std(ITF);
% ts4 = (mass - mean(mass)) / std(mass);
% 
% ts5 = (forcing  - mean(forcing)) / std(forcing);
% ts6 = (mixing  - mean(mixing)) / std(mixing);
% ts7 = (skew  - mean(skew)) / std(skew);
% ts8 = (neutral  - mean(neutral)) / std(neutral);
% ts9 = (numerical  - mean(numerical)) / std(numerical);
ts2 = horizontala;
ts3 = ITFa;
ts4 = massa;

ts5 = forcinga;
ts6 = mixinga;
ts7 = skewa;
ts8 = neutrala;
ts9 = numerical;

ts10 = dVa;

lags = -20:20;
regs = zeros(size(lags));
for ii=1:length(lags)
    lag = lags(ii);
    if (lag >= 0)
        ts1sh = ts1(1:(tL-lag));
        ts2sh = ts2(1+lag:tL);          % horizontal transport
        ts3sh = ts3(1+lag:tL);          % ITF
        ts4sh = ts4(1+lag:tL);          % surface volume

        ts5sh = ts5(1+lag:tL);          % surface forcing
        ts6sh = ts6(1+lag:tL);          % vertical mixing
        ts7sh = ts7(1+lag:tL);          % skew diffusion
        ts8sh = ts8(1+lag:tL);          % neutral diffusion
        ts9sh = ts9(1+lag:tL);          % numerical diffusion

        ts10sh = ts10(1+lag:tL);        % dWWVdt regression as well
    
    else
        ts1sh = ts1(1-lag:tL);
        ts2sh = ts2(1:tL+lag);          % horiztonal transport
        ts3sh = ts3(1:tL+lag);          % ITF
        ts4sh = ts4(1:tL+lag);          % surface volume
        
        ts5sh = ts5(1:tL+lag);          % surface forcing
        ts6sh = ts6(1:tL+lag);          % vertical mixing
        ts7sh = ts7(1:tL+lag);          % skew diffusion
        ts8sh = ts8(1:tL+lag);          % neutral diffusion
        ts9sh = ts9(1:tL+lag);          % numerical diffusion
        
        ts10sh = ts10(1:tL+lag);        % dWWVdt regression as well 
    end        
    tLlag = length(ts1sh);
%     ts2sh_normalized = ts2sh;
%     ts3sh_normalized = ts3sh;
%     ts4sh_normalized = ts4sh;
%     
%     ts5sh_normalized = ts5sh;
%     ts6sh_normalized = ts6sh;
%     ts7sh_normalized = ts7sh;
%     ts8sh_normalized = ts8sh;
%     ts9sh_normalized = ts9sh;
%     
    dVr(ii)           = sum(ts1sh.*ts10sh)/tLlag;
    horizontalr(ii)   = sum(ts1sh.*ts2sh)/tLlag;
    ITFr(ii)          = sum(ts1sh.*ts3sh)/tLlag;
    massr(ii)         = sum(ts1sh.*ts4sh)/tLlag;
    
    forcingr(ii)   = sum(ts1sh.*ts5sh)/tLlag;
    mixingr(ii)    = sum(ts1sh.*ts6sh)/tLlag;
    skewr(ii)      = sum(ts1sh.*ts7sh)/tLlag;
    neutralr(ii)   = sum(ts1sh.*ts8sh)/tLlag;
    numericalr(ii) = sum(ts1sh.*ts9sh)/tLlag;
end


%% ~~~~~~~~ plotting  routine for lag regression with N34 ~~~~~~~~ %% 
directory = '/home/z5180028/MSC_thesis/access_figures/';
% load([directory 'WMT_time_series_1979-1991_rev3.mat']);

boom; antarctica = nature;
figure('units', 'pixels', 'position', [0 0 1600 800]);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
subplot(1,2,1)
h1 = plot(lags, horizontalr, 'color', RdYlBu(10,:), 'linewidth', 2); hold on;
h2 = plot(lags, ITFr, 'color', RdYlGn(60,:), 'linewidth', 2); hold on;
h3 = plot(lags, massr, 'color', RdYlBu(21,:), 'linewidth', 2); hold on;
% h11 = plot(lags, (horizontalr + ITFr + massr), 'color', 'k', 'linewidth', 2); hold on;
h11 = plot(lags, dVr, 'color', 'k', 'linewidth', 2.5); hold on;
line([0 0],[-100, 100], 'color', [.83 .83 .83]) % horizontal and vertical lines

hXLabel = xlabel('Lag [months]');
hYLabel = ylabel('[Sv $^{\circ}$C$^{-1}$]','interpreter','latex', 'color', RdYlBu(60,:));
set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 25, 'color', RdYlBu(60,:));
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -50:1:30   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

  
ylim([-5 2]);
pbaspect([2 1 1]);                        % aspect ratios: x, y, z

clear hTitle hXLabel hYLabel;


text(-27, 3, 'a) N34 and adiabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);

h55 = legend([h11 h1 h2 h3], ...
    'Change in Warm Water Volume', ...
    '$\mathcal{T}_{5^{\circ}\mathrm{N}+5^{\circ}\mathrm{S}}$: Meridional transport', ...
    '$\mathcal{T}_{ITF}$: Indonesian Throughflow', ...
    '$\mathcal{J}$: Surface volume', ...
    'location', 'southwest', 'orientation', 'vertical');
set(h55, 'interpreter', 'latex', 'fontsize', 11,'EdgeColor', [.83 .83 .83]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
subplot(1,2,2)

h5 = plot(lags, forcingr, 'color', RdYlBu(60,:), 'linewidth', 2); hold on;
h6 = plot(lags, mixingr, 'color', antarctica(15,:), 'linewidth', 2); hold on;
h7 = plot(lags, skewr, 'color', antarctica(5,:), 'linewidth', 2); hold on;
h8 = plot(lags, neutralr, 'color', RdYlGn(10,:), 'linewidth', 2); hold on;
h9 = plot(lags, numericalr, 'color', [187,0,187]/255, 'linewidth', 2); hold on;
% h10 = plot(lags, (forcingr + mixingr + skewr + neutralr + numericalr), ...
%     'color', 'k', 'linewidth', 2); hold on;
line([0 0],[-100, 100], 'color', [.83 .83 .83]) % horizontal and vertical lines

hXLabel = xlabel('Lag [months]');
hYLabel = ylabel('[Sv $^{\circ}$C$^{-1}$]','interpreter','latex', 'color', RdYlBu(60,:));
set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 25, 'color', RdYlBu(60,:));
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -50:1:30   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

text(-27, 3, 'b) N34 and diabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);
  
 
ylim([-5 2]);
pbaspect([2 1 1]);                        % aspect ratios: x, y, z
clear hTitle hXLabel hYLabel h1 h2;

h55 = legend([h5 h6 h7 h8 h9], ...
    '$\mathcal{G_F}$: Surface forcing', ...
    '$\mathcal{G_M}$: Vertical mixing', ...
    '$\mathcal{G_E}$: Skew diffusion', ...
    '$\mathcal{G_N}$: Neutral diffusion', ...    
    '$\mathcal{G_I}$: Numerical mixing', ...  
    'location', 'southwest', 'orientation', 'vertical');
set(h55, 'interpreter', 'latex', 'fontsize', 11,'EdgeColor', [.83 .83 .83]);


% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = '/home/z5180028/MSC_thesis/access_figures/';
print('-dpng','-r300', [directory 'WMT_time_series_lag_regression_N34']);


%% consistency check for N34 regression
figure(2)
a = horizontalr + ITFr + massr + forcingr + mixingr + skewr + neutralr + numericalr;
plot(lags,a); hold on;
plot(lags,dVr, 'linewidth', 2.5)
legend('sum of all fluxes', 'dWWVdt', 'location', 'best')
title('WWV balance terms lag regression with N34')

% they are very similar but not 100% the same...hmm, strange
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = '/home/z5180028/MSC_thesis/access_figures/';
print('-dpng','-r300', [directory 'WMT_time_series_lag_regression_check']);


%% calculate lag regression on dWWVdt ~~~~~~~~
t = 1:456;
tL = length(t);
ts1 = dVa;
% ts2 = (horizontal - mean(horizontal)) / std(horizontal);
% ts3 = (ITF - mean(ITF)) / std(ITF);
% ts4 = (mass - mean(mass)) / std(mass);
% 
% ts5 = (forcing  - mean(forcing)) / std(forcing);
% ts6 = (mixing  - mean(mixing)) / std(mixing);
% ts7 = (skew  - mean(skew)) / std(skew);
% ts8 = (neutral  - mean(neutral)) / std(neutral);
% ts9 = (numerical  - mean(numerical)) / std(numerical);
ts2 = horizontala;
ts3 = ITFa;
ts4 = massa;

ts5 = forcinga;
ts6 = mixinga;
ts7 = skewa;
ts8 = neutrala;
ts9 = numerical;
ts10 = nino;

lags = -20:20;
regs = zeros(size(lags));
for ii=1:length(lags)
    lag = lags(ii);
    if (lag >= 0)
        ts1sh = ts1(1:(tL-lag));
        ts2sh = ts2(1+lag:tL);          % horizontal transport
        ts3sh = ts3(1+lag:tL);          % ITF
        ts4sh = ts4(1+lag:tL);          % surface volume

        ts5sh = ts5(1+lag:tL);          % surface forcing
        ts6sh = ts6(1+lag:tL);          % vertical mixing
        ts7sh = ts7(1+lag:tL);          % skew diffusion
        ts8sh = ts8(1+lag:tL);          % neutral diffusion
        ts9sh = ts9(1+lag:tL);          % numerical diffusion
        ts10sh = ts10(1+lag:tL);        % N34 time series as well

    
    else
        ts1sh = ts1(1-lag:tL);
        ts2sh = ts2(1:tL+lag);          % horiztonal transport
        ts3sh = ts3(1:tL+lag);          % ITF
        ts4sh = ts4(1:tL+lag);          % surface volume
        
        ts5sh = ts5(1:tL+lag);          % surface forcing
        ts6sh = ts6(1:tL+lag);          % vertical mixing
        ts7sh = ts7(1:tL+lag);          % skew diffusion
        ts8sh = ts8(1:tL+lag);          % neutral diffusion
        ts9sh = ts9(1:tL+lag);          % numerical diffusion
        ts10sh = ts10(1:tL+lag);          % N34
        
    end        
    tLlag = length(ts1sh);
%     ts2sh_normalized = ts2sh;
%     ts3sh_normalized = ts3sh;
%     ts4sh_normalized = ts4sh;
%     
%     ts5sh_normalized = ts5sh;
%     ts6sh_normalized = ts6sh;
%     ts7sh_normalized = ts7sh;
%     ts8sh_normalized = ts8sh;
%     ts9sh_normalized = ts9sh;
%     
    horizontalr(ii)   = sum(ts1sh.*ts2sh)/tLlag;
    ITFr(ii)          = sum(ts1sh.*ts3sh)/tLlag;
    massr(ii)         = sum(ts1sh.*ts4sh)/tLlag;
    
    forcingr(ii)   = sum(ts1sh.*ts5sh)/tLlag;
    mixingr(ii)    = sum(ts1sh.*ts6sh)/tLlag;
    skewr(ii)      = sum(ts1sh.*ts7sh)/tLlag;
    neutralr(ii)   = sum(ts1sh.*ts8sh)/tLlag;
    numericalr(ii) = sum(ts1sh.*ts9sh)/tLlag;
    ninor(ii)      = sum(ts1sh.*ts10sh)/tLlag;
end


%% ~~~~~~~~ plotting  routine for lag regression with dWWVdt ~~~~~~~~ %% 
directory = '/home/z5180028/MSC_thesis/access_figures/';
% load([directory 'WMT_time_series_1979-1991_rev3.mat']);

boom; antarctica = nature;
figure('units', 'pixels', 'position', [0 0 1600 800]);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
subplot(1,2,1)
h1 = plot(lags, horizontalr, 'color', RdYlBu(10,:), 'linewidth', 2); hold on;
h2 = plot(lags, ITFr, 'color', RdYlGn(60,:), 'linewidth', 2); hold on;
h3 = plot(lags, massr, 'color', RdYlBu(21,:), 'linewidth', 2); hold on;
% h11 = plot(lags, ninor, 'color', 'k', 'linewidth', 2); hold on;
line([0 0],[-100, 200], 'color', [.83 .83 .83]) % horizontal and vertical lines

hXLabel = xlabel('Lag [months]');
hYLabel = ylabel('Sv Sv$^{-1}$','interpreter','latex', 'color', RdYlBu(60,:));
set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 25, 'color', RdYlBu(60,:));
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -500:50:300   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

  
ylim([-50 150]);
pbaspect([.8 1 1]);                        % aspect ratios: x, y, z

clear hTitle hXLabel hYLabel;


text(-27, 162, 'a) dWWVdt and adiabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);

h55 = legend([h1 h2 h3], ...
    '$\mathcal{T}_{5^{\circ}\mathrm{N}+5^{\circ}\mathrm{S}}$: Meridional transport', ...
    '$\mathcal{T}_{ITF}$: Indonesian Throughflow', ...
    '$\mathcal{J}$: Surface volume', ...
    'location', 'southeast', 'orientation', 'vertical');
set(h55, 'interpreter', 'latex', 'fontsize', 11,'EdgeColor', [.83 .83 .83]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
subplot(1,2,2)

h5 = plot(lags, forcingr, 'color', RdYlBu(60,:), 'linewidth', 2); hold on;
h6 = plot(lags, mixingr, 'color', antarctica(15,:), 'linewidth', 2); hold on;
h7 = plot(lags, skewr, 'color', antarctica(5,:), 'linewidth', 2); hold on;
h8 = plot(lags, neutralr, 'color', RdYlGn(10,:), 'linewidth', 2); hold on;
h9 = plot(lags, numericalr, 'color', [187,0,187]/255, 'linewidth', 2); hold on;
h10 = plot(lags, (forcingr + mixingr + skewr + neutralr + numericalr), ...
    'color', 'k', 'linewidth', 2); hold on;
line([0 0],[-100, 200], 'color', [.83 .83 .83]) % horizontal and vertical lines

hXLabel = xlabel('Lag [months]');
hYLabel = ylabel('Sv Sv$^{-1}$','interpreter','latex', 'color', RdYlBu(60,:));
set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 25, 'color', RdYlBu(60,:));
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -500:50:300   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

text(-27, 162, 'b) dWWVdt and diabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);
  
 
ylim([-50 150]);
pbaspect([.8 1 1]);                        % aspect ratios: x, y, z
clear hTitle hXLabel hYLabel h1 h2;

h55 = legend([h10 h5 h6 h7 h8 h9 h10], ...
    'Total diabatic', ...
    '$\mathcal{G_F}$: Surface forcing', ...
    '$\mathcal{G_M}$: Vertical mixing', ...
    '$\mathcal{G_E}$: Skew diffusion', ...
    '$\mathcal{G_N}$: Neutral diffusion', ...    
    '$\mathcal{G_I}$: Numerical mixing', ...  
    'location', 'northeast', 'orientation', 'vertical');
set(h55, 'interpreter', 'latex', 'fontsize', 11,'EdgeColor', [.83 .83 .83]);


% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = '/home/z5180028/MSC_thesis/access_figures/';
print('-dpng','-r300', [directory 'WMT_time_series_lag_regression_dWWVdt']);





