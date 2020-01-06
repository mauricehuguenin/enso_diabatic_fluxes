%% EXP1: Plotting latitude integrated OHCa during December of the first year

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  %
%                                                                         %
%                     hmaurice, 23.11.2017, 10:19 AEST                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  %


% load workspace
f3 = 'EXP1_and_EXP2_timeseries_wwv_transport_5S_5N_clim';
srv = 'H:\Maurice_ENSO_Data\'
load('E:\2017 Fr¸hlingssemester\Master Seminar 1\workspace_EXP1_analysis_composite_ninos_rev3.mat');
antarctica = nature;
RdBu_short  = cbrewer('div', 'RdBu', 21, 'PCHIP');
RdYlGn  = cbrewer('div', 'RdYlGn', 60, 'PCHIP');

% prepare mask                  % importante!       % wwv region 20S - 20N
                                % wwv_mask(117:840, 417:580) = 20S - 20N
                                % wwv_mask(117:840, 466:531) = 8S - 8N
                                % wwv_mask(117:840, 478:518) = 5S - 5N
Reds = cbrewer('seq', 'Reds', 16, 'PCHIP');
Blues = cbrewer('seq', 'Blues', 16, 'PCHIP');


wwv_mask = wwv_mask(91:811,478:518);
wwv_mask_2S = wwv_mask_2S(91:811,478:518);
wwv_mask_2S_borneo = wwv_mask_2S_borneo(91:811,478:518);
lon = lon(91:811,478:518);
lat = lat(91:811,478:518);
testmap(lon, lat, wwv_mask_2S_borneo);

antarctica =    [150,   0,  17; 165,   0,  33; 200,   0,  40; ...
                 216,  21,  47; 247,  39,  53; 255,  61,  61; ...
                 255, 120,  86; 255, 172, 117; 255, 214, 153; ...
                 255, 241, 188; 255, 255, 255; 188, 249, 255; ...
                 153, 234, 255; 117, 211, 255;  86, 176, 255; ...
                  61, 135, 255;  40,  87, 255;  24,  28, 247; ...
                  30,   0, 230;  36,   0, 216;  45,   0, 200]*1./255;

%% preamble and data preparation for warm water volume

                                % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
                                % location_tahiti        = [523, 427] %
                                % location_center_nino34 = [540, 498] %
                                % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
                           
                                
            % importante!       % wwv region 20S - 20N
                                % wwv_mask(117:840, 417:580) = 20S - 20N
                                % wwv_mask(117:840, 466:531) = 8S - 8N
                                
                                
first = 'EXP1';   second = 'EXP2';  restart = 'restart000';   
string = [first '_and_' second '_' restart];

tic;
% load in potential temperature
% p1 = ['/srv/ccrc/data15/z5180028/MSC_thesis_mom_output/' first '_and_' second '_' restart '_windstress/output000/'];
pc = 'H:\Maurice_ENSO_Data\EXP0_control_run/';


neutral_rho = getnc([pc 'ocean_wmass_clim.nc'], 'neutral'); % temperature space array


% calculating WWV from the restart file as a mean of 5 climatological years

% loading in climatology (i.e. the five restart files and creating the
% mean)
temp_resc1 = squeeze(permute(getnc([pc 'restart000/ocean_temp_salt.res.nc'], 'temp', ...
         [1,1,478,91], [-1,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
temp_resc2 = squeeze(permute(getnc([pc 'restart001/ocean_temp_salt.res.nc'], 'temp', ...
         [1,1,478,91], [-1,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
temp_resc3 = squeeze(permute(getnc([pc 'restart002/ocean_temp_salt.res.nc'], 'temp', ...
         [1,1,478,91], [-1,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
temp_resc4 = squeeze(permute(getnc([pc 'restart003/ocean_temp_salt.res.nc'], 'temp', ...
         [1,1,478,91], [-1,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
temp_resc5 = squeeze(permute(getnc([pc 'restart004/ocean_temp_salt.res.nc'], 'temp', ...
         [1,1,478,91], [-1,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
temp_resc6 = squeeze(permute(getnc([pc 'restart004/ocean_temp_salt.res.nc'], 'temp', ...
         [1,1,478,91], [-1,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
temp_resc = squeeze(nanmean(cat(4, temp_resc1, temp_resc2, temp_resc3, temp_resc4, temp_resc5), 4)); 
clear temp_resc1 temp_resc2 temp_resc3 temp_resc4 temp_resc5 temp_resc6;
     


% now finding temperature within it which is higher than 20 degrees
maskcr = temp_resc;
indices = find(maskcr < 20.5);  maskcr(indices) = 0; clear indices;
indices = find(maskcr >= 20.5); maskcr(indices) = 1; clear indices;



% loading in climatology
dztc = squeeze(permute(getnc([pc 'ocean_snap_clim.nc'], 'dzt', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 

% mask is now my array containing all grid cells 8S-8N and 156E-95W which
% have temperatures higher than 20 degrees
for z = 1:50
        volumec_res(:,:,z) = maskcr(:,:,z) .* areacello(91:811,478:518) .* dztc(:,:,z) .* wwv_mask_2S_borneo;
end


dV(1) = squeeze(nansum(squeeze(nansum(squeeze(nansum(volumec_res, 1)), 1)), 2))


%% calculating the 12 WWV values from the snapshot files
                               
% loading in climatology
thetaoc = squeeze(permute(getnc([pc 'ocean_snap_clim.nc'], 'temp', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 


% now finding temperature within it which is higher than 20 degrees
maskc = thetaoc;
indices = find(maskc < 20);  maskc(indices) = 0; clear indices;
indices = find(maskc >= 20); maskc(indices) = 1; clear indices;


% loading in climatology
dztc = squeeze(permute(getnc([pc 'ocean_snap_clim.nc'], 'dzt', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 

% mask is now my array containing all grid cells 5S-5N and 100E-60W which
% have temperatures higher than 20 degrees
volumec = nan(721, 41, 50, 12);
for z = 1:50
    for t = 1:12
        volumec(:,:,z, t) = maskc(:,:,z, t) .* areacello(91:811,478:518) .* dztc(:,:,z, t) .* wwv_mask_2S_borneo;
    end
end

size(volumec)
dV(2:13) = nansum(reshape(volumec, [721*41*50 12]));

plot(dV); % this is now my WWV before taking the derivative, i.e. 
% it has 13 time steps which get reduced to 12 once I take the derivative

delta_t = 30*24*60*60;
dV = diff(dV) / delta_t / 1e6 % nochmals durch 1 Million um sch√∂n in Sverdrup zu plotten

plot(dV)
nanmean(dV)

dV = cat(2, dV, dV, dV, dV, dV, dV)


%% calculating the horizontal transport and vertical transport as residual


% loading in climatology
ty_transc = squeeze(nansum(permute(getnc([pc 'ocean_clim.nc'], 'ty_trans_nrho', ...
         [1,46,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)); 
ty_transc = cat(3, ty_transc, ty_transc, ty_transc, ty_transc, ty_transc, ty_transc);

% calculating the transport across the different transects
for i = 1:721           % transect to the north, +5 degrees
    northc(i,:) = squeeze(ty_transc(i,41,:)) .* wwv_mask_2S_borneo(i,41);
end
for i = 50:127           % transect to the south, -2 degrees
    ITFc(i,:) = squeeze(-ty_transc(i,12,:)) .* wwv_mask_2S_borneo(i,12);
end
for i =  1:721           % transect to the south, -5 degrees
    % in my ACCESS-OM2 run I calculate over the grid cells 128:721 as
    % the first part of the 5∞S transect is empty (i.e. the bottom
    % left corner of the wwv_mask_2S_borneo
    southc(i,:) = squeeze(-ty_transc(i,1,:)) .* wwv_mask_2S_borneo(i,1);
end

northc = nansum(northc);
southc = nansum(southc); 
ITFc = nansum(ITFc);

plot(northc); hold on; plot(southc); hold on; plot(ITFc); hold on;

% multiplying by (-1) in order to achieve

%           ^
%           |           transport in
%           |
% ------------------------------------------ zero line
%           |
%           |           transport out
%           v   

northc = -(northc); 
southc = -(southc);
ITFc = -(ITFc);
%
horizontal      = northc + southc;      % total meridional transport
vertical        = dV - (horizontal + ITFc);  % total vertical tranport

% yes, all calculations now consistent with the ones I did in my
% ACCESS-OM2 run

plot(dV, 'linewidth', 2.5); hold on; line([1 72], [0 0], 'color', 'k');
plot(horizontal); hold on; plot(ITFc, 'color', 'g'); hold on;
plot(vertical); hold on;

% safespot: 08. 03. 2018, 16:16


%% creating budget with water mass transformation framework 

% defining constants
rho_0 = 1035;       % reference air density [kg m^3]
cp = 3992.1;         % specific heat capacity of seawater [J kg^-1 K^-1]
dT = 0.5;           % the temperature difference of the binned temperature 


                % ~~~~~~~~~~~~~~~~~ GM ~~~~~~~~~~~~~~~~~ %                
% climatology
cbtc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'temp_vdiffuse_diff_cbt_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
cbtc = cat(3, cbtc, cbtc, cbtc, cbtc, cbtc, cbtc);

% climatology
nonc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'temp_nonlocal_KPP_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
nonc = cat(3, nonc, nonc, nonc, nonc, nonc, nonc);

% units are in [W m^2]

         % ~~~~~~~~~~~ the calculation for GM is here ~~~~~~~~~~~ %        
for t = 1:72
    GMc(:,:,t) = (cbtc(:,:,t) + nonc(:,:,t)) .* wwv_mask_2S_borneo;
% units are in = W m^-2
end
for t = 1:72
    GMc(:,:,t) = GMc(:,:,t) .* (1 ./ (rho_0 .* cp .* dT)) .* areacello(91:811, 478:518) ./ 1e6;
                          % durch 1e6 teilen f√ºr Einheit Sv [10^6 m^3 s^-1]
end
mixingc = squeeze(nansum(squeeze(nansum(GMc,1)),1));

%           ^
%           |           transport in
%           |
% ------------------------------------------ zero line
%           |
%           |           transport out
%           v   

hold on; plot(vertical); hold on;
plot(mixingc); hold on; legend('vertical mixing'); 


                % ~~~~~~~~~~~~~~~~~ GF ~~~~~~~~~~~~~~~~~ %                

% climatology
sbcc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'temp_vdiffuse_sbc_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1]));
sbcc = cat(3, sbcc, sbcc, sbcc, sbcc, sbcc, sbcc);

% climatology
swc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'sw_heat_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
swc = cat(3, swc, swc, swc, swc, swc, swc);

% climatology
frazilc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'frazil_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
frazilc = cat(3, frazilc, frazilc, frazilc, frazilc, frazilc, frazilc);


% climatology
etac = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'temp_eta_smooth_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
etac = cat(3, etac, etac, etac, etac, etac, etac);

% units are in [W m^2]

         % ~~~~~~~~~~~ the calculation for GF is here ~~~~~~~~~~~ %        
for t = 1:72
    GFc(:,:,t) = (sbcc(:,:,t) + swc(:,:,t) + frazilc(:,:,t) + etac(:,:,t)) .* wwv_mask_2S_borneo;
% units are in = W m^-2
end

for t = 1:72
    GFc(:,:,t) = GFc(:,:,t) .* (1 ./ (rho_0 .* cp .* dT)) .* areacello(91:811, 478:518) ./ 1e6;
                          % durch 1e6 teilen f√ºr Einheit Sv [10^6 m^3 s^-1]
end
mixingc = squeeze(nansum(squeeze(nansum(GMc,1)),1));

forcingc = squeeze(nansum(squeeze(nansum(GFc,1)),1));

%           ^
%           |           transport in
%           |
% ------------------------------------------ zero line
%           |
%           |           transport out
%           v   

hold on; plot(vertical); hold on;
plot(forcingc); hold on; legend('surface forcing'); 
plot(mixingc); hold on; legend('vertical mixing'); 


                % ~~~~~~~~~~~~~~~~~ GJ ~~~~~~~~~~~~~~~~~ %                
% climatology
massc = squeeze(nansum(permute(getnc([pc 'ocean_clim.nc'], 'mass_pmepr_on_nrho', ...
         [1,46,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]),3)); 
massc = cat(3, massc, massc, massc, massc, massc, massc);


         % ~~~~~~~~~~~ the calculation for GJ is here ~~~~~~~~~~~ %        
for t = 1:72
    GJc(:,:,t) = (massc(:,:,t)) .* wwv_mask_2S_borneo .* 1e-3 ./ 1e6; % umrechnen von [kg s‚?ª1] zu [m^3 s^-1]
                                 % ./ 1e6 Umrechnung auf Sv [10^6 m^3 s^-1]
% units are in = [m^3 s^-1]
end

massc = squeeze(nansum(squeeze(nansum(GJc,1)),1));

%           ^
%           |           transport in
%           |
% ------------------------------------------ zero line
%           |
%           |           transport out
%           v   



implicit = vertical - (mixingc + forcingc + massc);

plot(implicit, 'linewidth', 5);

% now since I have all the fluxes calculated I create the 5-month moving
% means of the vertical ones too

plot(vertical); hold on;
plot(forcingc); hold on; legend('surface forcing'); 
plot(mixingc); hold on; legend('vertical mixing'); 
plot(massc); hold on; legend('volume flux'); 
plot(implicit, 'linewidth', 5);


%% creating moving average for data
dV = movmean(dV(1:12),3);
horizontal = movmean(horizontal(1:12), 3);
ITFc = movmean(ITFc(1:12), 3);
massc = movmean(massc(1:12),3);
  
forcingc = movmean(forcingc(1:12),3);
mixingc = movmean(mixingc(1:12),3);
implicit = dV - horizontal - ITFc - forcingc - mixingc;


%% clean up  workspace before plotting
clear ans areacello bathymetry cto ctoc dzt dztc h7 h8 h9 lat lev_bnds lon;
clear mask maskc  so soc_3 soc_4 thetao thetaoc volcello; 
clear volume volumec z h3 h4 h5 h1 h2 net_sfcc;
clear soc i sstc tempc;
clear p1 p2 p3 p4 p5 pc tx_trans ty_trans wwv_region;
clear a neutral_rho;
clear cbt cbtc cp delta_t dT eta etac forcing  frazil frazilc GF;
clear GFc GM GMc GJ GJc maskcr maskr mass mixing non nonc;
clear rho_0 sbc sbcc sw swc wwv_mask wwv_mask_2S ty_transc t;


%% ~~~~~~~~~~~ plotting  routine for meridional figure ~~~~~~~~~~~ %% 

clc
% delta_t = 30*24*60*60;
% nanmean(dV(1:12)) * delta_t * 1e6
% nanmean(horizontal(1:12))* delta_t * 1e6
% nanmean(ITFc(1:12))* delta_t * 1e6
% nanmean(massc(1:12))* delta_t * 1e6
% nanmean(forcingc(1:12))* delta_t * 1e6
% nanmean(mixingc(1:12))* delta_t * 1e6
% nanmean(implicit(1:12))* delta_t * 1e6


figure('units', 'pixels', 'position', [0 0 1080 1080]);
subplot(2,1,1)
% figure('units', 'pixels', 'position', [0 0 1500 960]);
% h6 = plot(1:12, massc, 'color', RdYlBu(21,:), 'linewidth', 2.5); hold on;
% h2 = plot(1:12, horizontal, 'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;
% h8 = plot(1:12, ITFc, 'color', RdYlGn(60,:), 'linewidth', 2.5); hold on;
h1 = plot(1:12, dV, 'color', [0 0 0], 'linewidth', 3); hold on;

% calculate annual trend
% trend_data = movmean(dV,3);
% trend_data = trend_data(1:12); a = (trend_data(12) - trend_data(1)) / 12;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~ create fancy plot ~~~~~~~~~~~~~~~~~~~~~~~~~~ %
hYLabel = ylabel(['Sv'],'interpreter','latex');
  
%     % draw arrows to denote transport in and out
% 	arrow([13.5 4], [13.5, 9], 'BaseAngle', 20, 'length', 2, 'width', .2, 'ends', [1], 'color', RdYlBu(60,:));
% 	text([12.2 12.2],[11.5 11.5],'volume increase', 'interpreter', 'latex', 'fontsize', 21, 'color', RdYlBu(60,:));
% 	arrow([13.5, -4], [13.5, -9], 'BaseAngle', 20, 'length', 2, 'width', .2, 'ends', [1], 'color', RdYlBu(60,:));
% 	text([12.2 12.2],[-11.5 -11.5],'volume decrease', 'interpreter', 'latex', 'fontsize', 21, 'color', RdYlBu(60,:));
 

set([hYLabel], 'interpreter', 'latex', 'Fontsize', 25);
set(hYLabel, 'color', [0 0 0]);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -60:10:60   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...
   


set(gca, 'Xtick', [1:2:12]);

ylim([-30 28]);
xlim([1 12]);
pbaspect([2 1 1]);                        % aspect ratios: x, y, z
% h55 = legend([h1 h2 h3], 'location', 'northeast', 'orientation', 'vertical', ...
%     'Change in WWV anomaly', 'Horizontal Transport', ...
%     'Vertical Transport as residual');
h55 = legend([h1 ], 'location', 'northwest', 'orientation', 'vertical', ...
    'Change in WWV anomaly', ...
    '$\mathcal{T}$(5$^{\circ}$N + 5$^{\circ}$S): Meridional transport', ...
    '$\mathcal{T}$(ITF): Indonesian Throughflow', ...
    '$\mathcal{J}$: Surface volume', ...
    '$\mathcal{G_F}$: Surface forcing (2.8)', ...
    '$\mathcal{G_M}$: Vertical mixing  (6.8)', ...
    '$\mathcal{G_I}$: Numerical mixing  (3.1)');
set(h55, 'interpreter', 'latex', 'fontsize', 18, 'textcolor', RdYlBu(60,:));

clear hTitle hXLabel hYLabel h1 h2;

text(-1, 31, 'a) Adiabatic processes', 'interpreter', 'latex', 'Fontsize', 25, ...
     'color', [0 0 0]);


subplot(2,1,2)
% figure('units', 'pixels', 'position', [0 0 1500 960]);

h4 = plot(1:12, mixingc, 'color', antarctica(15,:), 'linewidth', 2.5); hold on;
h7 = plot(1:12, implicit, 'color', [187,0,187]/255, 'linewidth', 2.5); hold on;
h5 = plot(1:12, forcingc, 'color', RdYlBu(60,:), 'linewidth', 2.5); hold on;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~ create fancy plot ~~~~~~~~~~~~~~~~~~~~~~~~~~ %
hXLabel = xlabel('Month');
hYLabel = ylabel(['Sv'],'interpreter','latex');
  
    % draw arrows to denote transport in and out
% 	arrow([13.5 4], [13.5, 9], 'BaseAngle', 20, 'length', 2, 'width', .2, 'ends', [1], 'color', RdYlBu(60,:));
% 	text([12.2 12.2],[11.5 11.5],'volume increase', 'interpreter', 'latex', 'fontsize', 21, 'color', RdYlBu(60,:));
% 	arrow([13.5, -4], [13.5, -9], 'BaseAngle', 20, 'length', 2, 'width', .2, 'ends', [1], 'color', RdYlBu(60,:));
% 	text([12.2 12.2],[-11.5 -11.5],'volume decrease', 'interpreter', 'latex', 'fontsize', 21, 'color', RdYlBu(60,:));
 

set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 25);
set(hYLabel, 'color', [0 0 0]);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -60:10:60   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...
   


set(gca, 'Xtick', [1:2:12]);

ylim([-30 28]);
xlim([1 12]);
pbaspect([2 1 1]);                        % aspect ratios: x, y, z
% h55 = legend([h1 h2 h3], 'location', 'northeast', 'orientation', 'vertical', ...
%     'Change in WWV anomaly', 'Horizontal Transport', ...
%     'Vertical Transport as residual');
h55 = legend([h5 h4  h7], 'location', 'southwest', 'orientation', 'vertical', ...
    '$\mathcal{G_F}$: Surface forcing', ...
    '$\mathcal{G_M}$: Vertical mixing', ...
    '$\mathcal{G_I}$: Numerical mixing');
set(h55, 'interpreter', 'latex', 'fontsize', 18, 'textcolor', RdYlBu(60,:));

clear hTitle hXLabel hYLabel h1 h2;

text(-1, 31, 'b) Diabatic processes', 'interpreter', 'latex', 'Fontsize', 25, ...
     'color', [0 0 0]);

% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = 'C:/Users/Maurice Huguenin/Desktop/';
print('-dpng','-r300', [directory f3]);
 

%% ~~~~~~~~~~~ plotting  routine for complete figure with subplots ~~~~~~~~~~~ %% 

figure('units', 'pixels', 'position', [0 0 1920 1080]);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
subplot(2,2,1)

h2 = plot(linspace(1,12,12), horizontal, 'color', RdYlBu(10,:), 'linewidth', 2); hold on;
% h3 = plot(linspace(1,72,72), vertical, 'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;
h1 = plot(linspace(1,12,12), dV, 'color', [0 0 0], 'linewidth', 2.5); hold on;
% h11 = plot(linspace(1,72,72), movmean(test,5), 'color', [0 0 0], 'linewidth', 3, 'linestyle', '--'); hold on;
h6 = plot(linspace(1,12,12), massc, 'color', RdYlBu(21,:), 'linewidth', 2); hold on;

% h5 = plot(linspace(1,72,72), forcinga_EXP1, 'color', RdYlBu(60,:), 'linewidth', 2); hold on;
% h4 = plot(linspace(1,72,72), mixinga_EXP1, 'color', antarctica(15,:), 'linewidth', 2); hold on;
% h7 = plot(linspace(1,72,72), implicit_EXP1, 'color', [187,0,187]/255, 'linewidth', 2); hold on;
h1 = plot(linspace(1,12,12), dV, 'color', [0 0 0], 'linewidth', 2.5); hold on;
h8 = plot(linspace(1,12,12), ITFc, 'color', RdYlGn(60,:), 'linewidth', 2); hold on;

hXLabel = xlabel('Month');
hYLabel = ylabel(['Sv'],'interpreter','latex');
set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 25);
set(hYLabel, 'color', [0 0 0]);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -40:10:30   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);


set(gca, 'Xtick', [1, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60, 66, 72]);

h55 = legend([h1 h2 h8 h6], 'location', 'northwest', 'orientation', 'vertical', ...
    'Change in WWV anomaly', ...
    '$\mathcal{T}_{\mathrm{5^{\circ}N + 5^{\circ}S}}$: Meridional transport', ...
    '$\mathcal{T}_{\mathrm{ITF}}$: Indonesian Throughflow', ...
    '$\mathcal{J}$: Surface volume', ...
    '$\mathcal{G_F}$: Surface forcing', ...
    '$\mathcal{G_M}$: Vertical mixing', ...
    '$\mathcal{G_I}$: Numerical mixing');
set(h55, 'interpreter', 'latex', 'fontsize', 15.5, 'textcolor', RdYlBu(60,:));

ylim([-30 25]);
xlim([1 12]);
pbaspect([1.5 1 1]);                        % aspect ratios: x, y, z
clear hTitle hXLabel hYLabel h1 h2;

% draw arrows to denote El Nino and La Nina events
arrow([14, 5], [14, 20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(60,:));
arrow([14, -5], [14, -20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(1,:));
text([15 15],[5 5],'Recharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(60,:), 'rotation', 90);
text([15 15],[-5 -5],'Discharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(1,:), 'rotation', 90, 'horizontalalignment', 'right');

text(0, 31, 'a) Adiabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
subplot(2,2,2)


h5 = plot(linspace(1,12,12), forcingc, 'color', RdYlBu(60,:), 'linewidth', 2); hold on;
h4 = plot(linspace(1,12,12), mixingc, 'color', antarctica(15,:), 'linewidth', 2); hold on;
h7 = plot(linspace(1,12,12), implicit, 'color', [187,0,187]/255, 'linewidth', 2); hold on;
hXLabel = xlabel('Month');
hYLabel = ylabel(['Sv'],'interpreter','latex');    
set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 25);
set(hYLabel, 'color', [0 0 0]);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -30:10:25   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

set(gca, 'Xtick', [1, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60, 66, 72]);


ylim([-30 25]);
xlim([1 12]);
pbaspect([1.5 1 1]);                        % aspect ratios: x, y, z
clear hTitle hXLabel hYLabel h1 h2;

h55 = legend([h5 h4 h7], 'location', 'southwest', 'orientation', 'vertical', ...
    '$\mathcal{G_F}$: Surface forcing', ...
    '$\mathcal{G_M}$: Vertical mixing', ...
    '$\mathcal{G_I}$: Numerical mixing');
set(h55, 'interpreter', 'latex', 'fontsize', 15.5, 'textcolor', RdYlBu(60,:));

% draw arrows to denote El Nino and La Nina events
arrow([13, 5], [13, 20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(60,:));
arrow([13, -5], [13, -20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(1,:));
text([14 14],[5 5],'Recharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(60,:), 'rotation', 90);
text([14 14],[-5 -5],'Discharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(1,:), 'rotation', 90, 'horizontalalignment', 'right');


text(0, 31, 'b) Diabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = 'C:/Users/Maurice Huguenin/Desktop/';
print('-dpng','-r300', [directory f3]);

% 



