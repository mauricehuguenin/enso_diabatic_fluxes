%% Plotting time series of idealized El Nino and La Nina simulations,
% i.e. GMSST and OHC as well as N34 and WWV anomalous time series

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  %
%                     hmaurice, 23.11.2017, 10:19 AEST                    %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  %


f1 = 'pnEXP1_and_pnEXP2_timeseries_nino34_and_wwv_and_taux_and_nsfc';
f2 = 'pnEXP1_and_pnEXP2_timeseries_nino34_and_wwv_and_taux_and_nsfc_ensemble difference';
srv = 'E:\2017 Herbstsemester\Application Paper\AMS LaTeX Package v4.3.2\matlab_scripts_and_powerpoint_to_create_figures\'
% load workspace data
load([srv 'workspace_EXP1_and_EXP2_polynomial_PC_composites_symmetric.mat']);
load(['H:/Maurice_ENSO_Data/workspace_EXP1_analysis_composite_ninos_rev3'])
% colours
antarctica = nature; % antarctica colour bar
RdBu_short = cbrewer('div', 'RdBu', 21, 'PCHIP');
Reds = cbrewer('seq', 'Reds', 16, 'PCHIP');
Blues = cbrewer('seq', 'Blues', 16, 'PCHIP');

plot(linspace(1,24, 2920), pnnino_comp_highres(1,:)); hold on; 
plot(linspace(1,24, 2920), pnnino_comp_highres(2,:)); hold on;


%% preamble and data preparation


for schlaufe = 1:2 % loop over El Nino and La Nina events
tic 
    if schlaufe == 1
       first = 'pnEXP1_composite_nino_windstress'; 
       string = ['pnEXP1_restart000'];

    elseif schlaufe == 2
       first = 'pnEXP2_composite_nina_windstress'; 
       string = ['pnEXP2_restart000'];
    end

% pEXP1_nino97


% load in potential temperature
p1 = ['H:/Maurice_ENSO_Data/' first '/output000/'];
p2 = ['H:/Maurice_ENSO_Data/' first '/output001/'];
p3 = ['H:/Maurice_ENSO_Data/' first '/output002/'];
p4 = ['H:/Maurice_ENSO_Data/' first '/output003/'];
pc = 'H:/Maurice_ENSO_Data/EXP0_control_run/';

% --- import potential temperature
thetao_1 = squeeze(permute(getnc([p1 'ocean.nc'], 'temp', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
         % [month, depth, lat, lon], [end], [stride]
thetao_2 = squeeze(permute(getnc([p2 'ocean.nc'], 'temp', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
thetao_3 = squeeze(permute(getnc([p3 'ocean.nc'], 'temp', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
thetao_4 = squeeze(permute(getnc([p4 'ocean.nc'], 'temp', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
thetaoc = squeeze(permute(getnc([pc 'ocean_clim.nc'], 'temp', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
% concatenate across the time dimension and remove climatology to get anomalies     
thetao = cat(3, thetao_1, thetao_2, thetao_3, thetao_4) - ...
    cat(3, thetaoc, thetaoc, thetaoc, thetaoc);
clear thetao_1 thetao_2 thetao_3 thetao_4 thetao_5 thetao_6 thetaoc;
     
% --- import practical salinity
so_1 = squeeze(permute(getnc([p1 'ocean.nc'], 'salt', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
         % [month, depth, lat, lon], [end], [stride]
so_2 = squeeze(permute(getnc([p2 'ocean.nc'], 'salt', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
so_3 = squeeze(permute(getnc([p3 'ocean.nc'], 'salt', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
so_4 = squeeze(permute(getnc([p4 'ocean.nc'], 'salt', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
soc = squeeze(permute(getnc([pc 'ocean_clim.nc'], 'salt', ...
         [1,1,478,440], [12,1,518,640], [1,1,1,1]), [4 3 2 1])); 
so = cat(3, so_1, so_2, so_3, so_4) - ...
    cat(3, soc, soc, soc, soc);
clear so_1 so_2 so_3 so_4 so_5 so_6 soc;
    

cto = gsw_CT_from_pt(so,thetao); % convert potential to conservative temperature

% create nino3.4 index for comparison with WWV
if string == 'pnEXP1_restart000'
   mom_nino_pnEXP1(1,:) = squeeze(nanmean(squeeze(nanmean(cto, 1)),1));
elseif string == 'pnEXP2_restart000'
   mom_nino_pnEXP2(1,:) = squeeze(nanmean(squeeze(nanmean(cto, 1)),1));
end

%
clear thetao test so soc thetaoc;

% loading in warm water volume, i.e. volume of water above 20 degrees C in 
% the Pacific region 5S - 5N and 120E - 80W
thetao_1 = squeeze(permute(getnc([p1 'ocean.nc'], 'temp', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
         % [month, depth, lat, lon], [end], [stride]
thetao_2 = squeeze(permute(getnc([p2 'ocean.nc'], 'temp', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
thetao_3 = squeeze(permute(getnc([p3 'ocean.nc'], 'temp', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
         % [month, depth, lat, lon], [end], [stride]
thetao_4 = squeeze(permute(getnc([p4 'ocean.nc'], 'temp', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
thetao = cat(4, thetao_1, thetao_2, thetao_3, thetao_4); 
clear thetao_1 thetao_2 thetao_3 thetao_4 thetao_5 thetao_6;

so_1 = squeeze(permute(getnc([p1 'ocean.nc'], 'salt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
so_2 = squeeze(permute(getnc([p2 'ocean.nc'], 'salt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
so_3 = squeeze(permute(getnc([p3 'ocean.nc'], 'salt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
so_4 = squeeze(permute(getnc([p4 'ocean.nc'], 'salt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
so = cat(4, so_1, so_2, so_3, so_4); 
clear so_1 so_2 so_3 so_4 so_5 so_6;

% loading in climatology
thetaoc = squeeze(permute(getnc([pc 'ocean_clim.nc'], 'temp', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
         % [month, depth, lat, lon], [end], [stride]
thetaoc = cat(4, thetaoc, thetaoc, thetaoc, thetaoc);

soc = squeeze(permute(getnc([pc 'ocean_clim.nc'], 'salt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
         % [month, depth, lat, lon], [end], [stride]
soc = cat(4, soc, soc, soc, soc);

cto = gsw_CT_from_pt(so,thetao); 
ctoc = gsw_CT_from_pt(soc,thetaoc);

% now finding temperature within it which is higher than 20 degrees
mask = cto;
maskc = ctoc;
mask(mask < 20) = 0; 
maskc(maskc < 20) = 0; 
mask(mask >= 20) = 1; 
maskc(maskc >= 20) = 1; 


dzt_1 = squeeze(permute(getnc([p1 'ocean.nc'], 'dzt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
dzt_2 = squeeze(permute(getnc([p2 'ocean.nc'], 'dzt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
dzt_3 = squeeze(permute(getnc([p3 'ocean.nc'], 'dzt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
dzt_4 = squeeze(permute(getnc([p4 'ocean.nc'], 'dzt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
dzt = cat(4, dzt_1, dzt_2, dzt_3, dzt_4); 
clear dzt_1 dzt_2 dzt_3 dzt_4 dzt_5 dzt_6;

% loading in climatology
dztc = squeeze(permute(getnc([pc 'ocean_clim.nc'], 'dzt', ...
         [1,1,478,161], [12,-1,518,880], [1,1,1,1]), [4 3 2 1])); 
dztc = cat(4, dztc, dztc, dztc, dztc);

% mask is now my array containing all grid cells 5S-5N and 100E-60W which
% have temperatures higher than 20 degrees
for z = 1:50
    for t = 1:48
        volume(:,:,z,t)  = mask(:,:,z,t) .* areacello(161:880,478:518) .* dzt(:,:,z,t) .* wwv_mask_2S_borneo(161:880,478:518);
        volumec(:,:,z,t) = maskc(:,:,z,t) .* areacello(161:880,478:518) .* dzt(:,:,z,t) .* wwv_mask_2S_borneo(161:880,478:518);
    end
end
s = nansum(reshape(volume, [720*41*50 48]));
t = nansum(reshape(volumec, [720*41*50 48]));

% create anomalous WWV and save as variable
if string == 'pnEXP1_restart000'
    wwv_pnEXP1(1,:) = s-t;
elseif string == 'pnEXP2_restart000'
    wwv_pnEXP2(1,:) = s-t;
end

% load in zonal wind stress anomalies over 5S - 5N and 160E - 160W
tau_x1 = permute(getnc([p1 'ocean.nc'], 'tau_x', ...
         [1,478,320], [12,518,480], [1,1,1]), [3 2 1]); 
         % [month, lat, lon], [end], [stride]
tau_x2 = permute(getnc([p2 'ocean.nc'], 'tau_x', ...
         [1,478,320], [12,518,480], [1,1,1]), [3 2 1]); 
tau_x3 = permute(getnc([p3 'ocean.nc'], 'tau_x', ...
         [1,478,320], [12,518,480], [1,1,1]), [3 2 1]); 
         % [month, lat, lon], [end], [stride]
tau_x4 = permute(getnc([p4 'ocean.nc'], 'tau_x', ...
         [1,478,320], [12,518,480], [1,1,1]), [3 2 1]); 
tau_x = cat(3, tau_x1, tau_x2, tau_x3, tau_x4); 
clear tau_x1 tau_x2 tau_x3 tau_x4 tau_x5 tau_x6;
     
tauc_x = permute(getnc([pc 'ocean_clim.nc'], 'tau_x', ...
         [1,478,320], [12,518,480], [1,1,1]), [3 2 1]); 
tauc_x = cat(3, tauc_x, tauc_x, tauc_x, tauc_x);

% creating wind stress anomalies relative to climatology
tau_x = tau_x - tauc_x; clear tauc_x;          

% testmap(lon(320:480, 478:518), lat(320:480, 478:518), tau_x(:,:,12)); 
% great, the correct region

s = squeeze(nanmean(squeeze(nanmean(tau_x,1)),1));

% naming convention
if string == 'pnEXP1_restart000'
    taux_pnEXP1(1,:) = s;
elseif string == 'pnEXP1_restart000'
    taux_pnEXP2(1,:) = s;
end

net_sfc_1 = getnc([p1 'ocean_scalar.nc'], 'total_net_sfc_heating');
net_sfc_2 = getnc([p2 'ocean_scalar.nc'], 'total_net_sfc_heating');
net_sfc_3 = getnc([p3 'ocean_scalar.nc'], 'total_net_sfc_heating');
net_sfc_4 = getnc([p4 'ocean_scalar.nc'], 'total_net_sfc_heating');
net_sfc_global = cat(1, net_sfc_1, net_sfc_2, net_sfc_3, net_sfc_4);
clear net_sfc_1 net_sfc_2 net_sfc_3 net_sfc_4 net_sfc_5 net_sfc_6;
% units are in [10^15 W m^-2]

net_sfcc = getnc([pc 'ocean_scalar_clim.nc'], 'total_net_sfc_heating');
net_sfcc = cat(1, net_sfcc, net_sfcc, net_sfcc, net_sfcc);
% units are in [10^15 W m^-2]

% creating global anomaly
net_sfc_global = net_sfc_global - net_sfcc;
if string == 'pnEXP1_restart000'
    net_sfc_global_pnEXP1(1,:) = net_sfc_global;
elseif string == 'pnEXP2_restart000'
    net_sfc_global_pnEXP2(1,:) = net_sfc_global;
end

% GMSST with area-weighting
area = permute(getnc([p1 'ocean_grid.nc'], 'area_t'), [2 1]);
area_nan = area;
    
% calculate sst
thetao_1 = squeeze(permute(getnc([p1 'ocean.nc'], 'temp', ...
         [1,1,1,1], [12,1,-1,-1], [1,1,1,1]), [4 3 2 1])); 
thetao_2 = squeeze(permute(getnc([p2 'ocean.nc'], 'temp', ...
         [1,1,1,1], [12,1,-1,-1], [1,1,1,1]), [4 3 2 1])); 
thetao_3 = squeeze(permute(getnc([p3 'ocean.nc'], 'temp', ...
         [1,1,1,1], [12,1,-1,-1], [1,1,1,1]), [4 3 2 1])); 
thetao_4 = squeeze(permute(getnc([p4 'ocean.nc'], 'temp', ...
         [1,1,1,1], [12,1,-1,-1], [1,1,1,1]), [4 3 2 1])); 
thetao = cat(3, thetao_1, thetao_2, thetao_3, thetao_4);
clear thetao_1 thetao_2 thetao_3 thetao_4 thetao_5 thetao_6;

% loading in climatology
thetaoc = squeeze(permute(getnc([pc 'ocean_clim.nc'], 'temp', ...
         [1,1,1,1], [12,1,-1,-1], [1,1,1,1]), [4 3 2 1])); 
thetaoc = cat(3, thetaoc, thetaoc, thetaoc, thetaoc);
   
    
    
for t = 1:48
%     s(t) = squeeze(nansum(nansum(thetao(:,:,t)) .* area,1),2)) ./ nansum(nansum(area_nan(isnan(thetao(:,:,t))),1),2);
    global_temp(t) = squeeze(nansum(nansum(thetao(:,:,t).*area,1),2)./nansum(nansum(area_nan,1),2));
    global_tempc(t) = squeeze(nansum(nansum(thetaoc(:,:,t).*area,1),2)./nansum(nansum(area_nan,1),2));
end 
    sst = global_temp - global_tempc; clear global_tempc;

if string == 'pnEXP1_restart000'
    sst_pnEXP1(1,:) = sst;
elseif string == 'pnEXP2_restart000'
    sst_pnEXP2(1,:) = sst;
end

% clean up  workspace before plotting

clear ans bathymetry cto ctoc dzt dztc h7 h8 h9 lev_bnds;
clear mask maskc s so soc_3 soc_4 t thetao thetaoc volcello; 
clear volume volumec z h3 h4 h5 h1 h2 net_sfcc;
clear soc i sstc tempc p1 p2 p3 p4 pc antarctica first second restart tau_x
clear wwv_mask sstc;

toc;
end    % end of for loup

clear area area_nan  basin_mask global_temp lat lon; 
clear net_sfc_global schlaufe sst wwv_mask_2S wwv_mask_2S_borneo;


%% plotting time series
% initialize figure
figure('units', 'pixels', 'position', [0 0 1920 1080]);

% I have four subplots here, instead of a loop through 1:4 I instead opted
% to copy-paste the code four times (not really good coding but it did the
% job back then)

% sst anomaly
subplot(2,2,1)
yyaxis left % two time series in one plot with axis to the left and right   
% --- plot time series: GMSST
h5 = plot(linspace(1,48,48), sst_pnEXP1, 'color', RdYlBu(60,:), 'linewidth', 2.5); hold on;
ylim([-0.15 0.15]);
hYLabel = ylabel(['[$^{\circ}$C]'],'interpreter','latex');
set([hYLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(60,:));

% draw arrows to denote El Nino and La Nina events
arrow([1.1, 0.165], [4, 0.165], 'BaseAngle', 70, 'length', 5, 'width', .6, 'ends', [2], 'color', RdYlBu(1,:));
arrow([21, 0.165], [23.9, 0.165], 'BaseAngle', 70, 'length', 5, 'width', .6, 'ends', [1], 'color', RdYlBu(1,:));
text([4 4],[0.165 0.165],'El Ni\~no event', 'interpreter', 'latex', 'fontsize', 21, 'color', RdYlBu(1,:));

arrow([24.1, 0.165], [29, 0.165], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [2], 'color', 'k');
arrow([43, 0.165], [47.9, 0.165], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [1], 'color', 'k');
text([30 30],[0.165 0.165],'spin-down', 'interpreter', 'latex', 'fontsize', 22, 'color', 'k');

set(gca, ...
    'Box'             , 'off'         , ...
    'TickDir'         , 'out'         , ...
    'TickLength'      , [.01 .01]     , ...
    'XMinorTick'      , 'off'         , ...
    'YMinorTick'      , 'off'         , ...
    'YGrid'           , 'on'          , ...
    'YTick'           , -0.15:0.05:0.15  , ...     % define grid, every 1000 m a grid line
    'XColor'          , RdYlBu(60,:)  , ...
    'YColor'          , RdYlBu(60,:)  , ...
    'ticklabelinterpreter', 'latex'   , ...
    'LineWidth'       , 1.25);

yyaxis right % two time series in one plot with axis to the left and right
% --- plot time series: OHC as the cumulative integrated net surface heat
% flux
h3 = plot(linspace(1,48,48), cumtrapz(net_sfc_global_pnEXP1)* ...
    60*60*24*365.25*4*squeeze(nansum(squeeze(nansum(areacello, 1)), 2)), 'color', RdYlBu(20,:), 'linewidth', 2.5); hold on;
set(gca, 'xtick',[]);
ylim([-3e23 3e23]);
    hYLabel = ylabel(['[J]'],'interpreter','latex');

    vline = line([24 24], [-5e23 5e23]); % draw boundary line between El Niño and La Niña simulation
    set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);
    vline = line([48 48], [-5e23 5e23]); % draw boundary line between El Niño and La Niña simulation
    set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);
    
    xlim([1 48]);

    set(gca, ...
       'XColor'          , RdYlBu(60,:)  , ...
       'YColor'          , RdYlBu(60,:));
    set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 21);
    set([hYLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(20,:));

text(-7, 4e23, 'a) El Ni\~no event', 'interpreter', 'latex', 'Fontsize', 30, ...
     'color', [0 0 0]);

% add legend
h55 = legend([h5 h3 ], 'GMSST anomaly', 'OHC anonmaly', ...
    'OHC anomaly', 'OHC$_{\mathrm{non-mirrored}}$', ...
    'location', 'southeast', 'orientation', 'vertical');
set(h55, 'interpreter', 'latex', 'fontsize', 18, 'textcolor', RdYlBu(60,:));

% custom labels
hXLabel = xlabel('Month', 'interpreter', 'latex', 'color', [1 1 1]);
set(gca, 'Xtick', [1, 6, 12, 18, 24, 30, 36, 42, 48]);

subplot(2,2,3)
yyaxis left   
% --- time series, idealized N34 index and the one used for deriving the
% atmospheric input for the simulation
h2 = plot(linspace(1,48,48), mom_nino_pnEXP1, 'color', [RdYlBu(50,:)], 'linewidth', 2.5); hold on;
h22 = plot(linspace(1,24, 2920), pnnino_comp_highres(1,:), 'color', RdYlBu(50,:), 'linewidth', 2.5, 'linestyle', '--'); hold on;
h23 = plot(linspace(25,48, 2920), pnnino_comp_highres_spin(1,:), 'color', RdYlBu(50,:), 'linewidth', 2.5, 'linestyle', '--'); hold on;
ylim([-2.5 2.5]);
    hYLabel = ylabel(['[$^{\circ}$C]'],'interpreter','latex');

% draw arrows to denote El Nino and La Nina events
arrow([1.1, 2.6], [4, 2.6], 'BaseAngle', 70, 'length', 5, 'width', .6, 'ends', [2], 'color', [1 1 1]);
arrow([21, 2.6], [23.9, 2.6], 'BaseAngle', 70, 'length', 5, 'width', .6, 'ends', [1], 'color', [1 1 1]);
text([5 5],[2.6 2.6],'El Ni\~no event', 'interpreter', 'latex', 'fontsize', 21, 'color', [1 1 1]);

arrow([24.1, 2.6], [29, 2.6], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [2], 'color', [1 1 1]);
arrow([43, 2.6], [47.9, 2.6], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [1], 'color', [1 1 1]);
text([30 30],[2.6 2.6],'spin-down', 'interpreter', 'latex', 'fontsize', 22, 'color', [1 1 1]);

set(gca, ...
    'Box'             , 'off'         , ...
    'TickDir'         , 'out'         , ...
    'TickLength'      , [.01 .01]     , ...
    'XMinorTick'      , 'off'         , ...
    'YMinorTick'      , 'off'         , ...
    'YGrid'           , 'on'          , ...
    'YTick'           , -2:1:2   , ...     % define grid, every 1000 m a grid line
    'XColor'          , RdYlBu(60,:)  , ...
    'YColor'          , RdYlBu(60,:)  , ...
    'ticklabelinterpreter', 'latex'   , ...
    'LineWidth'       , 1.25);

set([hYLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(50,:));

yyaxis right
% --- plot time series: anomalous warm water volume anomlies
h4 = plot(linspace(1,48,48), wwv_pnEXP1, 'color', [RdYlBu(10,:)], 'linewidth', 2.5); hold on;
ylim([-2.5e14 2.5e14]);
hYLabel = ylabel(['[m$^{3}$]'],'interpreter','latex'); % add unit

vline = line([24 24], [-3e14 3e14]); % draw boundary line 
set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);
vline = line([48 48], [-3e14 3e14]); % draw boundary line 
set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);
    
% set custom colours
set(gca, ...
    'XColor'          , RdYlBu(60,:)  , ...
    'YColor'          , RdYlBu(60,:));
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 21);
set([hYLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(10,:));

h66 = legend([h22 h2 h4], 'N34$_{\mathrm{ideal.}}$','N34 anomaly', ...
    'WWV anomaly', 'WWV$_{\mathrm{non-mirrored}}$', ...
    'location', 'southeast', 'orientation', 'vertical');
set(h66, 'interpreter', 'latex', 'fontsize', 18, 'textcolor', RdYlBu(60,:));


% properties for full plot
hXLabel = xlabel('Month', 'interpreter', 'latex');
set(gca, 'Xtick', [1, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72]);

 
% correlation between simulated and idealized nino34
 a = interp1(1:48, mom_nino_pnEXP1, linspace(1, 48, 5840), 'linear'); 
 b(1:2920) = pnnino_comp_highres(1,:); b(2921:5840) = pnnino_comp_highres_spin(1,:);
 corrcoef(a, b)
 % add text labels
 text(-7, 3.25e14, 'c) El Ni\~no event', 'interpreter', 'latex', 'Fontsize', 30, ...
     'color', [0 0 0]);
 text(16, 1.75e14, '$r = 0.97$', 'interpreter', 'latex', 'Fontsize', 21, ...
     'color', RdYlBu(50,:));


    xlim([1 48]);


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% plot for La Nina event (subplots b and d)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% (this is just copy-paste from the two subplots above with different data)

% sst anomaly
subplot(2,2,2)
yyaxis left    
h5 = plot(linspace(1,48,48), sst_pnEXP2, 'color', RdYlBu(60,:), 'linewidth', 2.5); hold on;
ylim([-0.15 0.15]);
hYLabel = ylabel(['[$^{\circ}$C]'],'interpreter','latex');
set([hYLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(60,:));

% draw arrows to denote El Nino and La Nina events
arrow([1.1, 0.165], [4, 0.165], 'BaseAngle', 70, 'length', 5, 'width', .6, 'ends', [2], 'color', RdYlBu(60,:));
arrow([21, 0.165], [23.9, 0.165], 'BaseAngle', 70, 'length', 5, 'width', .6, 'ends', [1], 'color', RdYlBu(60,:));
text([4 4],[0.165 0.165],'La Ni\~na event', 'interpreter', 'latex', 'fontsize', 21, 'color', RdYlBu(60,:));

arrow([24.1, 0.165], [29, 0.165], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [2], 'color', 'k');
arrow([43, 0.165], [47.9, 0.165], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [1], 'color', 'k');
text([30 30],[0.165 0.165],'spin-down', 'interpreter', 'latex', 'fontsize', 22, 'color', 'k');

set(gca, ...
    'Box'             , 'off'         , ...
    'TickDir'         , 'out'         , ...
    'TickLength'      , [.01 .01]     , ...
    'XMinorTick'      , 'off'         , ...
    'YMinorTick'      , 'off'         , ...
    'YGrid'           , 'on'          , ...
    'YTick'           , -0.15:0.05:0.15   , ...     % define grid, every 1000 m a grid line
    'XColor'          , RdYlBu(60,:)  , ...
    'YColor'          , RdYlBu(60,:)  , ...
    'ticklabelinterpreter', 'latex'   , ...
    'LineWidth'       , 1.25);

yyaxis right
h3 = plot(linspace(1,48,48), cumtrapz(net_sfc_global_pnEXP2)* ...
    60*60*24*365.25*4*squeeze(nansum(squeeze(nansum(areacello, 1)), 2)), 'color', RdYlBu(20,:), 'linewidth', 2.5); hold on;
set(gca, 'xtick',[]);
ylim([-3e23 3e23]);
hYLabel = ylabel(['[J]'],'interpreter','latex');

vline = line([24 24], [-5e23 5e23]); % draw boundary line between El Niño and La Niña simulation
set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);
vline = line([48 48], [-5e23 5e23]); % draw boundary line between El Niño and La Niña simulation
set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);

xlim([1 48]);
set(gca, ...
    'XColor'          , RdYlBu(60,:)  , ...
    'YColor'          , RdYlBu(60,:));
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 21);
set([hYLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(20,:));

hXLabel = xlabel('Month', 'interpreter', 'latex', 'color', [1 1 1]);
set(gca, 'Xtick', [1, 6, 12, 18, 24, 30, 36, 42, 48]);

text(-7, 4e23, 'b) La Ni\~na event', 'interpreter', 'latex', 'Fontsize', 30, ...
     'color', [0 0 0]);


subplot(2,2,4)
yyaxis left    
h2 = plot(linspace(1,48,48), mom_nino_pnEXP2, 'color', [RdYlBu(50,:)], 'linewidth', 2.5); hold on;
h22 = plot(linspace(1,24, 2920), -pnnino_comp_highres(1,:), 'color', RdYlBu(50,:), 'linewidth', 2.5, 'linestyle', '--'); hold on;
h23 = plot(linspace(25,48, 2920), -pnnino_comp_highres_spin(1,:), 'color', RdYlBu(50,:), 'linewidth', 2.5, 'linestyle', '--'); hold on;
ylim([-2.5 2.5]);
    hYLabel = ylabel(['[$^{\circ}$C]'],'interpreter','latex');

% draw arrows with white colour to keep aspect ratio of subplots the same
arrow([1.1, 2.6], [4, 2.6], 'BaseAngle', 70, 'length', 5, 'width', .6, 'ends', [2], 'color', [1 1 1]);
arrow([21, 2.6], [23.9, 2.6], 'BaseAngle', 70, 'length', 5, 'width', .6, 'ends', [1], 'color', [1 1 1]);
text([4 4],[2.6 2.6],'El Ni\~no event', 'interpreter', 'latex', 'fontsize', 21, 'color', [1 1 1]);

arrow([24.1, 2.6], [29, 2.6], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [2], 'color', [1 1 1]);
arrow([43, 2.6], [47.9, 2.6], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [1], 'color', [1 1 1]);
text([30 30],[2.6 2.6],'spin-down', 'interpreter', 'latex', 'fontsize', 22, 'color', [1 1 1]);

set(gca, ...
    'Box'             , 'off'         , ...
    'TickDir'         , 'out'         , ...
    'TickLength'      , [.01 .01]     , ...
    'XMinorTick'      , 'off'         , ...
    'YMinorTick'      , 'off'         , ...
    'YGrid'           , 'on'          , ...
    'YTick'           , -2:1:2   , ...     % define grid, every 1000 m a grid line
    'XColor'          , RdYlBu(60,:)  , ...
    'YColor'          , RdYlBu(60,:)  , ...
    'ticklabelinterpreter', 'latex'   , ...
    'LineWidth'       , 1.25);
set([hYLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(50,:));

yyaxis right
h4 = plot(linspace(1,48,48), wwv_pnEXP2, 'color', [RdYlBu(10,:)], 'linewidth', 2.5); hold on;
ylim([-2.5e14 2.5e14]);
hYLabel = ylabel(['[m$^{3}$]'],'interpreter','latex');

vline = line([24 24], [-3e14 3e14]); % draw boundary line between El Niño and La Niña simulation
set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);
vline = line([48 48], [-3e14 3e14]); % draw boundary line between El Niño and La Niña simulation
set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);

set(gca, ...
    'XColor'          , RdYlBu(60,:)  , ...
    'YColor'          , RdYlBu(60,:));
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 21);
set([hYLabel], 'interpreter', 'latex', 'Fontsize', 21, 'color', RdYlBu(10,:));

% properties for full plot
hXLabel = xlabel('Month', 'interpreter', 'latex');
set(gca, 'Xtick', [1, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72]);

text(-7, 3.25e14, 'd) La Ni\~na event', 'interpreter', 'latex', 'Fontsize', 30, ...
     'color', [0 0 0]);


xlim([1 48]);



% % save plot
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
    directory = 'C:/Users/Maurice Huguenin/Desktop/';
    print('-dpng','-r300', [directory f1]);


