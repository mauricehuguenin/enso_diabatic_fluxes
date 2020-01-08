%% Purpose: Create Warm Water Volume (WWV) budget for the period 1979-2016
%           from the ACCESS-OM2 JRA55-do iaf simulation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %                         
%                        Maurice Huguenin-Virchaux                        %                         
%                     m.huguenin-virchaux@unsw.edu.au                     %                         
%                          29.07.2019, 08:48 AEST                         %                         
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
base = '/srv/ccrc/data67/z5180028/MSC_thesis_access_output/'; % new data from ACCESS stored in here
% load workspace, e.g. my warm water volume mask
load('/srv/ccrc/data15/z5180028/workspace_EXP1_analysis_composite_ninos_rev3.mat');

first = 'output020'; % first output in time period, output020 = year 1979
first_num = 20;
ts_len = 57; % last output entry, i.e. ts_len = 25 -> output025

%% [1.3 sec.] load in workspace
tic;
% colours for plotting
antarctica = nature;
RdBu_short  = cbrewer('div', 'RdBu', 21, 'PCHIP');
RdYlGn  = cbrewer('div', 'RdYlGn', 60, 'PCHIP');
Reds  = cbrewer('seq', 'Reds', 16, 'PCHIP');
Blues  = cbrewer('seq', 'Blues', 16, 'PCHIP');

wwv_mask_2S_borneo = wwv_mask_2S_borneo(91:811,478:518);

clear ans basin_mask bathymetry filename2017 wwv_mask 
clear wwv_mask_2S areacello volcello lon lat

% temperature space array -> the 74 bins
neutral_rho = getnc([base first '/ocean/ocean_wmass.nc'], 'neutral'); 
area_t = permute(getnc([base first '/ocean/ocean_grid.nc'], 'area_t'), [2 1]);
area_t = area_t(91:811,478:518); % only select WWV region area
lon = permute(getnc([base first '/ocean/ocean_grid.nc'], 'geolon_t'), [2 1]);
lat = permute(getnc([base first '/ocean/ocean_grid.nc'], 'geolat_t'), [2 1]);

% [lat lon] = meshgrid(lat, lon); % loading in lon and lat values for ACCESS
lon = lon(91:811,478:518); % in fact they are the same as in MOM025
lat = lat(91:811,478:518);

toc;         


%% [43/489 sec. for 11/27 yrs] preamble and data preparation for warm water volume
tic;

% loading in warm water volume, i.e. volume of water above 20 degrees C in 
% the Pacific region 5S - 5N and 120E - 80W, same as McGregor et al., 2014

warning off;
% load in potential temperature
% for the warm water volume I need one more extra data entry, i.e. I also
% load in data from December 1978 because later on I take the derivative
% and dWWV would otherwise be too short
thetao = squeeze(permute(getnc([base  'output019/ocean/ocean_snap.nc'], ...
            'temp', [12,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
%                   [month, depth, lat, lon], [end], [stride]
% load in depth levels
dzt = squeeze(permute(getnc([base 'output019/ocean/ocean_snap.nc'], ...
            'dzt', [12,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
   
for i = (first_num):ts_len % load in subsequent years and concatenate
    if i < 10
        a = ['0' num2str(i)];
    else
        a = num2str(i);
    end
    thetao = cat(4,thetao, squeeze(permute(getnc([base 'output0' a '/ocean/ocean_snap.nc'], ...
           'temp', [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1]))); 
    dzt = cat(4,dzt, squeeze(permute(getnc([base 'output0' a '/ocean/ocean_snap.nc'], ...
        'dzt', [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])));     
end
fprintf('~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ \n');
fprintf(['--- WWV --- :' first ' to output0' num2str(i)])
size(thetao) % ok, it gets concatenated correctly    



% now finding temperature within it which is higher than 20 degrees
mask = thetao - 273.15; % convert from degrees K to degrees C

indices = find(mask < neutral_rho(47));  mask(indices) = 0; clear indices;
indices = find(mask >= neutral_rho(47)); mask(indices) = 1; clear indices;

% mask is now my array containing all grid cells 5S-5N and 100E-60W which
% have temperatures higher than 20 degrees
volume = nan(721, 41, 50, (ts_len-first_num+1)*12);
for z = 1:50
    for t = 1:(ts_len-first_num+1)*12 + 1
        volume(:,:,z,t)  = mask(:,:,z,t) .* area_t .* dzt(:,:,z,t) .* wwv_mask_2S_borneo;
    end
end

size(volume)
% calculate anomaly
wwv = squeeze(nansum(squeeze(nansum(squeeze(nansum(volume, 1)), 1)), 1));
% plot(wwv); hold on;

% by taking the derivative, I need one more data entry -> I copy the first
% entry twice but this is not an issue as my final analysis will be
% covering the period within the time series 1979-2016
% dV(1) = wwv(1);               
% dV(2:length(wwv)+1) = wwv;

dV = diff(wwv) / (30*24*60*60) / 1e6; % nochmals durch 1 Million um in Sverdrup zu plotten

% % ok, now I have the derivative of the wwv time series but not as an anomaly
% % next step: create anomaly
% clear dVc;
% for t = 1:12 % loop over months
%     dVc(t) =  nanmean(dV(t:12:end)); % calculate mean over all individual months
% end
% % repeat climatology n-times to make it as long as the wwv data
% dVc = repmat(dVc,1,(ts_len-first_num+1));
% 
% % calculate anomaly
% dV = dV - dVc; clear dVc;
% figure(2);
% plot(dV); hold on;


% safety test with observations
srv2 = '/home/z5180028/MSC_thesis/access_matlab_scripts/'
load([srv2 'workspace_NOAA_PMEL_WWVa.mat']);

wwv_obs(1:12) = nan;
wwv_obs(13:(ts_len-first_num+1)*12) = NOAA_PMEL_WWVa(1:(ts_len-first_num+1)*12-12,2);


% % load in upper ocean temperature as well (T300)
% clear thetao;
% thetao = squeeze(permute(getnc([base first '/ocean/ocean.nc'], ...
%             'temp', [1,1,478,91], [12,20,518,811], [1,1,1,1]), [4 3 2 1])); 
% 
% for i = (first_num+1):ts_len % load in subsequent years and concatenate
%     if i < 10
%         a = ['0' num2str(i)];
%     else
%         a = num2str(i);
%     end
%     thetao = cat(4,thetao, squeeze(permute(getnc([base 'output0' a '/ocean/ocean_snap.nc'], ...
%            'temp', [1,1,478,91], [12,20,518,811], [1,1,1,1]), [4 3 2 1]))); 
% end
% 
% fprintf('~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ \n');
% fprintf(['T300:' first ' to output0' num2str(i)])
% 
% % depth-weighting here
% for z = 1:20
%     if z == 1
%         thetao(:,:,z,:) = thetao(:,:,z,:);
%     else
%         thetao(:,:,z,:) * (lev_bnds(z)-lev_bnds(z-1)) ./ lev_bnds(20);
%     end
% end
% 
% % calculate depth-averaged temperature over the WWV region
% thetao = squeeze(nanmean(thetao,3));
% for t = 1:((ts_len-first_num+1)*12)
%     thetao(:,:,t) = thetao(:,:,t) .* wwv_mask_2S_borneo;
% end
% 
% 
% 
% % take mean over lon and lat and convert to °C
% T300 = squeeze(nanmean(squeeze(nanmean(thetao - 273.15,1)),1));   
% 
% 
% 
% create anomaly for simulated wwv and T300
clear wwvc T300c;
for t = 1:12 % loop over months
    wwvc(t)         =  nanmean(wwv(t:12:end)); % calculate mean over all Jan/Feb/etc.
%     T300c(t)         =  nanmean(T300(t:12:end)); % calculate mean over all Jan/Feb/etc.
end

% repeat climatology n-times to make it as long as the wwv data
wwvc         =  repmat(wwvc,1,(ts_len-first_num+1));
% T300c         =  repmat(T300c,1,(ts_len-first_num+1));

wwva = wwv(2:(ts_len-first_num+1)*12+1) - wwvc; % create anomaly but only select Jan 1979-Dec 2016 period
% T300a = T300 - T300c; % calculate anomaly and convert to °C

% creating 5-month running mean of wwv as in Meinen and McPhaden 2000
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
wwva = filter(b,1,wwva);


% plot(wwva); hold on; plot(wwv_obs)


%% [superfast] plot observed, simulated WWVa 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
ticks = linspace(1,length(wwva),length(wwva));
figure('units', 'pixels', 'position', [0 0 2300 800]);
subplot(2,1,1)
% -------------------------------------------------------------------------
h6 = plot(ticks, wwv_obs, 'color', 'k', 'linewidth', 2); hold on;
h2 = plot(ticks, wwva, 'color', RdYlBu(10,:), 'linewidth', 2); hold on;
% h3 = plot(ticks, T300a, 'color', RdYlBu(50,:), 'linewidth', 2); hold on;

% hXLabel = xlabel('Year');
hYLabel = ylabel(['m$^3$'],'interpreter','latex');
set([hYLabel], 'interpreter', 'latex', 'Fontsize', 25, 'color', RdYlBu(60,:));
set(hYLabel, 'color', [0 0 0]);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -1e18:2e14:1e18   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

set(gca,'XLim',[1 length(wwv)],'XTick',[0:24:length(wwv)])
set(gca,'xticklabel',({'1979','1981','','1985',...
    '','1989','','1993', '','1997','','2001','',...
    '2005','','2009','','2013','', '2017'}))
   
ylim([-4e14 4e14]);
xlim([1 length(wwva)]);
pbaspect([5 1 1]);                        % aspect ratios: x, y, z

clear hTitle hXLabel hYLabel;


% text(-17, 30, 'a) ACCESS-OM2: adiabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
%     'color', [0 0 0]);

h55 = legend([h6 h2], ...
    'WWVa NOAA PMEL', 'WWVa ACCESS-OM2', ...
    'location', 'southwest', 'orientation', 'horizontal');
set(h55, 'interpreter', 'latex', 'fontsize', 11,'EdgeColor', [.83 .83 .83]);

r = corrcoef(wwva(13:end), wwv_obs(13:end));
text(230, 2.5e14,['r = ' num2str(r(2))], 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);


% % finished fancy plot
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
% directory = '/home/z5180028/MSC_thesis/access_figures/';
% print('-dpng','-r300', [directory 'WWV_and_OBS_time_series_1979-2016']);
% % 

% clear workspace
clear ans bathymetry cto ctoc dzt dztc maskc s so soc;
clear t test thetao thetaoc volcello z;
clear volume volume_res volumec volumec_res temp_res temp_resc;
clear a A b f1 f2 f3 i;
clear a ans filename2018 i mask maskc NOAA_PMEL_WWVa t thetao thetaoc;
clear z clear wwv_filter test h2 h55 h6 r;
toc;


%% [42/155 sec. for 11/27 years] calculating the horizontal transport 
tic; 
rho_0 = 1035;       % reference air density [kg m^3]
cp = 3992.1;        % specific heat capacity of seawater [J kg^-1 K^-1]
dT = 0.5;           % the temperature difference of the binned temperature 
nlev = 46;

% load in volume transport above 20°C (i.e. from temperature level 46 all
% the way up to 74)
ty_trans_nrho = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'ty_trans_nrho', [1,nlev,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)); 
% load in Gent-McWilliams skew diffusion and submesoscale diffusion as well
temp_yflux_gm = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'temp_yflux_gm_on_nrho', [1,nlev,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)); 
temp_yflux_submeso = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'temp_yflux_submeso_on_nrho', [1,nlev,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)); 

for i = (first_num+1):ts_len % load in subsequent years and concatenate
    if i < 10
        a = ['0' num2str(i)];
    else
        a = num2str(i);
    end
    ty_trans_nrho = cat(3,ty_trans_nrho, ...
            squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
            'ty_trans_nrho', [1,nlev,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)));
    temp_yflux_gm = cat(3,temp_yflux_gm, ...
            squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
            'temp_yflux_gm_on_nrho', [1,nlev,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)));
    temp_yflux_submeso = cat(3,temp_yflux_submeso, ...
            squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
            'temp_yflux_submeso_on_nrho', [1,nlev,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)));
end
fprintf('~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ \n');
fprintf(['horizontal: ' first ' to output0' num2str(i)])
size(ty_trans_nrho) % ok, it gets concatenated correctly    


% convert gm and submeso heat fluxes to volume fluxes
% [W] = [J / s]
% [J m3 kg K / s kg J K] = [m3 / s]
for t = 1:(ts_len-first_num+1)*12
    temp_yflux_gm(:,:,t) = (1 / (rho_0 * cp * dT)) .* temp_yflux_gm(:,:,t);
    temp_yflux_submeso(:,:,t) = (1 / (rho_0 * cp * dT)) .* temp_yflux_submeso(:,:,t);
end


% --- checking if all the meridional transport calculation is done
% --- correctly



% construction site
testmap(lon, lat, ty_trans_nrho(:,:,12))
test = ty_trans_nrho(:,:,12);
test(:,40:41) = 1e10;
testmap(lon, lat, test)

% check after holidays


for i = 1:721 % loop over all the longitude values
    % transport across 5°N to higher latitudes
    north(i,:) = squeeze(ty_trans_nrho(i,41,:)) * wwv_mask_2S_borneo(i,41);
        north_gm(i,:) = squeeze(temp_yflux_gm(i,41,:)) * wwv_mask_2S_borneo(i,41);
        north_submeso(i,:) = squeeze(temp_yflux_submeso(i,41,:)) * wwv_mask_2S_borneo(i,41);
end
% minus sign for transport in the Southern hemisphere as we want transport
% across 5°S bzw. 2°S to higher latitudes down south and ty_trans is
% defined as the South-North water mass transport
for i = 50:127
    % transport across the ITF to higher latitudes
    ITF(i,:) = squeeze(-ty_trans_nrho(i,12,:)) * wwv_mask_2S_borneo(i,41);
        ITF_gm(i,:) = squeeze(-temp_yflux_gm(i,12,:)) * wwv_mask_2S_borneo(i,41);
        ITF_submeso(i,:) = squeeze(-temp_yflux_submeso(i,12,:)) * wwv_mask_2S_borneo(i,41);
end
for i = 128:721 
    % transport across 5°S to higher latitudes
    south(i,:) = squeeze(-ty_trans_nrho(i,1,:)) * wwv_mask_2S_borneo(i,41);
        south_gm(i,:) = squeeze(-temp_yflux_gm(i,1,:)) * wwv_mask_2S_borneo(i,41);
        south_submeso(i,:) = squeeze(-temp_yflux_submeso(i,1,:)) * wwv_mask_2S_borneo(i,41);
end

north = nansum(north,1); 
    north_gm = nansum(north_gm,1); 
    north_submeso = nansum(north_submeso,1); 
ITF = nansum(ITF,1); 
    ITF_gm = nansum(ITF_gm,1); 
    ITF_submeso = nansum(ITF_submeso,1); 
south = nansum(south,1); 
    south_gm = nansum(south_gm,1); 
    south_submeso = nansum(south_submeso,1); 

horizontal = -(north + south) / 1e9; % total meridional transport in [Sv]
    horizontal_gm = -(north_gm + south_gm) / 1e9; % total meridional transport in [Sv]
    horizontal_submeso = -(north_submeso + south_submeso) / 1e9; % total meridional transport in [Sv]
ITF = -ITF / 1e9; 
    ITF_gm = -ITF_gm / 1e9; 
    ITF_submeso = -ITF_submeso / 1e9;
skew = horizontal_gm + horizontal_submeso + ITF_gm + ITF_submeso;

% this is important here: I define negative transport during EN, i.e. WWV
% transport OUT of the region, and the opposite during LN, i.e. WWV into
% the region

% -> that means: during discharge phases, the horizontal transport is
% negative!

% % subtracting the climatological time series
% clear horizontalc horizontal_gmc horizontal_submesoc ITFc ITF_gmc ITF_submesoc;
% for t = 1:12 % loop over months
%     horizontalc(t) =  nanmean(horizontal(t:12:end)); % calculate mean over all individual months
%         horizontal_gmc(t) =  nanmean(horizontal_gm(t:12:end)); 
%         horizontal_submesoc(t) =  nanmean(horizontal_submeso(t:12:end)); 
%     ITFc(t) =  nanmean(ITF(t:12:end));
%     ITF_gmc(t) =  nanmean(ITF_gm(t:12:end)); 
%     ITF_submesoc(t) =  nanmean(ITF_submeso(t:12:end)); 
% 
% end
% % repeat climatology n-times to make it as long as the wwv data
% horizontalc = repmat(horizontalc,1,(ts_len-first_num+1));
%     horizontal_gmc = repmat(horizontal_gmc,1,(ts_len-first_num+1));
%     horizontal_submesoc = repmat(horizontal_submesoc,1,(ts_len-first_num+1));
% ITFc = repmat(ITFc,1,(ts_len-first_num+1));
%     ITF_gmc = repmat(ITF_gmc,1,(ts_len-first_num+1));
%     ITF_submesoc = repmat(ITF_submesoc,1,(ts_len-first_num+1));
% 
% % create anomalies
% horizontal = horizontal - horizontalc; clear horizontalc;
%     horizontal_gm = horizontal_gm - horizontal_gmc; clear horizontal_gmc;
%     horizontal_submeso = horizontal_submeso - horizontal_submesoc; clear horizontal_submesoc;
% ITF = ITF - ITFc; clear ITFc;
%     ITF_gm = ITF_gm - ITF_gmc; clear ITF_gmc;
%     ITF_submeso = ITF_submeso - ITF_submesoc; clear ITF_submesoc;

% % calculate the vertical volume flux as the residual
% vertical = dV - horizontal - horizontal_gm - horizontal_submeso - ...
%     ITF - ITF_gm - ITF_submeso;


% safespot: 08. 03. 2018, 16:16 % wow, that was a while ago!
% safespot: 06. 09. 2019, 14:33
% safespot: 25. 09. 2019, 09:36
% safespot: 02. 10. 2019, 16:55
% safespot: 03. 10. 2019, 16:47 
% safespot: 08. 10. 2019, 15:38 I do not see where there could be a 
%                               problem until here
% safespot: 17. 10. 2019, 10:52

clear a ans b i ITFa northa southa t
clear north north_gm north_sub south south_gm south_sub;
clear ty_trans ty_trans_gm ty_trans_nrho_gm ty_trans_nrho_sub ty_trans_sub;
clear north_submeso south_submeso temp_yflux_gm 
clear temp_yflux_submeso ty_trans_nrho windowSize horizontal_gm ...
    horizontal_submeso ITF_gm ITF_submeso
toc;


%% [94/190 sec. for 11/27 years] creating budget with water mass transformation framework 
tic; % start clock timer
rho_0 = 1035;       % reference air density [kg m^3]
cp = 3992.1;        % specific heat capacity of seawater [J kg^-1 K^-1]
dT = 1;           % the temperature difference of the binned temperature 
nlev = 46; nlev2 = 47;

% diagnostics for GM - vertical mixing
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cbt = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'temp_vdiffuse_diff_cbt_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3));
non = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'temp_nonlocal_KPP_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3));

% diagnostics for GF - surface forcing
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sbc = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'temp_vdiffuse_sbc_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3));
sw = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'sw_heat_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3));
frazil = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'frazil_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3));
eta = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'temp_eta_smooth_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3));

% diagnostics for J - surface volume flux
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mass = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'mass_pmepr_on_nrho', [1,nlev,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]),3));

% diagnostics for GN - neutral surface volume flux
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
diff = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'neutral_diffusion_on_nrho_temp', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3));
k33 = squeeze(nansum(permute(getnc([base first '/ocean/ocean_wmass.nc'], ...
    'temp_vdiffuse_k33_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3));

for i = (first_num+1):ts_len % load in subsequent years and concatenate
    if i < 10
        a = ['0' num2str(i)];
    else
        a = num2str(i);
    end
    % GM - vertical mixing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cbt = cat(3,cbt, ...
        squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
        'temp_vdiffuse_diff_cbt_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3)));   
    non = cat(3,non, ...
        squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
        'temp_nonlocal_KPP_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3)));  
   
    % GF - surface forcing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sbc = cat(3,sbc, ...
        squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
        'temp_vdiffuse_sbc_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3)));   
    sw = cat(3,sw, ...
        squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
        'sw_heat_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3)));   
    frazil = cat(3,frazil, ...
        squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
        'frazil_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3))); 
    eta = cat(3,eta, ...
        squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
        'temp_eta_smooth_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3))); 
    
    % J - surface volume flux ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mass = cat(3,mass, ...
        squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
        'mass_pmepr_on_nrho', [1,47,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]),3))); 
    
    % N - neutral diffusion volume flux ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    diff = cat(3,diff, ...
        squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
        'neutral_diffusion_on_nrho_temp', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3))); 
    k33 = cat(3,k33, ...
        squeeze(nansum(permute(getnc([base 'output0' a '/ocean/ocean_wmass.nc'], ...
        'temp_vdiffuse_k33_on_nrho', [1,nlev,478,91], [12,nlev2,518,811], [1,1,1,1]), [4 3 2 1]),3)));  

end
fprintf('~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ ~~~~~~~~~~ \n');
fprintf(['WMT: ' first ' to output0' num2str(i)])
size(cbt) % ok, it gets concatenated correctly    


% ~~~~~~~~~~ the calculation for the WMT volume fluxes is here ~~~~~~~~~~ %
clear GM GF GJ GN GR
for t = 1:(ts_len-first_num+1)*12
    GM(:,:,t) = (cbt(:,:,t) + non(:,:,t)) .* ...
        wwv_mask_2S_borneo .* (1 ./ (rho_0 .* cp .* dT)) .* ...
        area_t ./ 1e6;
    
    GF(:,:,t) = (sbc(:,:,t) + sw(:,:,t) + frazil(:,:,t) + eta(:,:,t)) .* ...
        wwv_mask_2S_borneo .* (1 ./ (rho_0 .* cp .* dT)) .* ...
        area_t ./ 1e6;
	% durch 1e6 teilen für Einheit Sv [10^6 m^3 s^-1]
    % units are in = W m^-2
    
    GJ(:,:,t) = (mass(:,:,t)) .* wwv_mask_2S_borneo .* ...
        1e-3 ./ 1e6; % umrechnen von [kg s�?�1] zu [m^3 s^-1]
       
    GN(:,:,t) = (diff(:,:,t) + k33(:,:,t)) .* ...
    wwv_mask_2S_borneo .* (1 ./ (rho_0 .* cp .* dT)) .* ...
    area_t ./ 1e6;
    % umrechnen von [kg s�?�1] zu [m^3 s^-1]
    end

mass = squeeze(nansum(squeeze(nansum(GJ,1)),1));
forcing = squeeze(nansum(squeeze(nansum(GF,1)),1));
mixing = squeeze(nansum(squeeze(nansum(GM,1)),1));
neutral = squeeze(nansum(squeeze(nansum(GN,1)),1));
% rest = squeeze(nansum(squeeze(nansum(GR,1)),1));


% % subtracting the climatological time series
% clear mixingc forcingc massc neutralc;
% for t = 1:12 % loop over months
%     massc(t) =  nanmean(mass(t:12:end)); % calculate mean over all individual months
%     forcingc(t) =  nanmean(forcing(t:12:end)); 
%     mixingc(t) =  nanmean(mixing(t:12:end)); 
%     neutralc(t) =  nanmean(neutral(t:12:end)); 
% 
% end
% % repeat climatology n-times to make it as long as the wwv data
% massc = repmat(massc,1,(ts_len-first_num+1));
% forcingc = repmat(forcingc,1,(ts_len-first_num+1));
% mixingc = repmat(mixingc,1,(ts_len-first_num+1));
% neutralc = repmat(neutralc,1,(ts_len-first_num+1));
% 
% % create anomalies
% mass = mass - massc; clear massc;
% forcing = forcing - forcingc; clear forcingc;
% mixing = mixing - mixingc; clear mixingc;
% neutral = neutral - neutralc; clear neutralc;

% % calculate the numerical volume flux as the residual
% implicit = dV - horizontal - horizontal_gm - horizontal_submeso - ...
%     ITF - ITF_gm - ITF_submeso - mass - forcing - mixing - neutral;
% 
toc; % end clock timer


%% [very fast] create anomalies, close budget and calculate 5-month running mean

clearvars -except ...
    dV horizontal ITF mass forcing mixing neutral skew numerical ...
    dVa horizontala ITFa massa forcinga mixinga neutrala skewa ...
    dVaf horizontalaf ITFaf massaf forcingaf mixingaf neutralaf skewaf ...
    antarctica table_EN table_LN ts_len implicit first_num ticks wwv EN LN ...
    RdYlBu RdYlGn Reds Blues;

% create anomalies
for t = 1:12 % loop over months
    dVc(t)         =  nanmean(dV(t:12:end)); % calculate mean over all Jan/Feb/etc.
    horizontalc(t) =  nanmean(horizontal(t:12:end)); 
    ITFc(t)        =  nanmean(ITF(t:12:end)); 
    massc(t)       =  nanmean(mass(t:12:end));
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    forcingc(t)    =  nanmean(forcing(t:12:end)); 
    mixingc(t)     =  nanmean(mixing(t:12:end)); 
    neutralc(t)    =  nanmean(neutral(t:12:end));
    skewc(t)       =  nanmean(skew(t:12:end));
end

% repeat climatology n-times to make it as long as the wwv data
dVc         =  repmat(dVc,1,(ts_len-first_num+1));
horizontalc =  repmat(horizontalc,1,(ts_len-first_num+1));
ITFc        =  repmat(ITFc,1,(ts_len-first_num+1));
massc       =  repmat(massc,1,(ts_len-first_num+1));
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
forcingc    =  repmat(forcingc,1,(ts_len-first_num+1));
mixingc     =  repmat(mixingc,1,(ts_len-first_num+1));
neutralc    =  repmat(neutralc,1,(ts_len-first_num+1));
skewc       =  repmat(skewc,1,(ts_len-first_num+1));

% calculate anomalies
dVa = dV - dVc;
horizontala = horizontal - horizontalc;
ITFa = ITF - ITFc;
massa = mass - massc;
% -------------------------------------------------------------------------
forcinga = forcing - forcingc;
mixinga = mixing - mixingc;
neutrala = neutral - neutralc;
skewa = skew - skewc;


% numerical mixing as the residual
numerical = dVa - horizontala - ITFa - massa - ...
    forcinga - mixinga - neutrala - skewa;
plot(numerical)

% calculate 5-month running mean of only the anomalous time series
windowSize = 5; b = (1/windowSize)*ones(1,windowSize);
dVaf = filter(b, 1, dVa);
horizontalaf = filter(b,1,horizontala); 
ITFaf = filter(b, 1, ITFa);
massaf = filter(b,1,massa);
% -------------------------------------------------------------------------    
forcingaf = filter(b,1,forcinga);
mixingaf = filter(b,1,mixinga);
neutralaf = filter(b,1,neutrala);
skewaf = filter(b,1,skewa);
numerical = filter(b,1,numerical);

clearvars -except ...
    dV horizontal ITF mass forcing mixing neutral skew numerical ...
    dVc horizontalc ITFc massc forcingc mixingc neutralc skewc...
    dVa horizontala ITFa massa forcinga mixinga neutrala skewa ...
    dVaf horizontalaf ITFaf massaf forcingaf mixingaf neutralaf skewaf ...
    antarctica table_EN table_LN ts_len implicit first_num ticks wwv EN LN ...
    RdYlBu RdYlGn Reds Blues;


%% clean up  workspace and summarize in tables

% preparing data for table with percentages
if ts_len == 57;
    % El Niño events in 1982/83, 1997/98 and 2015/16
    EN = [48,61, 224, 238, 438, 452];
    % La Niña events in 1988/89, 2077/08 and 2010/11
    LN = [114, 129, 343, 356, 382, 390];
else
     EN = [48,61,48,61,48,61];       
     LN = [114,129,114,129,114,129];
end

% create time series for conversion
%        J,  F,  M,  A,  M,  J,  J,  A,  S,  O,  N,  D
time = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
time = repmat(time,1,(ts_len-first_num+1)); % repeat time series

seconds = 86400;        % so many seconds per day
time = time .* seconds; % so many seconds each month

table_EN = nan(10,8);
table_LN = nan(10,8);
% take sum over full discharge period during EN 192/1983
for i = 1:3                 % select always two points of the time series
                             % i.e. El Nino event No.1 = 37:60 in the time
                             % series
    if i == 1
        col = 1;
    elseif i == 2
        col = 4;
    elseif i == 3
        col = 7;
    end
    % create table for discharge periods during El Nino
    table_EN(1,col) = nansum(dVaf(1,EN(i*2-1):EN(i*2)) .* time(EN(i*2-1):EN(i*2)));

    table_EN(3,col) =nansum(horizontalaf(EN(i*2-1):EN(i*2)) .* time(EN(i*2-1):EN(i*2)));
    table_EN(4,col) =nansum(ITFaf(EN(i*2-1):EN(i*2)) .* time(EN(i*2-1):EN(i*2)));
    table_EN(5,col) =nansum(massaf(EN(i*2-1):EN(i*2)) .* time(EN(i*2-1):EN(i*2)));

    table_EN(7,col) = nansum(forcingaf(EN(i*2-1):EN(i*2)) .* time(EN(i*2-1):EN(i*2)));
    table_EN(8,col) = nansum(mixingaf(EN(i*2-1):EN(i*2)) .* time(EN(i*2-1):EN(i*2)));
    table_EN(9,col) = nansum(neutralaf(EN(i*2-1):EN(i*2)) .* time(EN(i*2-1):EN(i*2)));
    table_EN(10,col) = nansum(skewaf(EN(i*2-1):EN(i*2)) .* time(EN(i*2-1):EN(i*2)));
    table_EN(11,col) = nansum(numerical(EN(i*2-1):EN(i*2)) .* time(EN(i*2-1):EN(i*2)));

    % same procedure but for the recharge periods during La Nina
    table_LN(1,col) = nansum(dVaf(1,LN(i*2-1):LN(i*2)) .* time(LN(i*2-1):LN(i*2)));

    
    table_LN(3,col) =nansum(horizontalaf(LN(i*2-1):LN(i*2)) .* time(LN(i*2-1):LN(i*2)));
    table_LN(4,col) =nansum(ITFaf(LN(i*2-1):LN(i*2)) .* time(LN(i*2-1):LN(i*2)));
    table_LN(5,col) =nansum(massaf(LN(i*2-1):LN(i*2)) .* time(LN(i*2-1):LN(i*2)));

    table_LN(7,col) = nansum(forcingaf(LN(i*2-1):LN(i*2)) .* time(LN(i*2-1):LN(i*2)));
    table_LN(8,col) = nansum(mixingaf(LN(i*2-1):LN(i*2)) .* time(LN(i*2-1):LN(i*2)));
    table_LN(9,col) = nansum(neutralaf(LN(i*2-1):LN(i*2)) .* time(LN(i*2-1):LN(i*2)));
    table_LN(10,col) = nansum(skewaf(LN(i*2-1):LN(i*2)) .* time(LN(i*2-1):LN(i*2)));
    table_LN(11,col) = nansum(numerical(LN(i*2-1):LN(i*2)) .* time(LN(i*2-1):LN(i*2)));

end

% get same values but in percentage as a contribution to the overall change
% in warm water volume
for f = 1:11
    % first EN
    table_EN(f,2) = table_EN(f,1) ./ table_EN(1,1) .* 100;
    % second EN
    table_EN(f,5) = table_EN(f,4) ./ table_EN(1,4) .* 100;
    % third EN
    table_EN(f,8) = table_EN(f,7) ./ table_EN(1,7) .* 100;

    table_LN(f,2) = table_LN(f,1) ./ table_LN(1,1) .* 100;
    table_LN(f,5) = table_LN(f,4) ./ table_LN(1,4) .* 100;
    table_LN(f,8) = table_LN(f,7) ./ table_LN(1,7) .* 100;
end

% convert to sverdrup and divide by 1e14 to get small numbers
% and round values to one-digit accuracy
table_EN(:,[1,4,7]) = round(table_EN(:,[1,4,7]) * 1e6 / 1e14,1);
% all values now with unit [10^{14} m^3]
% same for La Ninas
table_LN(:,[1,4,7]) = round(table_LN(:,[1,4,7]) * 1e6 / 1e14,1);
table_LN(:,1) = round(table_LN(:,1),1); % rounding to one digit

open table_LN
open table_EN

% round percentage values also, zero digit accuracy here as it is only 
% important how much the percentage value is
table_EN(:,[2,5,8]) = round(table_EN(:,[2,5,8]),0); % no digit accuracy
table_LN(:,[2,5,8]) = round(table_LN(:,[2,5,8]),0);


%% safety tests
figure('units', 'pixels', 'position', [0 0 1600 900]);

antarctica = nature;

% adiabatic volume fluxes
subplot(4,4,1)
plot(dVa, 'color', [0 0 0], 'linewidth', 1.5); hold on; 
title('dWWV/dt')
ylim([-35 25]);
subplot(4,4,2)
plot(horizontala, 'color', RdYlBu(10,:), 'linewidth', 1.5); hold on;
title('ty\_trans')
subplot(4,4,3)
plot(ITFa, 'color', RdYlGn(60,:), 'linewidth', 1.5); hold on; 
title('ITF')
ylim([-35 25]);
subplot(4,4,4)
plot(massa, 'color', RdYlBu(21,:), 'linewidth', 1.5); hold on;
title('surface volume')
% -------------------------------------------------------------------------    
subplot(4,4,5)
plot(forcinga, 'color', RdYlBu(60,:), 'linewidth', 1.5); hold on;title('vertical residual')
title('surface forcing')
ylim([-35 25]);
subplot(4,4,6)
plot(mixinga, 'color', antarctica(15,:), 'linewidth', 1.5); hold on;
title('vertical mixing')
ylim([-35 25]);
subplot(4,4,7)
plot(neutrala, 'color', RdYlGn(10,:), 'linewidth', 1.5); hold on;title('vertical residual')
title('neutral diff.')
subplot(4,4,8)
plot(skewa, 'color', RdYlGn(10,:), 'linewidth', 1.5); hold on;title('vertical residual')
title('skew diff.')

% subplot(4,4,12)
% plot(rest, 'color', RdYlGn(10,:), 'linewidth', 1.5); hold on;title('vertical residual')
% title('submeso + rivermix + sfc\_hflux')
% ylim([-35 25]);

subplot(4,4,9)
plot(numerical, 'color', [187,0,187]/255, 'linewidth', 1.5); hold on;
title('numerical mixing')
ylim([-35 25]);
line([0 (ts_len-first_num+1)*12],[0,0])


% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = '/home/z5180028/MSC_thesis/access_figures/';
print('-dpng','-r300', [directory 'WMT_time_series_safety_test']);
 % when only considering 11 years

 
%% ~~~~~~~~~~~~~~ plotting  routine for meridional figure ~~~~~~~~~~~~~~ %% 
directory = '/home/z5180028/MSC_thesis/access_figures/';
load([directory 'WMT_time_series_1979-2016.mat']);

% general info of what I load with the .mat file:
% there are four time series for all budget terms
% (1) dV   -> original time series, here the change in Warm Water Volume
% (2) dVa  -> the anomalous time series
% (3) dVc  -> the climatological time series
% (4) dVaf -> the 5-month running mean filtered anomalous time series
% -------------------------------------------------------------------------

load('/home/z5180028/MSC_thesis/access_matlab_scripts/workspace_regression_patterns_PC1_equal_nino34_rev2.mat', 'PC1');
nino = PC1; clear PC1;

boom; antarctica = nature;
figure('units', 'pixels', 'position', [0 0 2000 2000]);
a = [3.1 1 1]; % aspect ratio of subplots

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
subplot(2,1,1)
ticks = linspace(1,length(dV),length(dV));

h2 = plot(ticks, horizontalaf, 'color', RdYlBu(10,:), 'linewidth', 1.5); hold on;
h6 = plot(ticks, massaf, 'color', RdYlBu(21,:), 'linewidth', 1.5); hold on;
h8 = plot(ticks, ITFaf, 'color', RdYlGn(60,:), 'linewidth', 1.5); hold on;
h1 = plot(ticks, dVaf, 'color', [0 0 0], 'linewidth', 1.5); hold on;

% hXLabel = xlabel('Year');
hYLabel = ylabel(['Sv'],'interpreter','latex');
set([hYLabel], 'interpreter', 'latex', 'Fontsize', 25, 'color', RdYlBu(60,:));
% set(hYLabel, 'color', [0 0 0]);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -50:10:30   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

set(gca,'XLim',[1 length(wwv)],'XTick',[0:24:length(wwv)])
set(gca,'xticklabel',({'','','','',...
    '','','','', '','','','','',...
    '','','','','','', ''}))
   
% shade EN and LN periods
xbars = [EN(1) EN(2)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [EN(3) EN(4)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [EN(5) EN(6)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;

xbars = [LN(1) LN(2)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [LN(3) LN(4)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [LN(5) LN(6)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;

ylim([-26.5 22]);
xlim([1 length(dV)]);
pbaspect([a]);                        % aspect ratios: x, y, z

clear hTitle hXLabel hYLabel;


text(-17, 26, 'a) Change in WWV anomaly and ad', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);

h55 = legend([h1 h2 h8 h6], ...
    'Change in WWV anomaly', ...
    '$\mathcal{T}_{5^{\circ}\mathrm{N}+5^{\circ}\mathrm{S}}$: Meridional transport', ...
    '$\mathcal{T}_{ITF}$: Indonesian Throughflow', ...
    '$\mathcal{J}$: Surface volume', ...
    'location', 'southwest', 'orientation', 'horizontal');
set(h55, 'interpreter', 'latex', 'fontsize', 11, 'Edgecolor', [.83 .83 .83]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
subplot(2,1,2)

h11 = plot(ticks, neutralaf, 'color', RdYlGn(10,:), 'linewidth', 1.5); hold on;
h5 = plot(ticks, forcingaf, 'color', RdYlBu(60,:), 'linewidth', 1.5); hold on;
h4 = plot(ticks, mixingaf, 'color', antarctica(15,:), 'linewidth', 1.5); hold on;
h6 = plot(ticks, skewaf, 'color', antarctica(5,:), 'linewidth', 1.5); hold on;
% h99 = plot(ticks, rest, 'color', antarctica(1,:), 'linewidth', 2); hold on;
h7 = plot(ticks, numerical, 'color', [187,0,187]/255, 'linewidth', 1.5); hold on;
h8 = plot(ticks, nino*10, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--'); hold on;

hXLabel = xlabel('Year');
hYLabel = ylabel(['Sv, $^{\circ}$C'],'interpreter','latex', 'color', RdYlBu(60,:));
set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 25, 'color', RdYlBu(60,:));
% set(hYLabel, 'color', [0 0 0]);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 25);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -50:10:30   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

text(-17, 26, 'b) N34 index and diabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);
  
 
set(gca,'XLim',[1 length(dV)],'XTick',[0:12:length(dV)])

% shade EN and LN periods
xbars = [EN(1) EN(2)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [EN(3) EN(4)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [EN(5) EN(6)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;

xbars = [LN(1) LN(2)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [LN(3) LN(4)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [LN(5) LN(6)];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-50 50 50 -50], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;

set(gca,'XLim',[1 length(wwv)],'XTick',[0:24:length(wwv)])
set(gca,'xticklabel',({'1979','1981','','1985',...
    '','1989','','1993', '','1997','','2001','',...
    '2005','','2009','','2013','', '2017'}))

ylim([-26.5 22]);
xlim([1 length(dV)]);
pbaspect([a]);                        % aspect ratios: x, y, z% 
clear hTitle hXLabel hYLabel h1 h2;

h55 = legend([h8 h5 h4 h6 h11 h7], ...
    'N34', ...    
    '$\mathcal{G_F}$: Surface forcing', ...
    '$\mathcal{G_M}$: Vertical mixing', ...
    '$\mathcal{G_E}$: Skew diffusion', ...
    '$\mathcal{G_N}$: Neutral diffusion', ...    
    '$\mathcal{G_I}$: Numerical mixing', ...
    'location', 'southwest', 'orientation', 'horizontal');
set(h55, 'interpreter', 'latex', 'fontsize', 11, 'Edgecolor', [.83 .83 .83]);

clearvars -except ...
    dV horizontal ITF mass forcing mixing neutral skew numerical ...
    dVc horizontalc ITFc massc forcingc mixingc neutralc skewc...
    dVa horizontala ITFa massa forcinga mixinga neutrala skewa ...
    dVaf horizontalaf ITFaf massaf forcingaf mixingaf neutralaf skewaf ...
    antarctica table_EN table_LN ts_len implicit first_num ticks wwv EN LN ...
    RdYlBu RdYlGn Reds Blues nino a;
% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = '/home/z5180028/MSC_thesis/access_figures/';
print('-dpng','-r300', [directory 'WMT_time_series_1979-2016']);

% save([directory 'WMT_time_series_1979-2016.mat']);





