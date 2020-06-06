%% EXP1: Plotting latitude integrated OHCa during December of the first year

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  %
%                                                                         %
%                     hmaurice, 23.11.2017, 10:19 AEST                    %
%                                                                         %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  %



                         first = 'pnEXP1_composite_nino_windstress';   
%                          first = 'pnEXP2_composite_nina_windstress'; 

antarctica =    [150,   0,  17; 165,   0,  33; 200,   0,  40; ...
                 216,  21,  47; 247,  39,  53; 255,  61,  61; ...
                 255, 120,  86; 255, 172, 117; 255, 214, 153; ...
                 255, 241, 188; 255, 255, 255; 188, 249, 255; ...
                 153, 234, 255; 117, 211, 255;  86, 176, 255; ...
                  61, 135, 255;  40,  87, 255;  24,  28, 247; ...
                  30,   0, 230;  36,   0, 216;  45,   0, 200]*1./255;

%% load in workspace
% load workspace
f1 = 'pnEXP1_and_pnEXP2_timeseries_wwv_transport_20S_20N';
f2 = 'pnEXP1_and_pnEXP2_timeseries_wwv_transport_8S_8N';
f3 = 'timeseries_wwv_transport_5S_5N';
srv = 'H:\Maurice_ENSO_Data\'
load('E:\2017 Fr¸hlingssemester\Master Seminar 1\workspace_EXP1_analysis_composite_ninos_rev3.mat');
antarctica = nature;
RdBu_short  = cbrewer('div', 'RdBu', 21, 'PCHIP');
RdYlGn  = cbrewer('div', 'RdYlGn', 60, 'PCHIP');
Reds  = cbrewer('seq', 'Reds', 16, 'PCHIP');
Blues  = cbrewer('seq', 'Blues', 16, 'PCHIP');

% prepare mask                  % importante!       % wwv region 20S - 20N
                                % wwv_mask(117:840, 417:580) = 20S - 20N
                                % wwv_mask(117:840, 466:531) = 8S - 8N
                                % wwv_mask(117:840, 478:518) = 5S - 5N


wwv_mask = wwv_mask(91:811,478:518);
wwv_mask_2S = wwv_mask_2S(91:811,478:518);
wwv_mask_2S_borneo = wwv_mask_2S_borneo(91:811,478:518);
lon = lon(91:811,478:518);
lat = lat(91:811,478:518);
testmap(lon, lat, wwv_mask_2S_borneo);

srv = 'H:\Maurice_ENSO_Data\'

%% preamble and data preparation for warm water volume

                                % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
                                % location_tahiti        = [523, 427] %
                                % location_center_nino34 = [540, 498] %
                                % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

                                
                                
            % importante!       % wwv region 20S - 20N
                                % wwv_mask(117:840, 417:580) = 20S - 20N
                                % wwv_mask(117:840, 466:531) = 8S - 8N
                                
                                
                             
             
                             

restart = 'restart000';
string = [first '_restart000'];

tic;
% load in potential temperature
p1 = [srv first '/output000/'];
p2 = [srv first '/output001/'];
pc = [srv 'EXP0_control_run/'];


neutral_rho = getnc([p1 'ocean_wmass.nc'], 'neutral'); % temperature space array


% loading in warm water volume, i.e. volume of water above 20 degrees C in 
% the Pacific region, same as McGregor et al., 2014
% the Pacific region 5S - 5N and 120E - 80W
thetao_1 = squeeze(permute(getnc([p1 'ocean_snap.nc'], 'temp', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
         % [month, depth, lat, lon], [end], [stride]
thetao_2 = squeeze(permute(getnc([p2 'ocean_snap.nc'], 'temp', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
thetao = cat(4, thetao_1, thetao_2); 
clear thetao_1 thetao_2 thetao_3 thetao_4 thetao_5 thetao_6;

% so_1 = squeeze(permute(getnc([p1 'ocean.nc'], 'salt', ...
%          [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
% so_2 = squeeze(permute(getnc([p2 'ocean.nc'], 'salt', ...
%          [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
% so = cat(4, so_1, so_2); 
% clear so_1 so_2 so_3 so_4 so_5 so_6;

% loading in climatology
thetaoc = squeeze(permute(getnc([pc 'ocean_snap_clim.nc'], 'temp', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
thetaoc = cat(4, thetaoc, thetaoc);

% soc = squeeze(permute(getnc([pc 'ocean_clim.nc'], 'salt', ...
%          [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
% soc = cat(4, soc, soc);

cto = thetao; 
ctoc = thetaoc;

% now finding temperature within it which is higher than 20 degrees
mask = cto;
maskc = ctoc;
indices = find(mask < 20);  mask(indices) = 0; clear indices;
indices = find(mask >= 20); mask(indices) = 1; clear indices;

indices = find(maskc < 20);  maskc(indices) = 0; clear indices;
indices = find(maskc >= 20); maskc(indices) = 1; clear indices;


dzt_1 = squeeze(permute(getnc([p1 'ocean_snap.nc'], 'dzt', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
dzt_2 = squeeze(permute(getnc([p2 'ocean_snap.nc'], 'dzt', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
dzt = cat(4, dzt_1, dzt_2); 
clear dzt_1 dzt_2 dzt_3 dzt_4 dzt_5 dzt_6;

% loading in climatology
dztc = squeeze(permute(getnc([pc 'ocean_snap_clim.nc'], 'dzt', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
dztc = cat(4, dztc, dztc);

% mask is now my array containing all grid cells 5S-5N and 100E-60W which
% have temperatures higher than 20 degrees
volume = nan(721, 41, 50, 72);
volumec = nan(721, 41, 50, 72);
for z = 1:50
    for t = 1:24
        volume(:,:,z,t)  = mask(:,:,z,t) .* areacello(91:811,478:518) .* dzt(:,:,z,t) .* wwv_mask_2S_borneo;
        volumec(:,:,z,t) = maskc(:,:,z,t) .* areacello(91:811,478:518) .* dzt(:,:,z,t) .* wwv_mask_2S_borneo;
    end
end

size(volume)
% s = nansum(reshape(volume, [640*41*50 48]));
% t = nansum(reshape(volumec, [640*41*50 48]));

wwv = squeeze(nansum(squeeze(nansum(squeeze(nansum(volume-volumec, 1)), 1)), 1));
% 
plot(wwv); hold on;
toc;

%% calculating the derivative of wwv
% this is a bit tricky since I have monthly mean output and the derivative
% only gives 11 values
% the 12th value for dV/dt is from the restart file


p1_new = p1(1:end-10);
temp_res = squeeze(permute(getnc([p1(1:end-10) 'restart000/ocean_temp_salt.res.nc'], 'temp', ...
         [1,1,478,91], [1,-1,518,811], [1,1,1,1]), [4 3 2 1])); 

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
maskr = temp_res;
maskcr = temp_resc;
maskr(maskr < 20) = 0; 
maskcr(maskcr < 20) = 0; 
maskr(maskr >= 20) = 1; 
maskcr(maskcr >= 20) = 1; 


dzt = squeeze(permute(getnc([p1 'ocean_snap.nc'], 'dzt', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 

% loading in climatology
dztc = squeeze(permute(getnc([pc 'ocean_snap_clim.nc'], 'dzt', ...
         [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 


% mask is now my array containing all grid cells 8S-8N and 156E-95W which
% have temperatures higher than 20 degrees
for z = 1:50
        volume_res (:,:,z)  = maskr(:,:,z) .* areacello(91:811,478:518) .* dzt(:,:,z) .* wwv_mask_2S_borneo;
        volumec_res(:,:,z) = maskcr(:,:,z) .* areacello(91:811,478:518) .* dztc(:,:,z) .* wwv_mask_2S_borneo;
end


dV(1) = squeeze(nansum(squeeze(nansum(squeeze(nansum(volume_res-volumec_res, 1)), 1)), 2));
dV(1) = dV(1) / (30*24*60*60) / 1e6; % nochmals durch 1 Million um sch√∂n in Sverdrup zu plotten
plot(dV); hold on;





% size(temp_res)
% dV(1) = squeeze(nansum(squeeze(nansum(squeeze(nansum(temp_res, 1)), 1)), 2));


% taking derivative and converting to Sverdrup
% dV/dt = (wwv_i - wwv_i-1)/dt
delta_t = 30*24*60*60;
for i = 2:24
    dV(i) = (wwv(i) - wwv(i-1));
end
dV = dV / delta_t / 1e6; % nochmals durch 1 Million um sch√∂n in Sverdrup zu plotten
plot(dV); hold on;

% creating 5-month running mean as in Meinen and McPhaden 2000
dV = movmean(dV, 5);

% clear workspace
clear ans bathymetry cto ctoc dzt dztc maskc s so soc;
clear t test thetao thetaoc volcello z;
clear volume volume_res volumec volumec_res temp_res temp_resc;

%% calculating the horizontal transport and vertical transport as residual


ty_trans_nrho_1 = squeeze(nansum(permute(getnc([p1 'ocean.nc'], 'ty_trans_nrho', ...
         [1,46,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)); 
ty_trans_nrho_2 = squeeze(nansum(permute(getnc([p2 'ocean.nc'], 'ty_trans_nrho', ...
         [1,46,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)); 
ty_trans = cat(3, ty_trans_nrho_1, ty_trans_nrho_2);
size(ty_trans)
clear ty_trans_nrho_1 ty_trans_nrho_2 ty_trans_nrho_3 ty_trans_nrho_4 ty_trans_nrho_5 ty_trans_nrho_6;

% loading in climatology
ty_transc = squeeze(nansum(permute(getnc([pc 'ocean_clim.nc'], 'ty_trans_nrho', ...
         [1,46,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]), 3)); 
ty_transc = cat(3, ty_transc, ty_transc);
size(ty_trans)

% calculating the transport across the different transects

for i = 1:721           % transect to the north, +5 degrees
    north(i,:) = squeeze(ty_trans(i,41,:)) .* wwv_mask_2S_borneo(i,41);
    northc(i,:) = squeeze(ty_transc(i,41,:)) .* wwv_mask_2S_borneo(i,41);
end
for i = 50:127           % transect to the south, -2 degrees
    ITF(i,:) = squeeze(-ty_trans(i,12,:)) .* wwv_mask_2S_borneo(i,12);
    ITFc(i,:) = squeeze(-ty_transc(i,12,:)) .* wwv_mask_2S_borneo(i,12);
end
for i =  1:721           % transect to the south, -5 degrees
    south(i,:) = squeeze(-ty_trans(i,1,:)) .* wwv_mask_2S_borneo(i,1);
    southc(i,:) = squeeze(-ty_transc(i,1,:)) .* wwv_mask_2S_borneo(i,1);
end

north = nansum(north); northc = nansum(northc);
south = nansum(south); southc = nansum(southc);
ITF = nansum(ITF); ITFc = nansum(ITFc);

plot(northc); hold on; plot(southc); hold on; plot(ITFc); hold on;
boom
plot(north); hold on; plot(south); hold on; plot(ITF); hold on;



% creating anomalies and creating the 5-month running mean as in Meinen and
% McPhaden 2000
% multiplying by (-1) in order to achieve

%           ^
%           |           transport in
%           |
% ------------------------------------------ zero line
%           |
%           |           transport out
%           v   

northa = -(north - northc); clear north northc;
southa = -(south - southc); clear south southc;
ITFa = -(ITF - ITFc); clear ITF ITFc;
%
horizontal      = northa + southa;      % total meridional transport

ITFa = movmean(ITFa, 5);
horizontal = movmean(horizontal, 5);    % taking the 5-month mean as in Meinen and McPhaden 2000
                                        % wait a bit until taking the 5-month mean of the vertical flux
vertical        = dV - (horizontal + ITFa);  % total vertical tranport


plot(dV, 'linewidth', 2.5); hold on; line([1 24], [0 0], 'color', 'k');
plot(horizontal); hold on; plot(ITFa, 'color', 'g'); hold on; plot(vertical); hold on;

% safespot: 08. 03. 2018, 16:16

%% creating budget with water mass transformation framework 

% defining constants
rho_0 = 1035;       % reference air density [kg m^3]
cp = 3992.1;         % specific heat capacity of seawater [J kg^-1 K^-1]
dT = 0.5;           % the temperature difference of the binned temperature 


                % ~~~~~~~~~~~~~~~~~ GM ~~~~~~~~~~~~~~~~~ %                
cbt1 = squeeze(permute(getnc([p1 'ocean_wmass.nc'], 'temp_vdiffuse_diff_cbt_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
cbt2 = squeeze(permute(getnc([p2 'ocean_wmass.nc'], 'temp_vdiffuse_diff_cbt_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
cbt = cat(3, cbt1, cbt2); 
clear cbt1 cbt2 cbt3 cbt4 cbt5 cbt6;
% climatology
cbtc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'temp_vdiffuse_diff_cbt_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
cbtc = cat(3, cbtc, cbtc);

non1 = squeeze(permute(getnc([p1 'ocean_wmass.nc'], 'temp_nonlocal_KPP_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
non2 = squeeze(permute(getnc([p2 'ocean_wmass.nc'], 'temp_nonlocal_KPP_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
non = cat(3, non1, non2);
clear non1 non2 non3 non4 non5 non6;
% climatology
nonc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'temp_nonlocal_KPP_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
nonc = cat(3, nonc, nonc);

% units are in [W m^2]

         % ~~~~~~~~~~~ the calculation for GM is here ~~~~~~~~~~~ %        
for t = 1:24
    GM(:,:,t) = (cbt(:,:,t) + non(:,:,t)) .* wwv_mask_2S_borneo;
    GMc(:,:,t) = (cbtc(:,:,t) + nonc(:,:,t)) .* wwv_mask_2S_borneo;
% units are in = W m^-2
end
for t = 1:24
GM(:,:,t) = GM(:,:,t) .* (1 ./ (rho_0 .* cp .* dT)) .* areacello(91:811, 478:518) ./ 1e6;
GMc(:,:,t) = GMc(:,:,t) .* (1 ./ (rho_0 .* cp .* dT)) .* areacello(91:811, 478:518) ./ 1e6;
                          % durch 1e6 teilen f√ºr Einheit Sv [10^6 m^3 s^-1]
end
                          
mixing = squeeze(nansum(squeeze(nansum(GM,1)),1));
mixingc = squeeze(nansum(squeeze(nansum(GMc,1)),1));

% creating anomaly last
mixinga = mixing - mixingc; 
%           ^
%           |           transport in
%           |
% ------------------------------------------ zero line
%           |
%           |           transport out
%           v   

hold on; plot(vertical); hold on;
plot(mixinga); hold on; legend('vertical mixing'); 

                % ~~~~~~~~~~~~~~~~~ GF ~~~~~~~~~~~~~~~~~ %                

sbc1 = squeeze(permute(getnc([p1 'ocean_wmass.nc'], 'temp_vdiffuse_sbc_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
sbc2 = squeeze(permute(getnc([p2 'ocean_wmass.nc'], 'temp_vdiffuse_sbc_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1]));
sbc  = cat(3, sbc1, sbc2); 
clear sbc1 sbc2 sbc3 sbc4 sbc5 sbc6;
% climatology
sbcc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'temp_vdiffuse_sbc_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1]));
sbcc = cat(3, sbcc, sbcc);

sw1 = squeeze(permute(getnc([p1 'ocean_wmass.nc'], 'sw_heat_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
sw2 = squeeze(permute(getnc([p2 'ocean_wmass.nc'], 'sw_heat_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
sw  = cat(3, sw1, sw2); 
clear sw1 sw2 sw3 sw4 sw5 sw6;
% climatology
swc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'sw_heat_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
swc = cat(3, swc, swc);

frazil1 = squeeze(permute(getnc([p1 'ocean_wmass.nc'], 'frazil_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
frazil2 = squeeze(permute(getnc([p2 'ocean_wmass.nc'], 'frazil_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
frazil  = cat(3, frazil1, frazil2); 
clear frazil1 frazil2 frazil3 frazil4 frazil5 frazil6;
% climatology
frazilc = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'frazil_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
frazilc = cat(3, frazilc, frazilc);


eta1 = squeeze(permute(getnc([p1 'ocean_wmass.nc'], 'temp_eta_smooth_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1]));
eta2 = squeeze(permute(getnc([p2 'ocean_wmass.nc'], 'temp_eta_smooth_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
eta  = cat(3, eta1, eta2); 
clear eta1 eta2 eta3 eta4 eta5 eta6;
% climatology
etac = squeeze(permute(getnc([pc 'ocean_wmass_clim.nc'], 'temp_eta_smooth_on_nrho', ...
         [1,46,478,91], [12,46,518,811], [1,1,1,1]), [4 3 2 1])); 
etac = cat(3, etac, etac);

% units are in [W m^2]

         % ~~~~~~~~~~~ the calculation for GF is here ~~~~~~~~~~~ %        
for t = 1:24
    GF(:,:,t) = (sbc(:,:,t) + sw(:,:,t) + frazil(:,:,t) + eta(:,:,t)) .* wwv_mask_2S_borneo;
    GFc(:,:,t) = (sbcc(:,:,t) + swc(:,:,t) + frazilc(:,:,t) + etac(:,:,t)) .* wwv_mask_2S_borneo;
% units are in = W m^-2
end

for t = 1:24
GF(:,:,t) = GF(:,:,t) .* (1 ./ (rho_0 .* cp .* dT)) .* areacello(91:811, 478:518) ./ 1e6;
GFc(:,:,t) = GFc(:,:,t) .* (1 ./ (rho_0 .* cp .* dT)) .* areacello(91:811, 478:518) ./ 1e6;
                          % durch 1e6 teilen f√ºr Einheit Sv [10^6 m^3 s^-1]
end

forcing = squeeze(nansum(squeeze(nansum(GF,1)),1));
forcingc = squeeze(nansum(squeeze(nansum(GFc,1)),1));

% creating anomaly last
forcinga = forcing - forcingc; 
%           ^
%           |           transport in
%           |
% ------------------------------------------ zero line
%           |
%           |           transport out
%           v   

hold on; plot(vertical); hold on;
plot(forcinga); hold on; legend('surface forcing'); 
plot(mixinga); hold on; legend('vertical mixing'); 


                % ~~~~~~~~~~~~~~~~~ GJ ~~~~~~~~~~~~~~~~~ %                
mass1 = squeeze(nansum(permute(getnc([p1 'ocean.nc'], 'mass_pmepr_on_nrho', ...
         [1,46,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]),3)); 
mass2 = squeeze(nansum(permute(getnc([p2 'ocean.nc'], 'mass_pmepr_on_nrho', ...
         [1,46,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]),3)); 
mass = cat(3, mass1, mass2); 
clear mass1 mass2 mass3 mass4 mass5 mass6 massc;
% climatology
massc = squeeze(nansum(permute(getnc([pc 'ocean_clim.nc'], 'mass_pmepr_on_nrho', ...
         [1,46,478,91], [12,74,518,811], [1,1,1,1]), [4 3 2 1]),3)); 
massc = cat(3, massc, massc);


         % ~~~~~~~~~~~ the calculation for GJ is here ~~~~~~~~~~~ %        
for t = 1:24
    GJ(:,:,t) = (mass(:,:,t)) .* wwv_mask_2S_borneo .* 1e-3 ./ 1e6; % umrechnen von [kg s‚?ª1] zu [m^3 s^-1]
    GJc(:,:,t) = (massc(:,:,t)) .* wwv_mask_2S_borneo .* 1e-3 ./ 1e6; % umrechnen von [kg s‚?ª1] zu [m^3 s^-1]
                                 % ./ 1e6 Umrechnung auf Sv [10^6 m^3 s^-1]
% units are in = [m^3 s^-1]
end

mass = squeeze(nansum(squeeze(nansum(GJ,1)),1));
massc = squeeze(nansum(squeeze(nansum(GJc,1)),1));

% creating anomaly last
massa = mass - massc; 
%           ^
%           |           transport in
%           |
% ------------------------------------------ zero line
%           |
%           |           transport out
%           v   

mixinga = movmean(mixinga, 5);
forcinga = movmean(forcinga, 5);
massa = movmean(massa, 5);



implicit = vertical(1:24) - (mixinga(1:24) + forcinga(1:24) + massa(1:24));

plot(implicit, 'linewidth', 5);

% now since I have all the fluxes calculated I create the 5-month moving
% means of the vertical ones too

plot(forcinga); hold on; legend('surface forcing'); 
plot(mixinga); hold on; legend('vertical mixing'); 
plot(massa); hold on; legend('volume flux'); 
plot(implicit, 'linewidth', 5);

%% clean up  workspace before plotting
clear ans areacello bathymetry cto ctoc dzt dztc h7 h8 h9 lat lev_bnds lon;
clear mask maskc  so soc_3 soc_4 thetao thetaoc volcello; 
clear volume volumec z h3 h4 h5 h1 h2 net_sfcc;
clear soc i sstc tempc;
clear p1 p2 p3 p4 p5 pc tx_trans ty_trans wwv_region;
clear a neutral_rho;
clear cbt cbtc cp delta_t dT eta etac forcing forcingc frazil frazilc GF;
clear GFc GM GMc GJ GJc maskcr maskr mass massc mixing mixingc non nonc;
clear rho_0 sbc sbcc sw swc wwv_mask wwv_mask_2S ty_transc t;

%% filling up time series with zero numbers for good picture ratio plot
horizontal(25:72) = 0;
ITFa(25:72) = 0;
vertical(25:72) = 0;
dV(25:72) = 0;
mixinga(25:72) = 0;
forcinga(25:72) = 0;
massa(25:72) = 0;
implicit(25:72) = 0;

%% save variables for plotting routine later
if first == 'pnEXP1_composite_nino_windstress'
    horizontal_EXP1 = horizontal;
    ITFa_EXP1 = ITFa;
    dV_EXP1 = dV;
    massa_EXP1 = massa;
    forcinga_EXP1 = forcinga;
    mixinga_EXP1 = mixinga;
    implicit_EXP1 = implicit;
elseif first == 'pnEXP2_composite_nina_windstress'
    horizontal_EXP2 = horizontal;
    ITFa_EXP2 = ITFa;
    dV_EXP2 = dV;
    massa_EXP2 = massa;
    forcinga_EXP2 = forcinga;
    mixinga_EXP2 = mixinga;
    implicit_EXP2 = implicit;
end

clear horizontal ITFa dV massa forcinga mixinga implicit;

%% preparing data for table with percentages

if first == 'pnEXP1_composite_nino_windstress'
                 d = [7:21];            % discharge phase
elseif first == 'pnEXP2_composite_nina_windstress' 
                 d = [7:21];           % recharge phase
end

%                 14 months
% creating array with time info since WWV = G_F * dt, etc.
seconds = 86400;        % so many seconds per day
time = [30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31];
%       J   J   A    S   O   N   D   J   F   M   A   M   J   J   A
time = time .* seconds;         % so many seconds each day

table_discharge = nan(7,2);
table_discharge(1) = nansum(dV_EXP1(d) .* time);
table_discharge(2) = nansum(horizontal_EXP1(d) .* time);
table_discharge(3) = nansum(ITFa_EXP1(d) .* time);

table_discharge(4) = nansum(massa_EXP1(d) .* time);
table_discharge(5) = nansum(forcinga_EXP1(d) .* time);
table_discharge(6) = nansum(mixinga_EXP1(d) .* time);
table_discharge(7) = nansum(implicit_EXP1(d) .* time);

% table_recharge = nan(7,2);
% table_recharge(1) = nansum(dV(r) .* time);
% table_recharge(2) = nansum(horizontal(r) .* time);
% table_recharge(3) = nansum(ITFa(r) .* time);
% 
% table_recharge(4) = nansum(forcinga(r) .* time);
% table_recharge(5) = nansum(mixinga(r) .* time);
% table_recharge(6) = nansum(massa(r) .* time);
% table_recharge(7) = nansum(implicit(r) .* time);

for i = 1:7
    table_discharge(i,2) = table_discharge(i,1) ./ table_discharge(1,1) .*100;
%     table_recharge(i,2) = table_recharge(i,1) ./ table_recharge(1,1);
end

% rounding
table_discharge(:,3) = round(table_discharge(:,2),0);
table_discharge(:,1) = round(table_discharge(:,1),1);
table_discharge(8,3) = sum(table_discharge(2:7,3));

% convert to sverdrup
table_discharge(:,1) = table_discharge(:,1) * 1e6;
table_discharge(:,1) = table_discharge(:,1) / 1e14; % nun sind es nur noch kleine zahlen
table_discharge(:,1) = round(table_discharge(:,1),1);

open table_discharge

%% clean up  workspace before plotting again
clear ans areacello bathymetry cto ctoc dzt dztc h7 h8 h9 lat lev_bnds lon;
clear mask maskc  so soc_3 soc_4 thetao thetaoc volcello; 
clear volume volumec z h3 h4 h5 h1 h2 net_sfcc;
clear soc i sstc tempc;
clear p1 p2 p3 p4 p5 pc tx_trans ty_trans wwv_region;
clear a neutral_rho;
clear cbt cbtc cp delta_t dT eta etac forcing forcingc frazil frazilc GF;
clear GFc GM GMc GJ GJc maskcr maskr mass massc mixing mixingc non nonc;
clear rho_0 sbc sbcc sw swc wwv_mask wwv_mask_2S ty_transc t;
clear horizontal vertical ITFa dV massa forcinga mixinga implicit dV;
clear d northa p1_new restart seconds southa time;


%% ~~~~~~~~~~~ plotting  routine for complete figure with subplots ~~~~~~~~~~~ %% 

figure('units', 'pixels', 'position', [0 0 1920 1080]);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
subplot(2,2,1)
% % % shading of period around El Nino
xbars = [7 21];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-30 30 30 -30], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [24 72];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-12 25 25 -12], ...
    [1 1 1], 'edgecolor', 'none'); hold on;

h1 = plot(linspace(1,72,72), dV_EXP1, 'color', [0 0 0], 'linewidth', 2.5); hold on;
h2 = plot(linspace(1,72,72), horizontal_EXP1, 'color', RdYlBu(10,:), 'linewidth', 2); hold on;
% h3 = plot(linspace(1,72,72), vertical, 'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;
h1 = plot(linspace(1,72,72), dV_EXP1, 'color', [0 0 0], 'linewidth', 2.5); hold on;
% h11 = plot(linspace(1,72,72), movmean(test,5), 'color', [0 0 0], 'linewidth', 3, 'linestyle', '--'); hold on;
h6 = plot(linspace(1,72,72), massa_EXP1, 'color', RdYlBu(21,:), 'linewidth', 2); hold on;

% h5 = plot(linspace(1,72,72), forcinga_EXP1, 'color', RdYlBu(60,:), 'linewidth', 2); hold on;
% h4 = plot(linspace(1,72,72), mixinga_EXP1, 'color', antarctica(15,:), 'linewidth', 2); hold on;
% h7 = plot(linspace(1,72,72), implicit_EXP1, 'color', [187,0,187]/255, 'linewidth', 2); hold on;
h8 = plot(linspace(1,72,72), ITFa_EXP1, 'color', RdYlGn(60,:), 'linewidth', 2); hold on;

hXLabel = xlabel('');
hYLabel = ylabel(['[Sv]'],'interpreter','latex');
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
  'YTick'           , -20:10:30   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);

% text to indicate blue and red discharge/recharge phases
text([7.5 7.5],[-15.5 -15.5],'discharge phase', 'interpreter', 'latex', 'fontsize', 21, 'color', RdYlBu(1,:));
 
% set(gca, 'Xtick', [1, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60, 66, 72]);

h55 = legend([h1 h2 h8 h6], 'location', 'northwest', 'orientation', 'vertical', ...
    'Rate of change in WWV anomaly', ...
    '$\mathcal{T}_{\mathrm{5^{\circ}N + 5^{\circ}S}}$: Meridional transport', ...
    '$\mathcal{T}_{\mathrm{ITF}}$: Indonesian Throughflow', ...
    '$\mathcal{J}$: Surface volume', ...
    '$\mathcal{G_F}$: Surface forcing', ...
    '$\mathcal{G_M}$: Vertical mixing', ...
    '$\mathcal{G_I}$: Numerical mixing');
set(h55, 'interpreter', 'latex', 'fontsize', 15, 'textcolor', RdYlBu(60,:),'Edgecolor', [.83 .83 .83]);

ylim([-20 27]);
xlim([1 24]);
pbaspect([1.5 1 1]);                        % aspect ratios: x, y, z
clear hTitle hXLabel hYLabel h1 h2;
set(gca, 'Xtick', [1, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60, 66, 72]);

% draw arrows to denote El Nino and La Nina events
arrow([26, 5], [26, 20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(60,:));
arrow([26, -5], [26, -20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(1,:));
text([28 28],[5 5],'Recharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(60,:), 'rotation', 90);
text([28 28],[-5 -5],'Discharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(1,:), 'rotation', 90, 'horizontalalignment', 'right');

text(-3, 31, 'a) Adiabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);
text(11, 34, 'El Ni\~no', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
subplot(2,2,2)
xbars = [7 21];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-30 30 30 -30], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;

h1 = plot(linspace(1,72,72), dV_EXP2, 'color', [0 0 0], 'linewidth', 2.5); hold on;
h2 = plot(linspace(1,72,72), horizontal_EXP2, 'color', RdYlBu(10,:), 'linewidth', 2); hold on;
% h3 = plot(linspace(1,72,72), vertical, 'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;
h1 = plot(linspace(1,72,72), dV_EXP2, 'color', [0 0 0], 'linewidth', 2.5); hold on;
% h11 = plot(linspace(1,72,72), movmean(test,5), 'color', [0 0 0], 'linewidth', 3, 'linestyle', '--'); hold on;
h6 = plot(linspace(1,72,72), massa_EXP2, 'color', RdYlBu(21,:), 'linewidth', 2); hold on;

% h5 = plot(linspace(1,72,72), forcinga_EXP2, 'color', RdYlBu(60,:), 'linewidth', 2); hold on;
% h4 = plot(linspace(1,72,72), mixinga_EXP2, 'color', antarctica(15,:), 'linewidth', 2); hold on;
% h7 = plot(linspace(1,72,72), implicit_EXP2, 'color', [187,0,187]/255, 'linewidth', 2); hold on;
h8 = plot(linspace(1,72,72), ITFa_EXP2, 'color', RdYlGn(60,:), 'linewidth', 2); hold on;

hXLabel = xlabel('');
hYLabel = ylabel(['[Sv]'],'interpreter','latex');    
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
  'YTick'           , -20:10:30   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

set(gca, 'Xtick', [1, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60, 66, 72]);

% text to indicate blue and red discharge/recharge phases
text([7.5 7.5],[22.5 22.5],'recharge phase', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(60,:));

ylim([-20 27]);
xlim([1 24]);
pbaspect([1.5 1 1]);                        % aspect ratios: x, y, z
clear hTitle hXLabel hYLabel h1 h2;

% draw arrows to denote El Nino and La Nina events
arrow([26, 5], [26, 20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(60,:));
arrow([26, -5], [26, -20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(1,:));
text([28 28],[5 5],'Recharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(60,:), 'rotation', 90);
text([28 28],[-5 -5],'Discharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(1,:), 'rotation', 90, 'horizontalalignment', 'right');

text(-3, 31, 'b) Adiabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);
text(11, 34, 'La Ni\~na', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
subplot(2,2,3)
% % % shading of period around El Nino
xbars = [7 21];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-30 30 30 -30], ...
    [Reds(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;
xbars = [24 72];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-12 25 25 -12], ...
    [1 1 1], 'edgecolor', 'none'); hold on;

% h2 = plot(linspace(1,72,72), horizontal_EXP1, 'color', RdYlBu(10,:), 'linewidth', 2); hold on;
% h3 = plot(linspace(1,72,72), vertical, 'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;
% h1 = plot(linspace(1,72,72), dV_EXP1, 'color', [0 0 0], 'linewidth', 2.5); hold on;
% h11 = plot(linspace(1,72,72), movmean(test,5), 'color', [0 0 0], 'linewidth', 3, 'linestyle', '--'); hold on;
% h6 = plot(linspace(1,72,72), massa_EXP1, 'color', RdYlBu(21,:), 'linewidth', 2); hold on;

h1 = plot(linspace(1,72,72), dV_EXP1, 'color', [0 0 0], 'linewidth', 2.5); hold on;
h5 = plot(linspace(1,72,72), forcinga_EXP1, 'color', RdYlBu(60,:), 'linewidth', 2); hold on;
h4 = plot(linspace(1,72,72), mixinga_EXP1, 'color', antarctica(15,:), 'linewidth', 2); hold on;
h7 = plot(linspace(1,72,72), implicit_EXP1, 'color', [187,0,187]/255, 'linewidth', 2); hold on;
% h8 = plot(linspace(1,72,72), ITFa_EXP1, 'color', RdYlGn(60,:), 'linewidth', 2); hold on;

hXLabel = xlabel('Month');
hYLabel = ylabel(['[Sv]'],'interpreter','latex');
set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 23);
set(hYLabel, 'color', [0 0 0]);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 23);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -20:10:30   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);

set(gca, 'Xtick', [1, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60, 66, 72]);

h55 = legend([h5 h4 h7], 'location', 'northwest', 'orientation', 'vertical', ...
    '$\mathcal{G_F}$: Surface forcing', ...
    '$\mathcal{G_M}$: Vertical mixing', ...
    '$\mathcal{G_I}$: Numerical mixing');
set(h55, 'interpreter', 'latex', 'fontsize', 15, 'textcolor', RdYlBu(60,:), 'Edgecolor', [.83 .83 .83]);

ylim([-20 27]);
xlim([1 24]);
pbaspect([1.5 1 1]);                        % aspect ratios: x, y, z
clear hTitle hXLabel hYLabel h1 h2;

% draw arrows to denote El Nino and La Nina events
arrow([26, 5], [26, 20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(60,:));
arrow([26, -5], [26, -20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(1,:));
text([28 28],[5 5],'Recharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(60,:), 'rotation', 90);
text([28 28],[-5 -5],'Discharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(1,:), 'rotation', 90, 'horizontalalignment', 'right');

text(-3, 31, 'c) Diabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
subplot(2,2,4)
xbars = [7 21];
h44 = patch([xbars(1) xbars(1), xbars(2), xbars(2)], [-30 30 30 -30], ...
    [Blues(10,:)], 'edgecolor', 'none', 'facealpha', 0.3); hold on;

% h2 = plot(linspace(1,72,72), horizontal_EXP2, 'color', RdYlBu(10,:), 'linewidth', 2); hold on;
% h3 = plot(linspace(1,72,72), vertical, 'color', RdYlBu(10,:), 'linewidth', 2.5); hold on;
h1 = plot(linspace(1,72,72), dV_EXP2, 'color', [0 0 0], 'linewidth', 2.5); hold on;
% h11 = plot(linspace(1,72,72), movmean(test,5), 'color', [0 0 0], 'linewidth', 3, 'linestyle', '--'); hold on;
% h6 = plot(linspace(1,72,72), massa_EXP2, 'color', RdYlBu(21,:), 'linewidth', 2); hold on;

h5 = plot(linspace(1,72,72), forcinga_EXP2, 'color', RdYlBu(60,:), 'linewidth', 2); hold on;
h4 = plot(linspace(1,72,72), mixinga_EXP2, 'color', antarctica(15,:), 'linewidth', 2); hold on;
h7 = plot(linspace(1,72,72), implicit_EXP2, 'color', [187,0,187]/255, 'linewidth', 2); hold on;
% h8 = plot(linspace(1,72,72), ITFa_EXP2, 'color', RdYlGn(60,:), 'linewidth', 2); hold on;

hXLabel = xlabel('Month');
hYLabel = ylabel(['[Sv]'],'interpreter','latex');    
set([hXLabel, hYLabel], 'interpreter', 'latex', 'Fontsize', 23);
set(hYLabel, 'color', [0 0 0]);
set(gca,'Fontname', 'Times New Roman', 'FontSize', 23);

set(gca, ...
  'Box'             , 'off'         , ...
  'TickDir'         , 'out'         , ...
  'TickLength'      , [.01 .01]     , ...
  'XMinorTick'      , 'off'         , ...
  'YMinorTick'      , 'off'         , ...
  'YGrid'           , 'on'          , ...
  'YTick'           , -20:10:30   , ...     % define grid, every 1000 m a grid line
  'XColor'          , RdYlBu(60,:)  , ...
  'YColor'          , RdYlBu(60,:)  , ...
  'ticklabelinterpreter', 'latex'   , ...
  'LineWidth'       , 1.25);
  % 'GridLineStyle'   , '--'          , ...

set(gca, 'Xtick', [1, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60, 66, 72]);

ylim([-20 27]);
xlim([1 24]);
pbaspect([1.5 1 1]);                        % aspect ratios: x, y, z
clear hTitle hXLabel hYLabel h1 h2;

% draw arrows to denote El Nino and La Nina events
arrow([26, 5], [26, 20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(60,:));
arrow([26, -5], [26, -20], 'BaseAngle', 80, 'length', 15, ...
    'width', 2.5, 'ends', [1], 'color', RdYlBu(1,:));
text([28 28],[5 5],'Recharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(60,:), 'rotation', 90);
text([28 28],[-5 -5],'Discharge', 'interpreter', 'latex', ...
    'fontsize', 21, 'color', RdYlBu(1,:), 'rotation', 90, 'horizontalalignment', 'right');

text(-3, 31, 'd) Diabatic fluxes', 'interpreter', 'latex', 'Fontsize', 25, ...
    'color', [0 0 0]);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = 'C:/Users/Maurice Huguenin/Desktop/';
%savefig([directory 'pnEXP1_and_pnEXP2_' f3 '_vertical_new_combined'])
print('-dpng','-r300', [directory 'pnEXP1_and_pnEXP2_' f3 '_vertical_new_combined_test']);

% 





