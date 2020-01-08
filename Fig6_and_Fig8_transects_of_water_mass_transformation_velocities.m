%% This script extracts an equatorial slice in the Pacific of data
%% from MOM025 runs


for s = 1
    if s == 1
           period = 'EN'
    baseD = 'G:/Maurice_ENSO_Data/pnEXP1_composite_nino_windstress/'; %Data Directory.
    outD = 'G:/Maurice_ENSO_Data/equatorial_slices_pnEXP1/'; %Data Directory.
    elseif s == 2
           period = 'LN'
    baseD = 'G:/Maurice_ENSO_Data/pnEXP2_composite_nina_windstress/'; %Data Directory.
    outD = 'G:/Maurice_ENSO_Data/equatorial_slices_pnEXP2/'; %Data Directory.
    elseif s == 3
        period = '1996-2001'
    baseD = 'G:/Maurice_ENSO_Data/pEXP9601_real_windstress/'; %Data Directory.
    outD = 'FG/Maurice_ENSO_Data/equatorial_slices_1996-2001/'; %Data Directory.
    end
    
    

model = 'MOM025';

for output= 0:5
tic

    base = [baseD sprintf('output%03d/',output)];
hname = [base 'ocean_heat.nc'];
% if (strfind(baseD,'01'))
%     fname = [base 'ocean_month.nc'];
     m3name = [base 'ocean.nc'];
% else
    fname = [base 'ocean.nc'];
% end
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
wname = [base 'ocean_wmass.nc'];
         
lon = ncread(hname,'geolon_t');lat = ncread(hname,'geolat_t');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonu = ncread(hname,'geolon_c');latu = ncread(hname,'geolat_c');

z = ncread(hname,'st_ocean');zL = length(z);

time = ncread(hname,'time');
tL = length(time);

Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

%%% Get equatorial slices of variables:
latsl = 0;
[tmp ltind] = min(abs(lat(1,:)-latsl));
[tmp ln1] = min(abs(lon(:,ltind)+240));
[tmp ln2] = min(abs(lon(:,ltind)+70));


temp = squeeze(ncread(fname,'temp',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
u = squeeze(ncread(fname,'u',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
v = squeeze(ncread(fname,'v',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
kappa = squeeze(ncread(fname,'diff_cbt_t',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
taux = squeeze(ncread(fname,'tau_x',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
tauy = squeeze(ncread(fname,'tau_y',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
mld = squeeze(ncread(fname,'mld',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
vdif = squeeze(ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
vnlc = squeeze(ncread(wname,'temp_nonlocal_KPP_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
pmer = squeeze(ncread(wname,'sfc_hflux_pme_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'temp_rivermix_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
sufc = squeeze(ncread(wname,'temp_vdiffuse_sbc_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'frazil_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'temp_eta_smooth_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'sw_heat_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
swrd = squeeze(ncread(wname,'sw_heat_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));

[Xt,Zt] = ndgrid(lon(ln1:ln2,ltind),z);
[Xu,Zu] = ndgrid(lonu(ln1:ln2,ltind),z);

name = [outD model sprintf('_output%03d',output) '_varsat_Eq.mat']
%save(name,'Xt','Zt','Xu','Zu','temp','u','v','kappa','taux','tauy','mld', ...
%     'vdif','vnlc','pmer','sufc','swrd', 'T');
toc;
end
end


%% Plot Equatorial Slices:


for s = 3
if s == 1
           WMT = 'vertical_mixing'
           period = 'EN'
           base = 'G:\Maurice_ENSO_Data\equatorial_slices_pnEXP1/';
    elseif s == 2
           WMT = 'vertical_mixing'
           period = 'LN'
           base = 'G:\Maurice_ENSO_Data\equatorial_slices_pnEXP2/';
    elseif s == 3
           WMT = 'surface_forcing'
           period = 'EN'
           base = 'G:\Maurice_ENSO_Data\equatorial_slices_pnEXP1/';
    elseif s == 4
           WMT = 'surface_forcing'
           period = 'LN'
           base = 'G:\Maurice_ENSO_Data\equatorial_slices_pnEXP2/';
    elseif s == 5
           WMT = 'vertical_mixing'
           period = 'CM'
           base = 'G:\Maurice_ENSO_Data\equatorial_slices_EXP0_control_run/';      
    elseif s == 6
           WMT = 'surface_forcing'
           period = 'CF'
           base = 'G:\Maurice_ENSO_Data\equatorial_slices_EXP0_control_run/'; 
    end
   
if s == 1 | s == 3
    first = 'pnEXP1_composite_nino_windstress';
    % load in potential temperature
    p1 = ['G:/Maurice_ENSO_Data/' first '/output000/'];

    % loading in warm water volume, i.e. volume of water above 20 degrees C in 
    % the Pacific region, same as McGregor et al., 2014
    % the Pacific region 5S - 5N and 120E - 80W
    thetao = squeeze(permute(getnc([p1 'ocean_snap.nc'], 'temp', ...
             [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
             % [month, depth, lat, lon], [end], [stride]
    thetao = squeeze(nanmean(thetao,2));
elseif s == 2 | s == 4
    first = 'pnEXP2_composite_nina_windstress';
        % load in potential temperature
    p1 = ['G:/Maurice_ENSO_Data/' first '/output000/'];

    % loading in warm water volume, i.e. volume of water above 20 degrees C in 
    % the Pacific region, same as McGregor et al., 2014
    % the Pacific region 5S - 5N and 120E - 80W
    thetao = squeeze(permute(getnc([p1 'ocean_snap.nc'], 'temp', ...
             [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
             % [month, depth, lat, lon], [end], [stride]
end

if s == 5 | s == 6
    pc = 'F:/Maurice_ENSO_Data/EXP0_control_run/';
    % loading in climatology
    thetaoc = squeeze(permute(getnc([pc 'ocean_snap_clim.nc'], 'temp', ...
             [1,1,478,91], [12,-1,518,811], [1,1,1,1]), [4 3 2 1])); 
    thetaoc = cat(4, thetaoc, thetaoc, thetaoc, thetaoc, thetaoc, thetaoc);
    thetaoc = squeeze(nanmean(thetaoc,2));
end

% Load Base Variables:
model = 'MOM025';
outputs = [0];

ndays = [31 28 31 30 31 30 31 31 30 31 30 31];

% Load Variable and calculate mean:
    load([base model sprintf(['_output%03d_varsat_Eq.mat'],outputs(1))]);
    
vars = {'temp','u','v','kappa','taux','tauy','mld','vdif','vnlc','pmer','sufc','swrd'};

for i=1:length(vars)
    eval([vars{i} 'a = ' vars{i} ';']);
end
for i=2:length(outputs)
    load([base model sprintf(['_output%03d_varsat_Eq.mat'],outputs(1))]);
    for i=1:length(vars)
        eval([vars{i} 'a = ' vars{i} 'a + ' vars{i} ';']);
    end
end
for i=1:length(vars)
    eval([vars{i} ' = ' vars{i} 'a/length(outputs);']);
    eval(['clear ' vars{i} 'a;']);
end
[xL,zL,tL] = size(temp);
TL = length(T);


% Depth of isotherms:
Zi = zeros(xL,TL,tL);
for ti=1:1
    for xi=500
        tvec = squeeze(temp(xi,:,ti));
        zvec = -Zt(xi,:);
        tvec(isnan(tvec)) = -1000;
        tvec = tvec - 0.01*(1:zL);
        Zi(xi,:,ti) = interp1(tvec,zvec,T,'linear');
        ind = find(~isnan(Zi(xi,:,ti)),1,'last');
        Zi(xi,(ind+1):end,ti) = max(zvec);%linspace(Zi(yi,ind;
    end
end
Xi = repmat(Xt(:,1),[1 TL]);

if WMT == 'surface_forcing'
var = cumsum(vdif,2,'reverse'); % Vertical Mixing Flux
clim = [-400 0];
sp = 10;
doWMT = 1;
end
Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3
var = (vdif+vnlc)/rho0/Cp*86400; % Vertical Mixing Transformation
if WMT == 'surface_forcing'
var = (pmer+sufc)/rho0/Cp*86400; % Surface Forcing Transformation
end
clim = [-1 1]*1e-5*86400; % FOR WMT
sp = 0.1*1e-5*86400;
doWMT = 1;

if s == 5 | s == 6
months = {[1:12],[3],[7],[11], [1:12], [1:12]}; % here he specifies the different subplots
                                % i.e. for annual mean, month 3, 7 and 11
else
    months = {[9:11],[1:12],[9:11],[1:12], [1:12], [1:12]}; % here he specifies the different subplots
                                % i.e. for annual mean, month 3, 7 and 11
end
monthsu01 = {[1:4],[3],[7],[11]};
labels = {'Annual','March','July','November'};

%Colormap:
cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% npts = length(cpts)
npts = 20;


if (doWMT)
    cmap = flipud(cbrewer('div', 'RdBu', 20, 'PCHIP'));
else
    cmap = flipud(cbrewer('div', 'RdBu', 20, 'PCHIP'));
end


RdBu = flipud(cbrewer('div', 'RdBu', 20, 'PCHIP'));
RdYlBu = cbrewer('div', 'RdBu', 60, 'PCHIP');
figure;
set(gcf,'Position',[1          36        1920         970]);
% for i=1:length(months)
for i=1:length(months)
subplot(3,2,i);
contourf(Xi,nanmonmean(Zi(:,:,months{i}),3,ndays(months{i})),nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
hold on;
Ztemp = Zt;
Ztemp(596:597,5:50) = nan;

% srv = 'F:/Maurice_ENSO_Data/'
% load([srv 'workspace_EXP1_analysis_composite_ninos_rev3.mat'], 'lon', 'lev_bnds');
% lon = lon(91:811,478:518);
% [Xt,Zt] = ndgrid(lon(:,1),lev_bnds);
if s == 5 | s == 6
    [c,h] = contour(Xt,-Ztemp,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[0:1:35],'-k');
    clabel(c,h,[0:3:35], 'interpreter', 'latex', 'fontsize', 12);
    [c1,h1] = contour(Xt,-Ztemp,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[20 20], 'color', [1 1 1],'linewidth',4.5); hold on;
    [c2,h2] = contour(Xt,-Ztemp,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[20 20],'-k','linewidth',3);
else
    [c,h] = contour(Xt,-Ztemp,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[0:1:35],'-k');
    clabel(c,h,[0:3:35], 'interpreter', 'latex', 'fontsize', 12);
    [c1,h1] = contour(Xt,-Ztemp,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[20 20], 'color', [1 1 1],'linewidth',4.5); hold on;
    [c2,h2] = contour(Xt,-Ztemp,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[20 20],'-k','linewidth',3);
end  
% plotting contour of mixed layer depth (mld)
if (strcmp(model,'MOM01'))
    mnu = monthsu01{i};
else
    mnu = months{i};
end
% ucol = [0.8706    0.4902         0];
% [c,h] = contour(Xu,-Zu,mean(u(:,:,mnu),3),[-2:0.2:-0.2],'--', ...
%                 'color',ucol);
% [c,h] = contour(Xu,-Zu,mean(u(:,:,mnu),3),[0.2:0.2:2],'-', ...
%                 'color',ucol);
plot(Xu(:,1),-monmean(mld(:,months{i}),2,ndays(months{i})),'--','color',RdYlBu(50,:),'linewidth',3);
% $$$ clabel(c,h,'color','w');
ylim([-300 0]);
xlim([-220 -80]);
cb = colorbar('color', [.1922 .2118 .5843], 'box', 'on', ...
    'tickdirection', 'in', 'ticklabelinterpreter', 'latex', 'Fontsize', 16);
limits = [-1.5 1.5];
set(cb, 'YTick', linspace(limits(1), limits(2),5));

if (doWMT)
    ylabel(cb,'[m day$^{-1}$]','interpreter', 'latex', 'FontSize', 16);
else
    ylabel(cb,'Wm$^{-2}$','interpreter', 'latex', 'FontSize', 16);
end
xlabel('Longitude','interpreter', 'latex', 'FontSize', 16);
ylabel('Depth [m]','interpreter', 'latex', 'FontSize', 16);
caxis(clim);
% text(-218,-288,labels{i},'Backgroundcolor','w','FontSize',25);

% LabelAxes(gca,i,20,0.008,0.95);
colormap(cmap);
set(gca, 'clim', [limits(1) limits(2)]); hold on;

set(gca,'XTickLabel',{'$140^{\circ}$W','$160^{\circ}$E','$180^{\circ}$W','$160^{\circ}$E','$140^{\circ}$W', ...
                         '$120^{\circ}$W', '$100^{\circ}$', '$80^{\circ}$W'}, ...
                         'ticklabelinterpreter', 'latex', 'FontSize', 16);

set(gca, ...
  'XColor'          , [.1922 .2118 .5843]  , ...
  'YColor'          , [.1922 .2118 .5843]);
set(gca,'layer','top'); % bring grid to the top layer
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 14);


if i == 1
    if s == 5
      text(-240, 32, 'a) $\mathcal{G_M}$: Annual Mean','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    elseif s == 6
      text(-240, 32, 'a) $\mathcal{G_F}$: Annual Mean','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    elseif s == 1
      text(-240, 32, 'c) $\mathcal{G_M}$: El Ni\~no months 9--11','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    elseif s == 2
      text(-240, 32, 'e) $\mathcal{G_M}$: La Ni\~na months 9--11','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    elseif s == 3
      text(-240, 32, 'd) $\mathcal{G_F}$: El Ni\~no months 9--11','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    elseif s == 4
      text(-240, 32, 'f) $\mathcal{G_F}$: La Ni\~na months 9--11','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    end
elseif i == 2 
    if s == 5
      text(-240, 32, 'b) $\mathcal{G_M}$: March','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    elseif s == 6
      text(-240, 32, 'b) $\mathcal{G_F}$: March','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    else
      text(-240, 32, 'b) $\mathcal{G_M}$: Climatology','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    end
    if s == 1 | s == 2 | s == 3 | s == 4
    h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
    end
elseif i == 3 
    if s == 5
      text(-240, 32, 'c) $\mathcal{G_M}$: July','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    elseif s == 6
      text(-240, 32, 'c) $\mathcal{G_F}$: July','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    else
      text(-240, 32, 'c) $\mathcal{G_M}$: Climatology','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    end
    if s == 1 | s == 2 | s == 3 | s == 4
    h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
    end
elseif i == 4
    if s == 5
      text(-240, 32, 'd) $\mathcal{G_M}$: November','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    elseif s == 6
      text(-240, 32, 'd) $\mathcal{G_F}$: November','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    else
      text(-240, 32, 'd) $\mathcal{G_M}$: Climatology','interpreter','latex','Fontsize',23, ...
    'color',[0 0 0]);
    end
    if s == 1 | s == 2 | s == 3 | s == 4
    h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
    end
elseif i == 5 
      text(-240, 31, 'e) $\mathcal{G_M}$: La Ni\~na months 9--11', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
elseif i == 6 
      text(-240, 31, 'f) $\mathcal{G_F}$: La Ni\~na months 9--11', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
end
end



set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = 'C:/Users/Maurice Huguenin/Desktop/';
f3 = 'EXP1_and_EXP2_equatorial_slices';
print('-dpng','-r100', [directory f3 '_' WMT '_' period]);
boom
end










