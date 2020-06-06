%% This script extracts an equatorial slice in the Pacific of data


%% Plot Equatorial Slices:

for s = 5:6
    if s == 1
           WMT = 'vertical_mixing'
           period = 'EN'
           base = 'H:\Maurice_ENSO_Data\equatorial_slices_EXP1_and_EXP2/';
    elseif s == 2
           WMT = 'vertical_mixing'
           period = 'LN'
           base = 'H:\Maurice_ENSO_Data\equatorial_slices_pnEXP2/';
    elseif s == 3
           WMT = 'surface_forcing'
           period = 'EN'
           base = 'H:\Maurice_ENSO_Data\equatorial_slices_EXP1_and_EXP2/';
    elseif s == 4
           WMT = 'surface_forcing'
           period = 'LN'
           base = 'H:\Maurice_ENSO_Data\equatorial_slices_EXP2_and_EXP1/';
    elseif s == 5
           WMT = 'vertical_mixing'
           period = 'CM'
           base = 'H:\Maurice_ENSO_Data\equatorial_slices_EXP0_control_run/'; 
    elseif s == 6
           WMT = 'surface_forcing'
           period = 'CM'
           base = 'H:\Maurice_ENSO_Data\equatorial_slices_EXP0_control_run/'; 
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
for ti=1:tL
    for xi=1:xL
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

months = {[1:12],[1:12],[9:11],[9:11], [1:12], [1:12]}; % here he specifies the different subplots
                                % i.e. for annual mean, month 3, 7 and 11
monthsu01 = {[1:4],[1],[3],[4]};
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
[c,h] = contour(Xt,-Ztemp,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[0:1:35],'-k');
% plot 20 degree isotherm bold
[c1,h1] = contour(Xt,-Zt,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[20 20], 'color', [1 1 1],'linewidth',4.5); hold on;
[c2,h2] = contour(Xt,-Zt,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[20 20],'-k','linewidth',3);


clabel(c,h,[12,15,23,25,27], 'interpreter', 'latex', 'fontsize', 20, 'color', 'g'); hold on;
clabel(c,h,[12,15,23,25,27], 'interpreter', 'latex', 'fontsize', 16, 'color', 'k');


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

% plot the mixed layer depth as a dashed line
plot(Xu(:,1),-monmean(mld(:,months{i}),2,ndays(months{i})),'--','color',RdYlBu(50,:),'linewidth',3);
% clabel(c,h,'color','w');
ylim([-300 0]);
xlim([-220 -80]);
cb = colorbar('color', [.1922 .2118 .5843], 'box', 'on', ...
    'tickdirection', 'in', 'ticklabelinterpreter', 'latex', 'Fontsize', 17);
limits = [-1.5 1.5];
set(cb, 'YTick', linspace(limits(1), limits(2),5));

if (doWMT)
    ylabel(cb,'(m day$^{-1}$)','interpreter', 'latex', 'FontSize', 17);
else
    ylabel(cb,'Wm$^{-2}$','interpreter', 'latex', 'FontSize', 17);
end
xlabel('Longitude','interpreter', 'latex', 'FontSize', 17);
ylabel('Depth [m]','interpreter', 'latex', 'FontSize', 17);
caxis(clim);
% text(-218,-288,labels{i},'Backgroundcolor','w','FontSize',25);

% LabelAxes(gca,i,20,0.008,0.95);
colormap(cmap);
set(gca, 'clim', [limits(1) limits(2)]); hold on;

set(gca,'XTickLabel',{'','$160^{\circ}$E','','$160^{\circ}$W','', ...
                         '$120^{\circ}$W', '', '$80^{\circ}$W'}, ...
                         'ticklabelinterpreter', 'latex', 'FontSize', 16);

set(gca, ...
  'XColor'          , [.1922 .2118 .5843]  , ...
  'YColor'          , [.1922 .2118 .5843]);
set(gca,'layer','top'); % bring grid to the top layer
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 17);


if i == 1                 
      text(-240, 32, 'a) Annual mean', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
elseif i == 2 
    text(-240, 32, 'b) -', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
elseif i == 3 
      text(-240, 31, 'c) SON', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
elseif i == 4 
      text(-240, 31, 'd) -', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
elseif i == 5 
      text(-240, 31, 'e) $\mathcal{G_M}$: Vertical mixing', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
elseif i == 6 
      text(-240, 31, '$\mathcal{G_F}$: Surface forcing', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
end
end



set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = 'C:/Users/Maurice Huguenin/Desktop/';
f3 = 'EXP1_and_EXP2_equatorial_slices';
print('-dpng','-r500', [directory f3 '_' WMT '_' period]);
boom
end










