
%% Plot Equatorial Slices:
tic;
WMT = 'vertical_mixing'
period = 'LN'
base = 'H:/Maurice_ENSO_Data/equatorial_slices_EXP1979-2016/';

% Load Base Variables:
model = 'ACCESS-OM2';
outputs = [29];

ndays = [31 28 31 30 31 30 31 31 30 31 30 31]; % number of days in each month

% Load Variable and calculate mean:
    load([base model sprintf(['_output%03d_varsat_Eq.mat'],outputs(1))]);
    
vars = {'temp','u','v','kappa','taux','tauy','mld','vdif','vnlc','pmer','sufc','swrd', 'T'};
temp = temp - 273.15; % convert temperature from Kelvin back to Celsius
                      % -> in MOM, temperatures were given in Kelvin
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

months = {[7],[7],[12],[12],[5],[5]}; % here he specifies the different subplots
                                % i.e. for annual mean, month 3, 7 and 11
monthsu01 = {[1:12],[3],[6],[11],[11],[11]};
% labels = {'Annual','March','July','November', 'test', 'test'};

%Colormap:
cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% npts = length(cpts)
npts = 20;


if (doWMT)
    cmap = flipud(cbrewer('div', 'RdBu', 20, 'PCHIP'));
else
    cmap = flipud(cbrewer('div', 'RdBu', 20, 'PCHIP'));
end

% my custom colourbars
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
    ylabel(cb,'[m day$^{-1}$]','interpreter', 'latex', 'FontSize', 17);
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
      text(-240, 32, 'a) July 1988', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
elseif i == 2 
h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
elseif i == 3 
      text(-240, 31, 'b) December 1988', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
elseif i == 4 
h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
elseif i == 5
text(-240, 31, 'c) May 1989', 'interpreter', 'latex', 'Fontsize', 23, ...
    'color', [0 0 0]);
elseif i == 6 
h5 = patch([200 -500 -500 200], [0, 0, -300, -300], ...
    [0.8875 0.8875 0.8875], 'edgecolor', 'none', 'facealpha', 1); hold on;
end
end



% finished fancy plot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = 'C:/Users/Maurice Huguenin/Desktop/';
print('-dpng','-r500', [directory 'Equatorial_transect_vertical_mixing_LN_' ...
    'output0' num2str(outputs)]);
toc;






