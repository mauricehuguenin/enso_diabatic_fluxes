%% Creating symmetric PC time series for EXP1: Composite of three strongest El Nino events
% hmaurice: 19.10.2017, 15:47

% Dimensions:
%            LAT  = 192
%            LON  = 94
%            TIME = 2920
      
            mirrored = 0;
%            mirrored = 1;


%% [works] preamble
load('workspace_regression_patterns_PC1_equal_nino34_rev2', 'PC');
load('C:\Users\Maurice Huguenin\Desktop\workspace_EXP1_and_EXP2_polynomial_PC_composites_symmetric');
PC(:,457) = 0;
RdBu_short = cbrewer('div', 'RdBu', 21, 'PCHIP');
RdYlBu = cbrewer('div', 'RdYlBu', 60, 'PCHIP');
RdYlBu_short = cbrewer('div', 'RdYlBu', 21, 'PCHIP');
time = 1:24;

f1 = 'pnEXP1_and_pnEXP2_idealized_symmetric_timeseries';
srv = '/srv/ccrc/data15/z5180028/MSC_thesis_workspaces/'

% load('workspace_EXP1_and_EXP2_polynomial_PC_composites_symmetric', 'pnnino_comp_highres');


if mirrored == 0
    pnnino_comp_highres = pnnino_comp_highres;
elseif mirrored == 1
    pnnino_comp_highres(:,1461:2920) = fliplr(pnnino_comp_highres(:,1:1460));

pnnino_comp_highres(1,1409:1511) = 1.746;
pnnino_comp_highres(1,1445:1485) = 1.748;
pnnino_comp_highres(2,1460:2920) = -pnnino_comp_highres(2,1460:2920);
pnnino_comp_highres(2,1461:2920) = -fliplr(pnnino_comp_highres(2,1:1460));


% creating a ~~~ 9th order polynomial ~~~ for the first PC of La Niña
mu = 1460.5
sigma = 843.08

% coefficients:
  p1 = -5.5362e-17
  p2 = -0.31692
  p3 = 8.0854e-17
  p4 = 1.0194
  p5 = -3.5378e-17
  
for x = 1:2920
    z = (x-mu)/sigma
    pnnino_comp_highres(2,x) = p1*z^4 + p2*z^3 + p3*z^2 + p4*z + p5; 
end

% slightly adjusting the start of my time series so it goes from zero
% amplitude
for i = 1:30
pnnino_comp_highres(2,i) = i* pnnino_comp_highres(2,30)/30;
end

pnnino_comp_highres(2,1461:2920) = -fliplr(pnnino_comp_highres(2,1:1460));

end

plot(pnnino_comp_highres(1,:)); hold on;
plot(pnnino_comp_highres(2,:)); hold on;

if mirrored == 1
pmnino_comp_highres = pnnino_comp_highres;
end





%% plotting routine
% load('workspace_EXP1_and_EXP2_polynomial_PC_composites_symmetric', 'nino_comp_highres', 'pnnino_comp_highres_spin');

figure('units', 'pixels', 'position', [0 0 1000 500]);

time = [1:24]
% EL NINO CASE
% original three PC1 time series
h1 = plot(time, PC(1,37:60), 'color', [0 0 0], 'linewidth', 1.); hold on;
h2 = plot(time, PC(1,216:239), 'color', [0.3 0.7 0.6], 'linewidth', 1.); hold on;
h3 = plot(time, PC(1,432:455), 'color', [0 0 0], 'linewidth', 1.); hold on;
% original three PC2 time series
h5 = plot(time, PC(2,37:60), 'color', [RdYlBu(1,1:3)], 'linewidth', 1); hold on;
h6 = plot(time, PC(2,216:239), 'color', [RdYlBu(1,1:3)], 'linewidth', 1); hold on;
h7 = plot(time, PC(2,432:455), 'color', [RdYlBu(1,1:3)], 'linewidth', 1); hold on;

vline = line([12 12], [-1 5]); % draw boundary line between El Niño and La Niña simulation
set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);

vline = line([24 24], [-5 5]); % draw boundary line between El Niño and La Niña simulation
set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);

vline = line([36 36], [-1 5]); % draw boundary line between El Niño and La Niña simulation
set(vline, 'color', [0.8784 0.8824 0.9373], 'linewidth', 1.25);

% idealized timeseries evolution
% h5 = plot(linspace(1, 24, 2920), nino_comp_highres(1,:), 'color', [0 0 0], ...
%     'linewidth', 3); hold on; 
% h9 = plot(linspace(1, 24, 2920), nino_comp_highres(2,:), 'color', [RdYlBu(1,1:3)], ...
%     'linewidth', 3); hold on;

a = nan(2,5840); 
a(1,1:2920) = pnnino_comp_highres(1,:); a(1,2921:5840) = pnnino_comp_highres_spin(1,:);
a(2,1:2920) = pnnino_comp_highres(2,:); a(2,2921:5840) = pnnino_comp_highres_spin(2,:);
h4 = plot(linspace(1, 24, 2920), a(1,1:2920), 'color', [0 0 0], ...
    'linewidth', 3); hold on;   
h8 = plot(linspace(1, 24, 2920), a(2,1:2920), 'color', [RdYlBu(1,1:3)], ...
    'linewidth', 3); hold on;   
h4 = plot(linspace(24, 48, 2920), a(1,2921:5840), 'color', [0 0 0], ...
    'linewidth', 3); hold on;   
h8 = plot(linspace(24, 48, 2920), a(2,2921:5840), 'color', [RdYlBu(1,1:3)], ...
    'linewidth', 3); hold on;   



    hYLabel = ylabel('Amplitude','interpreter','latex');
    hXLabel = xlabel('Month','interpreter','latex');
    set([hYLabel, hXLabel], 'interpreter', 'latex', 'Fontsize', 25, 'color', RdYlBu(60,:));

   
    
    set(gca, ...
      'Box'             , 'off'         , ...
      'TickDir'         , 'out'         , ...
      'TickLength'      , [.01 .01]     , ...
      'XMinorTick'      , 'off'         , ...
      'YMinorTick'      , 'off'         , ...
      'YGrid'           , 'on'          , ...
      'YTick'           , -1.5:1:3   , ...     % define grid, every 1000 m a grid line
      'XColor'          , RdYlBu(60,:)  , ...
      'YColor'          , RdYlBu(60,:)  , ...
      'ticklabelinterpreter', 'latex'   , ...
      'LineWidth'       , 1.25);
      'GridLineStyle'   , '--'          , ...

    
	set(gca, ...
        'XColor'          , RdYlBu(60,:)  , ...
        'YColor'          , RdYlBu(60,:));
    set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 25);
    

set(gca, 'Xtick', [3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48]);
pbaspect([2 1 1]);                        % aspect ratios: x, y, z

h27 = legend([h1 h5 h4 h8], ...
    'N34 of strong events', 'PC2 of strong events', ...
    'N34$_{\mathrm{ideal.}}$', 'PC2$_{\mathrm{ideal.}}$', ...
    'location', 'northeast', 'orientation', 'vertical');
set(h27, 'interpreter', 'latex', 'fontsize', 25, 'textcolor', RdYlBu(60,:));



ylim([-1.5 2.5]);
xlim([1 48]);

% pbaspect([6 1 1]);                        % aspect ratios: x, y, z

    arrow([1.1, -1.25], [6, -1.25], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [2], 'color', RdYlBu(1,:));
    arrow([18, -1.25], [23.9, -1.25], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [1], 'color', RdYlBu(1,:));
    text([6.75 6.75],[-1.25 -1.25],'El Ni\~no event', 'interpreter', 'latex', 'fontsize', 22, 'color', RdYlBu(1,:));

    arrow([24.1, -1.25], [31, -1.25], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [2], 'color', 'k');
    arrow([41, -1.25], [47.9, -1.25], 'BaseAngle', 70, 'length', 5, 'width', 1, 'ends', [1], 'color', 'k');
    text([32 32],[-1.25 -1.25],'spin-down', 'interpreter', 'latex', 'fontsize', 22, 'color', 'k');


% % finished fancy plot
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    set(gcf, 'color', 'w', 'PaperPositionMode', 'auto');
directory = 'E:\2017 Herbstsemester\Application Paper\AMS LaTeX Package v4.3.2\figures\patterns\';
%     print('-dpng','-r1000', [directory f1]);





