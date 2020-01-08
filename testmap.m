%% testmap script
function h = testmap(lon, lat, data)
RdYlBu              = cbrewer('div', 'RdYlBu', 60, 'PCHIP'); 
RdBu_short          = cbrewer('div', 'RdBu', 21, 'PCHIP'); 
RdYlBu_short        = cbrewer('div', 'RdYlBu', 21, 'PCHIP'); 
% myblue = the color I define for labels -> [.1922 .2118 .5843]

figure(27);
colormap(flipud(RdBu_short));

% tropical pacific ocean centered
% m_proj('miller','lat', [-30 30], 'lon',[100 300]);
% full world projection
% m_proj('miller','lat', [-60 60], 'lon',[18 378]);
m_proj('equidistant cylindrical','lat',[-90 90],'long',[22 382]);



h=m_pcolor(lon, lat, data); set(h,'linestyle','none'); hold on;
h=m_pcolor(lon+360, lat, data); set(h,'linestyle','none'); hold on;
set(h,'linestyle','none');
hold on;

m_coast('color',[.1 .1 .1]); %'patch', [.83 .83 .83]);
m_grid('box', 'on', 'xtick', 10, 'ytick', 10, 'tickdir', 'in', 'yaxislocation', 'left', ...
    'Fontname', 'Times New Roman', 'Fontsize', 12, 'color', RdYlBu(60,:));

title('testmap', 'units', 'normalized', 'position', [1 1], ...
    'HorizontalAlignment', 'right', 'Fontname', 'Times New Roman', 'Fontsize', 12);

h3 = colorbar('color', [.1922 .2118 .5843], 'box', 'on', 'tickdirection', 'in', ...
    'Fontname', 'Times New Roman', 'Fontsize', 12);
h = ylabel(h3, '[Testunit]', 'color', RdYlBu(60,:), ...
    'Fontname', 'Times New Roman', 'Fontsize', 12);


% set color limit %
% set(gca, 'clim', [-40 40]); hold on;













% draw rectangular area in Ni�o3.4 region: [-5 5], [170W, 120W]
% m_line([190, 190],[-5, 5],'color', [0 0 0], 'linewidth', .5); % m_line([lat],[lon])
% m_line([240, 240],[-5, 5], 'color', [0 0 0], 'linewidth', .5); % m_line([lat],[lon])
% m_line([190, 240],[-5, -5], 'color', [0 0 0], 'linewidth', .5); % m_line([lat],[lon])
% m_line([190, 240],[5, 5], 'color', [0 0 0], 'linewidth', .5); % m_line([lat],[lon])
% m_text(225.5, -2.3, 'Ni�o3.4','color', [0 0 0], 'fontsize', 10, ...
%     'Fontname', 'Times New Roman', 'Fontsize',6); hold on;

% set(gcf, 'PaperPositionMode', 'auto');
% directory = 'E:\2017 Fr�hlingssemester\Master Seminar 1\ProposalTemplate\figures\';
% print('-depsc','-r600', [directory 'testmap']);

