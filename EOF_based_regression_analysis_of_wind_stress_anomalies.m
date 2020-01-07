%% Regression Analysis: 
% Creating the two wind stress anomaly spatial patterns and their principle 
% component time series as in McGregor et al., 2014

% Dimensions of the Era-Interim data:
%            LAT  = 480                             % 0.75x0.75 degree grid
%            LON  = 241
%            TIME = 456  (UNLIMITED)                % monthly data from 
%                                                     Jan 1979 - Dec 2016

% loading in NOAA ERSSTv4 Nino3.4 index   
load('workspace_NOAA_ERSST_nino34.mat', 'nino'); 
  

%% [0.02s] preamble and data preparation
tic;
% custom colours
RdYlBu = cbrewer('div', 'RdYlBu', 60, 'PCHIP'); 
RdYlBu_short = cbrewer('div', 'RdYlBu', 21, 'PCHIP');
RdBu_short = cbrewer('div', 'RdBu', 21, 'PCHIP');

p1 = '/srv/ccrc/data15/z5180028/MSC_thesis_datasets/era_interim_anomaly_fields_CDO/';
% this is my file with my wind stresses
f1 = 'era_interim_tau_x_tau_y_monthly_anomaly_fields_1979_2016';
% load in instantaneous eastward/northward turbulent surface stress
tau_x                = getnc([p1 f1], 'iews');   
tau_y                = getnc([p1 f1], 'inss');   
tau_x = permute(tau_x(1:456,:,:), [3 2 1]);        % permute to get correct
tau_y = permute(tau_y(1:456,:,:), [3 2 1]);        % dimensions, i.e.
                                                   % [480, 241, 456]

lon = getnc([p1 f1], 'longitude'); % load in longitude and latitude values
lat = getnc([p1 f1], 'latitude'); 
[lat,lon]=meshgrid(lat,lon); % create meshgrid for plotting
time = linspace(1979, 2016, 456);

toc;    
  

%% [3.12s] EOF1: regressing NOAA ERSSTv4's Niño3.4 onto wind stress anomalies

% variables
[xL,yL] = size(lon);                      % dimensions of my map [480, 241]
tL = length(time);                        % [456]

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% regressing the first PC score
% X_{ideal.} = 1/n * sum(X(x,y,t) * PC(t)) / var(PC(t))

for i = 1:length(time)
     s(:,:,i) = (tau_x(:,:,i) .* nino(i) ./ (std(nino)^2));
     t(:,:,i) = (tau_y(:,:,i) .* nino(i) ./ (std(nino)^2));
end
clear i;
EOFs(:,:,1) = nanmean(s,3);     % take mean over time period
EOFs(:,:,2) = nanmean(t,3);

% checking with figures, great!
testmap(lon, lat, EOFs(:,:,1));

% EOFs(:,:,1) = first zonal (East-West) wind stress pattern
% EOFs(:,:,2) = first meridional (North-South) wind stress pattern


%% [13.41s] EOF2: removing data associated with the first EOF1 spatial
% pattern and calculating the first EOF of the residuals (which is then
% EOF2)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% EOF region as in McGregor et al., 2014: 10S - 10N & 100E - 60W
%                                         lat(135:401, 107:134);
% calculation only over this domain = [-100.5 300 -9.75 9.75];
EOFdomain = [100.5 300 -9.75 9.75];

% part from Ryan, first defining variables for manual EOF analysis later on
[tmp lt1] = min(abs(lat(1,:)-EOFdomain(3)));
[tmp lt2] = min(abs(lat(1,:)-EOFdomain(4)));
if (lt1>lt2); tmp=lt2; lt2=lt1;lt1=tmp; end;
[tmp ln1] = min(abs(lon(:,1)-EOFdomain(1)));
[tmp ln2] = min(abs(lon(:,1)-EOFdomain(2)));

lon_new = lon(ln1:ln2);
lat_new = lat(lt1:lt2);
Y1 = tau_x(135:401,107:134,:);
Y2 = tau_x(135:401,107:134,:);
[EOFxL,EOFyL,tmp] = size(Y1);

% subtracting the EOF1 data
for i = 1:length(time) 
    residualx(:,:,i) = tau_x(135:401,107:134,i) - ...
                                         EOFs(135:401,107:134,1) * nino(i);   
    residualy(:,:,i) = tau_y(135:401,107:134,i) - ...
                                         EOFs(135:401,107:134,2) * nino(i); 
end

% reshaping arrays to be able to better calculate the remaining variance
a = reshape(residualx, [267*28, 456])';
b = reshape(residualy, [267*28, 456])';

% ~~~~~~~~~~~~~~~~~~~~~~~~~~ EOF Analysis here ~~~~~~~~~~~~~~~~~~~~~~~~~~ %
A = cat(2, a, b);
C = A'*A;                       % A' = transposed matrix
nmod = 1;                       % what we want, i.e. the 1st EOF (mode)
ntot = 30;                      % in total 30 EOFs are calculated
tmp = eigs(C,ntot);
varfrac = tmp./sum(tmp);        % calculate variance fraction
[V,D] = eigs(C,nmod);           % calculate Eigenvalues
PC(2,:) = A*V;                  % reconstruct the principle component time 
                                % series

% scaling
std1 = std(nino);
std2 = std(PC(2,:));
PC(1,:) = nino;                 % first time series  = Nino3.4
PC(2,:) = PC(2,:) .* std2;      % second time series = PC2

% regressing the new-found time series onto my anomaly fields so I get the
% complete global map instead of the Pacific region only
clear a b;
for i = 1:length(time)
     a(:,:,i) = (tau_x(:,:,i) .* PC(2,i)) ./ (std(PC(2,:).^2));
     b(:,:,i) = (tau_y(:,:,i) .* PC(2,i)) ./ (std(PC(2,:).^2));
end
% ~~~~~~~~~ --------------------------------------------------- ~~~~~~~~~ %
EOFs(:,:,3) = nanmean(a,3);         % taking the time-mean as stated in
EOFs(:,:,4) = nanmean(b,3);         % the regression equation
toc;

% great: now I have two time series PC(1,:) = Nino3.4 and PC(2,:) = PC2
% and four EOF patterns: for each time series the zonal and meridional 
% components
% -> EOFs(:,:,1) = first EOF of zonal wind stress
% -> EOFs(:,:,2) = first EOF of meridional wind stress
% -> EOFs(:,:,3) = second "   "    zonal    "   "
% -> EOFs(:,:,4) = second "   "    meridional    "   "



%% [0.22s] Removing every 8th components so I do not have as many wind 
% stress arrows on my map -> clearer visualization
% so the map looks cleaner
tic;
clear s t;
s(:,:,1) = EOFs(:,:,1); 
s(:,:,2) = EOFs(:,:,2); 
t(:,:,1) = EOFs(:,:,3); 
t(:,:,2) = EOFs(:,:,4); 

% only take every fith column and row
for i = 1:480
    for l = 1:241
        for r = 1:2
        if mod(i,8) == 0                % if division by 8 ends up with no
            s(i,:,r) = s(i,:,r);        % rest, then good
            t(i,:,r) = t(i,:,r);
            if mod(l,8) == 0;
                s(:,l,r) = s(:,l,r);
                t(:,l,r) = t(:,l,r);
            else                        % else I just put an empty cell 
                s(:,l,r) = nan;         % there
                t(:,l,r) = nan;
            end
        else
            s(i,l,r) = nan;
            t(i,l,r) = nan;
        end
    end
    end
end
arrows(:,:,1:2) = s; 
arrows(:,:,3:4) = t; 
% clear workspace before saving 
clear s t i l r a A ans b C D f1 nmod ntot p1 pattern1 RdBu_short RdYlBu;
clear RdYlBu_short residual std1 std2 tau_x tau_y tmp u10 V v10 varfrac;
clear EOF3 EOF4 nino;
clear EOFdomain EOFxL EOFyL lat_new lon_new ln1 ln2 lt1 lt2 residualx;
clear residualy tL xL Y1 Y2 yL;
toc;
  

%% [super fast] save workspace

filename2017 = 'workspace_regression_patterns_PC1_equal_nino34_rev4';
save(filename2017);







