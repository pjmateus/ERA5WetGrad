%% Example: GNSS wet gradients from a minimal ERA5-like data structure
% This example is self-contained and uses a synthetic wet-refractivity
% field. Replace only the data-construction section with ERA5 arrays; the
% call to estimateERA5WetGradients remains unchanged.

clear
close all

%% Minimal ERA5-like input
era5Data.time = datetime(2025,1,1,[0;6],0,0,'TimeZone','UTC');
era5Data.lon_deg = (-8:0.25:-2).';
era5Data.lat_deg = (37:0.25:43).';
heightLevels_m = (0:1000:12000).';

[LON,LAT,Z] = ndgrid(era5Data.lon_deg,era5Data.lat_deg,heightLevels_m);
baseNw = 85*exp(-Z/2200);

% Two anisotropic atmospheric states. Array order must be
% longitude-by-latitude-by-level-by-time.
era5Data.z_m = repmat(Z,1,1,1,2);
era5Data.Nw = zeros(size(era5Data.z_m));
era5Data.Nw(:,:,:,1) = baseNw.* ...
    (1+0.020*(LON+5)+0.025*(LAT-40));
era5Data.Nw(:,:,:,2) = baseNw.* ...
    (1+0.030*(LON+5)-0.015*(LAT-40));

%% Minimal station input
stations.lon_deg = [-5.0;-4.5];
stations.lat_deg = [40.0;39.5];
stations.height_m = [100;250];

%% Options
opt = struct;
opt.rayStep_m = 500;
opt.az_deg = 0:45:315;
opt.el_deg = [10 15 20 30 45 55];
opt.gradientConstant = 0.0007;
opt.useParallel = false;
opt.verbose = true;

%% Estimate ZWD, north-south gradient, and east-west gradient
gradient = estimateERA5WetGradients(era5Data,stations,opt);

%% Compact diagnostic plot
figure('Color','w','Position',[100 100 1100 330]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile
imagesc(1000*gradient.Gns_m.')
colorbar
xlabel('ERA5 time index')
ylabel('Station')
title('G_{ns} [mm]')

nexttile
imagesc(1000*gradient.Gew_m.')
colorbar
xlabel('ERA5 time index')
ylabel('Station')
title('G_{ew} [mm]')

nexttile
imagesc(1000*gradient.fitRMS_m.')
colorbar
xlabel('ERA5 time index')
ylabel('Station')
title('Fit RMS [mm]')

set(findall(gcf,'Type','axes'),'FontSize',11)
