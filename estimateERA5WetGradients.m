function result = estimateERA5WetGradients(era5Data,stations,opt)
%ESTIMATEERA5WETGRADIENTS Estimate GNSS wet gradients from ERA5 3-D fields.
%
%   result = estimateERA5WetGradients(era5Data,stations)
%   result = estimateERA5WetGradients(era5Data,stations,opt)
%
% This standalone function traces synthetic straight rays through 3-D wet
% refractivity fields and fits the standard three-parameter GNSS model
%
%   SWD(az,el) = mw(el)*ZWD + mg(el)* ...
%                (Gns*cos(az) + Gew*sin(az)),
%
% where
%
%   mw(el) = 1/sin(el)
%   mg(el) = 1/(sin(el)*tan(el) + gradientConstant).
%
% No temporal interpolation is performed. Each selected ERA5 state
% produces one independent set of gradients at every selected station.
%
% MINIMUM ERA5 INPUT
% -------------------------------------------------------------------------
% era5Data is a scalar structure containing only these five fields:
%
%   time      [Nt x 1] datetime
%       UTC epoch of each 3-D atmospheric state. A timezone is recommended
%       but not required. The function preserves the supplied datetimes.
%
%   lon_deg   [Nlon x 1] numeric
%       ERA5 grid longitudes [degrees]. Either [-180,180] or [0,360] is
%       accepted. Values must be unique. Increasing or decreasing input is
%       accepted and is sorted internally.
%
%   lat_deg   [Nlat x 1] numeric
%       ERA5 grid latitudes [degrees]. Values must be unique. Increasing or
%       decreasing input is accepted and is sorted internally.
%
%   z_m       [Nlon x Nlat x Nlevel x Nt] numeric
%       Geopotential or geometric height of each ERA5 pressure-level point
%       [m above mean sea level]. A 3-D array is accepted when Nt=1.
%
%   Nw        [Nlon x Nlat x Nlevel x Nt] numeric
%       Wet refractivity [N-units] on exactly the same grid as z_m. The
%       refractive-index relation is n_w-1 = Nw*1e-6. A 3-D array is
%       accepted when Nt=1.
%
% Nw can be prepared from ERA5 pressure p [Pa], temperature T [K], and
% specific humidity q [kg kg^-1], for example with Rueger (2002) constants:
%
%   epsilon = 18.01528/28.9647;
%   e  = q.*p ./ (epsilon + (1-epsilon).*q);       % water vapour [Pa]
%   k1 = 77.6890; k2 = 71.2952; k3 = 375463.0;    % hPa convention
%   k2Prime_Pa = (k2-epsilon*k1)*1e-2;
%   k3_Pa      = k3*1e-2;
%   Nw = k2Prime_Pa.*e./T + k3_Pa.*e./T.^2;
%
% STATION INPUT
% -------------------------------------------------------------------------
% stations is a scalar structure containing only:
%
%   lon_deg   [Ns x 1] station longitude [degrees]
%   lat_deg   [Ns x 1] station latitude [degrees]
%   height_m  [Ns x 1] station ellipsoidal/orthometric height [m]
%
% The station height must use the same vertical datum as era5Data.z_m.
%
% OPTIONS
% -------------------------------------------------------------------------
% opt is optional. Supported fields and defaults are:
%
%   az_deg = 0:45:315
%       Synthetic-ray azimuths [degrees clockwise from north].
%
%   el_deg = [10 15 20 30 45 55]
%       Synthetic-ray elevation angles [degrees above the horizon].
%
%   rayStep_m = 500
%       Distance interval used to sample each synthetic ray [m].
%
%   maxRayDistance_m = 500e3
%       Maximum ray length [m]. A ray is rejected if it leaves the
%       horizontal ERA5 domain before reaching the top model level.
%
%   gradientConstant = 0.0007
%       Constant in the Chen-Herring gradient mapping function.
%
%   minValidRays = 6
%       Minimum number of successfully integrated rays required for a fit.
%       At least three geometrically independent rays are always required.
%
%   stationIndex = 1:Ns
%       Indices of stations to process.
%
%   timeIndex = 1:Nt
%       Indices of ERA5 states to process.
%
%   earthRadius_m = 6371000
%       Spherical-Earth radius used for local horizontal displacements [m].
%
%   useParallel = false
%       If true, use PARFOR over stations. If Parallel Computing Toolbox or
%       a pool is unavailable, the function warns and continues serially.
%
%   numWorkers = []
%       Requested number of process workers when a new pool is opened.
%       Empty uses the MATLAB default. Existing pools are never resized or
%       closed by this function.
%
%   verbose = true
%       Print configuration, progress, and output ranges.
%
% OUTPUT
% -------------------------------------------------------------------------
% result is a structure with:
%
%   package, version              package identification
%   time                         [NselectedTime x 1] datetime
%   stationIndex                 [1 x NselectedStation]
%   timeIndex                    [1 x NselectedTime]
%   ZWD_m                        [NselectedTime x NselectedStation]
%   Gns_m                        [NselectedTime x NselectedStation]
%   Gew_m                        [NselectedTime x NselectedStation]
%   gradientMagnitude_m          hypot(Gns_m,Gew_m)
%   gradientAzimuth_deg          direction of maximum gradient, clockwise
%                                from north: atan2(Gew,Gns)
%   fitRMS_m                     RMS of slant-delay fit residuals [m]
%   numberValidRays              number of rays used in each fit
%   conditionNumber              condition number of the design matrix
%   valid                        true where ZWD, Gns, and Gew are finite
%   az_deg, el_deg               angular geometry actually used
%   options                      validated options
%   parallel, numberWorkers      execution information
%   elapsed_s                    total elapsed wall time
%
% METHOD AND LIMITATIONS
% -------------------------------------------------------------------------
% - Rays are straight in a local east-north-up frame; atmospheric bending
%   is not modelled.
% - Local horizontal offsets are converted using a spherical Earth.
% - z_m and Nw are interpolated horizontally at every pressure level and
%   then vertically to each ray point.
% - Below the lowest model level, the lowest-level Nw is used. Integration
%   stops when the ray reaches the highest model level.
% - Regional grids that cross the longitude wrap discontinuity should be
%   unwrapped by the caller before use.
%
% The routine requires MATLAB R2019b or newer. Parallel Computing Toolbox
% is optional.

if nargin<3 || isempty(opt)
    opt = struct;
end

[era5,station,opt] = validateAndPrepareInputs(era5Data,stations,opt);

stationIndex = opt.stationIndex;
timeIndex = opt.timeIndex;
nStations = numel(stationIndex);
nTimes = numel(timeIndex);
lonGrid_deg = era5.lon_deg;
latGrid_deg = era5.lat_deg;
stationLon_deg = station.lon_deg;
stationLat_deg = station.lat_deg;
stationHeight_m = station.height_m;
selectedStationLon_deg = stationLon_deg(stationIndex);
selectedStationLat_deg = stationLat_deg(stationIndex);
selectedStationHeight_m = stationHeight_m(stationIndex);

[AZ,EL] = ndgrid(opt.az_deg(:),opt.el_deg(:));
rayAz_deg = AZ(:);
rayEl_deg = EL(:);
nSyntheticRays = numel(rayAz_deg);

[useParallel,numberWorkers] = configureParallelExecution(opt);

ZWD_m = nan(nTimes,nStations);
Gns_m = nan(nTimes,nStations);
Gew_m = nan(nTimes,nStations);
fitRMS_m = nan(nTimes,nStations);
numberValidRays = zeros(nTimes,nStations);
conditionNumber = nan(nTimes,nStations);

if opt.verbose
    printHeader(era5,opt,nTimes,nStations,nSyntheticRays, ...
        useParallel,numberWorkers);
end

totalTimer = tic;
for jt = 1:nTimes
    stateTimer = tic;
    it = timeIndex(jt);
    z_m = era5.z_m(:,:,:,it);
    Nw = era5.Nw(:,:,:,it);

    nLevels = size(z_m,3);
    heightInterpolant = cell(nLevels,1);
    refractivityInterpolant = cell(nLevels,1);
    for k = 1:nLevels
        heightInterpolant{k} = griddedInterpolant( ...
            {lonGrid_deg,latGrid_deg},z_m(:,:,k),'linear','none');
        refractivityInterpolant{k} = griddedInterpolant( ...
            {lonGrid_deg,latGrid_deg},Nw(:,:,k),'linear','none');
    end

    ZWD_thisTime = nan(1,nStations);
    Gns_thisTime = nan(1,nStations);
    Gew_thisTime = nan(1,nStations);
    RMS_thisTime = nan(1,nStations);
    raysThisTime = zeros(1,nStations);
    conditionThisTime = nan(1,nStations);

    if useParallel
        parfor js = 1:nStations
            [ZWD_thisTime(js),Gns_thisTime(js),Gew_thisTime(js), ...
                RMS_thisTime(js),raysThisTime(js),conditionThisTime(js)] = ...
                processStation(selectedStationLon_deg(js), ...
                selectedStationLat_deg(js),selectedStationHeight_m(js), ...
                lonGrid_deg,latGrid_deg, ...
                heightInterpolant,refractivityInterpolant,rayAz_deg, ...
                rayEl_deg,opt);
        end
    else
        for js = 1:nStations
            [ZWD_thisTime(js),Gns_thisTime(js),Gew_thisTime(js), ...
                RMS_thisTime(js),raysThisTime(js),conditionThisTime(js)] = ...
                processStation(selectedStationLon_deg(js), ...
                selectedStationLat_deg(js),selectedStationHeight_m(js), ...
                lonGrid_deg,latGrid_deg, ...
                heightInterpolant,refractivityInterpolant,rayAz_deg, ...
                rayEl_deg,opt);
        end
    end

    ZWD_m(jt,:) = ZWD_thisTime;
    Gns_m(jt,:) = Gns_thisTime;
    Gew_m(jt,:) = Gew_thisTime;
    fitRMS_m(jt,:) = RMS_thisTime;
    numberValidRays(jt,:) = raysThisTime;
    conditionNumber(jt,:) = conditionThisTime;

    if opt.verbose
        fprintf('ERA5 state %d/%d: %s (%.2f s)\n',jt,nTimes, ...
            char(string(era5.time(it),'yyyy-MM-dd HH:mm:ss')), ...
            toc(stateTimer));
    end
end
elapsed_s = toc(totalTimer);

result = struct;
result.package = 'ERA5WetGrad';
result.version = '0.1.0';
result.time = era5.time(timeIndex);
result.stationIndex = stationIndex;
result.timeIndex = timeIndex;
result.ZWD_m = ZWD_m;
result.Gns_m = Gns_m;
result.Gew_m = Gew_m;
result.gradientMagnitude_m = hypot(Gns_m,Gew_m);
result.gradientAzimuth_deg = mod(atan2d(Gew_m,Gns_m),360);
result.fitRMS_m = fitRMS_m;
result.numberValidRays = numberValidRays;
result.conditionNumber = conditionNumber;
result.valid = isfinite(ZWD_m) & isfinite(Gns_m) & isfinite(Gew_m);
result.az_deg = opt.az_deg;
result.el_deg = opt.el_deg;
result.options = opt;
result.parallel = useParallel;
result.numberWorkers = numberWorkers;
result.elapsed_s = elapsed_s;
result.method = ['Straight-ray integration of ERA5 wet refractivity and ' ...
    'least-squares fit of ZWD, north-south gradient, and east-west gradient.'];

if opt.verbose
    printSummary(result);
end
end


function [era5,station,opt] = validateAndPrepareInputs(era5Data,stations,opt)
if ~isstruct(era5Data) || ~isscalar(era5Data)
    error('era5Data must be a scalar structure.');
end
requiredERA5 = {'time','lon_deg','lat_deg','z_m','Nw'};
if ~all(isfield(era5Data,requiredERA5))
    error('era5Data must contain: %s.',strjoin(requiredERA5,', '));
end
if ~isdatetime(era5Data.time) || isempty(era5Data.time) || ...
        any(isnat(era5Data.time(:)))
    error('era5Data.time must be a nonempty datetime array without NaT.');
end
if numel(unique(era5Data.time(:)))~=numel(era5Data.time)
    error('era5Data.time must not contain duplicated epochs.');
end

era5 = struct;
era5.time = era5Data.time(:);
era5.lon_deg = double(era5Data.lon_deg(:));
era5.lat_deg = double(era5Data.lat_deg(:));
era5.z_m = double(era5Data.z_m);
era5.Nw = double(era5Data.Nw);

if numel(era5.lon_deg)<2 || numel(era5.lat_deg)<2 || ...
        any(~isfinite(era5.lon_deg)) || any(~isfinite(era5.lat_deg)) || ...
        numel(unique(era5.lon_deg))~=numel(era5.lon_deg) || ...
        numel(unique(era5.lat_deg))~=numel(era5.lat_deg)
    error('ERA5 longitude/latitude vectors must be finite and unique.');
end
if ~isreal(era5.z_m) || ~isreal(era5.Nw) || ...
        ~isequal(size(era5.z_m),size(era5.Nw))
    error('era5Data.z_m and era5Data.Nw must be real arrays of equal size.');
end
expectedSize = [numel(era5.lon_deg),numel(era5.lat_deg)];
if size(era5.z_m,1)~=expectedSize(1) || ...
        size(era5.z_m,2)~=expectedSize(2) || size(era5.z_m,3)<2 || ...
        size(era5.z_m,4)~=numel(era5.time)
    error(['z_m/Nw must have size Nlon-by-Nlat-by-Nlevel-by-Nt, ' ...
        'with at least two levels.']);
end
if ~any(isfinite(era5.z_m(:))) || ~any(isfinite(era5.Nw(:)))
    error('z_m and Nw must contain finite values.');
end
if any(isinf(era5.z_m(:))) || any(isinf(era5.Nw(:))) || ...
        any(era5.Nw(isfinite(era5.Nw))<0)
    error('z_m/Nw may contain NaN but not Inf; finite Nw must be nonnegative.');
end

[era5.lon_deg,longitudeOrder] = sort(era5.lon_deg,'ascend');
[era5.lat_deg,latitudeOrder] = sort(era5.lat_deg,'ascend');
longitudeChanged = ~isequal(longitudeOrder,(1:numel(longitudeOrder)).');
latitudeChanged = ~isequal(latitudeOrder,(1:numel(latitudeOrder)).');
if longitudeChanged || latitudeChanged
    era5.z_m = era5.z_m(longitudeOrder,latitudeOrder,:,:);
    era5.Nw = era5.Nw(longitudeOrder,latitudeOrder,:,:);
end

if ~isstruct(stations) || ~isscalar(stations)
    error('stations must be a scalar structure.');
end
requiredStations = {'lon_deg','lat_deg','height_m'};
if ~all(isfield(stations,requiredStations))
    error('stations must contain: %s.',strjoin(requiredStations,', '));
end
station = struct;
station.lon_deg = double(stations.lon_deg(:));
station.lat_deg = double(stations.lat_deg(:));
station.height_m = double(stations.height_m(:));
nStation = numel(station.lon_deg);
if nStation<1 || numel(station.lat_deg)~=nStation || ...
        numel(station.height_m)~=nStation
    error('Station longitude, latitude, and height must have equal lengths.');
end
if any(~isfinite(station.lon_deg)) || any(~isfinite(station.lat_deg)) || ...
        any(~isfinite(station.height_m)) || ...
        any(station.lat_deg < -90 | station.lat_deg > 90)
    error('Station coordinates/heights must be finite and latitude in [-90,90].');
end

opt = parseOptions(opt,nStation,numel(era5.time));
end


function opt = parseOptions(userOptions,nStations,nTimes)
defaults = struct;
defaults.az_deg = 0:45:315;
defaults.el_deg = [10 15 20 30 45 55];
defaults.rayStep_m = 500;
defaults.maxRayDistance_m = 500e3;
defaults.gradientConstant = 0.0007;
defaults.minValidRays = 6;
defaults.stationIndex = 1:nStations;
defaults.timeIndex = 1:nTimes;
defaults.earthRadius_m = 6371000;
defaults.useParallel = false;
defaults.numWorkers = [];
defaults.verbose = true;

if ~isstruct(userOptions) || ~isscalar(userOptions)
    error('opt must be a scalar structure.');
end
unknown = setdiff(fieldnames(userOptions),fieldnames(defaults));
if ~isempty(unknown)
    error('Unknown option field(s): %s.',strjoin(unknown,', '));
end
opt = defaults;
names = fieldnames(userOptions);
for k = 1:numel(names)
    if ~isempty(userOptions.(names{k}))
        opt.(names{k}) = userOptions.(names{k});
    end
end

opt.az_deg = mod(double(opt.az_deg(:).'),360);
opt.el_deg = double(opt.el_deg(:).');
opt.stationIndex = double(opt.stationIndex(:).');
opt.timeIndex = double(opt.timeIndex(:).');

if isempty(opt.az_deg) || any(~isfinite(opt.az_deg)) || ...
        isempty(opt.el_deg) || any(~isfinite(opt.el_deg)) || ...
        any(opt.el_deg<=0 | opt.el_deg>=90)
    error('az_deg must be finite; el_deg must lie strictly between 0 and 90.');
end
if ~isPositiveScalar(opt.rayStep_m) || ...
        ~isPositiveScalar(opt.maxRayDistance_m) || ...
        opt.rayStep_m>=opt.maxRayDistance_m
    error('rayStep_m and maxRayDistance_m must be positive, with step < max.');
end
if ~isnumeric(opt.gradientConstant) || ~isscalar(opt.gradientConstant) || ...
        ~isfinite(opt.gradientConstant) || opt.gradientConstant<0
    error('gradientConstant must be a finite nonnegative scalar.');
end
nSyntheticRays = numel(opt.az_deg)*numel(opt.el_deg);
if ~isIntegerScalar(opt.minValidRays) || opt.minValidRays<3 || ...
        opt.minValidRays>nSyntheticRays
    error('minValidRays must be an integer from 3 to the number of rays.');
end
if ~isIndexVector(opt.stationIndex,nStations)
    error('stationIndex must contain valid station indices.');
end
if ~isIndexVector(opt.timeIndex,nTimes)
    error('timeIndex must contain valid ERA5 time indices.');
end
if ~isPositiveScalar(opt.earthRadius_m)
    error('earthRadius_m must be a positive finite scalar.');
end
if ~islogical(opt.useParallel) || ~isscalar(opt.useParallel) || ...
        ~islogical(opt.verbose) || ~isscalar(opt.verbose)
    error('useParallel and verbose must be logical scalars.');
end
if ~isempty(opt.numWorkers) && ...
        (~isIntegerScalar(opt.numWorkers) || opt.numWorkers<1)
    error('numWorkers must be empty or a positive integer.');
end
end


function tf = isPositiveScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && value>0;
end


function tf = isIntegerScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value==fix(value);
end


function tf = isIndexVector(value,upperBound)
tf = isnumeric(value) && isvector(value) && ~isempty(value) && ...
    all(isfinite(value)) && all(value==fix(value)) && ...
    all(value>=1 & value<=upperBound) && numel(unique(value))==numel(value);
end


function [useParallel,numberWorkers] = configureParallelExecution(opt)
useParallel = opt.useParallel;
numberWorkers = 1;
if ~useParallel
    return
end
if isempty(ver('parallel'))
    warning(['Parallel Computing Toolbox is unavailable; continuing ' ...
        'with serial execution.']);
    useParallel = false;
    return
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(opt.numWorkers)
            pool = parpool('Processes');
        else
            pool = parpool('Processes',opt.numWorkers);
        end
    end
    numberWorkers = pool.NumWorkers;
catch exception
    warning('estimateERA5WetGradients:ParallelUnavailable','%s', ...
        ['Parallel pool unavailable (' exception.message ...
        '); continuing serially.']);
    useParallel = false;
    numberWorkers = 1;
end
end


function [ZWD_m,Gns_m,Gew_m,rms_m,nValid,conditionNumber] = ...
        processStation(stationLon_deg,stationLat_deg,stationHeight_m, ...
        lonGrid_deg,latGrid_deg,heightInterpolant, ...
        refractivityInterpolant,rayAz_deg,rayEl_deg,opt)
ZWD_m = NaN;
Gns_m = NaN;
Gew_m = NaN;
rms_m = NaN;
nValid = 0;
conditionNumber = NaN;

if ~all(isfinite([stationLon_deg stationLat_deg stationHeight_m]))
    return
end
stationLon_deg = normalizeLongitude(stationLon_deg,lonGrid_deg);

nRays = numel(rayAz_deg);
slantWetDelay_m = nan(nRays,1);
for r = 1:nRays
    slantWetDelay_m(r) = traceStraightRay( ...
        stationLon_deg,stationLat_deg,stationHeight_m, ...
        rayAz_deg(r),rayEl_deg(r),lonGrid_deg,latGrid_deg, ...
        heightInterpolant,refractivityInterpolant,opt);
end

valid = isfinite(slantWetDelay_m) & isfinite(rayAz_deg) & ...
    isfinite(rayEl_deg);
nValid = nnz(valid);
if nValid<opt.minValidRays
    return
end

azimuth_deg = rayAz_deg(valid);
elevation_deg = rayEl_deg(valid);
observed_m = slantWetDelay_m(valid);
wetMapping = 1./sind(elevation_deg);
gradientMapping = 1./(sind(elevation_deg).*tand(elevation_deg) + ...
    opt.gradientConstant);
design = [wetMapping,gradientMapping.*cosd(azimuth_deg), ...
    gradientMapping.*sind(azimuth_deg)];
if rank(design)<3
    return
end

parameters = design\observed_m;
ZWD_m = parameters(1);
Gns_m = parameters(2);
Gew_m = parameters(3);
residual_m = observed_m-design*parameters;
rms_m = sqrt(mean(residual_m.^2));
conditionNumber = cond(design);
end


function slantWetDelay_m = traceStraightRay(lon0_deg,lat0_deg,z0_m, ...
        azimuth_deg,elevation_deg,lonGrid_deg,latGrid_deg, ...
        heightInterpolant,refractivityInterpolant,opt)
slantWetDelay_m = NaN;

cosElevation = cosd(elevation_deg);
eastDirection = cosElevation*sind(azimuth_deg);
northDirection = cosElevation*cosd(azimuth_deg);
upDirection = sind(elevation_deg);
if upDirection<=0
    return
end

distance_m = (0:opt.rayStep_m:opt.maxRayDistance_m).';
east_m = distance_m*eastDirection;
north_m = distance_m*northDirection;
zRay_m = z0_m+distance_m*upDirection;
latRay_deg = lat0_deg+rad2deg(north_m/opt.earthRadius_m);
cosLatitude = cosd(lat0_deg);
if abs(cosLatitude)<1e-8
    return
end
lonRay_deg = lon0_deg+rad2deg(east_m/(opt.earthRadius_m*cosLatitude));

inside = lonRay_deg>=lonGrid_deg(1) & lonRay_deg<=lonGrid_deg(end) & ...
    latRay_deg>=latGrid_deg(1) & latRay_deg<=latGrid_deg(end);
firstOutside = find(~inside,1,'first');
if ~isempty(firstOutside)
    lastInside = firstOutside-1;
    distance_m = distance_m(1:lastInside);
    lonRay_deg = lonRay_deg(1:lastInside);
    latRay_deg = latRay_deg(1:lastInside);
    zRay_m = zRay_m(1:lastInside);
end
numberPoints = numel(distance_m);
if numberPoints<2
    return
end

numberLevels = numel(heightInterpolant);
profileHeight_m = nan(numberPoints,numberLevels);
profileNw = nan(numberPoints,numberLevels);
for k = 1:numberLevels
    profileHeight_m(:,k) = heightInterpolant{k}(lonRay_deg,latRay_deg);
    profileNw(:,k) = refractivityInterpolant{k}(lonRay_deg,latRay_deg);
end

rayNw = nan(numberPoints,1);
reachedTop = false;
topIndex = NaN;
for i = 1:numberPoints
    zProfile_m = profileHeight_m(i,:).';
    NwProfile = profileNw(i,:).';
    valid = isfinite(zProfile_m) & isfinite(NwProfile);
    zProfile_m = zProfile_m(valid);
    NwProfile = NwProfile(valid);
    if numel(zProfile_m)<2
        break
    end
    [zProfile_m,order] = sort(zProfile_m,'ascend');
    NwProfile = NwProfile(order);
    [zProfile_m,uniqueIndex] = unique(zProfile_m,'stable');
    NwProfile = NwProfile(uniqueIndex);
    if numel(zProfile_m)<2
        break
    end
    if zRay_m(i)>=zProfile_m(end)
        rayNw(i) = NwProfile(end);
        reachedTop = true;
        topIndex = i;
        break
    elseif zRay_m(i)<=zProfile_m(1)
        rayNw(i) = NwProfile(1);
    else
        rayNw(i) = interp1(zProfile_m,NwProfile,zRay_m(i),'linear');
    end
end
if ~reachedTop
    return
end

distance_m = distance_m(1:topIndex);
rayNw = rayNw(1:topIndex);
valid = isfinite(distance_m) & isfinite(rayNw);
distance_m = distance_m(valid);
rayNw = rayNw(valid);
if numel(distance_m)<2
    return
end
slantWetDelay_m = 1e-6*trapz(distance_m,rayNw);
end


function longitude_deg = normalizeLongitude(longitude_deg,gridLongitude_deg)
longitude_deg = double(longitude_deg);
if all(gridLongitude_deg>=0) && max(gridLongitude_deg)>180
    longitude_deg = mod(longitude_deg,360);
elseif min(gridLongitude_deg)<0 && max(gridLongitude_deg)<=180
    longitude_deg = mod(longitude_deg+180,360)-180;
end
end


function printHeader(era5,opt,nTimes,nStations,nRays,useParallel,nWorkers)
fprintf('\nERA5 -> GNSS wet-gradient estimation\n');
fprintf('ERA5 states              : %d selected / %d available\n', ...
    nTimes,numel(era5.time));
fprintf('Stations                 : %d\n',nStations);
fprintf('Synthetic rays per fit   : %d\n',nRays);
fprintf('Ray step / maximum       : %.1f m / %.1f km\n', ...
    opt.rayStep_m,opt.maxRayDistance_m/1000);
fprintf('Parallel execution       : %s',string(useParallel));
if useParallel
    fprintf(' (%d workers)',nWorkers);
end
fprintf('\n\n');
end


function printSummary(result)
fprintf('\nERA5 wet-gradient estimation finished in %.2f s.\n', ...
    result.elapsed_s);
fprintf('Valid station/time fits  : %d / %d\n', ...
    nnz(result.valid),numel(result.valid));
printFiniteRange('ZWD',result.ZWD_m,1000,'mm');
printFiniteRange('Gns',result.Gns_m,1000,'mm');
printFiniteRange('Gew',result.Gew_m,1000,'mm');
printFiniteRange('Fit RMS',result.fitRMS_m,1000,'mm');
end


function printFiniteRange(label,values,scale,unit)
values = scale*values(isfinite(values));
if ~isempty(values)
    fprintf('%-25s: %.6g to %.6g %s\n',label,min(values),max(values),unit);
end
end
