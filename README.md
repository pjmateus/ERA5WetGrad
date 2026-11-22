# ERA5WetGrad

**ERA5WetGrad** is a standalone MATLAB package for estimating GNSS zenith
wet delay and horizontal wet gradients from ERA5 three-dimensional wet
refractivity fields.

It traces synthetic straight rays through the atmosphere and estimates:

- zenith wet delay (`ZWD`);
- north–south wet gradient (`Gns`);
- east–west wet gradient (`Gew`).

The package is independent of any tomography workflow and uses a deliberately
small input interface. It does not require ERA5 PWV, mapping functions,
tomographic grids, file indices, or other application-specific products.

## Main features

- Minimal ERA5 input: time, longitude, latitude, height, and wet refractivity.
- Multiple ERA5 epochs and stations in one call.
- Configurable azimuths, elevations, ray step, and gradient mapping constant.
- Automatic sorting of increasing or decreasing ERA5 longitude/latitude axes.
- Support for both `[-180,180]` and `[0,360]` longitude conventions.
- Optional station-parallel processing with automatic serial fallback.
- Fit diagnostics: RMS, number of valid rays, and condition number.
- No dependencies on the original tomography codebase.

## Requirements

- MATLAB R2019b or newer.
- Parallel Computing Toolbox is optional.

Serial execution is the default and requires no additional toolbox.

## Repository contents

| File | Purpose |
|---|---|
| `estimateERA5WetGradients.m` | Standalone estimator and complete API documentation. |
| `example_estimateERA5WetGradients.m` | Runnable synthetic example and diagnostic figure. |
| `README.md` | Package overview and usage guide. |

## Installation

Clone or download the repository and add its folder to the MATLAB path:

```matlab
addpath('/path/to/ERA5WetGrad')
```

No compilation or installation step is required.

## Minimum input data

### ERA5 structure

`era5Data` contains exactly the atmospheric quantities needed by the ray
tracer:

| Field | Dimensions | Units | Description |
|---|---:|---:|---|
| `time` | `Nt x 1` | datetime | Epoch of each ERA5 state, normally UTC. |
| `lon_deg` | `Nlon x 1` | degree | ERA5 longitudes. |
| `lat_deg` | `Nlat x 1` | degree | ERA5 latitudes. |
| `z_m` | `Nlon x Nlat x Nlevel x Nt` | m | Geopotential or geometric height. |
| `Nw` | `Nlon x Nlat x Nlevel x Nt` | N-unit | Wet refractivity on the same grid as `z_m`. |

A three-dimensional `z_m`/`Nw` array is accepted when `Nt = 1`. Missing
atmospheric values may be `NaN`; infinite values and negative finite `Nw`
are rejected.

### Station structure

| Field | Dimensions | Units | Description |
|---|---:|---:|---|
| `lon_deg` | `Ns x 1` | degree | Station longitude. |
| `lat_deg` | `Ns x 1` | degree | Station latitude. |
| `height_m` | `Ns x 1` | m | Station height in the same vertical datum as `z_m`. |

## Basic usage

```matlab
era5Data.time    = timeUTC;
era5Data.lon_deg = longitude;
era5Data.lat_deg = latitude;
era5Data.z_m     = geopotentialHeight_m;
era5Data.Nw      = wetRefractivity_Nunits;

stations.lon_deg  = stationLongitude;
stations.lat_deg  = stationLatitude;
stations.height_m = stationHeight;

result = estimateERA5WetGradients(era5Data,stations);
```

The atmospheric arrays must be ordered as:

```text
longitude × latitude × pressure level × time
```

Every input ERA5 state is processed independently. ERA5WetGrad does not
perform temporal interpolation or temporal averaging.

## Configuring the calculation

Options are supplied in a scalar structure:

```matlab
opt = struct;
opt.rayStep_m = 500;
opt.az_deg = 0:45:315;
opt.el_deg = [10 15 20 30 45 55];
opt.gradientConstant = 0.0007;
opt.maxRayDistance_m = 500e3;
opt.minValidRays = 6;
opt.useParallel = true;
opt.numWorkers = 8;
opt.verbose = true;

result = estimateERA5WetGradients(era5Data,stations,opt);
```

| Option | Default | Description |
|---|---:|---|
| `az_deg` | `0:45:315` | Synthetic-ray azimuths, clockwise from north. |
| `el_deg` | `[10 15 20 30 45 55]` | Synthetic-ray elevations above the horizon. |
| `rayStep_m` | `500` | Sampling interval along each ray. |
| `maxRayDistance_m` | `500e3` | Maximum ray length. |
| `gradientConstant` | `0.0007` | Chen–Herring gradient mapping constant. |
| `minValidRays` | `6` | Minimum rays required for a three-parameter fit. |
| `stationIndex` | all | Station subset to process. |
| `timeIndex` | all | ERA5-state subset to process. |
| `earthRadius_m` | `6371000` | Spherical-Earth radius used by the local geometry. |
| `useParallel` | `false` | Use `parfor` over stations. |
| `numWorkers` | `[]` | Worker count when opening a new process pool. |
| `verbose` | `true` | Display progress and result ranges. |

The complete validation rules for every option are documented in the help
header of `estimateERA5WetGradients.m`.

## Output

The returned structure contains:

| Field | Units | Description |
|---|---:|---|
| `package`, `version` | — | Package identification (`ERA5WetGrad`, `0.1.0`). |
| `ZWD_m` | m | Zenith wet delay. |
| `Gns_m` | m | North–south wet gradient. |
| `Gew_m` | m | East–west wet gradient. |
| `gradientMagnitude_m` | m | `hypot(Gns_m,Gew_m)`. |
| `gradientAzimuth_deg` | degree | Direction of maximum gradient, clockwise from north. |
| `fitRMS_m` | m | RMS of slant-delay fit residuals. |
| `numberValidRays` | — | Rays used by each fit. |
| `conditionNumber` | — | Condition number of the design matrix. |
| `valid` | logical | True when ZWD, `Gns`, and `Gew` are finite. |
| `time`, `timeIndex` | — | Selected ERA5 epochs and original indices. |
| `stationIndex` | — | Original indices of the selected stations. |
| `options` | — | Validated options used in the calculation. |
| `elapsed_s` | s | Total elapsed wall time. |

All two-dimensional results have dimensions:

```text
selected ERA5 time × selected station
```

## Method

For each station and ERA5 state, the code generates rays at the requested
azimuths and elevations. Wet slant delay is obtained by integrating wet
refractivity along each valid ray:

$$
\mathrm{SWD}=10^{-6}\int N_w(s)\,ds.
$$

The GNSS model fitted by least squares is

$$
\mathrm{SWD}(a,e)=m_w(e)\,\mathrm{ZWD}
+m_g(e)\left(G_{ns}\cos a+G_{ew}\sin a\right),
$$

with

$$
m_w(e)=\frac{1}{\sin e},
\qquad
m_g(e)=\frac{1}{\sin e\tan e+C},
$$

where $C$ is `gradientConstant`.

Wet refractivity can be prepared from ERA5 pressure $p$ [Pa], temperature
$T$ [K], and specific humidity $q$ [kg kg\(^{-1}\)]. The function help gives
the complete conversion, including the refractivity constants used by the
example workflow.

## Runnable example

Run:

```matlab
example_estimateERA5WetGradients
```

The example constructs two synthetic anisotropic atmospheric states and two
stations, estimates the gradients, and plots `Gns`, `Gew`, and fit RMS.

## Adapting an existing atmospheric-state structure

If an existing workflow stores one structure per ERA5 state with the fields
`time`, `lon`, `lat`, `z`, and `Nw`, reduce it before calling ERA5WetGrad:

```matlab
state = fullERA5Data.field;

era5Data = struct;
era5Data.time    = vertcat(state.time);
era5Data.lon_deg = state(1).lon(:);
era5Data.lat_deg = state(1).lat(:);
era5Data.z_m     = cat(4,state.z);
era5Data.Nw      = cat(4,state.Nw);

gradientStations.lon_deg  = stations.lon_deg(:);
gradientStations.lat_deg  = stations.lat_deg(:);
gradientStations.height_m = stations.z_m(:);

clear fullERA5Data state
result = estimateERA5WetGradients(era5Data,gradientStations);
```

All states must use the same horizontal grid and array ordering.

## Assumptions and limitations

- Rays are straight; atmospheric bending is not modeled.
- Horizontal offsets use a local east–north–up frame and spherical Earth.
- Horizontal interpolation is linear at every pressure level.
- Vertical interpolation is linear in height.
- Below the lowest model level, the lowest-level wet refractivity is used.
- Integration stops at the highest model level.
- A ray is rejected if it leaves the ERA5 horizontal domain before reaching
  the atmospheric top.
- Regional grids crossing the longitude wrap discontinuity must be unwrapped
  before use.
- Station height and ERA5 height must share a consistent vertical datum.

## Validation

The supplied synthetic example produces valid solutions for every
station/time pair. The standalone implementation was also regression-tested
against the original optimized implementation using identical atmospheric
fields and angular geometry; maximum differences in ZWD, `Gns`, `Gew`, and
fit RMS were all exactly zero.

## Citation

Users should cite the ERA5 dataset and the gradient-mapping formulation
appropriate to their application. An archival citation or DOI can be added
when a release is deposited in a repository such as Zenodo.

## License

ERA5WetGrad is released under the MIT License. See [`LICENSE`](LICENSE) for
the full license text.
