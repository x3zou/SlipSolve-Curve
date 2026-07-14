# Legacy Inventory: InSAR Subsampling

This note records what the current Myanmar InSAR subsampling scripts actually do, before SlipSolve wraps or refactors them.

## Entry Scripts

### `sample_insar_data.m`

Bundled location: `external/legacy_matlab/geodetic_data/sample_insar_data.m`

Role: project-specific driver for initial quadtree sampling.

Hard-coded settings:

- Originally added its source checkout recursively to the MATLAB path; SlipSolve-curve now adds only the bundled dependency directories it needs.
- Uses fault trace `./trace.S2`.
- Uses local coordinate origin `lonc = 96.05`, `latc = 20.75`, `ref_lon = 96.05`.
- Uses `data_list = 'sample_list_alos_newmethod.txt'` by default.
- Loads `mesh_trace.txt`, divides it by `1e3`, and passes it as a fault-aware sampling trace in kilometers.
- Calls `make_insar_data(..., 'method', 'quadtree', ..., 'mesh_trace', mesh_trace)`.

The original driver expected to run from its own data directory because paths such as `./trace.S2`, `mesh_trace.txt`, and `./ALOS2/A152` were relative. SlipSolve-curve resolves configured inputs from the project layout instead.

### `sample_from_model.m`

Bundled location: `external/legacy_matlab/geodetic_data/sample_from_model.m`

Role: model-informed resampling after a preliminary inversion.

Hard-coded settings:

- Originally added its source checkout recursively to the MATLAB path; SlipSolve-curve now adds only the bundled dependency directories it needs.
- Uses fault trace `trace.ll`.
- Uses the same local origin as the initial sampler.
- Uses `iint = 3`.
- Loads `mesh_trace.txt`.
- Uses `data_list = 'sample_list_syn_noise_newgeometry.txt'`.
- Calls `make_insar_data_from_model(...)` with `nan_frac_max = 1`, `high = 1`, and `mesh_trace`.

## Core Initial Sampler

### `make_insar_data.m`

Bundled location: `external/legacy_matlab/geodetic_inversion-master/sampling/make_insar_data.m`

Role: reads one or more InSAR tracks from a sample-list text file, performs quadtree subsampling, saves a `.mat` file per track, and opens a comparison figure.

Important dependencies:

- `grdread2`
- `ll2xy`
- `expandGridNaN`
- `quadtree_unstructured2`
- `apply_unstructured_quadtree`
- `plot_insar_sample_new`
- `calc_insar_cov` and `get_insar_varigram` if covariance is enabled

Sample-list columns used by current initial sampler:

1. Track directory, relative to current working directory
2. Data grid file, usually `los_ll.grd`
3. Data type: `los`, `rng`, or azimuth/other
4. Nominal target point count, currently not used by the new method
5. Minimum longitude bound
6. Maximum longitude bound
7. Minimum latitude bound
8. Maximum latitude bound
9. Minimum quadtree cell width, in kilometers
10. Maximum quadtree cell width, in kilometers
11. RMS/variance threshold, in data units
12. Output `.mat` file name
13. X padding, in number of grid columns
14. Y padding, in number of grid rows

Example row from `sample_list_alos_newmethod.txt`:

```text
./ALOS2/A152 los_ll.grd los 1000 -1000 1000 18 23.5 1.2 30 1.3 los_samp0_1.mat 900 0
```

Processing behavior:

- Reads the LOS/range/azimuth grid with `grdread2`.
- Pads the grid with NaNs using `expandGridNaN` so quadtree cells remain closer to square.
- Reads matching look-vector grids:
  - LOS: `look_e.grd`, `look_n.grd`, `look_u.grd`
  - Low-resolution LOS: `look_e_low.grd`, `look_n_low.grd`, `look_u_low.grd`
  - High-resolution LOS: `look_e_high.grd`, `look_n_high.grd`, `look_u_high.grd`
  - Similar suffixes exist for `rng` and azimuth data.
- Crops by lon/lat bounds from the sample-list row.
- Converts lon/lat grid points to local Cartesian coordinates with `ll2xy`.
- Subtracts the origin coordinate from `lonc`, `latc`, then converts meters to kilometers for quadtree subdivision.
- Runs the fault-aware `quadtree_unstructured2` on the displacement values.
- Applies the same quadtree cells to look vectors with `apply_unstructured_quadtree`.
- Saves `sampled_insar_data = [x_m, y_m, displacement, look_e, look_n, look_u]`.
- Also saves `rms_out`; saves `covd` if covariance is enabled.
- Produces a figure from `plot_insar_sample_new` comparing raw cropped data and sampled points.

Output units:

- `sampled_insar_data(:,1:2)` are local x/y coordinates in meters.
- `sampled_insar_data(:,3)` is the sampled observed displacement in the original grid unit. The comments indicate centimeters.
- `sampled_insar_data(:,4:6)` are look-vector components.
- `rms_out` is per-cell spread/statistic from quadtree sampling.

## Fault-Aware Quadtree Helper

### `quadtree_unstructured2.m`

Bundled location: `external/legacy_matlab/geodetic_data/quadtree_unstructured2.m`

Role: recursive unstructured quadtree sampler with geometric side selection near a fault polyline.

Behavior:

- Splits a cell if RMS exceeds `rms_min` and the half-cell dimensions are larger than `width_min`, or if the half-cell dimensions exceed `width_max`.
- Otherwise finalizes the cell.
- If a final cell touches/crosses the fault trace, classifies points by geometric side of the polyline and averages only the side corresponding to the cell centroid, with fallback to the majority side.
- Supports `stat = 'mean'` or `stat = 'median'`.

Important options:

- `rms_min`
- `nan_frac_max`
- `width_min`
- `width_max`
- `fault_tol`
- `minPixSide`
- `stat`

## Visualization Contract For SlipSolve

Each completed sampling run should save at least:

- `figures/insar_subsampling/<track_name>_raw_vs_sampled.fig`
- `figures/insar_subsampling/<track_name>_raw_vs_sampled.png`
- Optional diagnostic figures:
  - quadtree cell extents over raw data
  - histogram of per-cell RMS/spread
  - map of sample density/cell size
  - fault-crossing cells highlighted

The interactive `.fig` files are important because users need to inspect sampling quality, adjust thresholds/cell sizes/padding, and rerun before continuing.

## Proposed SlipSolve `insarSub` Contract

The first wrapper should normalize the per-track legacy outputs into a single struct:

```matlab
insarSub = struct();
insarSub.tracks(k).name
insarSub.tracks(k).type
insarSub.tracks(k).sourceDirectory
insarSub.tracks(k).sourceGrid
insarSub.tracks(k).sampleFile
insarSub.tracks(k).x
insarSub.tracks(k).y
insarSub.tracks(k).los
insarSub.tracks(k).lookE
insarSub.tracks(k).lookN
insarSub.tracks(k).lookU
insarSub.tracks(k).rms
insarSub.tracks(k).covariance
insarSub.tracks(k).quadtree
insarSub.tracks(k).figures
insarSub.all.x
insarSub.all.y
insarSub.all.los
insarSub.all.lookE
insarSub.all.lookN
insarSub.all.lookU
insarSub.metadata.origin
insarSub.metadata.units
insarSub.metadata.sampleListFile
```

The wrapper should initially call the legacy code with minimal behavior changes, then load the produced `.mat` files and assemble this struct.
