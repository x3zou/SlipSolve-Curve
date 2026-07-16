# SlipSolve-curve

SlipSolve-curve is a curved-fault variant of [SlipSolve](https://github.com/ZeyuJin/geodetic_inversion) for reproducible MATLAB
coseismic-slip inversion on triangular fault meshes. It preserves the
numerical methods and conventions of the legacy workflow while making inputs,
stage boundaries, parameters, figures, and restart products explicit.

The implemented workflow includes InSAR quadtree sampling, curved geometry,
legacy triangular meshing, homogeneous TDE inversion, model-based resampling,
EDGRN conversion, full XYZ composite Greens, and the final altered-node
layered inversion. External inputs are treated as read-only, and generated
products stay inside the SlipSolve-curve project.

## Requirements

- MATLAB R2025b is the tested version.
- Optimization Toolbox is required for `lsqlin`.
- Statistics and Machine Learning Toolbox is required by helpers such as
  `knnsearch`.
- Required legacy MATLAB numerical helpers are internalized under
  `external/legacy_matlab`; no separate legacy-code checkout is required.
- EDGRN is required only when creating a new layered Green database. It is not
  needed for the bundled final-inversion example.

## Quick Start: Bundled Composite Example

The supplied example intentionally runs only the final composite inversion.
Large example grids and Green matrices are distributed separately through
Google Drive and are not stored in GitHub. Before running MATLAB:

1. Download the Google Drive folder named
   `SlipSolve-curve-example-data`. Link is [here](https://drive.google.com/file/d/1Jktskt3gbN6-XIxi6REiWo71u5H28UgI/view?usp=sharing).
2. Copy its `data` and `greens` folders into the current SlipSolve-curve directory.
3. In MATLAB, verify the installation:

```matlab
addpath examples/myanmar_demo
verify_example_data
```

The verifier checks all 31 downloadable files against the tracked size
manifest. For a full byte-level check on macOS or Linux, run
`shasum -a 256 -c examples/myanmar_demo/example_data_sha256.txt` from the
project root.

**Maintainer:** replace `GOOGLE_DRIVE_DOWNLOAD_LINK` here with the public
Google Drive share URL after uploading the prepared folder.

Then, from the project root, run:

```matlab
run examples/myanmar_demo/run_example.m
```

The example reads all inversion and plotting parameters from
`config/example_project.m`, but forces both workflow endpoints to
`final_inversion`. It reuses these project-local products:

- `geometry/triangular_mesh.mat`
- Seven copied `samp3` files under `data/raw/legacy_samp3`
- `greens/composite_full_xyz/G_e.mat`, `G_n.mat`, and `G_u.mat`
- `greens/composite_sampled_tracks/sampled_composite_greens.mat`
- `data/raw/legacy_reference/vel2.mat`

The Google Drive bundle also supplies the full-resolution source `.grd` files
used by initial sampling, optional final data/model/residual comparisons, and
the independent Sentinel-2 forward-model example.

Therefore it does **not** run quadtree sampling, geometry construction, TDE,
EDGRN conversion, or the long composite-Green calculation. The three full
component matrices occupy about 4.2 GB, but they are precomputed rather than
rebuilt. Main outputs are `inversions/final_result.mat`, the simple-triangular
model, and interactive figures under `figures/final_layered_inversion`.

The normal nodal-slip, triangular-slip, and sampled data/model/residual
figures are enabled. The much slower optional full-resolution residual plots
are disabled by default; set
`cfg.visualization.finalInversion.fullResolutionFitPlots.enabled=true` to
generate them as well. Figures pop up in MATLAB Desktop and render offscreen
when the same script is run through headless or batch MATLAB.

To change the example, edit `cfg.finalInversion` and
`cfg.visualization.finalInversion` in `config/example_project.m`, then rerun
the same example script.

## Run Modes

| Goal | Entry point | Stage range |
|---|---|---|
| Try the supplied precomputed example | `examples/myanmar_demo/run_example.m` | Final composite inversion only |
| Build a new project through TDE | `workflows/run_quick_inversion.m` | Any subset of stages 1-6 |
| Build or run layered/composite products | `workflows/run_full_inversion.m` | Any subset of stages 1-9 |
| Sweep inversion smoothness | `workflows/run_l_curve.m` | Optional post-processing |
| Compute shallow slip deficit | `workflows/run_shallow_slip_deficit.m` | Optional post-processing |
| Predict independent observations | `workflows/run_forward_model.m` | Optional post-processing |

`run_quick_inversion.m` cannot run the layered or final stages. Use
`run_full_inversion.m` whenever `startStep` or `stopStep` is
`layered_greens`, `composite_greens`, or `final_inversion`.

## Workflow Map

| # | Stage name | Principal required input | Principal output |
|---|---|---|---|
| 1 | `subsample_insar` | Raw displacement and look-vector grids | `data/processed/insar_quadtree.mat` |
| 2 | `fault_geometry` | Surface trace and depth controls | `geometry/fault_geometry.mat` |
| 3 | `triangular_mesh` | Fault geometry | `geometry/triangular_mesh.mat` |
| 4 | `quick_inversion` | Initial samples and mesh | `inversions/quick_result.mat` and `model.grd` |
| 5 | `model_based_sampling` | Quick model and original grids | Per-track model-based sample MAT files |
| 6 | `resampled_tde_inversion` | Mesh and sampled observations | `inversions/resampled_tde_result.mat` |
| 7 | `layered_greens` | Completed EDGRN text database | Converted EDGRN MAT database |
| 8 | `composite_greens` | Mesh and converted EDGRN database | Full `G_e`, `G_n`, and `G_u` matrices |
| 9 | `final_inversion` | Mesh, XYZ Greens, and sampled data | Final nodal and triangular slip models |

Every computational stage saves an interactive `.fig` and a PNG when its
visualization is enabled. Stage-specific controls are grouped under
`cfg.visualization`.

## Figure Controls And Coordinate Origins

Every workflow figure has its user-facing style controls in the corresponding
`cfg.visualization.<figureGroup>` block. Depending on the figure, these include:

- `showTitle` or `showTitles`, plus `title` or `panelTitles`
- `axisRange`, in the units documented beside that setting
- `colorbarRange`, or separate data/model/residual or component ranges
- `colormap`
- `view` for 3-D figures

An empty range (`[]`) keeps automatic limits. Multi-track full-resolution fit
figures can override title, map range, and colorbar range in each track entry,
so different satellites do not have to share plotting limits. Figures without
a scalar color field, such as a line-profile plot, expose line colors and axis
ranges instead of an inapplicable colorbar.

The geographic reference is fully user configurable:

```matlab
cfg.insar.origin.lon = 96.05;    % longitude of local x=0
cfg.insar.origin.lat = 20.75;    % latitude of local y=0
cfg.insar.origin.refLon = 96.05; % ll2xy reference longitude/central meridian
```

Replace these values for a new study area. All geographic grids, fault traces,
mesh products, sampled coordinates, and full composite Greens in one project
must be built with the same `cfg.insar.origin`; changing it after those products
exist requires rebuilding or consistently transforming them.

`cfg.visualization.finalInversion.plotOrigin` and the corresponding
`fullResolutionFitPlots.plotOrigin` are display-only origins. They recenter map
axes without changing the inversion coordinates or Green matrices. Set them
equal to `cfg.insar.origin` to display the native local coordinate system.

## Stage Selection And Restarts

Set the first and last stage near the beginning of
`config/example_project.m`:

```matlab
cfg.workflow.startStep = "subsample_insar";
cfg.workflow.stopStep = "model_based_sampling";
cfg.workflow.pauseAfterStage = true;
```

Valid stage names are the nine names in the workflow map. A restarted stage
loads its required saved predecessor product. If that file is missing, the
workflow reports which earlier stage must be run.

With `cfg.visualization.visible=true`, figures pop up during execution. With
`cfg.workflow.pauseAfterStage=true`, interactive MATLAB sessions pause after
each stage so the user can inspect the figure before continuing. Set the pause
to false for unattended runs. Saved `.fig` files remain interactive even when
`visible=false`.

## Stages 1-5: Data, Geometry, Quick Model, And Resampling

### Stage 1: Initial InSAR Sampling

Add one `cfg.insar.tracks` entry per AZO, LOS, or range data set. Each entry
points to a displacement grid and its east/north/up look-vector grids. Initial
sampling parameters are independent for every track:

```matlab
cfg.insar.tracks(1).initialSampling.boundsLonLat = [];
cfg.insar.tracks(1).initialSampling.minCellSizeKm = 0.5;
cfg.insar.tracks(1).initialSampling.maxCellSizeKm = 20;
cfg.insar.tracks(1).initialSampling.rmsThreshold = 0.0025;
cfg.insar.tracks(1).initialSampling.nanFractionMax = 1;
```

The implementation uses `quadtree_unstructured2.m`, including its fault-aware
near-trace handling. Smaller cells and lower RMS thresholds usually retain
more observations. Inspect each data/quadtree comparison figure before moving
on. The saved sample table uses local x/y in metres and displacement in the
configured inversion unit.

### Stage 2: Curved Fault Geometry

Provide one surface trace plus either explicit depth controls or surface
segments with dip angles. For multiple disconnected faults, create one
`cfg.geometry.faults` entry per fault. The legacy-exact method reproduces the
`fit_curved_planeMM.m` surface fitting sequence. Detailed computed/reference
geometry examples are in [Stage 2-3 Details](#stage-2-3-details-geometry-and-mesh).

Output `geometry/fault_geometry.mat` records the fitted surface, controls,
fault identities, units, and provenance. The figure shows traces, depth
controls, and the fitted surface.

### Stage 3: Triangular Mesh

The default `legacy_exact_main_interpolate` method follows
`main_interpolate.m`: depth-dependent spacing, Delaunay triangulation,
point-source triangle IDs, and triangle neighbors. Mesh density is controlled
by `cfg.mesh.biasL`, `biasW`, `depthRatio`, and `depthKm`.

Output `geometry/triangular_mesh.mat` is required by both TDE and composite
inversions. Inspect `figures/triangular_mesh/triangular_mesh_preview.fig`
before inversion.

### Stage 4: Quick Homogeneous TDE Inversion

The quick inversion uses the initially sampled data, homogeneous-half-space
TDE Greens, the configured smoothing and boundary conditions, and the mesh.
Its purpose is to produce a first-order model for adaptive resampling, not to
replace the later resampled TDE or composite result.

Important controls are under `cfg.quickInversion`: smoothing weight and
strike/dip ratio, zero-slip boundary strength, polarity, and per-data-type
weights. It saves `inversions/quick_result.mat` and predicts `model.grd` for
every configured track under `inversions/quick_model_grids`.

### Stage 5: Model-Based Sampling

This stage builds a new quadtree on each quick `model.grd`, then applies that
tree to the original data and look vectors. Its parameters are deliberately
separate from initial sampling and remain per track:

```matlab
cfg.insar.tracks(1).modelBasedSampling.minCellSizeKm = 2;
cfg.insar.tracks(1).modelBasedSampling.maxCellSizeKm = 2.5;
cfg.insar.tracks(1).modelBasedSampling.rmsThreshold = 50;
cfg.insar.tracks(1).modelBasedSampling.nanFractionMax = 1;
```

To tune only this step after the quick model exists:

```matlab
cfg.workflow.startStep = "model_based_sampling";
cfg.workflow.stopStep = "model_based_sampling";
run workflows/run_quick_inversion.m
```

The comparison figure shows dense observations before sampling and the values
retained by the model-based tree. Outputs are placed under
`data/processed/model_based_sampling/<track_name>`.

## Stage 6: Resampled TDE Inversion

Users who do not need layered Green's functions or the later composite-Green's
inversion can end the workflow at the homogeneous-half-space TDE inversion.
There are two supported ways to do this.

### Run A New Project Through TDE And Stop

This is the normal route when the user supplies new displacement grids and
wants SlipSolve-curve to sample them, construct the geometry and mesh, run the quick
model, resample from that model, then perform the TDE inversion:

```matlab
cfg.workflow.startStep = "subsample_insar";
cfg.workflow.stopStep = "resampled_tde_inversion";
cfg.workflow.pauseAfterStage = true;
cfg.resampledTdeInversion.inputMode = "model_based_sampling";
```

Run `workflows/run_quick_inversion.m` after setting these values. The final
TDE result is saved as `inversions/resampled_tde_result.mat`; its slip and
data/model/residual figures are written to `figures/resampled_tde_inversion`.
No layered-Green's or final-inversion stage is run.

### Run Only An Existing TDE Inversion

When geometry and mesh products already exist, the TDE stage can be run by
itself. `triangular_mesh.mat` must already be present under the project output
directory. For the supplied legacy-reference data, use:

```matlab
cfg.workflow.startStep = "resampled_tde_inversion";
cfg.workflow.stopStep = "resampled_tde_inversion";
cfg.resampledTdeInversion.inputMode = "legacy_samp3_reference";
```

For a previously completed new-workflow run, use the same start/stop values
but set `inputMode` to `"model_based_sampling"`. In that case
`inversions/quick_result.mat` must also exist because it records the
model-based sampled observations.

### Add, Remove, And Control TDE Observations

For a new project, the list `cfg.insar.tracks` in
`config/example_project.m` is the definitive list of observations that enter
the TDE inversion. One configured track produces one model-based sampled
dataset. Add a track by duplicating a complete track entry, give it a unique
`name`, and update its `dataType`, `dataFile`, and `lookEFile`, `lookNFile`,
and `lookUFile` paths. Remove a track by deleting its entry from this list.

Each grid must be co-registered with its three look-vector grids. The sampled
MAT file stores `[x y displacement lookE lookN lookU]`: x/y are metres,
look-vector components are dimensionless, and displacement must use the
configured inversion unit (centimetres in the legacy example). For an AZO
track, the directory containing `dataFile` must also
contain the legacy `sinF.dat` and `cosF.dat` heading files used by the TDE
Green's-function calculation.

The number of observations per track is set independently through that
track's `modelBasedSampling` block:

```matlab
cfg.insar.tracks(1).modelBasedSampling.boundsLonLat = [];
cfg.insar.tracks(1).modelBasedSampling.minCellSizeKm = 0.5;
cfg.insar.tracks(1).modelBasedSampling.maxCellSizeKm = 20;
cfg.insar.tracks(1).modelBasedSampling.rmsThreshold = 1;
```

Smaller cell sizes or a lower `rmsThreshold` generally retain more data;
larger cells, a higher threshold, or a smaller `boundsLonLat` window retain
less. Inspect the model-sampling comparison figure for each track before
continuing. The optional `cfg.resampledTdeInversion.trackWeights` controls the
relative influence of AZO, LOS, and range data types; it does not change the
number of observations.

For the copied legacy-reference mode, the list
`cfg.resampledTdeInversion.legacySamples` is the definitive input list. Remove
an entry to omit that data set, or add an entry with its `name`, `type`,
`greenOption`, and project-local `relativeFile`. AZO reference datasets must
keep their matching `sinF.dat` and `cosF.dat` beside the sample MAT file.

To reproduce the legacy `MMInversionTri.m` TDE inversion input set, use the
copied `samp3` reference files:

```matlab
cfg.workflow.startStep = "resampled_tde_inversion";
cfg.workflow.stopStep = "resampled_tde_inversion";
cfg.resampledTdeInversion.inputMode = "legacy_samp3_reference";
```

This mode uses the same seven already-resampled files as
`make_fault_from_insar_curve.m`: four Sentinel AZO `azo_samp3.mat` tracks
and three ALOS2 LOS `los_samp3.mat` tracks, copied under
`data/raw/legacy_samp3`. It also uses the same legacy smoothness
(`smoothness = 10e-1`), `smoothingMatrix_test(..., ss_ratio=1, ds_ratio=3)`,
bottom/left/right zero-slip constraints, `bounds_new`, and `lsqlin`.

To instead invert the model-based samples produced by the new workflow, switch
the input mode and run from model-based sampling through TDE:

```matlab
cfg.workflow.startStep = "model_based_sampling";
cfg.workflow.stopStep = "resampled_tde_inversion";
cfg.resampledTdeInversion.inputMode = "model_based_sampling";
```

Initial and model-based sampling parameters are intentionally separate. See
[Stage 1](#stage-1-initial-insar-sampling) and
[Stage 5](#stage-5-model-based-sampling) for the active per-track controls.

## Output Ownership And Units

Treat all external source-data and legacy-code directories as read-only.
Generated files are written only beneath the configured project output root.

Quick inversion products are saved under the SlipSolve-curve project:

- `inversions/quick_result.mat`
- `inversions/quick_model_grids/<track_name>/model.grd`
- `data/processed/model_based_sampling/<track_name>/*_samp<iteration>.mat`
- `figures/model_sampling/<track_name>_model_sampling_comparison.fig`
- `figures/model_sampling/<track_name>_model_sampling_comparison.png`
- `inversions/resampled_tde_result.mat`
- `figures/resampled_tde_inversion/*.fig`
- `figures/resampled_tde_inversion/*.png`

The model-based sampling follows the legacy method: build the quadtree on the
quick model grid, apply that tree to the original data and look vectors, and
save `[x y data lookE lookN lookU]`. Local x/y are metres; the legacy Myanmar
displacements are centimetres. The comparison figure shows dense data before
sampling versus data after model-based sampling.

## Stage 7: Layered Green's Functions With EDGRN

SlipSolve-curve does not install or execute EDGRN. Users install and run EDGRN
themselves, then SlipSolve-curve converts the completed fundamental Green's-function
text tables into project-local MAT files using the same fields and reading
method as the legacy `getedgrn.m`.

EDGRN/EDCMP is distributed by Rongjiang Wang's group at GFZ for static
co-seismic deformation in a layered elastic half-space. Obtain the source from
the [GFZ EDGRN/EDCMP tool page](https://www.gfz.de/en/section/physics-of-earthquakes-and-volcanoes/infrastructure/tool-development-lab/)
and follow the `READ.ME` packaged with that distribution. Its documented Unix
build sequence is:

```sh
cd edgrnf77-2.0
make
mv edgrn ..
cd ../edcmpf77-2.0
make
mv edcmp ..
```

Before compiling, review the array limits in the EDGRN include files and make
them large enough for the requested radial-distance and source-depth grids.
The EDGRN input template itself documents the required receiver depth, radial
grid, source-depth grid, output directory, and layered velocity model.

### User-Managed EDGRN Run

Keep the raw EDGRN input and text database inside the SlipSolve-curve project. The
recommended layout is:

```text
SlipSolve-curve/
  data/raw/edgrn/my_earth_model/
    edgrnMM
    edgrnfcts/
      ridge.ss
      ridge.ds
      ridge.cl
  greens/edgrn/my_earth_model/
    edgrn_ss_edgrnMM.mat
    edgrn_ds_edgrnMM.mat
    edgrn_cl_edgrnMM.mat
```

The `data/raw/edgrn` location holds the user-managed raw EDGRN product. The
`greens/edgrn` location is reserved for SlipSolve-curve's converted MAT files.
Do not place generated MAT products in an external source-data directory.

1. Copy an EDGRN input template such as `edgrn.inp` or `edgrnMM` to
   `data/raw/edgrn/my_earth_model`. Set the receiver depth, radial-distance
   range and count, source-depth range and count, output directory, and
   layered Earth model.
2. Create `data/raw/edgrn/my_earth_model/edgrnfcts`. In the supplied input
   format, set the fifth non-comment line to the relative output directory
   and three text database names:

```text
'./edgrnfcts/'  'ridge.ss'  'ridge.ds'  'ridge.cl'
```

   These are the strike-slip (`ss`), dip-slip (`ds`), and CL databases.
3. From the directory containing the input file, run EDGRN with the input on
   standard input:

```sh
/path/to/edgrn < edgrnMM
```

   EDGRN must finish successfully and create all three text files before
   returning to MATLAB. Their names and locations must remain consistent with
   the fifth non-comment input line; relative output paths are relative to the
   input file's directory.

### Convert EDGRN Output In SlipSolve-curve

Point the configuration at the exact input file used for the completed run:

```matlab
cfg.layeredGreens.edgrnInputFile = fullfile(cfg.project.outputRoot, ...
    "data", "raw", "edgrn", "my_earth_model", "edgrnMM");
cfg.layeredGreens.outputDirectory = fullfile(cfg.project.outputRoot, ...
    "greens", "edgrn", "my_earth_model");
cfg.layeredGreens.overwriteExistingMatFiles = false;

cfg.workflow.startStep = "layered_greens";
cfg.workflow.stopStep = "layered_greens";
```

Then run `workflows/run_full_inversion.m`. The workflow reads the existing
triangular mesh, converts the EDGRN SS, DS, and CL tables, and writes:

- `greens/edgrn/<database>/edgrn_ss_<input_name>.mat`
- `greens/edgrn/<database>/edgrn_ds_<input_name>.mat`
- `greens/edgrn/<database>/edgrn_cl_<input_name>.mat`
- `greens/layered_greens.mat` with input/output provenance
- `figures/layered_greens/layered_edgrn_uz_diagnostic.fig` and `.png`

The MAT files contain the legacy variable names: `nr`, `nz`, `r1`, `r2`, `z1`,
`z2`, `zrec0`, `lambda`, `mu`, `uz`, `ur`, `ut`, `ezz`, `err`, `ett`, `ezr`,
`ert`, `etz`, and `duz_dr`. This makes them compatible with the later legacy
layered-Green's calculation routines.

## Stage 8: Full-Resolution Composite XYZ Green's Functions

After `layered_greens` has converted the EDGRN database, the
`composite_greens` stage follows `Comb_Green_Layered_Nodes_xyz.m` and
`comb_green_nodes_xyz.m` to calculate unprojected east, north, and vertical
responses over an adaptive rectilinear XY grid. This stage is expected to take
a very long time and may require substantial RAM and disk space.

The principal grid controls are deliberately exposed in
`config/example_project.m`, in meters, with the current legacy values as
defaults:

```matlab
cfg.compositeGreens.xRange = [-240e3 240e3];
cfg.compositeGreens.yRange = [-350e3 290e3];
cfg.compositeGreens.D1 = 1.3e3;
cfg.compositeGreens.D2 = 2.3e3;
cfg.compositeGreens.M = 9e3;
```

`D1` is the fine spacing inside the fault-trace bounding box expanded by `M`.
`D2` is the spacing outside that box and must be larger than `D1`. Before the
expensive calculation begins, SlipSolve-curve opens and saves an interactive grid
preview. With `pauseAfterGridPreview=true`, an interactive MATLAB session waits
for the user to inspect that figure and press Enter.

To run only this stage after the mesh and layered EDGRN MAT products exist:

```matlab
cfg.workflow.startStep = "composite_greens";
cfg.workflow.stopStep = "composite_greens";
cfg.compositeGreens.pauseAfterGridPreview = true;
cfg.compositeGreens.overwriteExistingFiles = false;
```

Then run `workflows/run_full_inversion.m`. Outputs are written under
`greens/composite_full_xyz`:

- `full_xyz_grid.mat`: row coordinates and adaptive-grid parameters
- `G_e.mat`: east response, variable `G_e`
- `G_n.mat`: north response, variable `G_n`
- `G_u.mat`: vertical response, variable `G_u`
- `greens/composite_greens.mat`: provenance, dimensions, and column mapping
- `figures/composite_greens/full_resolution_grid_preview.fig` and `.png`

For one fault, the numerical sequence and parameter ordering reproduce the
legacy script: special hybrid treatment for the top two node-depth levels,
followed by deeper layered node responses, with all strike-slip columns before
all dip-slip columns. For multiple disconnected faults, each fault is processed
independently and the global order is all faults' strike-slip columns followed
by all faults' dip-slip columns.

## Stage 9: Final Layered Inversion From Full XYZ Greens

The `final_inversion` stage follows `MMInversion_altered.m`, but it no longer
needs a separately computed `layered_green_nodes_new.mat` for every satellite
track. SlipSolve-curve linearly interpolates the shared full-resolution `G_e`, `G_n`,
and `G_u` matrices at every sampled observation and projects them as follows:

- LOS/RNG: `lookE*G_e + lookN*G_n + lookU*G_u`, using columns 4:6 of
  `sampled_insar_data`.
- AZO: `sinF*G_e + cosF*G_n`, using `sinF.dat` and `cosF.dat` beside the
  sample MAT file, or explicit `sinF`/`cosF` track fields.

Any number of datasets may be used. Each `cfg.finalInversion.tracks` entry
needs `name`, `type`, and either `sampleFile` or a `relativeFile` beneath
`cfg.finalInversion.samplesRoot`. For example:

```matlab
cfg.finalInversion.samplesRoot = fullfile(cfg.project.outputRoot, "data", "processed");
cfg.finalInversion.tracks = [
    struct("name", "track_1", "type", "los", ...
           "sampleFile", "/path/to/track_1/los_samp3.mat")
    struct("name", "track_2", "type", "azo", ...
           "sampleFile", "/path/to/track_2/azo_samp3.mat", ...
           "sinF", -0.19, "cosF", 0.982)
];
```

The default configuration uses the four copied Sentinel-1 AZO files and three
copied ALOS-2 LOS files. It retains the altered top-layer geometry,
`smoothingMatrix_laplace2`, `smoothness=7e-1` with the internal `1e-4` scale,
the `6e-4/3e-4/3e-4` bottom/left/right boundaries, polarity `[-1 0 0]`, and
legacy 550 cm / 1000 cm component limits. Each track is divided by its sample
count, matching `calc_weight_insar_error(ones(N,1))`.

To run only this stage after `geometry/triangular_mesh.mat` and the three full
component matrices exist:

```matlab
cfg.workflow.startStep = "final_inversion";
cfg.workflow.stopStep = "final_inversion";
cfg.finalInversion.interpolation.reuseCached = true;
```

Run `workflows/run_full_inversion.m`. The first run creates
`greens/composite_sampled_tracks/sampled_composite_greens.mat`; later changes
to smoothness, bounds, boundaries, or plot styling reuse it. After changing
samples, geometry, grid parameters, or component Greens, set
`cfg.finalInversion.interpolation.overwriteCached=true` once to rebuild it.

Outputs are:

- `inversions/final_result.mat`: inversion diagnostics, predictions, and the
  reconstructed nodal and triangular models.
- `inversions/simple_triangular_model.mat` and `.txt`: ten columns
  `[element_id, longitude_deg, latitude_deg, depth_m, strike_deg, dip_deg,
  strike_slip_m, dip_slip_m, area_m2, shear_modulus_pa]`.
- `figures/final_layered_inversion/layered_nodal_slip.fig` and `.png`, using
  the legacy two-panel strike/dip view.
- `figures/final_layered_inversion/simple_triangular_slip.fig` and `.png`.
- One interactive data/model/residual `.fig` and `.png` per input track.

Plot controls are grouped under `cfg.visualization.finalInversion`, including
colormap, colorbar ranges, XYZ ranges, view, line width, and title visibility.
The two horizontal nodal-slip colorbar positions are independently adjustable
with `strikeColorbarPosition` and `dipColorbarPosition`, using normalized
`[left bottom width height]` figure coordinates.

### Output Slip Polarity

The saved TDE array `tdeResult.slipModel(:,2:3)`, final array
`finalResult.hybridSlipModel(:,4:5)`, and `simple_triangular_model` files all
preserve the legacy Green-function convention:

- Positive strike slip is **left-lateral**; negative strike slip is
  **right-lateral**.
- For the current Myanmar fault orientation, positive dip slip moves the
  west/right side up; negative dip slip moves the east/left side up.

The TDE polarity figure displays these stored signs directly. The legacy-style
final slip figure changes only the strike-slip display sign, so positive
strike-slip colors in that final figure mean right-lateral motion. Its dip-slip
panel is not sign-flipped. Always use the stored convention above when reading
MAT/TXT model products programmatically. For another fault orientation,
determine the dip-slip hanging-wall direction from that mesh's strike/dip
convention rather than assuming west/east.

### Optional Full-Resolution Data/Model/Residual Figures

Set the following option to reproduce the `plotGrdModelDifference2.m` style:

```matlab
cfg.visualization.finalInversion.fullResolutionFitPlots.enabled = true;
```

Each entry in `fullResolutionFitPlots.tracks` identifies the final inversion
track, full model grid, full observed-data grid, look-vector grids for LOS,
colorbar range, map range, and title. The four default AZO entries preserve the
legacy distinction between `azo_ll_low2.grd` (model evaluation) and
`res_deramped_azo_low.grd` (full observed data). The LOS entries use
`los_ll_low.grd` plus `look_e_low.grd`, `look_n_low.grd`, and `look_u_low.grd`.

The numerical sequence follows the legacy routine:

1. Evaluate the solved model on the full composite `G_e/G_n/G_u` grid.
2. Project AZO using `sinF/cosF`; for LOS, interpolate the three look vectors
   to the composite grid before projection.
3. Interpolate the projected model through a fault-side barrier so values are
   not smeared across the rupture.
4. Plot sampled data, sampled model, and the filled full-grid residual using a
   shared, zero-centered color scale.

Interactive and PNG figures are written to
`figures/final_layered_inversion/full_resolution_fits`. When
`saveResidualGrids=true`, raw full-resolution residual grids are written under
`inversions/full_resolution_residuals`, never into the legacy data folders.
This option reloads all three large component matrices after the inversion;
disable it during rapid smoothness/boundary experiments when those figures do
not need to be regenerated.

## Optional Features

The following helpers are separate entry points, not workflow stage names.
They consume saved inversion products and do not rebuild geometry or Greens.

| Optional feature | Entry point or switch | Required saved product |
|---|---|---|
| Full-resolution fit plots | `cfg.visualization.finalInversion.fullResolutionFitPlots.enabled` | Final model and original full grids |
| TDE or composite L-curve | `workflows/run_l_curve.m` | Corresponding saved inversion result |
| TDE or composite shallow slip deficit | `workflows/run_shallow_slip_deficit.m` | Corresponding saved inversion result |
| Independent-data forward model | `workflows/run_forward_model.m` | Final composite inversion result and full XYZ Greens |

### L-Curve Smoothness Helper

`workflows/run_l_curve.m` repeats either the resampled TDE inversion or the
final composite-Green inversion over user-selected smoothness values. It uses
the same Green matrices, weighting, smoothing operators, boundary conditions,
polarity constraints, component limits, and internal smoothness scaling as
the corresponding production inversion. It does not overwrite the normal TDE
or final inversion result.

First run the inversion of interest once so its saved products exist. Then edit
the clearly labeled L-curve section near the end of `config/example_project.m`:

```matlab
cfg.lCurve.inversionType = "composite"; % "composite" or "tde"
cfg.lCurve.smoothnessValues = logspace(-1, 1, 12);
cfg.lCurve.misfitMetric = "raw_sse";
cfg.lCurve.saveSolutions = true;
```

Run:

```matlab
run workflows/run_l_curve.m
```

The horizontal axis is the actual RMS roughness of the solved model,
`RMS(H*m)`, not the requested smoothness value. This corrects a plotting quirk
in `L_curve_altered.m`, which labeled its reciprocal input parameter as model
roughness even though the inversion function returned measured roughness.
The default vertical axis, `"raw_sse"`, reproduces the legacy quantity named
`RMS_misfit`; despite that name, it is the unweighted residual sum of squares.
It can instead be `"weighted_sse"`, `"weighted_l2"`, or `"raw_rms"`.

No preferred or optimal point is marked. Each point can instead be labeled
with its requested smoothness value. Plot scale, color, marker, label, and
title controls are under `cfg.visualization.lCurve`. To reuse the exact legacy
trial parameterization, assign `1./legacyRoughnessValues` to
`cfg.lCurve.smoothnessValues`, because the old loop used
`smoothness=1/roughness(i)`.

Outputs are saved as interactive `.fig` and `.png` files under
`figures/l_curve`, with all trial metrics in `.mat` and `.csv` files under
`inversions/l_curve`. The TDE helper reuses `resampled_tde_result.mat`; the
composite helper reuses `final_result.mat` and the sampled composite-Green
cache. Both can be computationally expensive, so a short list of two or three
values is useful for an initial test before running a dense sweep.

### Shallow Slip Deficit Helper

`workflows/run_shallow_slip_deficit.m` is optional post-processing for the
saved composite/nodal model, the resampled TDE model, or both. It does not
rerun or alter either inversion. Select the model products near the end of
`config/example_project.m`:

```matlab
cfg.shallowSlipDeficit.modelTypes = ["composite" "tde"];
cfg.shallowSlipDeficit.component = "magnitude";
cfg.shallowSlipDeficit.referenceDepthRangeKm = [0 Inf];
cfg.shallowSlipDeficit.shallowDepthMaxKm = 5;
```

Run:

```matlab
run workflows/run_shallow_slip_deficit.m
```

The calculation follows the active method in
`Shallow_Slip_Deficit_Segmented_New.m` while making its implied deficit
explicit. At each depth, total slip amplitude is
`sqrt(strike_slip^2 + dip_slip^2)`. Each segment is normalized by the maximum
mean profile slip within `referenceDepthRangeKm`, and the profile deficit is:

```text
deficit_fraction(depth) = 1 - mean_slip(depth) / reference_slip
```

Two scalar summaries are reported for every segment:

- `shallowest_deficit_fraction` uses the shallowest model depth available.
  This is 0 km for the current nodal model and the shallowest triangle-centroid
  depth for TDE.
- `mean_shallow_deficit_fraction` compares the support-weighted mean slip from
  the surface through `shallowDepthMaxKm` with the same reference slip.

The two model types use geometry-appropriate depth profiles:

- **Composite/nodal:** arithmetic mean at each native nodal layer, matching
  the attached legacy script.
- **TDE:** triangle-centroid depth levels with area-weighted mean slip. This
  avoids bias from the unequal patch areas in the depth-dependent triangular
  mesh. Set `depthGrouping="fixed_bins"` to use a custom `binEdgesKm` or
  `binWidthKm` instead.

`component` can also be `"strike_magnitude"` or `"dip_magnitude"`. These use
absolute component amplitudes, so the stored left/right-lateral polarity does
not cancel the spatial mean.

#### Segment Configuration

Each segment can select fault IDs and local model-coordinate bounds in km.
Ranges use `(lower, upper]`, matching the non-overlapping inequalities in the
legacy south/middle/north calculation:

```matlab
cfg.shallowSlipDeficit.segments = [
    struct("name","south", "faultIds",[], ...
           "xRangeKm",[-Inf Inf], "yRangeKm",[-Inf -100])
    struct("name","middle","faultIds",[], ...
           "xRangeKm",[-Inf Inf], "yRangeKm",[-100 50])
    struct("name","north", "faultIds",[], ...
           "xRangeKm",[-Inf Inf], "yRangeKm",[50 Inf])
];
```

An empty `faultIds` includes every fault. For disconnected faults, use values
such as `faultIds=1` or `faultIds=[1 3]` to calculate separate profiles. The
default configuration analyzes one `whole_fault` segment.

For each model type, the helper creates:

- An interactive normalized-slip/deficit `.fig` and PNG under
  `figures/shallow_slip_deficit`.
- A complete depth-profile CSV and a scalar-summary CSV under
  `inversions/shallow_slip_deficit`.
- `inversions/shallow_slip_deficit/shallow_slip_deficit.mat`, containing all
  profiles, definitions, provenance, and figure paths.

The first analysis of an older TDE result extracts a compact
`tde_model_for_postprocessing.mat` cache from the large inversion file.
Subsequent changes to segments, components, depth grouping, or plot styling
reuse that compact file. Future TDE inversions write it automatically.

### Independent-Data Forward Modeling

`workflows/run_forward_model.m` projects the saved final composite model onto
observations that were not used in the inversion. It replaces the legacy
data-specific `G_new*u` step in `forwardmodel.m` with the same procedure used
by the final inversion:

1. Bilinearly interpolate the shared full-resolution `G_e`, `G_n`, and `G_u`
   matrices at every independent sample coordinate.
2. Project those component Greens into AZO, LOS/range, east, north, or up.
3. Multiply the resulting sampled Green matrix by the saved final slip vector.
4. Evaluate the solved model over the complete XYZ-Green grid.
5. Interpolate that projected field through the fault-side barrier onto the
   original full-resolution independent-data grid.
6. Plot full-resolution data, forward prediction, and `data - model` residual.

The helper is forward calculation only: it never adds the independent data to
the inversion and never changes the saved slip model. It currently applies to
the final composite/nodal model because the full XYZ matrix columns correspond
to that model parameterization.

Configure any number of tracks under `cfg.forwardModel.tracks`. Each sample
MAT file must contain `sampled_insar_data` with x/y in metres and displacement
in centimetres. LOS/range samples also use look-vector columns 4:6. AZO tracks
use explicit `sinF`/`cosF` values or `sinF.dat` and `cosF.dat` beside the sample
file. Direct `east`, `north`, and `up` types need only columns 1:3.

For the full-resolution comparison, every track also needs `dataGridFile`.
LOS/range tracks additionally need `lookEFile`, `lookNFile`, and `lookUFile`.
These grids are used only for projection and plotting; they are not added to
the inversion. `cfg.forwardModel.fullResolution.faultTraceFile` must contain
local x/y fault coordinates in metres in the full-Green coordinate system.

The bundled configuration uses the copied independent Sentinel-2 AZO sample:

```matlab
cfg.forwardModel.tracks(1).name = "Sentinel2_independent_AZO";
cfg.forwardModel.tracks(1).type = "azo";
cfg.forwardModel.tracks(1).relativeFile = fullfile("SEN2", "azo_samp3.mat");
cfg.forwardModel.tracks(1).dataGridFile = fullfile("SEN2", "azo_ll.grd");
```

After `final_inversion` exists, run:

```matlab
run workflows/run_forward_model.m
```

The first run interpolates and caches one sampled composite Green matrix per
independent track under `greens/forward_model_tracks`. Later slip-model changes
reuse that cache because the Green matrix depends on the geometry, sample
coordinates, and projection, not on slip. Set
`cfg.forwardModel.interpolation.overwriteCached=true` after changing those
inputs. The helper also caches the solved full XYZ component response and each
track's interpolated full model/residual grids. Unchanged reruns therefore
reload the caches instead of rereading the multi-gigabyte Green matrices or
repeating the fault-side interpolation.

Outputs include:

- Per-track prediction MAT and CSV files under `inversions/forward_model`.
- Per-track full-resolution model and residual GRD files under
  `inversions/forward_model` when `fullResolution.saveGrids=true`.
- `inversions/forward_model/forward_model_result.mat` with provenance and
  diagnostics.
- Interactive `.fig` and PNG data/model/residual plots under
  `figures/forward_model`.

Color maps, zero-centered ranges, map extent, and titles are controlled under
`cfg.visualization.forwardModel`. Set
`cfg.forwardModel.fullResolution.enabled=false` only when a quicker
sampled-point comparison is preferred.


## Stage 2-3 Details: Geometry And Mesh

SlipSolve-curve can either compute curved geometry from user inputs or, for legacy QA/debugging, load already saved legacy geometry products.

### Compute From Trace And Segments

This is the normal user workflow. Provide a surface trace, deep segment files, and one dip per segment:

```matlab
cfg.geometry.mode = "legacy_fit_curved_planeMM";
cfg.geometry.surfaceTraceFiles = "/path/to/trace.txt";
cfg.geometry.segmentFiles = [
    "/path/to/Segment_001.txt"
    "/path/to/Segment_002.txt"
];
cfg.geometry.segmentDipDegrees = [75 80];

cfg.geometry.legacyExact.useReferenceProducts = false;
cfg.mesh.method = "legacy_exact_main_interpolate";
cfg.mesh.legacyExact.useReferenceProducts = false;
```

This follows the legacy `fit_curved_planeMM.m` and `main_interpolate.m` recipe inside the streamlined workflow: dense trace interpolation, projected deep control lines, `gridFitInterpolate`, depth-dependent vertex spacing, `delaunayTriangulation`, point-source triangle IDs, and triangle neighbors.

### Use Existing Legacy Products

If a user already has saved legacy geometry products and wants to reuse them exactly, enable reference-product mode:

```matlab
cfg.geometry.legacyExact.useReferenceProducts = true;
cfg.geometry.legacyExact.referenceMeshFile = "/path/to/CurveMesh1_dense_new.mat";
cfg.geometry.legacyExact.referencePointFile = "/path/to/CurvePoint1_dense_new.mat";

cfg.mesh.legacyExact.useReferenceProducts = true;
cfg.mesh.legacyExact.referenceGeometryFile = "/path/to/geometry1_dense_new.mat";
```

Expected variables:

- `CurveMesh*.mat`: `x1m`, `y1m`, `d1m`
- `CurvePoint*.mat`: `x1p`, `y1p`, `d1p`
- `geometry*.mat`: `ID1`, `n1`, `DT1`, `Vx1`

Reference-product mode is mainly for reproducing an existing legacy run or validating the new wrapper against old outputs. For new projects, prefer computing from traces and segments.

### Multiple Disconnected Faults

Disconnected faults should be configured as separate entries in `cfg.geometry.faults`. SlipSolve-curve fits each fault independently and then combines the active fault meshes for downstream inversion:

```matlab
cfg.geometry.faults(1).name = "fault_a";
cfg.geometry.faults(1).surfaceTraceFiles = "/path/to/fault_a_trace.txt";
cfg.geometry.faults(1).mode = "legacy_fit_curved_planeMM";
cfg.geometry.faults(1).segmentFiles = [
    "/path/to/fault_a_segment_001.txt"
    "/path/to/fault_a_segment_002.txt"
];
cfg.geometry.faults(1).segmentDipDegrees = [75 80];

cfg.geometry.faults(2).name = "fault_b";
cfg.geometry.faults(2).surfaceTraceFiles = "/path/to/fault_b_trace.txt";
cfg.geometry.faults(2).mode = "legacy_fit_curved_planeMM";
cfg.geometry.faults(2).segmentFiles = [
    "/path/to/fault_b_segment_001.txt"
    "/path/to/fault_b_segment_002.txt"
];
cfg.geometry.faults(2).segmentDipDegrees = [65 70];

cfg.mesh.activeFaults = "all";
```

To run only selected faults:

```matlab
cfg.mesh.activeFaults = ["fault_a" "fault_b"];
```

For multiple faults, the recommended path is to compute each fault from its own trace and segment controls. The current reference-product shortcut is intended for single legacy-product sets unless the project is extended with explicit per-fault reference geometry files.

## Layout

```text
config/              User-editable project configuration
data/raw/            Raw input data copied or linked by users
data/processed/      Subsampled/resampled observations
geometry/            Fault geometry and triangular meshes
greens/              EDGRN and composite Green's function products
inversions/          Quick and final inversion results
src/+slipsolve/      MATLAB package code
workflows/           User-facing run scripts
docs/                Design notes and implementation roadmap
external/            Internalized MATLAB dependencies and external-tool notes
examples/            Example projects
tests/               MATLAB test notes and future automated tests
```
