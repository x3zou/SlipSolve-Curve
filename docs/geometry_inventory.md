# Geometry Inventory

This note records the current Myanmar curved-geometry workflow and the first SlipSolve geometry contract.

## Legacy Scripts

### `fit_curved_planeMM.m`

Bundled location: `external/legacy_matlab/Inversion/curved_geometry/fit_curved_planeMM.m`

Role: builds a curved fault surface from a surface trace plus projected deep control lines.

Current hard-coded behavior:

- Originally added its source checkout recursively to the MATLAB path; SlipSolve-curve now adds only the bundled dependency directories it needs.
- Uses origin `lonc = 96.05`, `latc = 20.75`, `ref_lon = 96.05`.
- Reads the configured surface trace (the example trace is bundled under `data/raw/example_geometry/`).
- Uses seven segment files from `segment_inversion_new/Segment_00*.txt`.
- Uses one dip angle per segment.
- Projects each segment to 25 km depth with `projectSegment3D`.
- Uses `gridFitInterpolate` / `gridfit` to interpolate a smooth surface.
- Can save `CurveMesh1_dense_new.mat` and `CurvePoint1_dense_new.mat`.
- Produces 3D diagnostic plots.

### `main_interpolate.m`

Bundled location: `external/legacy_matlab/Inversion/curved_geometry/main_interpolate.m`

Role: takes the interpolated curved surface grids and creates a triangular tessellation.

Current hard-coded behavior:

- Loads `CurveMesh1_dense_new.mat` and `CurvePoint1_dense_new.mat`.
- Builds a depth-dependent vertex distribution with `depth_dependent_point`.
- Uses `delaunayTriangulation` in a projected coordinate plane.
- Saves `geometry1_dense_new.mat` when `save_output = 1`.
- Produces a 3D triangle preview.

## SlipSolve Geometry Contract

The current reusable stage is:

```matlab
faultGeometry = slipsolve.stages.build_fault_geometry(cfg);
```

It currently:

- Reads one or more user-provided independent faults through `cfg.geometry.faults`.
- Fits disconnected faults independently so no artificial surface is interpolated between them.
- Converts lon/lat traces to local coordinates using `cfg.insar.origin`.
- Builds depth controls using either:
  - per-segment files and per-segment dips, or
  - a uniform dip angle.
- Fits a smooth curved surface with `gridfit`, using `cfg.geometry.surfaceFitSmoothness`.
- Can reproduce the legacy `fit_curved_planeMM.m` recipe through `cfg.geometry.mode = "legacy_fit_curved_planeMM"`.
- Can optionally load existing legacy `CurveMesh*.mat` and `CurvePoint*.mat` products for validation or exact reuse.
- Saves `geometry/fault_geometry.mat`.
- Saves interactive and PNG previews under `figures/fault_geometry/`.

The current surface grid is the geometry QA product that should be inspected before triangular meshing.

## SlipSolve Triangular Mesh Contract

The reusable meshing stage is:

```matlab
mesh = slipsolve.stages.build_triangular_mesh(cfg, faultGeometry);
```

It currently:

- Uses `cfg.mesh.method = "legacy_exact_main_interpolate"` for the Myanmar legacy-compatible path.
- Reproduces the useful behavior of `main_interpolate.m`: depth-dependent point spacing, `delaunayTriangulation`, point-source triangle IDs, and triangle neighbors.
- Supports multiple disconnected faults by meshing each `faultGeometry.faults(k)` independently and combining active faults into one `mesh` struct.
- Can optionally load an existing `geometry*.mat` reference product for validation or exact reuse.
- Saves `geometry/triangular_mesh.mat`.
- Saves interactive and PNG previews under `figures/triangular_mesh/`.
