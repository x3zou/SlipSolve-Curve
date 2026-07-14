# SlipSolve Roadmap

## Phase 1: Contracts

- Document the inputs and outputs of each legacy script.
- Create stable config fields for paths, quadtree settings, geometry, mesh controls, Green's functions, and inversion parameters.
- Add validation before expensive computation starts.

## Phase 2: Wrappers

- Wrap `sample_insar_data.m` as `slipsolve.stages.subsample_insar`.
- Wrap curved geometry tools as `slipsolve.stages.build_fault_geometry`.
- Wrap triangular mesh generation as `slipsolve.stages.build_triangular_mesh`.
- Wrap `MMInversionTri.m` as `slipsolve.stages.quick_inversion`.

## Phase 3: Layered Green's Functions

- Wrap EDGRN input generation and execution.
- Standardize paths for EDGRN database outputs.
- Wrap composite Green's function assembly.
- Add indexing that maps each subsampled observation to the corresponding composite Green's function rows.

## Phase 4: Final Inversion

- Wrap the altered inversion code.
- Save solution, residuals, roughness, regularization metadata, and plots.
- Add rerun controls so expensive products are reused unless config changes.

## Phase 5: User Experience

- Add plots after every major stage.
- Add example projects.
- Add MATLAB tests for config validation and geometry conventions.
- Write a short user guide with a minimal end-to-end example.

