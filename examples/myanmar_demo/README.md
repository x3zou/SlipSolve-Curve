# Myanmar Demo

This is the bundled final composite-inversion example.

The script reads `config/example_project.m` but forces the workflow to start
and stop at `final_inversion`. It reuses the project-local triangular mesh,
seven copied `samp3` observation files, precomputed `G_e/G_n/G_u`, and the
sampled composite-Green cache. It does not run quadtree sampling, geometry,
TDE, EDGRN conversion, or composite-Green construction.

The full Green matrices, sampled-Green cache, and source `.grd` files are
intentionally excluded from Git history. Download the Google Drive folder
`SlipSolve-curve-example-data`, then merge its `data` and `greens`
folders into the repository root without changing any names or subfolders.
Run `verify_example_data` before the inversion to confirm all 31 files and
their expected sizes are present.

From the SlipSolve-curve project root, run:

```matlab
addpath examples/myanmar_demo
verify_example_data
run examples/myanmar_demo/run_example.m
```

Generated results are written under `inversions` and `figures`. Edit the
`finalInversion` and `visualization.finalInversion` sections of
`config/example_project.m` to change smoothness, constraints, data tracks, or
plot styling. Normal sampled-data and slip figures are enabled. The optional
full-resolution fit figures are disabled by default because they reload and
evaluate the three large component matrices; enable
`cfg.visualization.finalInversion.fullResolutionFitPlots.enabled` when those
figures are needed. MATLAB Desktop opens the figures interactively; headless
or batch runs render them offscreen.
