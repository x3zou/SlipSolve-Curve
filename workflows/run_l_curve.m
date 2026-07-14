%RUN_L_CURVE Sweep smoothness for the TDE or composite inversion.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(repoRoot, "src"));

cfg = slipsolve.config.load(fullfile(repoRoot, "config", "example_project.m"));
cfg = slipsolve.config.apply_defaults(cfg);
slipsolve.config.validate(cfg, "full");
lCurveResult = slipsolve.analysis.run_l_curve(cfg);

disp("L-curve smoothness sweep completed.");
