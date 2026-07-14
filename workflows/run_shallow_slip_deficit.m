%RUN_SHALLOW_SLIP_DEFICIT Analyze shallow slip for TDE and/or composite models.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(repoRoot, "src"));

cfg = slipsolve.config.load(fullfile(repoRoot, "config", "example_project.m"));
cfg = slipsolve.config.apply_defaults(cfg);
slipsolve.config.validate(cfg, "quick");
shallowSlipDeficitResult = slipsolve.analysis.run_shallow_slip_deficit(cfg);

disp("Shallow slip deficit analysis completed.");
