%RUN_FULL_INVERSION Run the complete SlipSolve inversion workflow.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(repoRoot, "src"));

cfg = slipsolve.config.load(fullfile(repoRoot, "config", "example_project.m"));
cfg = slipsolve.config.apply_defaults(cfg);
slipsolve.config.validate(cfg, "full");
results = slipsolve.workflow.run(cfg, "full");

disp("Full inversion workflow completed.");
