%RUN_QUICK_INVERSION Run the SlipSolve quick inversion workflow.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(repoRoot, "src"));

cfg = slipsolve.config.load(fullfile(repoRoot, "config", "example_project.m"));
cfg = slipsolve.config.apply_defaults(cfg);
slipsolve.config.validate(cfg, "quick");
results = slipsolve.workflow.run(cfg, "quick");

disp("Quick inversion workflow completed.");
