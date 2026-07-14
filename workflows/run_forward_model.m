%RUN_FORWARD_MODEL Project the final composite model onto independent data.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(repoRoot, "src"));

cfg = slipsolve.config.load(fullfile(repoRoot, "config", "example_project.m"));
cfg = slipsolve.config.apply_defaults(cfg);
slipsolve.config.validate(cfg, "full");
forwardModelResult = slipsolve.analysis.run_forward_model(cfg);

disp("Independent-data forward modeling completed.");
