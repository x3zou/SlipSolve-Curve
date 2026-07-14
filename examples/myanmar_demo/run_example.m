%RUN_EXAMPLE Run only the bundled final composite inversion example.

exampleDirectory = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(exampleDirectory));
addpath(fullfile(repoRoot, "src"));
addpath(exampleDirectory);

% Stop early with explicit download/copy instructions when bundle files are
% absent or incomplete.
verify_example_data(repoRoot);

cfg = slipsolve.config.load(fullfile(repoRoot, "config", "example_project.m"));
cfg = slipsolve.config.apply_defaults(cfg);

% Keep the example fast: all expensive upstream products are already present.
cfg.workflow.startStep = "final_inversion";
cfg.workflow.stopStep = "final_inversion";

% Desktop runs open interactive figures; headless/batch runs render offscreen.
if ~usejava("desktop")
    cfg.visualization.visible = false;
end

slipsolve.config.validate(cfg, "full");
results = slipsolve.workflow.run(cfg, "full");

disp("Bundled final composite inversion example completed.");
