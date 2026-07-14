function validate(cfg, mode)
%VALIDATE Validate a SlipSolve configuration before running workflows.

if nargin < 2
    mode = "full";
end

requiredTopLevel = ["project", "paths", "insar", "geometry", "mesh", "quickInversion"];
for i = 1:numel(requiredTopLevel)
    if ~isfield(cfg, requiredTopLevel(i))
        error("SlipSolve:MissingConfigField", "Missing required config section: %s", requiredTopLevel(i));
    end
end

if strlength(string(cfg.project.name)) == 0
    error("SlipSolve:InvalidConfig", "cfg.project.name must be set.");
end

if strlength(string(cfg.paths.sampleListFile)) == 0
    warning("SlipSolve:UnsetInput", "cfg.paths.sampleListFile is not set yet.");
elseif exist(cfg.paths.sampleListFile, "file") ~= 2
    warning("SlipSolve:MissingInput", "Sample-list file does not exist: %s", cfg.paths.sampleListFile);
end

if strlength(string(cfg.paths.faultTrace)) == 0 && isempty(cfg.geometry.surfaceTraceFiles)
    warning("SlipSolve:UnsetInput", "No fault trace path is set yet.");
end

if mode == "full"
    if ~isfield(cfg, "layeredGreens") || ~isfield(cfg, "compositeGreens") || ~isfield(cfg, "finalInversion")
        error("SlipSolve:MissingConfigField", "Full mode requires layeredGreens, compositeGreens, and finalInversion sections.");
    end

    if strlength(string(cfg.paths.edgrnExecutable)) == 0
        warning("SlipSolve:UnsetInput", "cfg.paths.edgrnExecutable is not set yet.");
    end
end

end
