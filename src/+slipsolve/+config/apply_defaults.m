function cfg = apply_defaults(cfg)
%APPLY_DEFAULTS Fill optional configuration fields.

defaults = struct();
defaults.project.root = pwd;
defaults.project.outputRoot = defaults.project.root;
defaults.runtime.overwrite = false;
defaults.runtime.saveIntermediate = true;
defaults.runtime.verbose = true;
defaults.workflow.startStep = "subsample_insar";
defaults.workflow.stopStep = "model_based_sampling";
defaults.workflow.pauseAfterStage = true;

defaults.paths.legacyRoot = "";
defaults.paths.legacyGeodeticFunctionsRoot = "";
defaults.paths.rawInsar = "";
defaults.paths.faultTrace = "";
defaults.paths.edgrnExecutable = "";

cfg = merge_structs(defaults, cfg);

end

function out = merge_structs(base, override)
out = base;
names = fieldnames(override);

for i = 1:numel(names)
    name = names{i};
    if isfield(out, name) && isstruct(out.(name)) && isstruct(override.(name))
        out.(name) = merge_structs(out.(name), override.(name));
    else
        out.(name) = override.(name);
    end
end

end
