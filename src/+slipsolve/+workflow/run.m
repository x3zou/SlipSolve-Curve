function results = run(cfg, mode)
%RUN Execute SlipSolve stages from cfg.workflow.startStep to stopStep.

if nargin < 2 || strlength(string(mode)) == 0
    mode = "full";
end
mode = string(mode);

p = slipsolve.project.paths(cfg);
slipsolve.project.ensure_directories(cfg);

steps = workflow_steps(mode);
startStep = get_workflow_field(cfg, "startStep", steps(1));
stopStep = get_workflow_field(cfg, "stopStep", steps(end));
startIndex = find_step_index(steps, startStep);
stopIndex = find_step_index(steps, stopStep);
if startIndex > stopIndex
    error("SlipSolve:InvalidWorkflow", ...
        "cfg.workflow.startStep (%s) must come before cfg.workflow.stopStep (%s).", ...
        startStep, stopStep);
end

results = struct();
results.cfg = cfg;
results.steps = steps(startIndex:stopIndex);

if startIndex > 1
    results = load_required_inputs(results, p, steps(startIndex));
end

for i = startIndex:stopIndex
    step = steps(i);
    fprintf("\n=== SlipSolve step: %s ===\n", step);
    [results, cfg] = run_one_step(results, cfg, step);
    pause_after_stage(cfg, step);
end

end

function steps = workflow_steps(mode)
quickSteps = ["subsample_insar", "fault_geometry", "triangular_mesh", ...
    "quick_inversion", "model_based_sampling", "resampled_tde_inversion"];
fullSteps = [quickSteps, "layered_greens", "composite_greens", "final_inversion"];

switch mode
    case "quick"
        steps = quickSteps;
    case "full"
        steps = fullSteps;
    otherwise
        error("SlipSolve:InvalidWorkflow", "Unknown workflow mode: %s", mode);
end
end

function [results, cfg] = run_one_step(results, cfg, step)
switch step
    case "subsample_insar"
        results.insarSub = slipsolve.stages.subsample_insar(cfg);

    case "fault_geometry"
        if ~isfield(results, "insarSub")
            results = load_insar_sub(results, slipsolve.project.paths(cfg));
        end
        results.faultGeometry = slipsolve.stages.build_fault_geometry(cfg, results.insarSub);

    case "triangular_mesh"
        if ~isfield(results, "faultGeometry")
            results = load_fault_geometry(results, slipsolve.project.paths(cfg));
        end
        results.mesh = slipsolve.stages.build_triangular_mesh(cfg, results.faultGeometry);

    case "quick_inversion"
        if ~isfield(results, "insarSub")
            results = load_insar_sub(results, slipsolve.project.paths(cfg));
        end
        if ~isfield(results, "mesh")
            results = load_mesh(results, slipsolve.project.paths(cfg));
        end
        cfgQuick = cfg;
        if isfield(cfgQuick.quickInversion, "modelSampling")
            cfgQuick.quickInversion.modelSampling.enabled = false;
        end
        results.quickResult = slipsolve.stages.quick_inversion(cfgQuick, results.mesh, results.insarSub);

    case "model_based_sampling"
        if ~isfield(results, "quickResult")
            results = load_quick_result(results, slipsolve.project.paths(cfg));
        end
        [results.modelSampling, results.quickResult] = ...
            slipsolve.stages.model_based_sampling(cfg, results.quickResult);

    case "resampled_tde_inversion"
        if ~isfield(results, "mesh")
            results = load_mesh(results, slipsolve.project.paths(cfg));
        end
        if needs_model_based_sampling(cfg) && ~isfield(results, "quickResult")
            results = load_quick_result(results, slipsolve.project.paths(cfg));
        end
        if isfield(results, "quickResult") && isfield(results.quickResult, "modelSampling")
            sampling = results.quickResult.modelSampling;
        else
            sampling = [];
        end
        results.resampledTdeResult = slipsolve.stages.resampled_tde_inversion(cfg, results.mesh, sampling);

    case "layered_greens"
        if ~isfield(results, "mesh")
            results = load_mesh(results, slipsolve.project.paths(cfg));
        end
        if ~isfield(results, "resampledObs")
            results.resampledObs = struct();
        end
        results.layeredGreens = slipsolve.stages.build_layered_greens(cfg, results.mesh, results.resampledObs);

    case "composite_greens"
        if ~isfield(results, "layeredGreens")
            loaded = load_required_file(slipsolve.project.paths(cfg).layeredGreens, ...
                "layeredGreens", "layered_greens");
            results.layeredGreens = loaded.layeredGreens;
        end
        if ~isfield(results, "mesh")
            results = load_mesh(results, slipsolve.project.paths(cfg));
        end
        if ~isfield(results, "resampledObs")
            results.resampledObs = struct();
        end
        results.compositeGreens = slipsolve.stages.build_composite_greens( ...
            cfg, results.layeredGreens, results.mesh, results.resampledObs);

    case "final_inversion"
        if ~isfield(results, "mesh")
            results = load_mesh(results, slipsolve.project.paths(cfg));
        end
        results.finalResult = slipsolve.stages.final_inversion(cfg, results.mesh);

    otherwise
        error("SlipSolve:InvalidWorkflow", "Unknown workflow step: %s", step);
end
end

function results = load_required_inputs(results, p, startStep)
switch startStep
    case "fault_geometry"
        results = load_insar_sub(results, p);
    case "triangular_mesh"
        results = load_fault_geometry(results, p);
    case "quick_inversion"
        results = load_insar_sub(results, p);
        results = load_mesh(results, p);
    case "model_based_sampling"
        results = load_quick_result(results, p);
    case "resampled_tde_inversion"
        results = load_mesh(results, p);
        if needs_model_based_sampling(results.cfg)
            results = load_quick_result(results, p);
        end
    case "layered_greens"
        results = load_mesh(results, p);
    case "composite_greens"
        loaded = load_required_file(p.layeredGreens, "layeredGreens", "layered_greens");
        results.layeredGreens = loaded.layeredGreens;
        results = load_mesh(results, p);
    case "final_inversion"
        results = load_mesh(results, p);
end
end

function tf = needs_model_based_sampling(cfg)
tf = true;
if isfield(cfg, "resampledTdeInversion") && isfield(cfg.resampledTdeInversion, "inputMode")
    tf = string(cfg.resampledTdeInversion.inputMode) ~= "legacy_samp3_reference";
end
end

function results = load_insar_sub(results, p)
if ~isfield(results, "insarSub")
    loaded = load_required_file(p.insarSubsampled, "insarSub", "subsample_insar");
    results.insarSub = loaded.insarSub;
end
end

function results = load_fault_geometry(results, p)
if ~isfield(results, "faultGeometry")
    loaded = load_required_file(p.faultGeometry, "faultGeometry", "fault_geometry");
    results.faultGeometry = loaded.faultGeometry;
end
end

function results = load_mesh(results, p)
if ~isfield(results, "mesh")
    loaded = load_required_file(p.triangularMesh, "mesh", "triangular_mesh");
    results.mesh = loaded.mesh;
end
end

function results = load_quick_result(results, p)
if ~isfield(results, "quickResult")
    loaded = load_required_file(p.quickResult, "quickResult", "quick_inversion");
    results.quickResult = loaded.quickResult;
end
end

function loaded = load_required_file(filePath, variableName, stageName)
if exist(filePath, "file") ~= 2
    error("SlipSolve:MissingStageProduct", ...
        "Cannot start here because %s is missing. Run stage '%s' first.", ...
        filePath, stageName);
end
loaded = load(filePath, variableName);
end

function index = find_step_index(steps, name)
name = string(name);
index = find(steps == name, 1);
if isempty(index)
    error("SlipSolve:InvalidWorkflow", ...
        "Unknown workflow step '%s'. Valid steps are: %s", name, strjoin(steps, ", "));
end
end

function value = get_workflow_field(cfg, fieldName, defaultValue)
if isfield(cfg, "workflow") && isfield(cfg.workflow, fieldName) && ~isempty(cfg.workflow.(fieldName))
    value = string(cfg.workflow.(fieldName));
else
    value = string(defaultValue);
end
end

function pause_after_stage(cfg, step)
if ~isfield(cfg, "workflow") || ~isfield(cfg.workflow, "pauseAfterStage") || ~cfg.workflow.pauseAfterStage
    drawnow;
    return
end

drawnow;
if is_interactive_matlab()
    fprintf("Finished '%s'. Inspect the figures and edit config if needed.\n", step);
    input("Press Enter to continue to the next configured step...", "s");
end
end

function tf = is_interactive_matlab()
tf = usejava("desktop");
try
    tf = tf && desktop("-inuse");
catch
end
end
