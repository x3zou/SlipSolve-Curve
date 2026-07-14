function lCurveResult = run_l_curve(cfg)
%RUN_L_CURVE Sweep user-facing smoothness for TDE or composite inversion.

if ~isfield(cfg, "lCurve") || ~get_field(cfg.lCurve, "enabled", false)
    error("SlipSolve:DisabledHelper", "Set cfg.lCurve.enabled = true before running the L-curve helper.");
end
slipsolve.project.ensure_directories(cfg);
add_legacy_paths(cfg);
p = slipsolve.project.paths(cfg);
settings = cfg.lCurve;
smoothness = double(settings.smoothnessValues(:));
if numel(smoothness) < 2 || any(~isfinite(smoothness) | smoothness <= 0)
    error("SlipSolve:InvalidLCurveConfig", ...
        "cfg.lCurve.smoothnessValues must contain at least two positive finite values.");
end

inversionType = normalize_inversion_type(settings.inversionType);
switch inversionType
    case "composite"
        problem = prepare_composite_problem(cfg, p);
    case "tde"
        problem = prepare_tde_problem(cfg, p);
end

nTrial = numel(smoothness);
nParameter = size(problem.greenWeighted, 2);
solutions = zeros(nParameter, nTrial);
weightedSse = zeros(nTrial, 1);
weightedL2 = zeros(nTrial, 1);
rawSse = zeros(nTrial, 1);
rawRms = zeros(nTrial, 1);
varianceReduction = zeros(nTrial, 1);
modelRoughness = zeros(nTrial, 1);
modelRoughnessL2 = zeros(nTrial, 1);
resnorm = zeros(nTrial, 1);
exitflag = zeros(nTrial, 1);

nData = size(problem.greenWeighted, 1);
nSmooth = problem.smoothingRows;
smoothRows = nData+(1:nSmooth);
systemMatrix = [problem.greenWeighted; problem.smoothingMatrix; problem.boundaryMatrix];
systemData = [problem.dataWeighted; zeros(nSmooth, 1); problem.boundaryData];
problem.greenWeighted = [];

options = optimset("LargeScale", "on", "DiffMaxChange", 1e-1, ...
    "DiffMinChange", 1e-12, "TolCon", 1e-12, "TolFun", 1e-12, ...
    "TolPCG", 1e-12, "TolX", 1e-12, "MaxIter", 1e9, ...
    "MaxPCGIter", 1e9, "Display", "off");
startingModel = [];
for k = 1:nTrial
    lambda = smoothness(k)*problem.internalSmoothnessScale;
    systemMatrix(smoothRows, :) = problem.smoothingMatrix.*lambda./nSmooth;
    fprintf("L-curve %s trial %d/%d: smoothness = %.8g\n", ...
        inversionType, k, nTrial, smoothness(k));
    [solution, resnorm(k), ~, exitflag(k)] = lsqlin( ...
        systemMatrix, double(systemData), [], [], [], [], ...
        problem.lowerBounds, problem.upperBounds, startingModel, options);
    if get_field(settings, "useWarmStart", false)
        startingModel = solution;
    end
    solutions(:, k) = solution;
    weightedResidual = systemMatrix(1:nData, :)*solution-problem.dataWeighted;
    weightedSse(k) = sum(weightedResidual.^2);
    weightedL2(k) = norm(weightedResidual);
    [rawSse(k), observationCount] = raw_data_sse(problem, solution);
    rawRms(k) = sqrt(rawSse(k)/observationCount);
    varianceReduction(k) = 100*(problem.rawDataEnergy-rawSse(k))/problem.rawDataEnergy;
    roughnessVector = problem.smoothingMatrix*solution;
    modelRoughness(k) = sqrt(sum(roughnessVector.^2)/numel(roughnessVector));
    modelRoughnessL2(k) = norm(roughnessVector);
    fprintf("  roughness = %.8e, weighted SSE = %.8e, raw RMS = %.8e\n", ...
        modelRoughness(k), weightedSse(k), rawRms(k));
end

[plottedMisfit, misfitLabel] = select_misfit(settings.misfitMetric, ...
    weightedSse, weightedL2, rawSse, rawRms);
baseName = inversionType+"_l_curve";
matFile = fullfile(p.lCurveDirectory, baseName+".mat");
tableFile = fullfile(p.lCurveDirectory, baseName+".csv");

lCurveResult = struct();
lCurveResult.inversionType = inversionType;
lCurveResult.legacyReference = "L_curve_altered.m";
lCurveResult.smoothnessValues = smoothness;
lCurveResult.internalSmoothnessValues = smoothness*problem.internalSmoothnessScale;
lCurveResult.modelRoughness = modelRoughness;
lCurveResult.modelRoughnessL2 = modelRoughnessL2;
lCurveResult.weightedSse = weightedSse;
lCurveResult.weightedL2 = weightedL2;
lCurveResult.rawSse = rawSse;
lCurveResult.rawRms = rawRms;
lCurveResult.varianceReduction = varianceReduction;
lCurveResult.plottedMisfit = plottedMisfit;
lCurveResult.misfitMetric = string(settings.misfitMetric);
lCurveResult.misfitLabel = misfitLabel;
lCurveResult.resnorm = resnorm;
lCurveResult.exitflag = exitflag;
lCurveResult.parameterCount = nParameter;
lCurveResult.observationCount = problem.observationCount;
if get_field(settings, "saveSolutions", true)
    lCurveResult.solutions = solutions;
else
    lCurveResult.solutions = [];
end
lCurveResult.sourceProduct = problem.sourceProduct;
lCurveResult.figureFiles = slipsolve.plotting.plot_l_curve(cfg, lCurveResult);
lCurveResult.matFile = matFile;
lCurveResult.tableFile = tableFile;
lCurveResult.createdAt = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");

metrics = table(smoothness, lCurveResult.internalSmoothnessValues, ...
    modelRoughness, modelRoughnessL2, weightedSse, weightedL2, ...
    rawSse, rawRms, varianceReduction, resnorm, exitflag, ...
    'VariableNames', {'smoothness', 'internal_smoothness', ...
    'model_roughness_rms', 'model_roughness_l2', 'weighted_sse', ...
    'weighted_l2', 'raw_sse', 'raw_rms', 'variance_reduction_percent', ...
    'solver_resnorm', 'exitflag'});
writetable(metrics, tableFile);
save(matFile, "lCurveResult", "-v7.3");
fprintf("Saved %s L-curve result: %s\n", inversionType, matFile);
end

function problem = prepare_composite_problem(cfg, p)
if exist(p.finalResult, "file") ~= 2
    error("SlipSolve:MissingStageProduct", ...
        "Run final_inversion once before the composite L-curve helper: %s", p.finalResult);
end
if exist(p.sampledCompositeGreens, "file") ~= 2
    error("SlipSolve:MissingStageProduct", ...
        "Composite sampled-Green cache is missing: %s", p.sampledCompositeGreens);
end
loadedResult = load(p.finalResult, "finalResult");
result = loadedResult.finalResult;
loadedGreen = load(p.sampledCompositeGreens, "G_raw");
greenRaw = loadedGreen.G_raw;
model = result.hybridSlipModel;
tracks = result.tracks(:);
dataRaw = result.data(:);
greenWeighted = greenRaw;
dataWeighted = dataRaw;
for k = 1:numel(tracks)
    rows = tracks(k).rowStart:tracks(k).rowEnd;
    greenWeighted(rows, :) = greenWeighted(rows, :)*tracks(k).weight;
    dataWeighted(rows) = dataWeighted(rows)*tracks(k).weight;
end
[H, hRows] = smoothingMatrix_laplace2(model(:, 9:15), model(:, 6:8), ...
    "ss_ratio", cfg.finalInversion.smoothingStrikeRatio, ...
    "ds_ratio", cfg.finalInversion.smoothingDipRatio);
[boundaryMatrix, boundaryData] = composite_boundaries(cfg.finalInversion, model);
[lowerBounds, upperBounds] = composite_bounds(cfg.finalInversion, model);

problem = struct();
problem.greenWeighted = greenWeighted;
problem.dataWeighted = dataWeighted;
problem.greenRaw = greenRaw;
problem.dataRaw = dataRaw;
problem.rawTracks = struct([]);
problem.smoothingMatrix = H;
problem.smoothingRows = hRows;
problem.boundaryMatrix = boundaryMatrix;
problem.boundaryData = boundaryData;
problem.lowerBounds = lowerBounds;
problem.upperBounds = upperBounds;
problem.internalSmoothnessScale = cfg.finalInversion.smoothingInternalScale;
problem.rawDataEnergy = sum(dataRaw.^2);
problem.observationCount = numel(dataRaw);
problem.sourceProduct = string(p.finalResult);
end

function problem = prepare_tde_problem(cfg, p)
if exist(p.resampledTdeResult, "file") ~= 2
    error("SlipSolve:MissingStageProduct", ...
        "Run resampled_tde_inversion once before the TDE L-curve helper: %s", ...
        p.resampledTdeResult);
end
loaded = load(p.resampledTdeResult, "tdeResult");
result = loaded.tdeResult;
tracks = result.tracks(:);
greenWeighted = vertcat(tracks.greenWeighted);
dataWeighted = vertcat(tracks.dataWeighted);
for k = 1:numel(tracks)
    tracks(k).greenWeighted = [];
    tracks(k).dataWeighted = [];
    tracks(k).prediction = [];
    tracks(k).residual = [];
end
model = result.slipModel;
[H, hRows] = smoothingMatrix_test(model(:, 7:9), ...
    "ss_ratio", cfg.resampledTdeInversion.smoothingStrikeRatio, ...
    "ds_ratio", cfg.resampledTdeInversion.smoothingDipRatio);
[boundaryMatrix, boundaryData] = tde_boundaries( ...
    cfg.resampledTdeInversion, model, result.points);
[lowerBounds, upperBounds] = tde_bounds(cfg.resampledTdeInversion, model);
dataRaw = vertcat(tracks.data);

problem = struct();
problem.greenWeighted = greenWeighted;
problem.dataWeighted = dataWeighted;
problem.greenRaw = [];
problem.dataRaw = dataRaw;
problem.rawTracks = tracks;
problem.smoothingMatrix = H;
problem.smoothingRows = hRows;
problem.boundaryMatrix = boundaryMatrix;
problem.boundaryData = boundaryData;
problem.lowerBounds = lowerBounds;
problem.upperBounds = upperBounds;
problem.internalSmoothnessScale = 1;
problem.rawDataEnergy = sum(dataRaw.^2);
problem.observationCount = numel(dataRaw);
problem.sourceProduct = string(p.resampledTdeResult);
end

function [sse, count] = raw_data_sse(problem, solution)
if ~isempty(problem.greenRaw)
    residual = problem.greenRaw*solution-problem.dataRaw;
    count = numel(residual);
    sse = sum(residual.^2);
    return
else
    sse = 0;
    count = 0;
    for k = 1:numel(problem.rawTracks)
        track = problem.rawTracks(k);
        residual = track.greenRaw*solution-track.data;
        sse = sse+sum(residual.^2);
        count = count+numel(residual);
    end
    return
end
end

function [matrix, data] = composite_boundaries(invCfg, model)
matrix = sparse(0, 2*size(model, 1));
data = zeros(0, 1);
faultIds = unique(model(:, 1)).';
names = ["bottom", "top", "left", "right"];
for name = names
    settings = invCfg.boundaryConditions.(name);
    if ~get_field(settings, "enabled", false), continue; end
    for faultId = faultIds
        if name == "bottom"
            option = max(model(model(:, 1) == faultId, 3));
        elseif name == "top"
            option = min(model(model(:, 1) == faultId, 3));
        else
            option = char(name);
        end
        [block, target] = zero_slip_boundary(model, faultId, option, settings.ratio);
        matrix = [matrix; sparse(block)]; %#ok<AGROW>
        data = [data; target]; %#ok<AGROW>
    end
end
end

function [lower, upper] = composite_bounds(invCfg, model)
n = size(model, 1);
constraints = invCfg.constraints;
lower = [-constraints.maxStrikeSlipCm*ones(n, 1); ...
    -constraints.maxDipSlipCm*ones(n, 1)];
upper = [constraints.maxStrikeSlipCm*ones(n, 1); ...
    constraints.maxDipSlipCm*ones(n, 1)];
polarity = constraints.polarity;
faultIds = unique(model(:, 1)).';
for index = 1:numel(faultIds)
    con = polarity;
    if index > 1, con = [polarity(1), -polarity(2), polarity(3)]; end
    rows = find(model(:, 1) == faultIds(index));
    if con(1) > 0, lower(rows) = 0; end
    if con(1) < 0, upper(rows) = 0; end
    if con(2) > 0, lower(rows+n) = 0; end
    if con(2) < 0, upper(rows+n) = 0; end
end
if ~isempty(constraints.customLowerBounds), lower = constraints.customLowerBounds(:); end
if ~isempty(constraints.customUpperBounds), upper = constraints.customUpperBounds(:); end
end

function [matrix, data] = tde_boundaries(invCfg, model, points)
matrix = [];
data = [];
for name = ["bottom", "left", "right"]
    settings = invCfg.boundaryConditions.(name);
    if ~get_field(settings, "enabled", false), continue; end
    [block, target] = zero_slip_boundary_curve(model, points, char(name), settings.ratio);
    matrix = [matrix; block]; %#ok<AGROW>
    data = [data; target]; %#ok<AGROW>
end
if get_field(invCfg.boundaryConditions.top, "enabled", false)
    error("SlipSolve:UnsupportedBoundaryCondition", ...
        "Legacy TDE L-curve does not support a top zero-slip boundary.");
end
end

function [lower, upper] = tde_bounds(invCfg, model)
nFault = max(model(:, 1));
counts = zeros(1, nFault+1);
for k = 1:nFault, counts(k+1) = nnz(model(:, 1) == k); end
if get_field(invCfg.constraints, "nonNegativeSlip", false)
    polarity = [1 1 0];
else
    polarity = invCfg.constraints.polarity;
end
[lower, upper] = bounds_new(nFault, 2, counts, 0, polarity, model);
if ~isempty(invCfg.constraints.customLowerBounds), lower = invCfg.constraints.customLowerBounds(:); end
if ~isempty(invCfg.constraints.customUpperBounds), upper = invCfg.constraints.customUpperBounds(:); end
end

function [values, label] = select_misfit(metric, weightedSse, weightedL2, rawSse, rawRms)
switch lower(string(metric))
    case "weighted_sse"
        values = weightedSse;
        label = "Weighted data misfit (sum of squares)";
    case "weighted_l2"
        values = weightedL2;
        label = "Weighted data misfit (L2 norm)";
    case "raw_sse"
        values = rawSse;
        label = "Unweighted data misfit (sum of squares)";
    case "raw_rms"
        values = rawRms;
        label = "Unweighted data RMS";
    otherwise
        error("SlipSolve:InvalidLCurveConfig", ...
            "Unknown cfg.lCurve.misfitMetric: %s", metric);
end
end

function type = normalize_inversion_type(value)
switch lower(string(value))
    case {"composite", "final", "final_inversion", "layered"}
        type = "composite";
    case {"tde", "resampled_tde", "resampled_tde_inversion"}
        type = "tde";
    otherwise
        error("SlipSolve:InvalidLCurveConfig", ...
            "cfg.lCurve.inversionType must be 'composite' or 'tde'.");
end
end

function add_legacy_paths(cfg)
root = string(cfg.paths.legacyRoot);
projectRoot = fileparts(root);
paths = [root; fullfile(projectRoot, "resolution_test"); ...
    fullfile(projectRoot, "geodetic_inversion-master", "inversion"); ...
    fullfile(projectRoot, "geodetic_inversion-master", "geometry"); ...
    fullfile(projectRoot, "geodetic_inversion-master", "OtherFunc")];
for k = 1:numel(paths)
    if exist(paths(k), "dir") == 7, addpath(char(paths(k))); end
end
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
