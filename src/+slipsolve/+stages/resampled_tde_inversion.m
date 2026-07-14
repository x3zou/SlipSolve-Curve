function tdeResult = resampled_tde_inversion(cfg, mesh, sampling)
%RESAMPLED_TDE_INVERSION Run legacy-style TDE inversion on resampled data.

p = slipsolve.project.paths(cfg);
slipsolve.project.ensure_directories(cfg);
add_legacy_paths(cfg);
invCfg = cfg.resampledTdeInversion;

if nargin < 2 || isempty(mesh)
    loaded = load(p.triangularMesh, "mesh");
    mesh = loaded.mesh;
end
if uses_legacy_samp3(invCfg)
    sampling = legacy_samp3_records(cfg);
elseif nargin < 3 || isempty(sampling)
    loaded = load(p.quickResult, "quickResult");
    if ~isfield(loaded.quickResult, "modelSampling") || isempty(loaded.quickResult.modelSampling)
        error("SlipSolve:MissingInput", ...
            "quickResult.modelSampling is missing. Run model_based_sampling before resampled_tde_inversion.");
    end
    sampling = loaded.quickResult.modelSampling;
end

[slipModel, points] = legacy_slip_model_from_mesh(mesh);
[greenRaw, greenWeighted, dataRaw, dataWeighted, trackResults] = build_resampled_greens(cfg, slipModel, points, sampling);

[smoothMatrix, smoothRows] = smoothingMatrix_test(slipModel(:, 7:9), ...
    "ss_ratio", invCfg.smoothingStrikeRatio, "ds_ratio", invCfg.smoothingDipRatio);

[boundaryMatrix, boundaryData, boundaryInfo] = build_boundary_constraints(invCfg, slipModel, points);

lambda = invCfg.smoothingWeight;
systemMatrix = [greenWeighted; smoothMatrix .* lambda ./ smoothRows; boundaryMatrix];
systemData = [dataWeighted; zeros(smoothRows, 1); boundaryData];

[lb, ub] = legacy_bounds(invCfg, slipModel);
if isfield(invCfg.constraints, "customLowerBounds") && ~isempty(invCfg.constraints.customLowerBounds)
    lb = invCfg.constraints.customLowerBounds(:);
end
if isfield(invCfg.constraints, "customUpperBounds") && ~isempty(invCfg.constraints.customUpperBounds)
    ub = invCfg.constraints.customUpperBounds(:);
end

options = optimset("LargeScale", "on", "DiffMaxChange", 1e-1, "DiffMinChange", 1e-12, ...
    "TolCon", 1e-12, "TolFun", 1e-12, "TolPCG", 1e-12, "TolX", 1e-12, ...
    "MaxIter", 1e9, "MaxPCGIter", 1e9, "Display", "off");
fprintf("Solving resampled TDE system: %d equations x %d parameters\n", ...
    size(systemMatrix, 1), size(systemMatrix, 2));
[slipVector, resnorm, residual, exitflag] = lsqlin(systemMatrix, double(systemData), ...
    [], [], [], [], lb, ub, [], options);

nPatch = size(slipModel, 1);
slipModel(:, 2) = slipVector(1:nPatch);
slipModel(:, 3) = slipVector(nPatch+1:2*nPatch);
prediction = greenRaw * slipVector;
roughnessVector = smoothMatrix * slipVector;

tdeResult = struct();
tdeResult.method = "legacy_mm_inversion_tri_resampled_tde";
tdeResult.legacyReference = "MMInversionTri.m / make_fault_from_insar_curve.m";
tdeResult.slipModel = slipModel;
tdeResult.points = points;
tdeResult.slipVector = slipVector;
tdeResult.prediction = prediction;
tdeResult.data = dataRaw;
tdeResult.residual = prediction - dataRaw;
tdeResult.resnorm = resnorm;
tdeResult.lsqlinResidual = residual;
tdeResult.exitflag = exitflag;
tdeResult.rms0 = sum(dataRaw .^ 2);
tdeResult.rms = sum((prediction - dataRaw) .^ 2);
tdeResult.varianceReduction = 100 .* (tdeResult.rms0 - tdeResult.rms) ./ tdeResult.rms0;
tdeResult.modelRoughness = sqrt(sum(roughnessVector .^ 2) ./ numel(roughnessVector));
tdeResult.smoothingRows = smoothRows;
tdeResult.boundaryInfo = boundaryInfo;
tdeResult.tracks = add_track_predictions(trackResults, prediction);
tdeResult.config = invCfg;
tdeResult.figures = slipsolve.plotting.plot_resampled_tde_inversion(cfg, tdeResult, mesh);

save(p.resampledTdeResult, "tdeResult", "-v7.3");
sourceInfo = dir(p.resampledTdeResult);
tdeModel = struct();
tdeModel.slipModel = tdeResult.slipModel;
tdeModel.points = tdeResult.points;
tdeModel.method = tdeResult.method;
tdeModel.sourceProduct = string(p.resampledTdeResult);
tdeModel.sourceBytes = sourceInfo.bytes;
tdeModel.sourceDatenum = sourceInfo.datenum;
tdeModel.createdAt = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");
save(p.tdePostprocessingModel, "tdeModel", "-v7.3");

fprintf("resampled TDE rms misfit (dat., res.) = %e %e (%f%%)\n", ...
    tdeResult.rms0, tdeResult.rms, tdeResult.varianceReduction);
fprintf("resampled TDE resnorm, resid. = %e %e\n", sqrt(resnorm), mean(residual));
fprintf("resampled TDE exitflag is %d\n", exitflag);

end

function add_legacy_paths(cfg)
pathsToAdd = strings(0, 1);
legacyProjectRoot = fileparts(string(cfg.paths.legacyRoot));
if strlength(legacyProjectRoot) > 0 && exist(legacyProjectRoot, "dir") == 7
    addpath(genpath(char(legacyProjectRoot)));
end
pathsToAdd(end+1) = string(cfg.paths.legacyRoot);
pathsToAdd(end+1) = fullfile(string(cfg.paths.legacyRoot), "TDE_Greens");
pathsToAdd(end+1) = fullfile(fileparts(string(cfg.paths.legacySamplingRoot)), "Greens");
pathsToAdd(end+1) = fullfile(fileparts(string(cfg.paths.legacySamplingRoot)), "inversion");
pathsToAdd(end+1) = string(cfg.paths.legacySamplingRoot);
pathsToAdd(end+1) = string(cfg.paths.legacyOtherFuncRoot);
pathsToAdd(end+1) = string(get_field(cfg.paths, "legacyGeodeticFunctionsRoot", ""));
pathsToAdd(end+1) = string(cfg.paths.geodeticRoot);

for k = 1:numel(pathsToAdd)
    if strlength(pathsToAdd(k)) > 0 && exist(pathsToAdd(k), "dir") == 7
        addpath(char(pathsToAdd(k)));
    end
end
end

function [slipModel, points] = legacy_slip_model_from_mesh(mesh)
nPatch = size(mesh.triangles, 1);
slipModel = zeros(nPatch, 9);
slipModel(:, 1) = mesh.faultId(:);
slipModel(:, 4:6) = mesh.triangles;
slipModel(:, 7:9) = legacy_neighbors(mesh);

nodeFaultId = node_fault_ids(mesh);
points = [mesh.nodes .* 1e3, nodeFaultId];
end

function neighborsOut = legacy_neighbors(mesh)
neighborsOut = nan(size(mesh.triangles, 1), 3);
triOffset = 0;
for k = 1:numel(mesh.sections)
    section = mesh.sections(k);
    nTri = size(section.triangles, 1);
    localNeighbors = section.neighbors;
    localNeighbors(localNeighbors == 0) = NaN;
    valid = ~isnan(localNeighbors);
    localNeighbors(valid) = localNeighbors(valid) + triOffset;
    neighborsOut(triOffset + (1:nTri), :) = localNeighbors;
    triOffset = triOffset + nTri;
end
end

function nodeFaultId = node_fault_ids(mesh)
nodeFaultId = zeros(size(mesh.nodes, 1), 1);
nodeOffset = 0;
for k = 1:numel(mesh.sections)
    section = mesh.sections(k);
    nNode = size(section.nodes, 1);
    nodeFaultId(nodeOffset + (1:nNode)) = section.faultIndex;
    nodeOffset = nodeOffset + nNode;
end
end

function [greenRaw, greenWeighted, dataRaw, dataWeighted, trackResults] = build_resampled_greens(cfg, slipModel, points, sampling)
sampling = sampling(:);
trackCells = cell(numel(sampling), 1);
for k = 1:numel(sampling)
    fprintf("Building resampled TDE Green %d/%d: %s\n", ...
        k, numel(sampling), char(string(sampling(k).name)));
    if isfield(sampling, "greenOption")
        trackCells{k} = build_one_legacy_samp3_green(cfg, slipModel, points, sampling(k));
    else
        trackDef = find_track_definition(cfg, sampling(k));
        trackCells{k} = build_one_resampled_green(cfg, slipModel, points, sampling(k), trackDef);
    end
end
trackResults = vertcat(trackCells{:});

greenRaw = vertcat(trackResults.greenRaw);
dataRaw = vertcat(trackResults.data);
greenWeighted = vertcat(trackResults.greenWeighted);
dataWeighted = vertcat(trackResults.dataWeighted);
end

function trackDef = find_track_definition(cfg, sampleRecord)
trackDefs = cfg.insar.tracks(:);
sampleName = string(sampleRecord.name);
for k = 1:numel(trackDefs)
    if string(trackDefs(k).name) == sampleName || string(trackDefs(k).dataFile) == string(sampleRecord.sourceGrid)
        trackDef = trackDefs(k);
        return
    end
end
error("SlipSolve:MissingInput", "No cfg.insar.tracks entry matches sample %s.", sampleName);
end

function out = build_one_resampled_green(cfg, slipModel, points, sampleRecord, trackDef)
loaded = load(sampleRecord.sampleFile, "sampled_insar_data");
sampledData = double(loaded.sampled_insar_data);
trackType = lower(string(trackDef.dataType));
if trackType == "azo"
    [sinF, cosF] = read_azo_heading(trackDef.dataFile);
    greenRaw = calc_green_insar_curve_AZO(slipModel, points, sampledData, "sinf", sinF, "cosf", cosF);
else
    greenRaw = calc_green_insar_curve(slipModel, points, sampledData);
end

out = struct();
out.name = string(trackDef.name);
out.type = string(trackDef.dataType);
out.sourceGrid = string(trackDef.dataFile);
out.sampleFile = string(sampleRecord.sampleFile);
out.sampleCount = size(sampledData, 1);
out.greenRaw = greenRaw;
out.data = sampledData(:, 3);
out.x = sampledData(:, 1);
out.y = sampledData(:, 2);
out.weight = track_weight(cfg, trackType) .* ones(size(sampledData, 1), 1);
out.greenWeighted = greenRaw .* out.weight;
out.dataWeighted = out.data .* out.weight;
out.rowStart = 0;
out.rowEnd = 0;
end

function out = build_one_legacy_samp3_green(cfg, slipModel, points, sampleRecord)
sampleFile = char(sampleRecord.sampleFile);
[greenRaw, greenWeightedBase, dataRaw, dataWeightedBase] = build_green_function( ...
    slipModel, sampleFile, char(sampleRecord.greenOption), "noramp", ...
    char(cfg.resampledTdeInversion.modelType), "pts", points);

loaded = load(sampleFile, "sampled_insar_data");
sampledData = double(loaded.sampled_insar_data);
trackType = lower(string(sampleRecord.type));
weight = track_weight(cfg, trackType);

out = struct();
out.name = string(sampleRecord.name);
out.type = string(sampleRecord.type);
out.sourceGrid = string(sampleRecord.sourceGrid);
out.sampleFile = string(sampleRecord.sampleFile);
out.sampleCount = size(sampledData, 1);
out.greenRaw = greenRaw;
out.greenWeighted = greenWeightedBase .* weight;
out.data = dataRaw;
out.dataWeighted = dataWeightedBase .* weight;
out.x = sampledData(:, 1);
out.y = sampledData(:, 2);
out.weight = weight .* ones(size(sampledData, 1), 1);
out.rowStart = 0;
out.rowEnd = 0;
end

function [sinF, cosF] = read_azo_heading(sourceGrid)
trackDir = fileparts(char(sourceGrid));
sinFile = fullfile(trackDir, "sinF.dat");
cosFile = fullfile(trackDir, "cosF.dat");
if exist(sinFile, "file") ~= 2 || exist(cosFile, "file") ~= 2
    error("SlipSolve:MissingInput", "AZO track requires sinF.dat and cosF.dat in %s", trackDir);
end
sinF = load(sinFile);
cosF = load(cosFile);
end

function weight = track_weight(cfg, trackType)
weight = 1;
if isfield(cfg.resampledTdeInversion, "trackWeights")
    weights = cfg.resampledTdeInversion.trackWeights;
    if trackType == "azo" && isfield(weights, "azo")
        weight = weights.azo;
    elseif trackType == "los" && isfield(weights, "los")
        weight = weights.los;
    elseif trackType == "rng" && isfield(weights, "rng")
        weight = weights.rng;
    end
end
end

function tf = uses_legacy_samp3(invCfg)
tf = isfield(invCfg, "inputMode") && string(invCfg.inputMode) == "legacy_samp3_reference";
end

function records = legacy_samp3_records(cfg)
root = string(cfg.resampledTdeInversion.legacySamplesRoot);
samples = cfg.resampledTdeInversion.legacySamples(:);
records = repmat(struct("name", "", "type", "", "greenOption", "", ...
    "sourceGrid", "", "sampleFile", ""), numel(samples), 1);
for k = 1:numel(samples)
    sampleFile = fullfile(root, string(samples(k).relativeFile));
    if exist(sampleFile, "file") ~= 2
        error("SlipSolve:MissingInput", "Legacy reference sample file is missing: %s", sampleFile);
    end
    records(k).name = string(samples(k).name);
    records(k).type = string(samples(k).type);
    records(k).greenOption = string(samples(k).greenOption);
    records(k).sourceGrid = sampleFile;
    records(k).sampleFile = sampleFile;
end
end

function [boundaryMatrix, boundaryData, boundaryInfo] = build_boundary_constraints(invCfg, slipModel, points)
boundaryMatrix = [];
boundaryData = [];
boundaryInfo = struct("name", strings(0, 1), "enabled", false(0, 1), "ratio", zeros(0, 1), "rows", zeros(0, 1));

defs = [
    struct("field", "bottom", "legacyName", "bottom")
    struct("field", "top", "legacyName", "top")
    struct("field", "left", "legacyName", "left")
    struct("field", "right", "legacyName", "right")
];

for k = 1:numel(defs)
    bc = get_field(invCfg.boundaryConditions, defs(k).field, struct());
    enabled = get_field(bc, "enabled", false);
    ratio = get_field(bc, "ratio", 5e-4);
    rows = 0;
    if enabled
        if defs(k).field == "top"
            error("SlipSolve:UnsupportedBoundaryCondition", ...
                "Top/surface zero-slip is not supported by the legacy zero_slip_boundary_curve helper. " + ...
                "MMInversionTri constrains bottom, left, and right while leaving the surface free.");
        end
        [W, d] = zero_slip_boundary_curve(slipModel, points, defs(k).legacyName, ratio);
        boundaryMatrix = [boundaryMatrix; W]; %#ok<AGROW>
        boundaryData = [boundaryData; d]; %#ok<AGROW>
        rows = size(W, 1);
    end
    boundaryInfo.name(end+1, 1) = defs(k).field;
    boundaryInfo.enabled(end+1, 1) = enabled;
    boundaryInfo.ratio(end+1, 1) = ratio;
    boundaryInfo.rows(end+1, 1) = rows;
end
end

function [lb, ub] = legacy_bounds(invCfg, slipModel)
nFault = max(slipModel(:, 1));
tSm = zeros(1, nFault + 1);
for k = 1:nFault
    tSm(k + 1) = nnz(slipModel(:, 1) == k);
end
if get_field(invCfg.constraints, "nonNegativeSlip", false)
    con = [1 1 0];
else
    con = get_field(invCfg.constraints, "polarity", [-1 0 0]);
end
[lb, ub] = bounds_new(nFault, 2, tSm, 0, con, slipModel);
end

function trackResults = add_track_predictions(trackResults, prediction)
rowStart = 1;
for k = 1:numel(trackResults)
    rowEnd = rowStart + trackResults(k).sampleCount - 1;
    trackResults(k).rowStart = rowStart;
    trackResults(k).rowEnd = rowEnd;
    trackResults(k).prediction = prediction(rowStart:rowEnd);
    trackResults(k).residual = trackResults(k).prediction - trackResults(k).data;
    rowStart = rowEnd + 1;
end
end

function value = get_field(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
