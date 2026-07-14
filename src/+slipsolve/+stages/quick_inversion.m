function quickResult = quick_inversion(cfg, mesh, insarSub)
%QUICK_INVERSION Run quick homogeneous TDE inversion.

p = slipsolve.project.paths(cfg);
slipsolve.project.ensure_directories(cfg);
add_legacy_paths(cfg);

if nargin < 2 || isempty(mesh)
    loaded = load(p.triangularMesh, "mesh");
    mesh = loaded.mesh;
end
if nargin < 3 || isempty(insarSub)
    loaded = load(p.insarSubsampled, "insarSub");
    insarSub = loaded.insarSub;
end

[slipModel, points] = legacy_slip_model_from_mesh(mesh);
[greenRaw, greenWeighted, dataRaw, dataWeighted, trackResults] = build_quick_greens(cfg, slipModel, points, insarSub);
[smoothMatrix, smoothRows] = smoothingMatrix_test(slipModel(:, 7:9), ...
    "ss_ratio", 1, "ds_ratio", cfg.quickInversion.smoothingDipRatio);

boundaryRatio = cfg.quickInversion.boundaryConstraintRatio;
[Wb, db] = zero_slip_boundary_curve(slipModel, points, "bottom", boundaryRatio);
[Wl, dl] = zero_slip_boundary_curve(slipModel, points, "left", boundaryRatio);
[Wr, dr] = zero_slip_boundary_curve(slipModel, points, "right", boundaryRatio);

lambda = cfg.quickInversion.smoothingWeight;
systemMatrix = [greenWeighted; smoothMatrix .* lambda ./ smoothRows; Wb; Wl; Wr];
systemData = [dataWeighted; zeros(smoothRows, 1); db; dl; dr];

[lb, ub] = legacy_bounds(cfg, slipModel);
options = optimset("LargeScale", "on", "DiffMaxChange", 1e-1, "DiffMinChange", 1e-12, ...
    "TolCon", 1e-12, "TolFun", 1e-12, "TolPCG", 1e-12, "TolX", 1e-12, ...
    "MaxIter", 1e9, "MaxPCGIter", 1e9, "Display", "off");
[slipVector, resnorm, residual, exitflag] = lsqlin(systemMatrix, double(systemData), ...
    [], [], [], [], lb, ub, [], options);

nPatch = size(slipModel, 1);
slipModel(:, 2) = slipVector(1:nPatch);
slipModel(:, 3) = slipVector(nPatch+1:2*nPatch);
prediction = greenRaw * slipVector;
roughnessVector = smoothMatrix * slipVector;

quickResult = struct();
quickResult.method = "legacy_homogeneous_tde";
quickResult.slipModel = slipModel;
quickResult.points = points;
quickResult.slipVector = slipVector;
quickResult.prediction = prediction;
quickResult.data = dataRaw;
quickResult.residual = prediction - dataRaw;
quickResult.resnorm = resnorm;
quickResult.lsqlinResidual = residual;
quickResult.exitflag = exitflag;
quickResult.rms0 = sum(dataRaw .^ 2);
quickResult.rms = sum((prediction - dataRaw) .^ 2);
quickResult.varianceReduction = 100 .* (quickResult.rms0 - quickResult.rms) ./ quickResult.rms0;
quickResult.modelRoughness = sqrt(sum(roughnessVector .^ 2) ./ numel(roughnessVector));
quickResult.tracks = add_track_predictions(trackResults, prediction);
quickResult.figures = slipsolve.plotting.plot_quick_inversion(cfg, quickResult, mesh);

if isfield(cfg.quickInversion, "modelGrid") && cfg.quickInversion.modelGrid.enabled
    quickResult.modelGrids = save_quick_model_grids(cfg, quickResult);
else
    quickResult.modelGrids = struct([]);
end

if isfield(cfg.quickInversion, "modelSampling") && cfg.quickInversion.modelSampling.enabled
    [quickResult.modelSampling, quickResult] = slipsolve.stages.model_based_sampling(cfg, quickResult);
else
    quickResult.modelSampling = struct();
end

save(p.quickResult, "quickResult", "-v7.3");

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
pathsToAdd(end+1) = string(cfg.paths.geodeticRoot);
pathsToAdd(end+1) = string(get_field(cfg.paths, "legacyGeodeticFunctionsRoot", ""));

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
    faultIndex = section.faultIndex;
    nodeFaultId(nodeOffset + (1:nNode)) = faultIndex;
    nodeOffset = nodeOffset + nNode;
end
end

function [greenRaw, greenWeighted, dataRaw, dataWeighted, trackResults] = build_quick_greens(cfg, slipModel, points, insarSub)
tracks = insarSub.tracks(:);
trackCells = cell(numel(tracks), 1);
for k = 1:numel(tracks)
    trackCells{k} = build_one_track_green(cfg, slipModel, points, tracks(k));
end
trackResults = vertcat(trackCells{:});

greenRaw = vertcat(trackResults.greenRaw);
dataRaw = vertcat(trackResults.data);
weights = vertcat(trackResults.weight);
greenWeighted = greenRaw .* weights;
dataWeighted = dataRaw .* weights;
end

function out = build_one_track_green(cfg, slipModel, points, track)
sampledData = double([track.x, track.y, track.los, track.lookE, track.lookN, track.lookU]);
trackType = lower(string(track.type));
if trackType == "azo"
    [sinF, cosF] = read_azo_heading(track.sourceGrid);
    greenRaw = calc_green_insar_curve_AZO(slipModel, points, sampledData, "sinf", sinF, "cosf", cosF);
else
    greenRaw = calc_green_insar_curve(slipModel, points, sampledData);
end

out = struct();
out.name = string(track.name);
out.type = string(track.type);
out.sourceGrid = string(track.sourceGrid);
out.sampleCount = size(sampledData, 1);
out.greenRaw = greenRaw;
out.data = sampledData(:, 3);
out.weight = track_weight(cfg, trackType) .* ones(size(sampledData, 1), 1);
out.rowStart = 0;
out.rowEnd = 0;
end

function weight = track_weight(cfg, trackType)
weight = 1;
if trackType == "azo" && isfield(cfg.quickInversion.trackWeights, "azo")
    weight = cfg.quickInversion.trackWeights.azo;
elseif trackType == "los" && isfield(cfg.quickInversion.trackWeights, "los")
    weight = cfg.quickInversion.trackWeights.los;
elseif trackType == "rng" && isfield(cfg.quickInversion.trackWeights, "rng")
    weight = cfg.quickInversion.trackWeights.rng;
end
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

function [lb, ub] = legacy_bounds(cfg, slipModel)
nFault = max(slipModel(:, 1));
tSm = zeros(1, nFault + 1);
for k = 1:nFault
    tSm(k + 1) = nnz(slipModel(:, 1) == k);
end
if cfg.quickInversion.nonNegativeSlip
    con = [1 1 0];
else
    con = cfg.quickInversion.signConstraint;
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

function modelGrids = save_quick_model_grids(cfg, quickResult)
trackDefs = model_grid_tracks(cfg);
modelGrids = repmat(struct("name", "", "sourceGrid", "", "modelFile", "", ...
    "outputDirectory", "", "pointCount", 0), numel(trackDefs), 1);
for k = 1:numel(trackDefs)
    modelGrids(k) = save_one_model_grid(cfg, quickResult, trackDefs(k));
end
end

function trackDefs = model_grid_tracks(cfg)
source = string(cfg.quickInversion.modelGrid.source);
switch source
    case "insar_tracks"
        trackDefs = cfg.insar.tracks(:);
    otherwise
        error("SlipSolve:UnsupportedOption", "Unsupported cfg.quickInversion.modelGrid.source: %s", source);
end
end

function out = save_one_model_grid(cfg, quickResult, trackDef)
p = slipsolve.project.paths(cfg);
[lon, lat, dataGrid] = grdread2(char(trackDef.dataFile));
[lookE, lookN, lookU] = read_look_grids(trackDef);
[lonMesh, latMesh] = meshgrid(lon, lat);
[originX, originY] = ll2xy(cfg.insar.origin.lon, cfg.insar.origin.lat, cfg.insar.origin.refLon);
[xMeters, yMeters] = ll2xy(lonMesh(:), latMesh(:), cfg.insar.origin.refLon);
xMeters = xMeters - originX;
yMeters = yMeters - originY;
dataInsar = double([xMeters, yMeters, dataGrid(:), lookE(:), lookN(:), lookU(:)]);

valid = isfinite(dataInsar(:, 3)) & all(isfinite(dataInsar(:, 4:6)), 2);
modelValues = nan(size(dataInsar, 1), 1);
chunkSize = cfg.quickInversion.modelGrid.chunkSize;
validIndex = find(valid);
for first = 1:chunkSize:numel(validIndex)
    last = min(first + chunkSize - 1, numel(validIndex));
    rows = validIndex(first:last);
    if isfield(cfg.runtime, "verbose") && cfg.runtime.verbose
        fprintf("... model.grd %s pixels %d-%d of %d\n", ...
            string(get_field(trackDef, "name", "track")), first, last, numel(validIndex));
    end
    modelValues(rows) = quick_model_prediction(trackDef, quickResult.slipModel, quickResult.points, dataInsar(rows, :));
end

modelGrid = reshape(modelValues, size(dataGrid));
trackName = string(get_field(trackDef, "name", "track"));
outputDir = fullfile(p.quickModelGrids, safe_file_name(trackName));
if exist(outputDir, "dir") ~= 7
    mkdir(outputDir);
end
modelFile = fullfile(outputDir, char(cfg.quickInversion.modelGrid.fileName));
grdwrite2(lon, lat, modelGrid, modelFile);

out = struct();
out.name = trackName;
out.sourceGrid = string(trackDef.dataFile);
out.modelFile = string(modelFile);
out.outputDirectory = string(outputDir);
out.pointCount = nnz(valid);
end

function model = quick_model_prediction(trackDef, slipModel, points, dataInsar)
trackType = lower(string(trackDef.dataType));
[east, north, up] = slip2xyz_tri(dataInsar(:, 1), dataInsar(:, 2), slipModel, points);
if trackType == "azo"
    [sinF, cosF] = read_azo_heading(trackDef.dataFile);
    model = north .* cosF + east .* sinF;
else
    model = east .* dataInsar(:, 4) + north .* dataInsar(:, 5) + up .* dataInsar(:, 6);
end
end

function [lookE, lookN, lookU] = read_look_grids(trackDef)
[~, ~, lookE] = grdread2(char(trackDef.lookEFile));
[~, ~, lookN] = grdread2(char(trackDef.lookNFile));
[~, ~, lookU] = grdread2(char(trackDef.lookUFile));
end

function name = safe_file_name(value)
name = regexprep(char(string(value)), "[^A-Za-z0-9_.-]", "_");
name = regexprep(name, "_+$", "");
if isempty(name)
    name = "track";
end
end

function value = get_field(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
