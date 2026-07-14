function [sampling, quickResult] = model_based_sampling(cfg, quickResult)
%MODEL_BASED_SAMPLING Resample data using the quick model-derived quadtree.

p = slipsolve.project.paths(cfg);
slipsolve.project.ensure_directories(cfg);
add_legacy_paths(cfg);

if nargin < 2 || isempty(quickResult)
    loaded = load(p.quickResult, "quickResult");
    quickResult = loaded.quickResult;
end

if ~isfield(quickResult, "modelGrids") || isempty(quickResult.modelGrids)
    error("SlipSolve:MissingInput", ...
        "quickResult.modelGrids is missing. Run the quick_inversion stage before model_based_sampling.");
end

source = string(cfg.quickInversion.modelSampling.source);
switch source
    case "insar_tracks"
        trackDefs = cfg.insar.tracks(:);
    otherwise
        error("SlipSolve:UnsupportedOption", "Unsupported cfg.quickInversion.modelSampling.source: %s", source);
end

meshTrace = [];
if isfield(cfg.paths, "meshTrace") && exist(cfg.paths.meshTrace, "file") == 2
    meshTrace = load(cfg.paths.meshTrace);
end

sampling = repmat(empty_sampling_record(), numel(trackDefs), 1);
for k = 1:numel(trackDefs)
    modelGrid = find_model_grid(quickResult.modelGrids, trackDefs(k));
    sampling(k) = sample_one_track(cfg, trackDefs(k), modelGrid, meshTrace);
end

write_project_sample_list(sampling, p.quickModelSampleList);

for k = 1:numel(sampling)
    sampling(k).sampleListFile = string(p.quickModelSampleList);
    sampling(k).iterationStep = cfg.quickInversion.modelSampling.iterationStep;
    sampling(k).source = source;
end

quickResult.modelSampling = sampling;
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
pathsToAdd(end+1) = string(get_field(cfg.paths, "legacyGeodeticFunctionsRoot", ""));
pathsToAdd(end+1) = string(cfg.paths.geodeticRoot);

for k = 1:numel(pathsToAdd)
    if strlength(pathsToAdd(k)) > 0 && exist(pathsToAdd(k), "dir") == 7
        addpath(char(pathsToAdd(k)));
    end
end
end

function record = empty_sampling_record()
record = struct("name", "", "sourceGrid", "", "modelFile", "", "sampleFile", "", ...
    "sampleCount", 0, "rawCount", 0, "modelTreeCount", 0, ...
    "sampleListFile", "", "iterationStep", 0, "source", "", ...
    "parameters", struct(), "figureFig", "", "figurePng", "");
end

function modelGrid = find_model_grid(modelGrids, trackDef)
trackName = string(get_field(trackDef, "name", "track"));
for k = 1:numel(modelGrids)
    sameName = string(modelGrids(k).name) == trackName;
    sameSource = string(modelGrids(k).sourceGrid) == string(trackDef.dataFile);
    if sameName || sameSource
        modelGrid = modelGrids(k);
        return
    end
end
error("SlipSolve:MissingInput", "No project model.grd was produced for track %s", trackName);
end

function record = sample_one_track(cfg, trackDef, modelGridInfo, meshTrace)
p = slipsolve.project.paths(cfg);
trackName = string(get_field(trackDef, "name", "track"));
sampling = model_based_sampling_options(cfg, trackDef);
[lon, lat, dataGrid] = grdread2(char(trackDef.dataFile));
[~, ~, modelGrid] = grdread2(char(modelGridInfo.modelFile));
[lookE, lookN, lookU] = read_sampling_look_grids(cfg, trackDef);

bounds = sampling.boundsLonLat;
if isempty(bounds)
    bounds = [min(lon), max(lon), min(lat), max(lat)];
end

lonIndex = find(lon >= bounds(1) & lon <= bounds(2));
latIndex = find(lat >= bounds(3) & lat <= bounds(4));
if isempty(lonIndex) || isempty(latIndex)
    error("SlipSolve:EmptySamplingRegion", "Track %s has no pixels inside boundsLonLat.", trackName);
end

lonCrop = lon(lonIndex);
latCrop = lat(latIndex);
dataCrop = dataGrid(latIndex, lonIndex);
modelCrop = modelGrid(latIndex, lonIndex);
lookECrop = lookE(latIndex, lonIndex);
lookNCrop = lookN(latIndex, lonIndex);
lookUCrop = lookU(latIndex, lonIndex);

padX = sampling.padX;
padY = sampling.padY;
[lonPad, latPad, dataCrop] = expandGridNaN(lonCrop, latCrop, dataCrop, padX, padY);
[~, ~, modelCrop] = expandGridNaN(lonCrop, latCrop, modelCrop, padX, padY);
[~, ~, lookECrop] = expandGridNaN(lonCrop, latCrop, lookECrop, padX, padY);
[~, ~, lookNCrop] = expandGridNaN(lonCrop, latCrop, lookNCrop, padX, padY);
[~, ~, lookUCrop] = expandGridNaN(lonCrop, latCrop, lookUCrop, padX, padY);

[lonMesh, latMesh] = meshgrid(lonPad, latPad);
[originX, originY] = ll2xy(cfg.insar.origin.lon, cfg.insar.origin.lat, cfg.insar.origin.refLon);
[xMeters, yMeters] = ll2xy(lonMesh(:), latMesh(:), cfg.insar.origin.refLon);
xKm = (xMeters - originX) ./ 1e3;
yKm = (yMeters - originY) ./ 1e3;

dataValues = dataCrop(:);
modelValues = modelCrop(:);
lookEValues = lookECrop(:);
lookNValues = lookNCrop(:);
lookUValues = lookUCrop(:);

dataIndex = 1:numel(dataValues);
dataExtent = [min(xKm), max(xKm), min(yKm), max(yKm)];
nanFractionMax = sampling.nanFractionMax;
widthMin = sampling.minCellSizeKm;
widthMax = sampling.maxCellSizeKm;
rmsMin = sampling.rmsThreshold;

[xModelSamp, yModelSamp, modelSamp, ~, dataTree] = quadtree_unstructured2( ...
    xKm, yKm, modelValues, dataIndex, dataExtent, meshTrace, ...
    "rms_min", rmsMin, "nan_frac_max", nanFractionMax, ...
    "width_min", widthMin, "width_max", widthMax, ...
    "fault_tol", sampling.faultToleranceKm, ...
    "minPixSide", sampling.minPixelsPerSide, ...
    "stat", sampling.statistic);

if ~isempty(meshTrace)
    [xDataSamp, yDataSamp, dataSamp, rmsOut, ~, dataTreeOut] = apply_unstructured_quadtree2( ...
        xKm, yKm, dataValues, dataTree, nanFractionMax, meshTrace, ...
        "sign_clean", true, ...
        "fault_tol", sampling.faultToleranceKm, ...
        "minPixSide", sampling.minPixelsPerSide, ...
        "stat", sampling.statistic);
else
    [xDataSamp, yDataSamp, dataSamp, rmsOut] = apply_unstructured_quadtree( ...
        xKm, yKm, dataValues, dataTree, nanFractionMax);
    dataTreeOut = dataTree;
end

[~, ~, lookESamp] = apply_unstructured_quadtree(xKm, yKm, lookEValues, dataTreeOut, nanFractionMax);
[~, ~, lookNSamp] = apply_unstructured_quadtree(xKm, yKm, lookNValues, dataTreeOut, nanFractionMax);
[~, ~, lookUSamp] = apply_unstructured_quadtree(xKm, yKm, lookUValues, dataTreeOut, nanFractionMax);
[~, ~, modelSampClean] = apply_unstructured_quadtree(xKm, yKm, modelValues, dataTreeOut, nanFractionMax);

sampledInsarData = double([xDataSamp .* 1e3, yDataSamp .* 1e3, dataSamp, lookESamp, lookNSamp, lookUSamp]);
validRows = all(isfinite(sampledInsarData), 2);
sampledInsarData = sampledInsarData(validRows, :);
if numel(rmsOut) == numel(validRows)
    rmsOut = rmsOut(validRows);
end
if numel(modelSampClean) == numel(validRows)
    modelSampClean = modelSampClean(validRows);
end

rawData = struct("x", xKm(:) .* 1e3, "y", yKm(:) .* 1e3, "los", dataValues(:));
sampledData = struct("x", sampledInsarData(:, 1), "y", sampledInsarData(:, 2), ...
    "los", sampledInsarData(:, 3), "model", modelSampClean);

outputDir = fullfile(p.modelBasedSampling, safe_file_name(trackName));
if exist(outputDir, "dir") ~= 7
    mkdir(outputDir);
end
sampleFile = fullfile(outputDir, sample_file_name(trackDef, cfg.quickInversion.modelSampling.iterationStep));
covd = [];
rms_out = rmsOut;
model_samp = modelSampClean;
sampled_insar_data = sampledInsarData;
source_grid = string(trackDef.dataFile);
model_grid = string(modelGridInfo.modelFile);
model_tree = struct("x", xModelSamp, "y", yModelSamp, "model", modelSamp, "cells", {dataTree});
save(sampleFile, "sampled_insar_data", "rms_out", "covd", "model_samp", ...
    "source_grid", "model_grid", "model_tree", "-v7.3");

figures = slipsolve.plotting.plot_model_sampling_comparison(cfg, trackName, rawData, sampledData);

record = empty_sampling_record();
record.name = trackName;
record.sourceGrid = string(trackDef.dataFile);
record.modelFile = string(modelGridInfo.modelFile);
record.sampleFile = string(sampleFile);
record.sampleCount = size(sampledInsarData, 1);
record.rawCount = nnz(isfinite(dataValues));
record.modelTreeCount = numel(dataTree);
record.parameters = sampling;
record.figureFig = string(figures.fig);
record.figurePng = string(figures.png);
end

function sampling = model_based_sampling_options(cfg, trackDef)
if isfield(trackDef, "modelBasedSampling") && ~isempty(trackDef.modelBasedSampling)
    source = trackDef.modelBasedSampling;
else
    source = struct();
end

defaultNanFraction = cfg.insar.quadtree.nanFractionMax;
if isfield(cfg.quickInversion, "modelSampling")
    defaultNanFraction = get_field(cfg.quickInversion.modelSampling, "nanFractionMax", defaultNanFraction);
end

sampling = struct();
sampling.boundsLonLat = get_field(source, "boundsLonLat", get_field(trackDef, "boundsLonLat", []));
sampling.minCellSizeKm = get_field(source, "minCellSizeKm", ...
    get_field(trackDef, "minCellSizeKm", cfg.insar.quadtree.minBlockSize ./ 1e3));
sampling.maxCellSizeKm = get_field(source, "maxCellSizeKm", ...
    get_field(trackDef, "maxCellSizeKm", cfg.insar.quadtree.maxBlockSize ./ 1e3));
sampling.rmsThreshold = get_field(source, "rmsThreshold", ...
    get_field(trackDef, "rmsThreshold", cfg.insar.quadtree.varianceThreshold));
sampling.nanFractionMax = get_field(source, "nanFractionMax", defaultNanFraction);
sampling.padX = get_field(source, "padX", get_field(trackDef, "padX", 0));
sampling.padY = get_field(source, "padY", get_field(trackDef, "padY", 0));
sampling.faultToleranceKm = get_field(source, "faultToleranceKm", 0.01);
sampling.minPixelsPerSide = get_field(source, "minPixelsPerSide", 2);
sampling.statistic = get_field(source, "statistic", "mean");
end

function [lookE, lookN, lookU] = read_sampling_look_grids(cfg, trackDef)
useHigh = isfield(cfg.quickInversion.modelSampling, "highResolutionLooks") && ...
    cfg.quickInversion.modelSampling.highResolutionLooks;
if ~useHigh
    [lookE, lookN, lookU] = read_look_grids(trackDef);
    return
end

trackDir = fileparts(char(trackDef.dataFile));
trackType = lower(string(trackDef.dataType));
switch trackType
    case "los"
        files = ["look_e_high.grd", "look_n_high.grd", "look_u_high.grd"];
    case "rng"
        files = ["look_e_rng_high.grd", "look_n_rng_high.grd", "look_u_rng_high.grd"];
    otherwise
        files = ["look_e_azo_high.grd", "look_n_azo_high.grd", "look_u_azo_high.grd"];
end

if all(isfile(fullfile(trackDir, files)))
    [~, ~, lookE] = grdread2(char(fullfile(trackDir, files(1))));
    [~, ~, lookN] = grdread2(char(fullfile(trackDir, files(2))));
    [~, ~, lookU] = grdread2(char(fullfile(trackDir, files(3))));
else
    [lookE, lookN, lookU] = read_look_grids(trackDef);
end
end

function [lookE, lookN, lookU] = read_look_grids(trackDef)
[~, ~, lookE] = grdread2(char(trackDef.lookEFile));
[~, ~, lookN] = grdread2(char(trackDef.lookNFile));
[~, ~, lookU] = grdread2(char(trackDef.lookUFile));
end

function write_project_sample_list(sampling, sampleList)
fid = fopen(sampleList, "w");
if fid < 0
    error("SlipSolve:FileWriteError", "Cannot write quick model sample list: %s", sampleList);
end
cleaner = onCleanup(@() fclose(fid));
for k = 1:numel(sampling)
    fprintf(fid, "%s %s %s\n", sampling(k).name, sampling(k).modelFile, sampling(k).sampleFile);
end
end

function name = sample_file_name(trackDef, iterationStep)
trackType = lower(string(trackDef.dataType));
switch trackType
    case "los"
        prefix = "los_samp";
    case "rng"
        prefix = "rng_samp";
    otherwise
        prefix = "azo_samp";
end
name = sprintf("%s%d.mat", prefix, iterationStep);
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
