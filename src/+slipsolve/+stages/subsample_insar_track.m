function trackOut = subsample_insar_track(cfg, trackDef)
%SUBSAMPLE_INSAR_TRACK Subsample one InSAR/range/azimuth grid.

add_sampling_paths(cfg);
sampling = initial_sampling_options(cfg, trackDef);

[lon, lat, dataGrid] = grdread2(char(trackDef.dataFile));
[lookE, lookN, lookU] = read_look_grids(trackDef);

padX = sampling.padX;
padY = sampling.padY;
if padX ~= 0 || padY ~= 0
    lonOriginal = lon;
    latOriginal = lat;
    [lon, lat, dataGrid] = expandGridNaN(lonOriginal, latOriginal, dataGrid, padX, padY);
    [~, ~, lookE] = expandGridNaN(lonOriginal, latOriginal, lookE, padX, padY);
    [~, ~, lookN] = expandGridNaN(lonOriginal, latOriginal, lookN, padX, padY);
    [~, ~, lookU] = expandGridNaN(lonOriginal, latOriginal, lookU, padX, padY);
end

boundsLonLat = sampling.boundsLonLat;
if isempty(boundsLonLat)
    ix = true(size(lon));
    iy = true(size(lat));
else
    ix = lon >= boundsLonLat(1) & lon <= boundsLonLat(2);
    iy = lat >= boundsLonLat(3) & lat <= boundsLonLat(4);
end

lonCrop = lon(ix);
latCrop = lat(iy);
dataCrop = dataGrid(iy, ix);
lookECrop = lookE(iy, ix);
lookNCrop = lookN(iy, ix);
lookUCrop = lookU(iy, ix);

[lonMesh, latMesh] = meshgrid(lonCrop, latCrop);
[originX, originY] = ll2xy(cfg.insar.origin.lon, cfg.insar.origin.lat, cfg.insar.origin.refLon);
[xMeters, yMeters] = ll2xy(lonMesh(:), latMesh(:), cfg.insar.origin.refLon);
xMeters = xMeters - originX;
yMeters = yMeters - originY;
xKm = xMeters ./ 1e3;
yKm = yMeters ./ 1e3;

dataValues = dataCrop(:);
lookEValues = lookECrop(:);
lookNValues = lookNCrop(:);
lookUValues = lookUCrop(:);
dataIndex = 1:numel(dataValues);
dataExtent = [min(xKm), max(xKm), min(yKm), max(yKm)];
faultTraceKm = read_fault_trace_km(cfg);

[xSampKm, ySampKm, dataSamp, rmsOut, dataTree] = quadtree_unstructured2( ...
    xKm, yKm, dataValues, dataIndex, dataExtent, faultTraceKm, ...
    "rms_min", sampling.rmsThreshold, ...
    "nan_frac_max", sampling.nanFractionMax, ...
    "width_min", sampling.minCellSizeKm, ...
    "width_max", sampling.maxCellSizeKm, ...
    "fault_tol", sampling.faultToleranceKm, ...
    "minPixSide", sampling.minPixelsPerSide, ...
    "stat", sampling.statistic);

[~, ~, lookESamp] = apply_unstructured_quadtree(xKm, yKm, lookEValues, dataTree, sampling.nanFractionMax);
[~, ~, lookNSamp] = apply_unstructured_quadtree(xKm, yKm, lookNValues, dataTree, sampling.nanFractionMax);
[~, ~, lookUSamp] = apply_unstructured_quadtree(xKm, yKm, lookUValues, dataTree, sampling.nanFractionMax);

sampled = double([xSampKm .* 1e3, ySampKm .* 1e3, dataSamp, lookESamp, lookNSamp, lookUSamp]);
validRows = all(isfinite(sampled), 2);
sampled = sampled(validRows, :);
rmsOut = rmsOut(validRows);

trackOut = struct();
trackOut.name = string(get_field(trackDef, "name", "insar_track"));
trackOut.type = string(get_field(trackDef, "dataType", "insar"));
trackOut.sourceGrid = string(trackDef.dataFile);
trackOut.x = sampled(:, 1);
trackOut.y = sampled(:, 2);
trackOut.los = sampled(:, 3);
trackOut.lookE = sampled(:, 4);
trackOut.lookN = sampled(:, 5);
trackOut.lookU = sampled(:, 6);
trackOut.rms = rmsOut;
trackOut.sampleCount = size(sampled, 1);
trackOut.rawFiniteCount = nnz(isfinite(dataValues));
trackOut.boundsLonLat = [min(lonCrop), max(lonCrop), min(latCrop), max(latCrop)];
trackOut.quadtree.minCellSizeKm = sampling.minCellSizeKm;
trackOut.quadtree.maxCellSizeKm = sampling.maxCellSizeKm;
trackOut.quadtree.rmsThreshold = sampling.rmsThreshold;
trackOut.quadtree.nanFractionMax = sampling.nanFractionMax;
trackOut.quadtree.faultToleranceKm = sampling.faultToleranceKm;
trackOut.quadtree.minPixelsPerSide = sampling.minPixelsPerSide;
trackOut.quadtree.statistic = sampling.statistic;
trackOut.quadtree.function = "quadtree_unstructured2";

rawData = struct();
rawData.x = xMeters;
rawData.y = yMeters;
rawData.los = dataValues;
trackOut.figures = slipsolve.plotting.plot_insar_subsampling(cfg, trackOut, rawData);

end

function add_sampling_paths(cfg)
pathsToAdd = [
    string(cfg.paths.geodeticRoot)
    string(cfg.paths.legacySamplingRoot)
    string(cfg.paths.legacyOtherFuncRoot)
    string(get_field(cfg.paths, "legacyGeodeticFunctionsRoot", ""))
];

for k = 1:numel(pathsToAdd)
    if strlength(pathsToAdd(k)) > 0 && exist(pathsToAdd(k), "dir") == 7
        addpath(char(pathsToAdd(k)));
    end
end
end

function [lookE, lookN, lookU] = read_look_grids(trackDef)
[~, ~, lookE] = grdread2(char(trackDef.lookEFile));
[~, ~, lookN] = grdread2(char(trackDef.lookNFile));
[~, ~, lookU] = grdread2(char(trackDef.lookUFile));
end

function sampling = initial_sampling_options(cfg, trackDef)
if isfield(trackDef, "initialSampling") && ~isempty(trackDef.initialSampling)
    source = trackDef.initialSampling;
else
    source = struct();
end

sampling = struct();
sampling.boundsLonLat = get_field(source, "boundsLonLat", get_field(trackDef, "boundsLonLat", []));
sampling.minCellSizeKm = get_field(source, "minCellSizeKm", ...
    get_field(trackDef, "minCellSizeKm", cfg.insar.quadtree.minBlockSize ./ 1e3));
sampling.maxCellSizeKm = get_field(source, "maxCellSizeKm", ...
    get_field(trackDef, "maxCellSizeKm", cfg.insar.quadtree.maxBlockSize ./ 1e3));
sampling.rmsThreshold = get_field(source, "rmsThreshold", ...
    get_field(trackDef, "rmsThreshold", cfg.insar.quadtree.varianceThreshold));
sampling.nanFractionMax = get_field(source, "nanFractionMax", cfg.insar.quadtree.nanFractionMax);
sampling.padX = get_field(source, "padX", get_field(trackDef, "padX", 0));
sampling.padY = get_field(source, "padY", get_field(trackDef, "padY", 0));
sampling.faultToleranceKm = get_field(source, "faultToleranceKm", cfg.insar.quadtree.faultTolerance ./ 1e3);
sampling.minPixelsPerSide = get_field(source, "minPixelsPerSide", cfg.insar.quadtree.minPixelsPerSide);
sampling.statistic = get_field(source, "statistic", cfg.insar.quadtree.statistic);
end

function faultTraceKm = read_fault_trace_km(cfg)
faultTraceKm = [];
if ~isfield(cfg.paths, "meshTrace") || strlength(string(cfg.paths.meshTrace)) == 0
    return
end

if exist(cfg.paths.meshTrace, "file") ~= 2
    return
end

faultTrace = load(cfg.paths.meshTrace);
if size(faultTrace, 2) >= 2
    faultTraceKm = faultTrace(:, 1:2) ./ 1e3;
end
end

function value = get_field(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
