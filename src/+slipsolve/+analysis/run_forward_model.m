function forwardModelResult = run_forward_model(cfg)
%RUN_FORWARD_MODEL Project the final composite model onto independent data.

if ~isfield(cfg, "forwardModel") || ...
        ~get_field(cfg.forwardModel, "enabled", false)
    error("SlipSolve:DisabledHelper", ...
        "Set cfg.forwardModel.enabled = true before running this helper.");
end
slipsolve.project.ensure_directories(cfg);
p = slipsolve.project.paths(cfg);
settings = cfg.forwardModel;
add_legacy_paths(cfg);
modelFile = string(get_field(settings, "modelResultFile", p.finalResult));
if exist(modelFile, "file") ~= 2
    error("SlipSolve:MissingStageProduct", ...
        "Run final_inversion before forward modeling: %s", modelFile);
end
loaded = load(modelFile, "finalResult");
if ~isfield(loaded, "finalResult") || ...
        ~isfield(loaded.finalResult, "slipVector")
    error("SlipSolve:InvalidStageProduct", ...
        "%s does not contain finalResult.slipVector.", modelFile);
end
slipVector = double(loaded.finalResult.slipVector(:));
nColumns = numel(slipVector);
[xVector, yVector] = full_grid_vectors(cfg);
definitions = settings.tracks(:);
if isempty(definitions)
    error("SlipSolve:InvalidForwardModelConfig", ...
        "cfg.forwardModel.tracks must contain at least one independent dataset.");
end
fullResponse = build_full_component_response(settings, definitions, ...
    slipVector, xVector, yVector, p, modelFile);

trackCells = cell(numel(definitions), 1);
for k = 1:numel(definitions)
    definition = definitions(k);
    [sampleFile, sampled] = load_sample(settings, definition);
    [projection, type] = observation_projection(definition, sampleFile, sampled);
    signature = cache_signature(settings, sampleFile, projection, ...
        xVector, yVector, nColumns);
    cacheFile = fullfile(p.forwardModelGreensDirectory, ...
        safe_name(definition.name)+"_sampled_composite_green.mat");
    [green, reused] = load_or_build_green(settings, sampleFile, ...
        sampled, projection, xVector, yVector, nColumns, cacheFile, signature);

    prediction = green*slipVector;
    data = sampled(:, 3);
    residual = data-prediction;
    track = struct();
    track.name = string(definition.name);
    track.type = type;
    track.sampleFile = sampleFile;
    track.x = sampled(:, 1);
    track.y = sampled(:, 2);
    track.data = data;
    track.prediction = prediction;
    track.residualDataMinusModel = residual;
    track.projectionEast = projection(:, 1);
    track.projectionNorth = projection(:, 2);
    track.projectionUp = projection(:, 3);
    track.sampledGreenCache = cacheFile;
    track.sampledGreenReused = reused;
    track.rmsData = sqrt(mean(data.^2, "omitnan"));
    track.rmsResidual = sqrt(mean(residual.^2, "omitnan"));
    track.varianceReduction = 100*(sum(data.^2)-sum(residual.^2))/sum(data.^2);
    track.outputMatFile = fullfile(p.forwardModelDirectory, ...
        safe_name(definition.name)+"_forward_prediction.mat");
    track.outputCsvFile = fullfile(p.forwardModelDirectory, ...
        safe_name(definition.name)+"_forward_prediction.csv");
    predictionTable = table(track.x, track.y, track.data, track.prediction, ...
        track.residualDataMinusModel, track.projectionEast, ...
        track.projectionNorth, track.projectionUp, ...
        VariableNames=["x_m", "y_m", "data_cm", "model_cm", ...
        "data_minus_model_cm", "projection_e", "projection_n", "projection_u"]);
    writetable(predictionTable, track.outputCsvFile);
    fullGrid = struct();
    if get_field(settings.fullResolution, "enabled", true)
        fullGrid = full_resolution_prediction(cfg, settings, definition, ...
            track, fullResponse, p);
        track.fullResolution = fullGrid.metadata;
    else
        track.fullResolution = struct();
    end
    track.figureFiles = slipsolve.plotting.plot_forward_model(cfg, track, fullGrid);
    forwardTrack = track;
    save(track.outputMatFile, "forwardTrack", "predictionTable", "-v7.3");
    trackCells{k} = track;
    fprintf("Forward model %s: RMS data %.4f cm, residual %.4f cm, VR %.3f%%\n", ...
        track.name, track.rmsData, track.rmsResidual, track.varianceReduction);
end

forwardModelResult = struct();
forwardModelResult.method = ...
    "full_xyz_interpolation_with_full_resolution_independent_data_fit";
forwardModelResult.legacyReference = "forwardmodel.m and plotForwardModeling2.m";
forwardModelResult.modelResultFile = modelFile;
forwardModelResult.componentFiles = settings.componentFiles;
forwardModelResult.residualConvention = "data_minus_model";
forwardModelResult.tracks = vertcat(trackCells{:});
forwardModelResult.createdAt = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");
forwardModelResult.resultFile = p.forwardModelResult;
save(p.forwardModelResult, "forwardModelResult", "-v7.3");
fprintf("Saved independent-data forward model: %s\n", p.forwardModelResult);
end

function response = build_full_component_response(settings, definitions, ...
    slipVector, xVector, yVector, p, modelFile)
response = struct();
if ~get_field(settings.fullResolution, "enabled", true)
    return
end
types = lower(string({definitions.type}));
cacheFile = fullfile(p.forwardModelGreensDirectory, ...
    "full_component_forward_response.mat");
signature = struct();
signature.modelResult = file_signature(modelFile);
signature.componentFiles = [file_signature(settings.componentFiles.east); ...
    file_signature(settings.componentFiles.north); ...
    file_signature(settings.componentFiles.vertical)];
signature.types = types;
signature.gridX = [numel(xVector), min(xVector), max(xVector), sum(xVector)];
signature.gridY = [numel(yVector), min(yVector), max(yVector), sum(yVector)];
if get_field(settings.fullResolution, "reuseComponentCache", true) && ...
        exist(cacheFile, "file") == 2
    cachedMetadata = load(cacheFile, "fullResponseMetadata");
    if isfield(cachedMetadata, "fullResponseMetadata") && ...
            isequaln(cachedMetadata.fullResponseMetadata.signature, signature)
        cached = load(cacheFile, "fullComponentResponse");
        response = cached.fullComponentResponse;
        fprintf("Loading cached full XYZ forward response: %s\n", cacheFile);
        return
    end
end
needEast = any(ismember(types, ["azo", "los", "rng", "east"]));
needNorth = any(ismember(types, ["azo", "los", "rng", "north"]));
needUp = any(ismember(types, ["los", "rng", "up", "vertical"]));
[X, Y] = meshgrid(xVector, yVector);
response.x = X(:);
response.y = Y(:);
response.gridSize = size(X);
if needEast
    response.east = component_prediction( ...
        settings.componentFiles.east, "G_e", slipVector);
end
if needNorth
    response.north = component_prediction( ...
        settings.componentFiles.north, "G_n", slipVector);
end
if needUp
    response.up = component_prediction( ...
        settings.componentFiles.vertical, "G_u", slipVector);
end
fullComponentResponse = response;
fullResponseMetadata = struct("signature", signature, ...
    "createdAt", datetime("now", "Format", "yyyy-MM-dd HH:mm:ss"));
save(cacheFile, "fullComponentResponse", "fullResponseMetadata", "-v7.3");
end

function prediction = component_prediction(filePath, variableName, slipVector)
fprintf("Evaluating full-resolution forward component %s.\n", variableName);
loaded = load(filePath, variableName);
if ~isfield(loaded, variableName)
    error("SlipSolve:MissingInput", "%s lacks %s.", filePath, variableName);
end
prediction = double(loaded.(variableName))*slipVector;
end

function fullGrid = full_resolution_prediction(cfg, settings, definition, ...
    track, response, p)
dataGridFile = resolve_track_file(settings, definition, "dataGridFile");
if exist(dataGridFile, "file") ~= 2
    error("SlipSolve:MissingInput", ...
        "Full-resolution independent data grid is missing: %s", dataGridFile);
end
name = safe_name(definition.name);
modelGridFile = fullfile(p.forwardModelDirectory, name+"_full_model.grd");
residualGridFile = fullfile(p.forwardModelDirectory, name+"_full_residual.grd");
cacheFile = fullfile(p.forwardModelDirectory, name+"_full_grid_cache.mat");
signature = struct();
signature.modelResult = file_signature(settings.modelResultFile);
signature.dataGrid = file_signature(dataGridFile);
signature.componentFiles = [file_signature(settings.componentFiles.east); ...
    file_signature(settings.componentFiles.north); ...
    file_signature(settings.componentFiles.vertical)];
signature.faultTrace = file_signature(settings.fullResolution.faultTraceFile);
signature.projection = [track.projectionEast(1), ...
    track.projectionNorth(1), track.projectionUp(1)];
signature.method = string(settings.fullResolution.interpolationMethod);
signature.extrapolation = string(settings.fullResolution.extrapolation);
if get_field(settings.fullResolution, "reuseGridCache", true) && ...
        exist(cacheFile, "file") == 2 && exist(modelGridFile, "file") == 2 && ...
        exist(residualGridFile, "file") == 2
    cached = load(cacheFile, "fullGridMetadata");
    if isfield(cached, "fullGridMetadata") && ...
            isequaln(cached.fullGridMetadata.signature, signature)
        [longitude, latitude, dataGrid] = grdread2(char(dataGridFile));
        [~, ~, modelGrid] = grdread2(char(modelGridFile));
        [~, ~, residualGrid] = grdread2(char(residualGridFile));
        [queryX, queryY] = geographic_grid_to_local( ...
            longitude, latitude, cfg.insar.origin);
        fault = read_fault_trace(settings.fullResolution.faultTraceFile);
        fullGrid = struct("x", queryX, "y", queryY, "data", dataGrid, ...
            "model", modelGrid, "residual", residualGrid, "fault", fault, ...
            "metadata", cached.fullGridMetadata);
        fprintf("Loading cached full-resolution forward grids: %s\n", cacheFile);
        return
    end
end
[longitude, latitude, dataGrid] = grdread2(char(dataGridFile));
[queryX, queryY] = geographic_grid_to_local( ...
    longitude, latitude, cfg.insar.origin);
projected = project_full_response(cfg, settings, definition, track, response);
fault = read_fault_trace(settings.fullResolution.faultTraceFile);
if size(fault, 1) < 2
    error("SlipSolve:InvalidForwardModelConfig", ...
        "forwardModel.fullResolution.faultTraceFile needs at least two local x/y rows.");
end
modelGrid = interpXYZ_barrier( ...
    [response.x, response.y, projected], queryX, queryY, fault, ...
    "Method", char(settings.fullResolution.interpolationMethod), ...
    "Extrap", char(settings.fullResolution.extrapolation));
modelGrid(~isfinite(dataGrid)) = NaN;
residualGrid = dataGrid-modelGrid;

if get_field(settings.fullResolution, "saveGrids", true)
    grdwrite2(longitude, latitude, modelGrid, char(modelGridFile));
    grdwrite2(longitude, latitude, residualGrid, char(residualGridFile));
else
    modelGridFile = "";
    residualGridFile = "";
end

valid = isfinite(dataGrid) & isfinite(modelGrid);
metadata = struct();
metadata.dataGridFile = dataGridFile;
metadata.modelGridFile = modelGridFile;
metadata.residualGridFile = residualGridFile;
metadata.gridSize = size(dataGrid);
metadata.validPixelCount = nnz(valid);
metadata.rmsData = sqrt(mean(dataGrid(valid).^2));
metadata.rmsResidual = sqrt(mean(residualGrid(valid).^2));
metadata.varianceReduction = 100*(sum(dataGrid(valid).^2)- ...
    sum(residualGrid(valid).^2))/sum(dataGrid(valid).^2);
metadata.cacheFile = cacheFile;
metadata.signature = signature;
fullGridMetadata = metadata;
save(cacheFile, "fullGridMetadata", "-v7.3");
fullGrid = struct("x", queryX, "y", queryY, "data", dataGrid, ...
    "model", modelGrid, "residual", residualGrid, "fault", fault, ...
    "metadata", metadata);
end

function projected = project_full_response(cfg, settings, definition, track, response)
type = lower(string(definition.type));
switch type
    case "azo"
        projected = track.projectionEast(1)*response.east + ...
            track.projectionNorth(1)*response.north;
    case {"los", "rng"}
        lookEFile = resolve_track_file(settings, definition, "lookEFile");
        lookNFile = resolve_track_file(settings, definition, "lookNFile");
        lookUFile = resolve_track_file(settings, definition, "lookUFile");
        [lookLon, lookLat, lookE] = grdread2(char(lookEFile));
        [~, ~, lookN] = grdread2(char(lookNFile));
        [~, ~, lookU] = grdread2(char(lookUFile));
        [lookX, lookY] = geographic_grid_to_local(lookLon, lookLat, cfg.insar.origin);
        east = interpXYZ([lookX(:), lookY(:), lookE(:)], ...
            response.x, response.y, "natural", "none");
        north = interpXYZ([lookX(:), lookY(:), lookN(:)], ...
            response.x, response.y, "natural", "none");
        up = interpXYZ([lookX(:), lookY(:), lookU(:)], ...
            response.x, response.y, "natural", "none");
        projected = response.east.*east + response.north.*north + response.up.*up;
    case "east"
        projected = response.east;
    case "north"
        projected = response.north;
    case {"up", "vertical"}
        projected = response.up;
    otherwise
        error("SlipSolve:InvalidForwardModelConfig", ...
            "Unknown independent-data type: %s", type);
end
end

function filePath = resolve_track_file(settings, definition, fieldName)
if ~isfield(definition, fieldName) || ...
        strlength(string(definition.(fieldName))) == 0
    error("SlipSolve:InvalidForwardModelConfig", ...
        "Forward track %s requires %s for full-resolution plotting.", ...
        definition.name, fieldName);
end
filePath = string(definition.(fieldName));
if ~isfile(filePath)
    filePath = fullfile(string(settings.samplesRoot), filePath);
end
end

function [X, Y] = geographic_grid_to_local(longitude, latitude, origin)
[longitudeGrid, latitudeGrid] = meshgrid(longitude, latitude);
[X, Y] = ll2xy(longitudeGrid(:), latitudeGrid(:), origin.refLon);
[originX, originY] = ll2xy(origin.lon, origin.lat, origin.refLon);
X = reshape(X-originX, size(longitudeGrid));
Y = reshape(Y-originY, size(latitudeGrid));
end

function fault = read_fault_trace(filePath)
if exist(filePath, "file") ~= 2
    error("SlipSolve:MissingInput", "Forward-model fault trace is missing: %s", filePath);
end
fault = readmatrix(filePath);
fault = double(fault(:, 1:2));
fault = fault(all(isfinite(fault), 2), :);
end

function [green, reused] = load_or_build_green(settings, sampleFile, ...
    sampled, projection, xVector, yVector, nColumns, cacheFile, signature)
reuse = get_field(settings.interpolation, "reuseCached", true);
overwrite = get_field(settings.interpolation, "overwriteCached", false);
if reuse && exist(cacheFile, "file") == 2
    metadataLoaded = load(cacheFile, "forwardGreenMetadata");
    if isfield(metadataLoaded, "forwardGreenMetadata") && ...
            isequaln(metadataLoaded.forwardGreenMetadata.signature, signature)
        loaded = load(cacheFile, "G_forward");
        green = loaded.G_forward;
        reused = true;
        fprintf("Loading cached independent sampled Green: %s\n", cacheFile);
        return
    end
end
if exist(cacheFile, "file") == 2 && ~overwrite
    error("SlipSolve:ProtectedOutput", ...
        "Forward Green cache no longer matches %s. Set overwriteCached=true: %s", ...
        sampleFile, cacheFile);
end
green = slipsolve.composite.interpolate_projected_greens( ...
    settings.componentFiles, xVector, yVector, sampled(:, 1), sampled(:, 2), ...
    projection, nColumns, settings.interpolation);
reused = false;
if get_field(settings.interpolation, "saveCached", true)
    forwardGreenMetadata = struct();
    forwardGreenMetadata.signature = signature;
    forwardGreenMetadata.sampleFile = sampleFile;
    forwardGreenMetadata.matrixSize = size(green);
    forwardGreenMetadata.columnOrder = ...
        "all_faults_strike_slip_then_all_faults_dip_slip";
    forwardGreenMetadata.createdAt = datetime("now", ...
        "Format", "yyyy-MM-dd HH:mm:ss");
    G_forward = green;
    save(cacheFile, "G_forward", "forwardGreenMetadata", "-v7.3");
end
end

function [sampleFile, sampled] = load_sample(settings, definition)
if isfield(definition, "sampleFile") && strlength(string(definition.sampleFile)) > 0
    sampleFile = string(definition.sampleFile);
elseif isfield(definition, "relativeFile")
    sampleFile = fullfile(string(settings.samplesRoot), string(definition.relativeFile));
else
    error("SlipSolve:InvalidForwardModelConfig", ...
        "Each forward-model track needs sampleFile or relativeFile.");
end
if exist(sampleFile, "file") ~= 2
    error("SlipSolve:MissingInput", "Independent sample is missing: %s", sampleFile);
end
loaded = load(sampleFile, "sampled_insar_data");
if ~isfield(loaded, "sampled_insar_data")
    error("SlipSolve:InvalidInput", "%s lacks sampled_insar_data.", sampleFile);
end
sampled = double(loaded.sampled_insar_data);
if size(sampled, 2) < 3 || any(~isfinite(sampled(:, 1:3)), "all")
    error("SlipSolve:InvalidInput", ...
        "%s must contain finite [x_m y_m displacement_cm] columns.", sampleFile);
end
end

function [projection, type] = observation_projection(definition, sampleFile, sampled)
type = lower(string(definition.type));
n = size(sampled, 1);
switch type
    case "azo"
        [sinF, cosF] = azo_projection(definition, sampleFile);
        projection = repmat([sinF cosF 0], n, 1);
    case {"los", "rng"}
        if size(sampled, 2) < 6
            error("SlipSolve:InvalidInput", ...
                "%s requires look-vector columns 4:6 for type %s.", sampleFile, type);
        end
        projection = sampled(:, 4:6);
    case "east"
        projection = repmat([1 0 0], n, 1);
    case "north"
        projection = repmat([0 1 0], n, 1);
    case {"up", "vertical"}
        projection = repmat([0 0 1], n, 1);
    otherwise
        error("SlipSolve:InvalidForwardModelConfig", ...
            "Unknown independent-data type: %s", type);
end
end

function [sinF, cosF] = azo_projection(definition, sampleFile)
if isfield(definition, "sinF") && isfield(definition, "cosF") && ...
        ~isempty(definition.sinF) && ~isempty(definition.cosF)
    sinF = double(definition.sinF);
    cosF = double(definition.cosF);
else
    folder = string(fileparts(sampleFile));
    sinFile = fullfile(folder, "sinF.dat");
    cosFile = fullfile(folder, "cosF.dat");
    if isfield(definition, "sinFFile") && strlength(string(definition.sinFFile)) > 0
        sinFile = string(definition.sinFFile);
    end
    if isfield(definition, "cosFFile") && strlength(string(definition.cosFFile)) > 0
        cosFile = string(definition.cosFFile);
    end
    if exist(sinFile, "file") ~= 2 || exist(cosFile, "file") ~= 2
        error("SlipSolve:MissingInput", ...
            "AZO forward data requires sinF.dat/cosF.dat or explicit sinF/cosF.");
    end
    sinF = readmatrix(sinFile);
    cosF = readmatrix(cosFile);
    sinF = double(sinF(1));
    cosF = double(cosF(1));
end
end

function [xVector, yVector] = full_grid_vectors(cfg)
legacyDirectory = fullfile(string(cfg.paths.legacyRoot), "layered_green");
if exist(legacyDirectory, "dir") == 7
    addpath(legacyDirectory);
end
if exist("adaptive_meshgrid", "file") ~= 2
    error("SlipSolve:MissingLegacyHelper", ...
        "adaptive_meshgrid.m is required to recover the full-Green grid vectors.");
end
trace = readmatrix(cfg.compositeGreens.meshTraceFile);
trace = trace(:, 1:2);
trace = trace(all(isfinite(trace), 2), :);
[~, ~, xVector, yVector] = adaptive_meshgrid(trace, ...
    cfg.compositeGreens.xRange, cfg.compositeGreens.yRange, ...
    cfg.compositeGreens.D1, cfg.compositeGreens.D2, cfg.compositeGreens.M);
xVector = double(xVector(:).');
yVector = double(yVector(:).');
end

function signature = cache_signature(settings, sampleFile, projection, ...
    xVector, yVector, nColumns)
signature = struct();
signature.sampleFile = file_signature(sampleFile);
signature.componentFiles = [file_signature(settings.componentFiles.east); ...
    file_signature(settings.componentFiles.north); ...
    file_signature(settings.componentFiles.vertical)];
signature.projection = [size(projection, 1), sum(projection, 1), ...
    sum(projection.^2, 1)];
signature.gridX = [numel(xVector), min(xVector), max(xVector), ...
    sum(xVector), sum(xVector.^2)];
signature.gridY = [numel(yVector), min(yVector), max(yVector), ...
    sum(yVector), sum(yVector.^2)];
signature.nColumns = nColumns;
signature.method = string(settings.interpolation.method);
end

function value = file_signature(filePath)
details = dir(filePath);
if isempty(details)
    error("SlipSolve:MissingInput", "Input is missing: %s", filePath);
end
value = struct("path", string(filePath), "bytes", details.bytes, ...
    "datenum", details.datenum);
end

function name = safe_name(value)
name = regexprep(string(value), "[^A-Za-z0-9_.-]", "_");
end

function add_legacy_paths(cfg)
paths = [string(cfg.paths.legacyRoot); ...
    string(cfg.paths.legacyOtherFuncRoot); ...
    fullfile(string(cfg.paths.legacyRoot), "layered_green")];
for path = paths(:).'
    if exist(path, "dir") == 7
        addpath(char(path));
    end
end
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
