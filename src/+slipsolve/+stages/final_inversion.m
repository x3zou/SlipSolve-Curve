function finalResult = final_inversion(cfg, mesh)
%FINAL_INVERSION Legacy altered-node inversion from full XYZ layered Greens.

p = slipsolve.project.paths(cfg);
slipsolve.project.ensure_directories(cfg);
invCfg = cfg.finalInversion;
if ~get_field(invCfg, "enabled", true)
    error("SlipSolve:DisabledStage", ...
        "Set cfg.finalInversion.enabled = true before running final_inversion.");
end
add_legacy_paths(cfg);

tracks = load_tracks(invCfg);
[xVector, yVector] = full_grid_vectors(cfg);
model = prepare_hybrid_model(mesh);
nParameter = size(model.hybrid, 1);

componentFiles = invCfg.componentFiles;
validate_component_file(componentFiles.east, "G_e", numel(xVector) * numel(yVector), 2*nParameter);
validate_component_file(componentFiles.north, "G_n", numel(xVector) * numel(yVector), 2*nParameter);
validate_component_file(componentFiles.vertical, "G_u", numel(xVector) * numel(yVector), 2*nParameter);

[greenRaw, cacheInfo] = sampled_composite_green(cfg, tracks, xVector, yVector, 2*nParameter);
[greenWeighted, dataRaw, dataWeighted, tracks] = weight_tracks(invCfg, greenRaw, tracks);

[smoothMatrix, smoothRows] = smoothingMatrix_laplace2( ...
    model.hybrid(:, 9:15), model.hybrid(:, 6:8), ...
    "ss_ratio", invCfg.smoothingStrikeRatio, ...
    "ds_ratio", invCfg.smoothingDipRatio);
lambda = invCfg.smoothingWeight * invCfg.smoothingInternalScale;
[boundaryMatrix, boundaryData, boundaryInfo] = boundary_constraints(invCfg, model.hybrid);

systemMatrix = [greenWeighted; smoothMatrix .* lambda ./ smoothRows; boundaryMatrix];
systemData = [dataWeighted; zeros(smoothRows, 1); boundaryData];
[lowerBounds, upperBounds] = inversion_bounds(invCfg, model.hybrid);

options = optimset("LargeScale", "on", "DiffMaxChange", 1e-1, ...
    "DiffMinChange", 1e-12, "TolCon", 1e-12, "TolFun", 1e-12, ...
    "TolPCG", 1e-12, "TolX", 1e-12, "MaxIter", 1e9, ...
    "MaxPCGIter", 1e9);
fprintf("Solving final layered system: %d equations x %d parameters\n", ...
    size(systemMatrix, 1), size(systemMatrix, 2));
[slipVector, resnorm, solverResidual, exitflag] = lsqlin( ...
    systemMatrix, double(systemData), [], [], [], [], ...
    lowerBounds, upperBounds, [], options);

model.hybrid(:, 4) = slipVector(1:nParameter);
model.hybrid(:, 5) = slipVector(nParameter+1:end);
prediction = greenRaw * slipVector;
residual = prediction - dataRaw;
tracks = attach_predictions(tracks, prediction);
[completeModel, triangleModel, simpleTriangular, exportInfo] = ...
    export_simple_triangular(cfg, invCfg, model);

roughnessVector = smoothMatrix * slipVector;
finalResult = struct();
finalResult.method = "legacy_MMInversion_altered_interpolated_full_xyz";
finalResult.legacyReference = "MMInversion_altered.m";
finalResult.interpolationReference = "full G_e/G_n/G_u projected at sampled observations";
finalResult.hybridSlipModel = model.hybrid;
finalResult.completeNodalSlipModel = completeModel;
finalResult.triangleModel = triangleModel;
finalResult.simpleTriangular = simpleTriangular;
finalResult.simpleTriangularColumns = exportInfo.columns;
finalResult.simpleTriangularMatFile = exportInfo.matFile;
finalResult.simpleTriangularTextFile = exportInfo.textFile;
finalResult.M0 = exportInfo.M0;
finalResult.Mw = exportInfo.Mw;
finalResult.maximumSlipCm = max(hypot(completeModel(:, 4), completeModel(:, 5)));
finalResult.slipVector = slipVector;
finalResult.data = dataRaw;
finalResult.prediction = prediction;
finalResult.residual = residual;
finalResult.tracks = tracks;
finalResult.sampledCompositeGreenCache = cacheInfo;
finalResult.resnorm = resnorm;
finalResult.solverResidual = solverResidual;
finalResult.exitflag = exitflag;
finalResult.rms0 = sum(dataRaw.^2);
finalResult.rms = sum(residual.^2);
finalResult.varianceReduction = 100 * (finalResult.rms0-finalResult.rms) / finalResult.rms0;
finalResult.modelRoughness = sqrt(sum(roughnessVector.^2) / numel(roughnessVector));
finalResult.boundaryInfo = boundaryInfo;
finalResult.config = invCfg;
fullGridResponse = build_full_grid_response(cfg, slipVector, xVector, yVector);
finalResult.figures = slipsolve.plotting.plot_final_layered_inversion( ...
    cfg, finalResult, fullGridResponse);
finalResult.createdAt = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");

save(p.finalResult, "finalResult", "-v7.3");
fprintf("Final layered rms misfit (data, residual) = %e %e (%.3f%%)\n", ...
    finalResult.rms0, finalResult.rms, finalResult.varianceReduction);
fprintf("Final layered roughness = %e; exitflag = %d\n", ...
    finalResult.modelRoughness, exitflag);
fprintf("Maximum slip = %.3f cm; Mw = %.4f\n", ...
    finalResult.maximumSlipCm, finalResult.Mw);
end

function response = build_full_grid_response(cfg, slipVector, xVector, yVector)
response = struct();
if ~isfield(cfg.visualization, "finalInversion") || ...
        ~isfield(cfg.visualization.finalInversion, "fullResolutionFitPlots") || ...
        ~get_field(cfg.visualization.finalInversion.fullResolutionFitPlots, "enabled", false)
    return
end
files = cfg.finalInversion.componentFiles;
fprintf("Evaluating solved model on the full composite XYZ grid.\n");
response.east = component_prediction(files.east, "G_e", slipVector);
response.north = component_prediction(files.north, "G_n", slipVector);
response.up = component_prediction(files.vertical, "G_u", slipVector);
[X, Y] = meshgrid(xVector, yVector);
response.x = X(:);
response.y = Y(:);
response.gridSize = size(X);
end

function prediction = component_prediction(filePath, variableName, slipVector)
fprintf("  Loading and evaluating %s\n", variableName);
loaded = load(filePath, variableName);
prediction = loaded.(variableName)*slipVector;
end

function add_legacy_paths(cfg)
root = string(cfg.paths.legacyRoot);
projectRoot = fileparts(root);
paths = [root; ...
    fullfile(root, "layered_green"); ...
    fullfile(root, "layered_green", "otherfunc"); ...
    fullfile(projectRoot, "resolution_test"); ...
    fullfile(projectRoot, "geodetic_inversion-master", "inversion"); ...
    fullfile(projectRoot, "geodetic_inversion-master", "geometry"); ...
    fullfile(projectRoot, "geodetic_inversion-master", "OtherFunc")];
for k = 1:numel(paths)
    if exist(paths(k), "dir") == 7
        addpath(char(paths(k)));
    end
end
end

function tracks = load_tracks(invCfg)
definitions = invCfg.tracks(:);
root = string(get_field(invCfg, "samplesRoot", ""));
tracks = repmat(struct("name", "", "type", "", "sampleFile", "", ...
    "sampledData", [], "x", [], "y", [], "data", [], ...
    "projectionEast", [], "projectionNorth", [], "projectionUp", [], ...
    "rowStart", 0, "rowEnd", 0, "sampleCount", 0, ...
    "prediction", [], "residual", []), numel(definitions), 1);
rowStart = 1;
for k = 1:numel(definitions)
    definition = definitions(k);
    if isfield(definition, "sampleFile") && strlength(string(definition.sampleFile)) > 0
        sampleFile = string(definition.sampleFile);
    elseif isfield(definition, "relativeFile")
        sampleFile = fullfile(root, string(definition.relativeFile));
    else
        error("SlipSolve:InvalidFinalConfig", ...
            "Each cfg.finalInversion.tracks entry needs sampleFile or relativeFile.");
    end
    if exist(sampleFile, "file") ~= 2
        error("SlipSolve:MissingInput", "Final inversion sample is missing: %s", sampleFile);
    end
    loaded = load(sampleFile, "sampled_insar_data");
    if ~isfield(loaded, "sampled_insar_data")
        error("SlipSolve:MissingInput", "%s lacks sampled_insar_data.", sampleFile);
    end
    sampled = double(loaded.sampled_insar_data);
    if size(sampled, 2) < 3
        error("SlipSolve:InvalidInput", "%s must contain at least [x y displacement].", sampleFile);
    end
    type = lower(string(definition.type));
    if type == "azo"
        [east, north] = azo_projection(definition, sampleFile);
        east = repmat(east, size(sampled, 1), 1);
        north = repmat(north, size(sampled, 1), 1);
        up = zeros(size(sampled, 1), 1);
    elseif type == "los" || type == "rng"
        if size(sampled, 2) < 6
            error("SlipSolve:InvalidInput", ...
                "%s needs look-vector columns 4:6 for %s projection.", sampleFile, type);
        end
        east = sampled(:, 4);
        north = sampled(:, 5);
        up = sampled(:, 6);
    else
        error("SlipSolve:InvalidFinalConfig", "Unknown final track type: %s", type);
    end
    n = size(sampled, 1);
    tracks(k).name = string(definition.name);
    tracks(k).type = type;
    tracks(k).sampleFile = sampleFile;
    tracks(k).sampledData = sampled;
    tracks(k).x = sampled(:, 1);
    tracks(k).y = sampled(:, 2);
    tracks(k).data = sampled(:, 3);
    tracks(k).projectionEast = east;
    tracks(k).projectionNorth = north;
    tracks(k).projectionUp = up;
    tracks(k).rowStart = rowStart;
    tracks(k).rowEnd = rowStart+n-1;
    tracks(k).sampleCount = n;
    rowStart = rowStart+n;
end
end

function [sinF, cosF] = azo_projection(definition, sampleFile)
if isfield(definition, "sinF") && isfield(definition, "cosF")
    sinF = double(definition.sinF);
    cosF = double(definition.cosF);
    return
end
folder = fileparts(sampleFile);
sinFile = fullfile(folder, "sinF.dat");
cosFile = fullfile(folder, "cosF.dat");
if isfield(definition, "sinFFile"), sinFile = string(definition.sinFFile); end
if isfield(definition, "cosFFile"), cosFile = string(definition.cosFFile); end
if exist(sinFile, "file") ~= 2 || exist(cosFile, "file") ~= 2
    error("SlipSolve:MissingInput", ...
        "AZO sample %s requires sinF.dat and cosF.dat or explicit sinF/cosF.", sampleFile);
end
sinF = readmatrix(sinFile);
cosF = readmatrix(cosFile);
sinF = sinF(1);
cosF = cosF(1);
end

function [xVector, yVector] = full_grid_vectors(cfg)
trace = readmatrix(cfg.compositeGreens.meshTraceFile);
trace = trace(:, 1:2);
trace = trace(all(isfinite(trace), 2), :);
[~, ~, xVector, yVector] = adaptive_meshgrid(trace, ...
    cfg.compositeGreens.xRange, cfg.compositeGreens.yRange, ...
    cfg.compositeGreens.D1, cfg.compositeGreens.D2, cfg.compositeGreens.M);
xVector = double(xVector(:).');
yVector = double(yVector(:).');
end

function validate_component_file(filePath, variableName, expectedRows, expectedColumns)
if exist(filePath, "file") ~= 2
    error("SlipSolve:MissingInput", "Composite component is missing: %s", filePath);
end
info = whos("-file", filePath, variableName);
if isempty(info)
    error("SlipSolve:MissingInput", "%s lacks variable %s.", filePath, variableName);
end
if ~isequal(info.size, [expectedRows expectedColumns])
    error("SlipSolve:GreenSizeMismatch", ...
        "%s is %s; expected [%d %d] for this grid/model.", ...
        variableName, mat2str(info.size), expectedRows, expectedColumns);
end
end

function model = prepare_hybrid_model(mesh)
sections = mesh.sections(:);
contexts = repmat(struct("name", "", "original", [], "mixed", [], ...
    "hybrid", [], "patches", {{}}, "rowStart", 0, "rowEnd", 0, ...
    "strikeInterpolant", [], "dipInterpolant", []), numel(sections), 1);
hybrid = [];
rowOffset = 0;
for k = 1:numel(sections)
    section = sections(k);
    grid = section.sourceSurfaceGrid;
    [strike, dip] = compute_geometry_pointSource(grid.x, grid.y, grid.z);
    strikeF = scatteredInterpolant(grid.pointX(:), grid.pointY(:), grid.pointZ(:), strike(:));
    dipF = scatteredInterpolant(grid.pointX(:), grid.pointY(:), grid.pointZ(:), dip(:));
    pointsM = double(section.nodes) .* 1e3;
    nodeStrike = strikeF(section.nodes(:, 1), section.nodes(:, 2), section.nodes(:, 3));
    nodeDip = dipF(section.nodes(:, 1), section.nodes(:, 2), section.nodes(:, 3));
    original = make_node_geometry(k, pointsM, double(section.triangles));
    original = [original, nodeStrike, nodeDip]; %#ok<AGROW>
    [mixed, patches] = alter_nodes_geometry_tri2(original);
    mixed = pad_legacy_mixed_model(mixed);
    localHybrid = removeFirstLayer(mixed);
    localHybrid(:, 1) = k;
    valid = isfinite(localHybrid(:, 9:15));
    localHybrid(:, 9:15) = localHybrid(:, 9:15) + rowOffset .* valid;
    localHybrid(:, 2) = rowOffset + (1:size(localHybrid, 1)).';
    contexts(k).name = string(section.name);
    contexts(k).original = original;
    contexts(k).mixed = mixed;
    contexts(k).hybrid = localHybrid;
    contexts(k).patches = patches;
    contexts(k).rowStart = rowOffset+1;
    contexts(k).rowEnd = rowOffset+size(localHybrid, 1);
    contexts(k).strikeInterpolant = strikeF;
    contexts(k).dipInterpolant = dipF;
    hybrid = [hybrid; localHybrid]; %#ok<AGROW>
    rowOffset = size(hybrid, 1);
end
model = struct("hybrid", hybrid, "sections", contexts);
end

function mixed = pad_legacy_mixed_model(mixed)
if size(mixed, 2) < 27
    difference = 27-size(mixed, 2);
    supplement = nan(size(mixed, 1), difference);
    mixed = [mixed(:, 1:8), mixed(:, 9:16-difference), ...
        supplement, mixed(:, 16-difference+1:end)];
end
if size(mixed, 2) ~= 27
    error("SlipSolve:LegacyGeometryMismatch", ...
        "Altered model has %d columns; legacy final inversion requires 27.", size(mixed, 2));
end
end

function [greenRaw, cacheInfo] = sampled_composite_green(cfg, tracks, xVector, yVector, nColumns)
p = slipsolve.project.paths(cfg);
settings = cfg.finalInversion.interpolation;
cacheFile = p.sampledCompositeGreens;
sampleCounts = [tracks.sampleCount].';
trackNames = string({tracks.name}).';
expectedRows = sum(sampleCounts);
signature = cache_signature(cfg, tracks, xVector, yVector, nColumns);
if get_field(settings, "reuseCached", true) && exist(cacheFile, "file") == 2
    cached = load(cacheFile, "sampledComposite");
    if isfield(cached, "sampledComposite") && ...
            isequal(cached.sampledComposite.matrixSize, [expectedRows nColumns]) && ...
            isequal(cached.sampledComposite.trackNames, trackNames) && ...
            isequal(cached.sampledComposite.sampleCounts, sampleCounts) && ...
            string(cached.sampledComposite.method) == string(settings.method) && ...
            isfield(cached.sampledComposite, "signature") && ...
            cache_signatures_match(cached.sampledComposite.signature, signature)
        fprintf("Loading cached sampled composite Green's function: %s\n", cacheFile);
        loaded = load(cacheFile, "G_raw");
        greenRaw = loaded.G_raw;
        cacheInfo = cached.sampledComposite;
        cacheInfo.reused = true;
        return
    end
end
if exist(cacheFile, "file") == 2 && ~get_field(settings, "overwriteCached", false)
    error("SlipSolve:ProtectedOutput", ...
        "Sampled Green cache does not match current inputs: %s. Set interpolation.overwriteCached=true.", cacheFile);
end

x = vertcat(tracks.x);
y = vertcat(tracks.y);
coefficients = [vertcat(tracks.projectionEast), ...
    vertcat(tracks.projectionNorth), vertcat(tracks.projectionUp)];
files = cfg.finalInversion.componentFiles;
greenRaw = slipsolve.composite.interpolate_projected_greens( ...
    files, xVector, yVector, x, y, coefficients, nColumns, settings);
sampledComposite = struct();
sampledComposite.method = string(settings.method);
sampledComposite.matrixSize = size(greenRaw);
sampledComposite.trackNames = trackNames;
sampledComposite.sampleCounts = sampleCounts;
sampledComposite.componentFiles = files;
sampledComposite.gridSize = [numel(yVector), numel(xVector)];
sampledComposite.columnOrder = "all_faults_strike_slip_then_all_faults_dip_slip";
sampledComposite.signature = signature;
sampledComposite.createdAt = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");
sampledComposite.reused = false;
G_raw = greenRaw;
save(cacheFile, "G_raw", "sampledComposite", "-v7.3");
cacheInfo = sampledComposite;
end

function signature = cache_signature(cfg, tracks, xVector, yVector, nColumns)
sampleFiles = string({tracks.sampleFile}).';
componentFiles = [string(cfg.finalInversion.componentFiles.east); ...
    string(cfg.finalInversion.componentFiles.north); ...
    string(cfg.finalInversion.componentFiles.vertical)];
signature = struct();
signature.sampleFiles = file_signatures(sampleFiles);
signature.componentFiles = file_signatures(componentFiles);
signature.gridX = [numel(xVector), min(xVector), max(xVector), sum(xVector), sum(xVector.^2)];
signature.gridY = [numel(yVector), min(yVector), max(yVector), sum(yVector), sum(yVector.^2)];
signature.projectionSums = [sum(vertcat(tracks.projectionEast)), ...
    sum(vertcat(tracks.projectionNorth)), sum(vertcat(tracks.projectionUp))];
signature.nColumns = nColumns;
end

function signatures = file_signatures(files)
signatures = repmat(struct("name", "", "bytes", 0), numel(files), 1);
for k = 1:numel(files)
    details = dir(files(k));
    [~, name, extension] = fileparts(files(k));
    signatures(k).name = name+extension;
    signatures(k).bytes = details.bytes;
end
end

function matches = cache_signatures_match(cached, current)
required = ["sampleFiles" "componentFiles" "gridX" "gridY" ...
    "projectionSums" "nColumns"];
if ~all(isfield(cached, required))
    matches = false;
    return
end
matches = isequaln(portable_file_signatures(cached.sampleFiles), ...
        current.sampleFiles) && ...
    isequaln(portable_file_signatures(cached.componentFiles), ...
        current.componentFiles) && ...
    isequaln(cached.gridX, current.gridX) && ...
    isequaln(cached.gridY, current.gridY) && ...
    isequaln(cached.projectionSums, current.projectionSums) && ...
    isequaln(cached.nColumns, current.nColumns);
end

function signatures = portable_file_signatures(input)
signatures = repmat(struct("name", "", "bytes", 0), numel(input), 1);
for k = 1:numel(input)
    if isfield(input, "name")
        signatures(k).name = string(input(k).name);
    else
        [~, name, extension] = fileparts(string(input(k).path));
        signatures(k).name = name+extension;
    end
    signatures(k).bytes = input(k).bytes;
end
end

function [greenWeighted, dataRaw, dataWeighted, tracks] = weight_tracks(invCfg, greenRaw, tracks)
greenWeighted = greenRaw;
dataRaw = vertcat(tracks.data);
dataWeighted = dataRaw;
for k = 1:numel(tracks)
    rows = tracks(k).rowStart:tracks(k).rowEnd;
    typeWeight = get_field(invCfg.trackWeights, char(tracks(k).type), 1);
    normalization = 1;
    if get_field(invCfg, "normalizeEachTrack", true)
        normalization = 1/tracks(k).sampleCount;
    end
    weight = typeWeight * normalization;
    greenWeighted(rows, :) = greenWeighted(rows, :) * weight;
    dataWeighted(rows) = dataWeighted(rows) * weight;
    tracks(k).weight = weight;
end
end

function [matrix, data, info] = boundary_constraints(invCfg, slipModel)
matrix = sparse(0, 2*size(slipModel, 1));
data = zeros(0, 1);
info = struct("name", {}, "faultId", {}, "ratio", {}, "rows", {});
faultIds = unique(slipModel(:, 1)).';
names = ["bottom", "top", "left", "right"];
for name = names
    settings = invCfg.boundaryConditions.(name);
    if ~get_field(settings, "enabled", false)
        continue
    end
    for faultId = faultIds
        if name == "bottom"
            option = max(slipModel(slipModel(:, 1) == faultId, 3));
        elseif name == "top"
            option = min(slipModel(slipModel(:, 1) == faultId, 3));
        else
            option = char(name);
        end
        [block, target] = zero_slip_boundary(slipModel, faultId, option, settings.ratio);
        matrix = [matrix; sparse(block)]; %#ok<AGROW>
        data = [data; target]; %#ok<AGROW>
        info(end+1) = struct("name", name, "faultId", faultId, ...
            "ratio", settings.ratio, "rows", size(block, 1)); %#ok<AGROW>
    end
end
end

function [lower, upper] = inversion_bounds(invCfg, slipModel)
n = size(slipModel, 1);
maxSS = invCfg.constraints.maxStrikeSlipCm;
maxDS = invCfg.constraints.maxDipSlipCm;
lower = [-maxSS*ones(n, 1); -maxDS*ones(n, 1)];
upper = [ maxSS*ones(n, 1);  maxDS*ones(n, 1)];
polarity = invCfg.constraints.polarity;
faultIds = unique(slipModel(:, 1)).';
for index = 1:numel(faultIds)
    constraint = polarity;
    if index > 1
        constraint = [polarity(1), -polarity(2), polarity(3)];
    end
    rows = find(slipModel(:, 1) == faultIds(index));
    if constraint(1) > 0, lower(rows) = 0; end
    if constraint(1) < 0, upper(rows) = 0; end
    if constraint(2) > 0, lower(rows+n) = 0; end
    if constraint(2) < 0, upper(rows+n) = 0; end
end
if ~isempty(invCfg.constraints.customLowerBounds)
    lower = invCfg.constraints.customLowerBounds(:);
end
if ~isempty(invCfg.constraints.customUpperBounds)
    upper = invCfg.constraints.customUpperBounds(:);
end
end

function tracks = attach_predictions(tracks, prediction)
for k = 1:numel(tracks)
    rows = tracks(k).rowStart:tracks(k).rowEnd;
    tracks(k).prediction = prediction(rows);
    tracks(k).residual = prediction(rows)-tracks(k).data;
    tracks(k).sampledData = [];
    tracks(k).projectionEast = [];
    tracks(k).projectionNorth = [];
    tracks(k).projectionUp = [];
end
end

function [complete, triangles, simple, info] = export_simple_triangular(cfg, invCfg, model)
complete = zeros(0, 17);
triangles = zeros(0, 7);
strikeAll = zeros(0, 1);
dipAll = zeros(0, 1);
for k = 1:numel(model.sections)
    context = model.sections(k);
    rows = context.rowStart:context.rowEnd;
    localHybrid = model.hybrid(rows, :);
    localHybrid(:, 2) = (1:size(localHybrid, 1)).';
    valid = isfinite(localHybrid(:, 9:15));
    localHybrid(:, 9:15) = localHybrid(:, 9:15)-context.rowStart+1 .* valid;
    completeLocal = updateNodes(context.original, context.mixed, localHybrid, context.patches);
    mixedUpdated = updateTriangleSlip(context.mixed, localHybrid);
    mixedUpdated(mixedUpdated(:, 3) >= 2, 4:5) = localHybrid(:, 4:5);
    triangleLocal = node_to_triangle(mixedUpdated);
    strike = context.strikeInterpolant(triangleLocal(:, 1)/1e3, ...
        triangleLocal(:, 2)/1e3, triangleLocal(:, 3)/1e3);
    dip = context.dipInterpolant(triangleLocal(:, 1)/1e3, ...
        triangleLocal(:, 2)/1e3, triangleLocal(:, 3)/1e3);
    completeOffset = size(complete, 1);
    completeLocal(:, 2) = completeOffset+(1:size(completeLocal, 1)).';
    validComplete = isfinite(completeLocal(:, 9:15));
    completeLocal(:, 9:15) = completeLocal(:, 9:15) + completeOffset.*validComplete;
    complete = [complete; completeLocal]; %#ok<AGROW>
    triangles = [triangles; triangleLocal]; %#ok<AGROW>
    strikeAll = [strikeAll; strike]; %#ok<AGROW>
    dipAll = [dipAll; dip]; %#ok<AGROW>
end

velocity = load(invCfg.velocityModelFile, "vel");
vel = velocity.vel;
vel(:, 1:2) = -vel(:, 1:2);
if isempty(invCfg.uniformShearModulusPa)
    [M0, Mw, mu] = triangle_moment_nodes_tri(triangles, vel);
else
    [M0, Mw, mu] = triangle_moment_nodes_tri( ...
        triangles, vel, "uniform_shear", invCfg.uniformShearModulusPa);
end
origin = cfg.insar.origin;
[originX, originY] = ll2xy(origin.lon, origin.lat, origin.refLon);
[longitude, latitude] = xy2ll(triangles(:, 1)+originX, ...
    triangles(:, 2)+originY, origin.refLon);
elementNumber = (1:size(triangles, 1)).';
simple = [elementNumber, longitude, latitude, triangles(:, 3), ...
    strikeAll, dipAll, triangles(:, 6)/100, triangles(:, 7)/100, ...
    triangles(:, 5), mu];
columns = ["element_id", "longitude_deg", "latitude_deg", "depth_m", ...
    "strike_deg", "dip_deg", "strike_slip_m", "dip_slip_m", ...
    "area_m2", "shear_modulus_pa"];
tri_simple = simple;
column_names = columns;
save(invCfg.simpleTriangularMatFile, "tri_simple", "column_names", "M0", "Mw", "-v7");
fid = fopen(invCfg.simpleTriangularTextFile, "w");
if fid < 0
    error("SlipSolve:OutputFailed", "Cannot write %s.", invCfg.simpleTriangularTextFile);
end
fprintf(fid, "%s\n", strjoin(columns, sprintf("\t")));
fclose(fid);
writematrix(simple, invCfg.simpleTriangularTextFile, ...
    "Delimiter", "tab", "WriteMode", "append");
info = struct("columns", columns, "matFile", string(invCfg.simpleTriangularMatFile), ...
    "textFile", string(invCfg.simpleTriangularTextFile), "M0", M0, "Mw", Mw);
end

function out = node_to_triangle(slipData)
results = zeros(0, 7);
triangleRows = find(slipData(:, 3) == 1);
for row = triangleRows.'
    points = reshape(slipData(row, 19:27), 3, 3).';
    center = mean(points, 1);
    strikeSlip = slipData(row, 4);
    dipSlip = slipData(row, 5);
    area = 0.5*norm(cross(points(2, :)-points(1, :), points(3, :)-points(1, :)));
    results(end+1, :) = [center, hypot(strikeSlip, dipSlip), area, strikeSlip, dipSlip]; %#ok<AGROW>
end
nodeRows = find(slipData(:, 3) >= 2);
nNode = numel(nodeRows);
globalToLocal = zeros(size(slipData, 1), 1);
globalToLocal(nodeRows) = 1:nNode;
adjacency = false(nNode);
for k = 1:nNode
    row = nodeRows(k);
    connected = slipData(row, 9:16);
    connected = connected(isfinite(connected) & connected >= 1 & connected <= size(slipData, 1));
    connected = connected(slipData(connected, 3) >= 2);
    adjacency(k, globalToLocal(connected)) = true;
end
adjacency = adjacency | adjacency.';
for i = 1:nNode-2
    neighborsI = find(adjacency(i, i+1:end))+i;
    for j = neighborsI
        neighborsJ = find(adjacency(j, j+1:end))+j;
        common = intersect(neighborsI, neighborsJ);
        for k = common
            rows = nodeRows([i j k]);
            points = slipData(rows, 6:8);
            strikeSlip = mean(slipData(rows, 4));
            dipSlip = mean(slipData(rows, 5));
            area = 0.5*norm(cross(points(2, :)-points(1, :), points(3, :)-points(1, :)));
            results(end+1, :) = [mean(points, 1), hypot(strikeSlip, dipSlip), ...
                area, strikeSlip, dipSlip]; %#ok<AGROW>
        end
    end
end
out = results;
end

function value = get_field(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
