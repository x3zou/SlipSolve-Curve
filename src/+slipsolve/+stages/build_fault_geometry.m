function faultGeometry = build_fault_geometry(cfg, insarSub)
%BUILD_FAULT_GEOMETRY Build independent smooth fault surfaces.

if nargin < 2
    insarSub = [];
end

p = slipsolve.project.paths(cfg);
slipsolve.project.ensure_directories(cfg);

faultDefs = normalize_fault_definitions(cfg);
faultCells = cell(numel(faultDefs), 1);
for k = 1:numel(faultDefs)
    faultCells{k} = build_one_fault(cfg, faultDefs(k), insarSub, k);
end
faults = vertcat(faultCells{:});

faultGeometry = struct();
faultGeometry.faults = faults;
faultGeometry.metadata.origin = cfg.insar.origin;
faultGeometry.metadata.units.horizontal = "km";
faultGeometry.metadata.units.depth = "km";
faultGeometry.metadata.faultCount = numel(faults);
faultGeometry.metadata.maxDepthKm = max([faults.maxDepthKm]);

% Backward-compatible aliases for single-fault consumers while downstream
% stages are being migrated to faultGeometry.faults(k).
faultGeometry.surfaceTraces = faults(1).surfaceTraces;
faultGeometry.depthControls = faults(1).depthControls;
faultGeometry.surfaceGrid = faults(1).surfaceGrid;
faultGeometry.figures = slipsolve.plotting.plot_fault_geometry(cfg, faultGeometry);

save(p.faultGeometry, "faultGeometry", "-v7.3");

end

function faultDefs = normalize_fault_definitions(cfg)
if isfield(cfg.geometry, "faults") && ~isempty(cfg.geometry.faults)
    rawFaults = cfg.geometry.faults(:);
else
    rawFaults = struct();
end

if isempty(rawFaults) || (~isfield(rawFaults, "surfaceTraceFiles") && ~isfield(rawFaults, "surfaceTraceFile"))
    rawFaults = struct();
    rawFaults.name = "fault_1";
    rawFaults.surfaceTraceFiles = cfg.geometry.surfaceTraceFiles;
end

faultCells = cell(numel(rawFaults), 1);
for k = 1:numel(rawFaults)
    raw = rawFaults(k);
    fault = struct();
    fault.name = string(get_field(raw, "name", "fault_" + k));
    fault.mode = string(get_field(raw, "mode", cfg.geometry.mode));

    if isfield(raw, "surfaceTraceFile")
        fault.surfaceTraceFiles = string(raw.surfaceTraceFile);
    else
        fault.surfaceTraceFiles = string(get_field(raw, "surfaceTraceFiles", cfg.geometry.surfaceTraceFiles));
    end

    fault.depthControlFiles = string(get_field(raw, "depthControlFiles", cfg.geometry.depthControlFiles));
    fault.uniformDipDegrees = get_field(raw, "uniformDipDegrees", cfg.geometry.uniformDipDegrees);
    fault.segmentFiles = string(get_field(raw, "segmentFiles", cfg.geometry.segmentFiles));
    fault.segmentDipDegrees = get_field(raw, "segmentDipDegrees", cfg.geometry.segmentDipDegrees);
    fault.maxDepth = get_field(raw, "maxDepth", cfg.geometry.maxDepth);
    fault.traceResampleSpacing = get_field(raw, "traceResampleSpacing", cfg.geometry.traceResampleSpacing);
    fault.downDipLevels = get_field(raw, "downDipLevels", cfg.geometry.downDipLevels);
    fault.surfaceFitSmoothness = get_field(raw, "surfaceFitSmoothness", cfg.geometry.surfaceFitSmoothness);
    fault.legacyExact = get_field(raw, "legacyExact", get_field(cfg.geometry, "legacyExact", struct()));

    validate_fault_definition(fault, k);
    faultCells{k} = fault;
end

faultDefs = vertcat(faultCells{:});
end

function validate_fault_definition(fault, index)
if isempty(fault.surfaceTraceFiles) || strlength(fault.surfaceTraceFiles(1)) == 0
    error("SlipSolve:MissingInput", "Fault %d must define at least one surface trace file.", index);
end

if fault.mode == "segment_dips" && numel(fault.segmentFiles) ~= numel(fault.segmentDipDegrees)
    error("SlipSolve:InvalidGeometry", ...
        "Fault %d has %d segment files but %d segment dip values.", ...
        index, numel(fault.segmentFiles), numel(fault.segmentDipDegrees));
end
end

function fault = build_one_fault(cfg, faultDef, insarSub, index)
surfaceTraceFiles = string(faultDef.surfaceTraceFiles);
surfaceTraces = cell(numel(surfaceTraceFiles), 1);
for k = 1:numel(surfaceTraceFiles)
    surfaceTraces{k} = read_surface_trace(surfaceTraceFiles(k), cfg, faultDef);
end

surfaceTrace = surfaceTraces{1};
depthControls = build_depth_controls(surfaceTrace, faultDef);
if faultDef.mode == "legacy_fit_curved_planeMM"
    surfaceGrid = build_legacy_exact_surface_grid(surfaceTraceFiles(1), cfg, faultDef);
else
    surfaceGrid = build_smooth_surface_grid(surfaceTrace, depthControls, cfg, faultDef);
end

fault = struct();
fault.name = faultDef.name;
fault.index = index;
fault.surfaceTraces = surfaceTraces;
fault.depthControls = depthControls;
fault.surfaceGrid = surfaceGrid;
fault.mode = faultDef.mode;
fault.maxDepthKm = faultDef.maxDepth ./ 1e3;
fault.uniformDipDegrees = faultDef.uniformDipDegrees;
fault.segmentDipDegrees = faultDef.segmentDipDegrees;
fault.surfaceFitSmoothness = faultDef.surfaceFitSmoothness;
fault.insarExtentAvailable = ~isempty(insarSub);
end

function surfaceGrid = build_legacy_exact_surface_grid(traceFile, cfg, faultDef)
add_legacy_curved_geometry_paths(cfg);

[xf1, yf1, ~] = read_data_from_file(char(traceFile), ...
    cfg.insar.origin.lon, cfg.insar.origin.lat, cfg.insar.origin.refLon);

legacy = faultDef.legacyExact;
miniKm = get_field(legacy, "miniKm", 5);
depthKm = get_field(legacy, "depthKm", faultDef.maxDepth ./ 1e3);
depthLevels = get_field(legacy, "depthLevels", 30);
surfaceFactor = get_field(legacy, "surfaceInterpolationFactor", 1e5);
bottomFactor = get_field(legacy, "bottomInterpolationFactor", 1e2);

if get_field(legacy, "useReferenceProducts", false)
    surfaceGrid = load_legacy_reference_surface_grid(legacy, faultDef);
    return
end

yq1 = linspace(min(yf1), max(yf1), length(yf1) .* surfaceFactor);
xq1 = interp1(yf1, xf1, yq1, "spline");
xOne = xq1(:);
yOne = yq1(:);
depthOne = zeros(length(yOne), 1);

for i = 1:numel(faultDef.segmentDipDegrees)
    segment = load(faultDef.segmentFiles(i));
    xBegin = segment(1) ./ 1e3;
    yBegin = segment(2) ./ 1e3;
    xEnd = segment(3) ./ 1e3;
    yEnd = segment(4) ./ 1e3;

    projected = project_segment_3d(xBegin, yBegin, xEnd, yEnd, depthKm, faultDef.segmentDipDegrees(i));
    yBottom = linspace(min(projected(:, 2)), max(projected(:, 2)), length(yf1) .* bottomFactor);
    xBottom = interp1(projected(:, 2), projected(:, 1), yBottom);
    dBottom = repmat(-depthKm, length(xBottom), 1);

    xOne = [xOne; xBottom(:)]; %#ok<AGROW>
    yOne = [yOne; yBottom(:)]; %#ok<AGROW>
    depthOne = [depthOne; dBottom]; %#ok<AGROW>
end

depthGrid = linspace(-depthKm, 0, depthLevels);
ny1 = ceil((max(yf1) - min(yf1)) ./ miniKm);
y1 = linspace(min(yf1), max(yf1), ny1);

[y1m, d1m, x1m, y1p, d1p, x1p] = gridFitInterpolate(yOne, depthOne, xOne, ...
    y1, depthGrid, faultDef.surfaceFitSmoothness);

surfaceGrid = struct();
surfaceGrid.x = x1m;
surfaceGrid.y = y1m;
surfaceGrid.z = d1m;
surfaceGrid.pointX = x1p;
surfaceGrid.pointY = y1p;
surfaceGrid.pointZ = d1p;
surfaceGrid.fitMethod = "legacy_fit_curved_planeMM";
surfaceGrid.fitCrossCoordinate = "x";
surfaceGrid.smoothness = faultDef.surfaceFitSmoothness;
surfaceGrid.controlPoints = zeros(0, 3);
surfaceGrid.previewOnly = false;
surfaceGrid.legacyExact = struct();
surfaceGrid.legacyExact.x1 = linspace(min(xf1), max(xf1), ceil((max(xf1) - min(xf1)) ./ miniKm));
surfaceGrid.legacyExact.y1 = y1;
surfaceGrid.legacyExact.d = depthGrid;
surfaceGrid.legacyExact.miniKm = miniKm;
surfaceGrid.legacyExact.depthKm = depthKm;
surfaceGrid.legacyExact.depthLevels = depthLevels;
surfaceGrid.legacyExact.surfaceInterpolationFactor = surfaceFactor;
surfaceGrid.legacyExact.bottomInterpolationFactor = bottomFactor;
surfaceGrid.legacyExact.controlPointCount = numel(xOne);
end

function surfaceGrid = load_legacy_reference_surface_grid(legacy, faultDef)
meshFile = string(get_field(legacy, "referenceMeshFile", ""));
pointFile = string(get_field(legacy, "referencePointFile", ""));
if strlength(meshFile) == 0 || exist(meshFile, "file") ~= 2
    error("SlipSolve:MissingInput", "Missing legacy reference mesh file: %s", meshFile);
end
if strlength(pointFile) == 0 || exist(pointFile, "file") ~= 2
    error("SlipSolve:MissingInput", "Missing legacy reference point file: %s", pointFile);
end

meshData = load(meshFile, "x1m", "y1m", "d1m");
pointData = load(pointFile, "x1p", "y1p", "d1p");

surfaceGrid = struct();
surfaceGrid.x = meshData.x1m;
surfaceGrid.y = meshData.y1m;
surfaceGrid.z = meshData.d1m;
surfaceGrid.pointX = pointData.x1p;
surfaceGrid.pointY = pointData.y1p;
surfaceGrid.pointZ = pointData.d1p;
surfaceGrid.fitMethod = "legacy_reference_fit_curved_planeMM";
surfaceGrid.fitCrossCoordinate = "x";
surfaceGrid.smoothness = faultDef.surfaceFitSmoothness;
surfaceGrid.controlPoints = zeros(0, 3);
surfaceGrid.previewOnly = false;
surfaceGrid.legacyExact = struct();
surfaceGrid.legacyExact.source = "reference_products";
surfaceGrid.legacyExact.referenceMeshFile = meshFile;
surfaceGrid.legacyExact.referencePointFile = pointFile;
surfaceGrid.legacyExact.depthKm = abs(min(surfaceGrid.z(:)));
surfaceGrid.legacyExact.depthLevels = size(surfaceGrid.z, 1);
end

function trace = read_surface_trace(traceFile, cfg, faultDef)
addpath(char(cfg.paths.legacyOtherFuncRoot));
data = load(traceFile);
lon = data(:, 1);
lat = data(:, 2);
[originX, originY] = ll2xy(cfg.insar.origin.lon, cfg.insar.origin.lat, cfg.insar.origin.refLon);
[x, y] = ll2xy(lon, lat, cfg.insar.origin.refLon);

trace = struct();
trace.file = string(traceFile);
trace.lon = lon;
trace.lat = lat;
trace.x = (x - originX) ./ 1e3;
trace.y = (y - originY) ./ 1e3;
trace.z = zeros(size(trace.x));
trace.resampled = resample_polyline(trace.x, trace.y, faultDef.traceResampleSpacing ./ 1e3);
end

function depthControls = build_depth_controls(surfaceTrace, faultDef)
if faultDef.mode == "segment_dips" && ~isempty(faultDef.segmentFiles)
    depthControls = build_segment_depth_controls(faultDef);
else
    depthControls = build_uniform_dip_control(surfaceTrace, faultDef);
end
end

function depthControls = build_segment_depth_controls(faultDef)
segmentFiles = string(faultDef.segmentFiles);
dips = faultDef.segmentDipDegrees;
depthControls = repmat(empty_depth_control(), numel(segmentFiles), 1);
depthKm = faultDef.maxDepth ./ 1e3;

for k = 1:numel(segmentFiles)
    segment = load(segmentFiles(k)) ./ 1e3;
    projected = project_segment_3d(segment(1), segment(2), segment(3), segment(4), depthKm, dips(k));

    control = empty_depth_control();
    control.file = segmentFiles(k);
    control.dipDegrees = dips(k);
    control.surfaceX = [segment(1); segment(3)];
    control.surfaceY = [segment(2); segment(4)];
    control.surfaceZ = [0; 0];
    control.bottomX = projected(:, 1);
    control.bottomY = projected(:, 2);
    control.bottomZ = projected(:, 3);
    depthControls(k) = control;
end
end

function depthControls = build_uniform_dip_control(surfaceTrace, faultDef)
xy = surfaceTrace.resampled.xy;
depthKm = faultDef.maxDepth ./ 1e3;
dipDegrees = faultDef.uniformDipDegrees;

bottom = nan(size(xy, 1), 3);
for k = 1:size(xy, 1)-1
    projected = project_segment_3d(xy(k, 1), xy(k, 2), xy(k+1, 1), xy(k+1, 2), depthKm, dipDegrees);
    bottom(k, :) = projected(1, :);
    if k == size(xy, 1)-1
        bottom(k+1, :) = projected(2, :);
    end
end

control = empty_depth_control();
control.dipDegrees = dipDegrees;
control.surfaceX = xy(:, 1);
control.surfaceY = xy(:, 2);
control.surfaceZ = zeros(size(xy, 1), 1);
control.bottomX = bottom(:, 1);
control.bottomY = bottom(:, 2);
control.bottomZ = bottom(:, 3);
depthControls = control;
end

function control = empty_depth_control()
control = struct();
control.file = "";
control.dipDegrees = NaN;
control.surfaceX = [];
control.surfaceY = [];
control.surfaceZ = [];
control.bottomX = [];
control.bottomY = [];
control.bottomZ = [];
end

function surfaceGrid = build_smooth_surface_grid(surfaceTrace, depthControls, cfg, faultDef)
add_gridfit_path(cfg);

top = surfaceTrace.resampled.xy;
bottom = sample_depth_controls(depthControls, faultDef.traceResampleSpacing ./ 1e3);

controlX = [top(:, 1); bottom(:, 1)];
controlY = [top(:, 2); bottom(:, 2)];
controlZ = [zeros(size(top, 1), 1); bottom(:, 3)];

xRange = range(top(:, 1));
yRange = range(top(:, 2));
fitCrossX = yRange >= xRange;

if fitCrossX
    fitAlong = controlY;
    fitCross = controlX;
    alongGrid = linspace(min(fitAlong), max(fitAlong), size(top, 1));
    crossLabel = "x";
else
    fitAlong = controlX;
    fitCross = controlY;
    alongGrid = linspace(min(fitAlong), max(fitAlong), size(top, 1));
    crossLabel = "y";
end

depthGrid = linspace(-faultDef.maxDepth ./ 1e3, 0, faultDef.downDipLevels);
crossGrid = gridfit(fitAlong, controlZ, fitCross, alongGrid, depthGrid, ...
    'smoothness', faultDef.surfaceFitSmoothness);
[alongMesh, depthMesh] = meshgrid(alongGrid, depthGrid);

if fitCrossX
    surfaceGrid.x = crossGrid;
    surfaceGrid.y = alongMesh;
else
    surfaceGrid.x = alongMesh;
    surfaceGrid.y = crossGrid;
end
surfaceGrid.z = depthMesh;
surfaceGrid.fitMethod = "gridfit";
surfaceGrid.fitCrossCoordinate = crossLabel;
surfaceGrid.smoothness = faultDef.surfaceFitSmoothness;
surfaceGrid.controlPoints = [controlX, controlY, controlZ];
surfaceGrid.previewOnly = false;

[pointAlong, pointDepth] = center_grid(alongGrid, depthGrid);
[pointAlongMesh, pointDepthMesh] = meshgrid(pointAlong, pointDepth);
pointCross = interp2(alongMesh, depthMesh, crossGrid, pointAlongMesh, pointDepthMesh, "linear");
if fitCrossX
    surfaceGrid.pointX = pointCross;
    surfaceGrid.pointY = pointAlongMesh;
else
    surfaceGrid.pointX = pointAlongMesh;
    surfaceGrid.pointY = pointCross;
end
surfaceGrid.pointZ = pointDepthMesh;
end

function add_gridfit_path(cfg)
gridfitPath = fullfile(string(cfg.paths.legacyCurvedGeometryRoot), "gridfitdir", "gridfitdir");
if exist(gridfitPath, "dir") == 7
    addpath(char(gridfitPath));
end
end

function add_legacy_curved_geometry_paths(cfg)
if isfield(cfg.paths, "legacyRoot") && strlength(string(cfg.paths.legacyRoot)) > 0 ...
        && exist(cfg.paths.legacyRoot, "dir") == 7
    addpath(genpath(char(cfg.paths.legacyRoot)));
end
if isfield(cfg.paths, "legacyCurvedGeometryRoot") && exist(cfg.paths.legacyCurvedGeometryRoot, "dir") == 7
    addpath(char(cfg.paths.legacyCurvedGeometryRoot));
end
add_gridfit_path(cfg);
end

function bottom = sample_depth_controls(depthControls, spacingKm)
bottom = [];
for k = 1:numel(depthControls)
    xy = resample_polyline(depthControls(k).bottomX(:), depthControls(k).bottomY(:), spacingKm).xy;
    z = depthControls(k).bottomZ(1) .* ones(size(xy, 1), 1);
    bottom = [bottom; xy, z]; %#ok<AGROW>
end
[~, keep] = unique(bottom(:, 1:2), "rows", "stable");
bottom = bottom(keep, :);
end

function [xCenter, yCenter] = center_grid(x, y)
xCenter = movmean(x, 2, "Endpoints", "discard");
yCenter = movmean(y, 2, "Endpoints", "discard");
end

function projected = project_segment_3d(x1, y1, x2, y2, depthKm, dipDegrees)
strike = [x2 - x1, y2 - y1];
strikeLength = norm(strike);
if strikeLength == 0
    error("SlipSolve:InvalidGeometry", "Cannot project a zero-length segment.");
end

s = strike ./ strikeLength;
dipDirection = [-s(2), s(1)];
horizontalOffset = depthKm ./ tan(deg2rad(dipDegrees));
projected = [
    x1 + horizontalOffset .* dipDirection(1), y1 + horizontalOffset .* dipDirection(2), -depthKm
    x2 + horizontalOffset .* dipDirection(1), y2 + horizontalOffset .* dipDirection(2), -depthKm
];
end

function out = resample_polyline(x, y, spacingKm)
xy = [x(:), y(:)];
dist = [0; cumsum(sqrt(sum(diff(xy, 1, 1).^2, 2)))];
if dist(end) == 0
    sampleDist = 0;
else
    sampleDist = (0:spacingKm:dist(end)).';
    if sampleDist(end) < dist(end)
        sampleDist = [sampleDist; dist(end)];
    end
end

out = struct();
out.distance = sampleDist;
out.xy = [interp1(dist, xy(:, 1), sampleDist), interp1(dist, xy(:, 2), sampleDist)];
end

function value = get_field(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
