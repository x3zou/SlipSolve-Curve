function mesh = build_triangular_mesh(cfg, faultGeometry)
%BUILD_TRIANGULAR_MESH Build triangular patches from fitted fault surfaces.

p = slipsolve.project.paths(cfg);
slipsolve.project.ensure_directories(cfg);

if nargin < 2 || isempty(faultGeometry)
    loaded = load(p.faultGeometry, "faultGeometry");
    faultGeometry = loaded.faultGeometry;
end

if isfield(faultGeometry, "faults")
    faults = faultGeometry.faults;
else
    faults = faultGeometry;
end

activeFaults = resolve_active_faults(cfg, faults);
sections = cell(numel(activeFaults), 1);
nodes = [];
triangles = [];
faultId = [];
sectionId = [];
nodeOffset = 0;

for i = 1:numel(activeFaults)
    k = activeFaults(i);
    section = mesh_one_fault(faults(k), k, i, cfg);
    sections{i} = section;

    nodes = [nodes; section.nodes]; %#ok<AGROW>
    triangles = [triangles; section.triangles + nodeOffset]; %#ok<AGROW>
    faultId = [faultId; repmat(k, size(section.triangles, 1), 1)]; %#ok<AGROW>
    sectionId = [sectionId; repmat(i, size(section.triangles, 1), 1)]; %#ok<AGROW>
    nodeOffset = nodeOffset + size(section.nodes, 1);
end

mesh = struct();
mesh.method = string(cfg.mesh.method);
mesh.units.horizontal = "km";
mesh.units.depth = "km";
mesh.nodes = nodes;
mesh.triangles = triangles;
mesh.faultId = faultId;
mesh.sectionId = sectionId;
mesh.sections = vertcat(sections{:});
mesh = add_triangle_metrics(mesh);
mesh.neighbors = triangle_neighbors(mesh.triangles);
mesh.metadata.faultCount = numel(faults);
mesh.metadata.activeFaults = activeFaults(:);
mesh.metadata.nodeCount = size(mesh.nodes, 1);
mesh.metadata.triangleCount = size(mesh.triangles, 1);
mesh.metadata.createdBy = "slipsolve.stages.build_triangular_mesh";
mesh.figures = slipsolve.plotting.plot_triangular_mesh(cfg, mesh, faultGeometry);

save(p.triangularMesh, "mesh", "-v7.3");

end

function activeFaults = resolve_active_faults(cfg, faults)
if ~isfield(cfg.mesh, "activeFaults") || isempty(cfg.mesh.activeFaults)
    activeFaults = 1:numel(faults);
    return
end

requested = cfg.mesh.activeFaults;
if (isstring(requested) || ischar(requested)) && isscalar(string(requested)) && string(requested) == "all"
    activeFaults = 1:numel(faults);
elseif isnumeric(requested)
    activeFaults = requested(:).';
elseif isstring(requested) || ischar(requested)
    names = string({faults.name});
    requested = string(requested);
    activeFaults = zeros(1, numel(requested));
    for i = 1:numel(requested)
        match = find(names == requested(i), 1);
        if isempty(match)
            error("SlipSolve:InvalidMeshConfig", "cfg.mesh.activeFaults names an unknown fault: %s", requested(i));
        end
        activeFaults(i) = match;
    end
else
    error("SlipSolve:InvalidMeshConfig", "cfg.mesh.activeFaults must be ""all"", numeric indices, or fault names.");
end

if any(activeFaults < 1) || any(activeFaults > numel(faults))
    error("SlipSolve:InvalidMeshConfig", "cfg.mesh.activeFaults contains an out-of-range fault index.");
end
end

function section = mesh_one_fault(fault, faultIndex, sectionIndex, cfg)
method = string(cfg.mesh.method);
switch method
    case "legacy_exact_main_interpolate"
        section = mesh_one_fault_legacy_exact(fault, faultIndex, sectionIndex, cfg);
    case "legacy_depth_dependent_delaunay"
        section = mesh_one_fault_legacy(fault, faultIndex, sectionIndex, cfg);
    case "structured_surface_grid"
        section = mesh_one_fault_structured(fault, faultIndex, sectionIndex, cfg);
    otherwise
        error("SlipSolve:InvalidMeshConfig", "Unknown cfg.mesh.method: %s", method);
end
end

function section = mesh_one_fault_legacy_exact(fault, faultIndex, sectionIndex, cfg)
gridData = fault.surfaceGrid;
[meshGrid, pointGrid, parameterName] = legacy_parameter_grids(gridData);

if isfield(cfg.mesh, "legacyExact") && get_field(cfg.mesh.legacyExact, "useReferenceProducts", false)
    section = mesh_one_fault_legacy_reference(fault, faultIndex, sectionIndex, cfg, pointGrid, parameterName);
    return
end

depthKm = get_mesh_depth_km(cfg, fault);
[along, cross, depth] = legacy_depth_dependent_points_exact(meshGrid, pointGrid, cfg, depthKm);
dt = delaunayTriangulation(along, depth);

if parameterName == "y"
    nodes = [cross, along, depth];
else
    nodes = [along, cross, depth];
end

pointTriangleId = pointLocation(dt, pointGrid(:, 1:2));
dtNeighbors = neighbors(dt, (1:size(dt.ConnectivityList, 1)).');

section = struct();
section.name = string(fault.name);
section.faultIndex = faultIndex;
section.sectionIndex = sectionIndex;
section.nodes = nodes;
section.triangles = dt.ConnectivityList;
section.gridSize = size(gridData.x);
section.sourceSurfaceGrid = gridData;
section.parameterName = parameterName;
section.parameterCoordinates = [along, depth];
section.pointSources = struct();
section.pointSources.parameterCoordinates = pointGrid(:, 1:2);
section.pointSources.crossCoordinate = pointGrid(:, 3);
section.pointSources.triangleId = pointTriangleId;
section.pointSources.coordinates = point_source_coordinates(pointGrid, parameterName);
section.legacy = struct();
section.legacy.biasL = cfg.mesh.biasL;
section.legacy.biasW = cfg.mesh.biasW;
section.legacy.depthRatio = cfg.mesh.depthRatio;
section.legacy.depthKm = depthKm;
section.legacy.interpolationMethod = string(cfg.mesh.interpolationMethod);
section.legacy.DTPoints = dt.Points;
section.legacy.DTConnectivityList = dt.ConnectivityList;
section.legacy.ID1 = pointTriangleId;
section.legacy.n1 = dtNeighbors;
section.legacy.Vx1 = nodes(:, 1);
section = add_triangle_metrics(section);
section.neighbors = dtNeighbors;
end

function section = mesh_one_fault_legacy_reference(fault, faultIndex, sectionIndex, cfg, pointGrid, parameterName)
referenceFile = string(get_field(cfg.mesh.legacyExact, "referenceGeometryFile", ""));
if strlength(referenceFile) == 0 || exist(referenceFile, "file") ~= 2
    error("SlipSolve:MissingInput", "Missing legacy reference geometry file: %s", referenceFile);
end

legacyData = load(referenceFile, "ID1", "n1", "DT1", "Vx1");
nodes = [legacyData.Vx1, legacyData.DT1.Points(:, 1), legacyData.DT1.Points(:, 2)];

section = struct();
section.name = string(fault.name);
section.faultIndex = faultIndex;
section.sectionIndex = sectionIndex;
section.nodes = nodes;
section.triangles = legacyData.DT1.ConnectivityList;
section.gridSize = size(fault.surfaceGrid.x);
section.sourceSurfaceGrid = fault.surfaceGrid;
section.parameterName = parameterName;
section.parameterCoordinates = legacyData.DT1.Points;
section.pointSources = struct();
section.pointSources.parameterCoordinates = pointGrid(:, 1:2);
section.pointSources.crossCoordinate = pointGrid(:, 3);
section.pointSources.triangleId = legacyData.ID1;
section.pointSources.coordinates = point_source_coordinates(pointGrid, parameterName);
section.legacy = struct();
section.legacy.source = "reference_products";
section.legacy.referenceGeometryFile = referenceFile;
section.legacy.DTPoints = legacyData.DT1.Points;
section.legacy.DTConnectivityList = legacyData.DT1.ConnectivityList;
section.legacy.ID1 = legacyData.ID1;
section.legacy.n1 = legacyData.n1;
section.legacy.Vx1 = legacyData.Vx1;
section = add_triangle_metrics(section);
section.neighbors = legacyData.n1;
end

function section = mesh_one_fault_legacy(fault, faultIndex, sectionIndex, cfg)
gridData = fault.surfaceGrid;
[meshGrid, pointGrid, parameterName] = legacy_parameter_grids(gridData);

depthKm = get_mesh_depth_km(cfg, fault);
[along, cross, depth] = legacy_depth_dependent_points(meshGrid, pointGrid, cfg, depthKm);
dt = delaunayTriangulation(along, depth);

if parameterName == "y"
    nodes = [cross, along, depth];
else
    nodes = [along, cross, depth];
end

pointAlong = pointGrid(:, 1);
pointDepth = pointGrid(:, 2);
pointCross = pointGrid(:, 3);
pointTriangleId = pointLocation(dt, pointAlong, pointDepth);
dtNeighbors = neighbors(dt, (1:size(dt.ConnectivityList, 1)).');

section = struct();
section.name = string(fault.name);
section.faultIndex = faultIndex;
section.sectionIndex = sectionIndex;
section.nodes = nodes;
section.triangles = dt.ConnectivityList;
section.gridSize = size(gridData.x);
section.sourceSurfaceGrid = gridData;
section.parameterName = parameterName;
section.parameterCoordinates = [along, depth];
section.pointSources = struct();
section.pointSources.parameterCoordinates = [pointAlong, pointDepth];
section.pointSources.crossCoordinate = pointCross;
section.pointSources.triangleId = pointTriangleId;
section.pointSources.coordinates = point_source_coordinates(pointGrid, parameterName);
section.legacy = struct();
section.legacy.biasL = cfg.mesh.biasL;
section.legacy.biasW = cfg.mesh.biasW;
section.legacy.depthRatio = cfg.mesh.depthRatio;
section.legacy.depthKm = depthKm;
section.legacy.interpolationMethod = string(cfg.mesh.interpolationMethod);
section = add_triangle_metrics(section);
section.neighbors = dtNeighbors;
end

function section = mesh_one_fault_structured(fault, faultIndex, sectionIndex, cfg)
gridData = fault.surfaceGrid;
nodes = [gridData.x(:), gridData.y(:), gridData.z(:)];

[nDepth, nStrike] = size(gridData.x);
triangles = zeros((nDepth - 1) * (nStrike - 1) * 2, 3);
row = 0;
for r = 1:nDepth-1
    for c = 1:nStrike-1
        v11 = sub2ind([nDepth, nStrike], r, c);
        v12 = sub2ind([nDepth, nStrike], r, c + 1);
        v21 = sub2ind([nDepth, nStrike], r + 1, c);
        v22 = sub2ind([nDepth, nStrike], r + 1, c + 1);

        if use_forward_diagonal(cfg, r, c)
            cellTriangles = [v11, v12, v22; v11, v22, v21];
        else
            cellTriangles = [v11, v12, v21; v12, v22, v21];
        end

        triangles(row + (1:2), :) = cellTriangles;
        row = row + 2;
    end
end

valid = all(isfinite(reshape(nodes(triangles, :), size(triangles, 1), [], 3)), [2 3]);
triangles = triangles(valid, :);
[nodes, triangles] = compact_nodes(nodes, triangles);

section = struct();
section.name = string(fault.name);
section.faultIndex = faultIndex;
section.sectionIndex = sectionIndex;
section.nodes = nodes;
section.triangles = triangles;
section.gridSize = [nDepth, nStrike];
section.sourceSurfaceGrid = gridData;
section = add_triangle_metrics(section);
section.neighbors = triangle_neighbors(section.triangles);
end

function [meshGrid, pointGrid, parameterName] = legacy_parameter_grids(gridData)
parameterName = "x";
if isfield(gridData, "fitCrossCoordinate") && string(gridData.fitCrossCoordinate) == "x"
    parameterName = "y";
elseif isfield(gridData, "fitCrossCoordinate") && string(gridData.fitCrossCoordinate) == "y"
    parameterName = "x";
elseif range(gridData.y(:), "omitnan") >= range(gridData.x(:), "omitnan")
    parameterName = "y";
end

if parameterName == "y"
    alongMesh = gridData.y;
    crossMesh = gridData.x;
    alongPoint = gridData.pointY;
    crossPoint = gridData.pointX;
else
    alongMesh = gridData.x;
    crossMesh = gridData.y;
    alongPoint = gridData.pointX;
    crossPoint = gridData.pointY;
end

meshGrid = [row_major(alongMesh), row_major(gridData.z), row_major(crossMesh)];
pointGrid = [row_major(alongPoint), row_major(gridData.pointZ), row_major(crossPoint)];
end

function values = row_major(gridValues)
values = reshape(gridValues.', [], 1);
end

function depthKm = get_mesh_depth_km(cfg, fault)
if isfield(cfg.mesh, "depthKm") && ~isempty(cfg.mesh.depthKm)
    depthKm = cfg.mesh.depthKm;
elseif isfield(fault, "maxDepthKm")
    depthKm = fault.maxDepthKm;
else
    depthKm = abs(min(fault.surfaceGrid.z(:)));
end
end

function [along, cross, depth] = legacy_depth_dependent_points(meshGrid, pointGrid, cfg, depthKm)
meshGrid = meshGrid(all(isfinite(meshGrid), 2), :);
pointGrid = pointGrid(all(isfinite(pointGrid), 2), :);

biasL = cfg.mesh.biasL;
biasW = cfg.mesh.biasW;
ratio = cfg.mesh.depthRatio;

topN = numel(unique(pointGrid(:, 1)));
pointDepths = unique(pointGrid(:, 2));
if numel(pointDepths) < 2
    error("SlipSolve:InvalidMesh", "Legacy meshing requires at least two point-source depth levels.");
end

topWidth = abs(pointDepths(2) - pointDepths(1)) .* ratio;
if abs(biasW - 1) < eps
    layerCount = max(1, ceil(depthKm ./ topWidth));
    layerThickness = [0, repmat(depthKm ./ layerCount, 1, layerCount)];
else
    layerCount = ceil(log(1 + depthKm ./ topWidth .* (biasW - 1)) ./ log(biasW));
    widthFactors = biasW .^ (0:layerCount-1);
    layerThickness = [0, depthKm ./ sum(widthFactors) .* widthFactors];
end

pointLayer = max(2, ceil(topN ./ (biasL .^ (1:layerCount+1))));
minAlong = min(meshGrid(:, 1));
maxAlong = max(meshGrid(:, 1));

along = [];
depth = [];
for i = 1:layerCount+1
    layerAlong = linspace(minAlong, maxAlong, pointLayer(i));
    if i == layerCount + 1
        layerDepth = -depthKm - eps;
    else
        layerDepth = -sum(layerThickness(1:i));
    end

    along = [along; layerAlong(:)]; %#ok<AGROW>
    depth = [depth; layerDepth .* ones(numel(layerAlong), 1)]; %#ok<AGROW>
end

method = char(string(cfg.mesh.interpolationMethod));
cross = griddata([meshGrid(:, 1); pointGrid(:, 1)], [meshGrid(:, 2); pointGrid(:, 2)], ...
    [meshGrid(:, 3); pointGrid(:, 3)], along, depth, method);

invalid = ~isfinite(cross);
if any(invalid)
    crossNearest = griddata([meshGrid(:, 1); pointGrid(:, 1)], [meshGrid(:, 2); pointGrid(:, 2)], ...
        [meshGrid(:, 3); pointGrid(:, 3)], along(invalid), depth(invalid), "nearest");
    cross(invalid) = crossNearest;
end
end

function [along, cross, depth] = legacy_depth_dependent_points_exact(meshGrid, pointGrid, cfg, depthKm)
meshGrid = meshGrid(all(isfinite(meshGrid), 2), :);
pointGrid = pointGrid(all(isfinite(pointGrid), 2), :);

biasL = cfg.mesh.biasL;
biasW = cfg.mesh.biasW;
ratio = cfg.mesh.depthRatio;

topN = numel(unique(pointGrid(:, 1)));
pointDepths = unique(pointGrid(:, 2));
topWidth = (pointDepths(2) - pointDepths(1)) .* ratio;
layerCount = ceil(log(1 + depthKm ./ topWidth .* (biasW - 1)) ./ log(biasW));
widthFactors = biasW .^ (0:layerCount-1);
layerThickness = [0, depthKm ./ sum(widthFactors) .* widthFactors];
pointLayer = ceil(topN ./ (biasL .^ (1:layerCount+1)));

along = [];
depth = [];
minAlong = min(meshGrid(:, 1));
maxAlong = max(meshGrid(:, 1));
for i = 1:layerCount+1
    layerAlong = linspace(minAlong, maxAlong, pointLayer(i));
    if i == layerCount + 1
        layerDepth = -depthKm - eps;
    else
        layerDepth = -sum(layerThickness(1:i));
    end

    along = [along; layerAlong(:)]; %#ok<AGROW>
    depth = [depth; layerDepth .* ones(numel(layerAlong), 1)]; %#ok<AGROW>
end

cross = griddata([meshGrid(:, 1); pointGrid(:, 1)], [meshGrid(:, 2); pointGrid(:, 2)], ...
    [meshGrid(:, 3); pointGrid(:, 3)], along, depth, "cubic");
end

function coordinates = point_source_coordinates(pointGrid, parameterName)
if parameterName == "y"
    coordinates = [pointGrid(:, 3), pointGrid(:, 1), pointGrid(:, 2)];
else
    coordinates = [pointGrid(:, 1), pointGrid(:, 3), pointGrid(:, 2)];
end
end

function tf = use_forward_diagonal(cfg, r, c)
pattern = "alternating";
if isfield(cfg.mesh, "diagonalPattern")
    pattern = string(cfg.mesh.diagonalPattern);
end

switch pattern
    case "forward"
        tf = true;
    case "backward"
        tf = false;
    case "alternating"
        tf = mod(r + c, 2) == 0;
    otherwise
        error("SlipSolve:InvalidMeshConfig", "Unknown cfg.mesh.diagonalPattern: %s", pattern);
end
end

function [nodesOut, trianglesOut] = compact_nodes(nodes, triangles)
used = unique(triangles(:));
map = zeros(size(nodes, 1), 1);
map(used) = 1:numel(used);
nodesOut = nodes(used, :);
trianglesOut = map(triangles);
end

function out = add_triangle_metrics(out)
tri = out.triangles;
v1 = out.nodes(tri(:, 1), :);
v2 = out.nodes(tri(:, 2), :);
v3 = out.nodes(tri(:, 3), :);

edge12 = v2 - v1;
edge13 = v3 - v1;
areaVector = cross(edge12, edge13, 2);
areaVectorNorm = vecnorm(areaVector, 2, 2);
normal = areaVector ./ areaVectorNorm;

strikeVector = cross(repmat([0 0 1], size(normal, 1), 1), normal, 2);
nearHorizontal = vecnorm(strikeVector, 2, 2) < eps;
strikeVector(nearHorizontal, :) = repmat([1 0 0], sum(nearHorizontal), 1);
strikeVector = normalize_rows(strikeVector);

dipVector = cross(normal, strikeVector, 2);
dipVector(dipVector(:, 3) > 0, :) = -dipVector(dipVector(:, 3) > 0, :);
dipVector = normalize_rows(dipVector);

out.centroids = (v1 + v2 + v3) ./ 3;
out.areaKm2 = 0.5 .* areaVectorNorm;
out.areaM2 = out.areaKm2 .* 1e6;
out.normalVector = normal;
out.strikeVector = strikeVector;
out.dipVector = dipVector;
out.strikeAzimuth = mod(atan2d(strikeVector(:, 1), strikeVector(:, 2)), 360);
out.dipDegrees = atan2d(abs(dipVector(:, 3)), hypot(dipVector(:, 1), dipVector(:, 2)));
end

function rows = normalize_rows(rows)
rowNorm = vecnorm(rows, 2, 2);
rows = rows ./ rowNorm;
end

function neighbors = triangle_neighbors(triangles)
nTri = size(triangles, 1);
neighbors = zeros(nTri, 3);
edges = [
    triangles(:, [1 2])
    triangles(:, [2 3])
    triangles(:, [3 1])
];
edges = sort(edges, 2);
triId = repmat((1:nTri).', 3, 1);
edgeId = [ones(nTri, 1); 2 .* ones(nTri, 1); 3 .* ones(nTri, 1)];

[uniqueEdges, ~, groupId] = unique(edges, "rows");
for g = 1:size(uniqueEdges, 1)
    members = find(groupId == g);
    if numel(members) == 2
        a = members(1);
        b = members(2);
        neighbors(triId(a), edgeId(a)) = triId(b);
        neighbors(triId(b), edgeId(b)) = triId(a);
    end
end
end

function value = get_field(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end
