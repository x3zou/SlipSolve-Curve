function buildInfo = build_full_xyz(cfg, mesh, dataInsar, greenName, matrixFiles)
%BUILD_FULL_XYZ Legacy node-based full-resolution G_e, G_n, and G_u.
%
% The numerical sequence follows Comb_Green_Layered_Nodes_xyz.m: construct
% the nodal model, evaluate deep nodes with comb_green_nodes_xyz, replace the
% top two levels with alter_top_green_nodes_xyz, and assemble SS before DS.

sections = mesh.sections(:);
if double(cfg.compositeGreens.pointSourcesPerTriangle) ~= 15
    error("SlipSolve:LegacyCompatibility", ...
        "The legacy comb_green_nodes_xyz and alter_top_green_nodes_xyz routines use exactly 15 point sources per triangle.");
end
contexts = cell(numel(sections), 1);
parameterCounts = zeros(numel(sections), 1);
for k = 1:numel(sections)
    contexts{k} = prepare_section(sections(k), cfg, k);
    parameterCounts(k) = contexts{k}.parameterCount;
end

nobs = size(dataInsar, 1);
nparam = sum(parameterCounts);
fprintf("Allocating full composite Green's matrices: %d observations x %d columns each.\n", ...
    nobs, 2 * nparam);
G_e = zeros(nobs, 2 * nparam);
G_n = zeros(nobs, 2 * nparam);
G_u = zeros(nobs, 2 * nparam);

sectionInfo = repmat(struct( ...
    "name", "", "faultIndex", 0, "originalNodeCount", 0, ...
    "topParameterCount", 0, "deepParameterCount", 0, "parameterCount", 0, ...
    "strikeColumns", [0 0], "dipColumns", [0 0], ...
    "parameterCoordinatesKm", zeros(0, 3), "parameterKind", strings(0, 1)), ...
    numel(sections), 1);

parameterOffset = 0;
for k = 1:numel(sections)
    context = contexts{k};
    fprintf("Building full XYZ composite Green's functions for fault %d/%d: %s\n", ...
        k, numel(sections), char(context.name));

    localCount = context.parameterCount;
    ssColumns = parameterOffset + (1:localCount);
    dsColumns = nparam + parameterOffset + (1:localCount);

    [topEast, topNorth, topUp] = alter_top_green_nodes_xyz( ...
        context.slipModel, dataInsar, char(greenName), context.strikeInterpolant, context.dipInterpolant);
    topCount = size(topEast, 2) / 2;
    if topCount ~= context.topParameterCount
        error("SlipSolve:LegacyGeometryMismatch", ...
            "Top Green's column count changed from %d to %d for fault %s.", ...
            context.topParameterCount, topCount, context.name);
    end
    topLocal = 1:topCount;

    G_e(:, ssColumns(topLocal)) = topEast(:, 1:topCount);
    G_e(:, dsColumns(topLocal)) = topEast(:, topCount+1:end);
    G_n(:, ssColumns(topLocal)) = topNorth(:, 1:topCount);
    G_n(:, dsColumns(topLocal)) = topNorth(:, topCount+1:end);
    G_u(:, ssColumns(topLocal)) = topUp(:, 1:topCount);
    G_u(:, dsColumns(topLocal)) = topUp(:, topCount+1:end);
    clear topEast topNorth topUp

    ndeep = numel(context.deepNodeIndices);
    for j = 1:ndeep
        nodeIndex = context.deepNodeIndices(j);
        if mod(j, 100) == 0 || j == 1 || j == ndeep
            fprintf("Fault %d/%d deep-node progress: [%d/%d], original node %d\n", ...
                k, numel(sections), j, ndeep, nodeIndex);
        end

        [Ge, Gn, Gu] = comb_green_nodes_xyz( ...
            nodeIndex, context.connectivity, context.pointsKm, dataInsar, ...
            context.strikeInterpolant, context.dipInterpolant, char(greenName));
        npatch = size(Ge, 2) / 2;
        localColumn = topCount + j;
        G_e(:, ssColumns(localColumn)) = sum(Ge(:, 1:npatch), 2);
        G_e(:, dsColumns(localColumn)) = sum(Ge(:, npatch+1:end), 2);
        G_n(:, ssColumns(localColumn)) = sum(Gn(:, 1:npatch), 2);
        G_n(:, dsColumns(localColumn)) = sum(Gn(:, npatch+1:end), 2);
        G_u(:, ssColumns(localColumn)) = sum(Gu(:, 1:npatch), 2);
        G_u(:, dsColumns(localColumn)) = sum(Gu(:, npatch+1:end), 2);
        clear Ge Gn Gu
    end

    sectionInfo(k).name = context.name;
    sectionInfo(k).faultIndex = context.faultIndex;
    sectionInfo(k).originalNodeCount = size(context.pointsKm, 1);
    sectionInfo(k).topParameterCount = topCount;
    sectionInfo(k).deepParameterCount = numel(context.deepNodeIndices);
    sectionInfo(k).parameterCount = localCount;
    sectionInfo(k).strikeColumns = [ssColumns(1), ssColumns(end)];
    sectionInfo(k).dipColumns = [dsColumns(1), dsColumns(end)];
    sectionInfo(k).parameterCoordinatesKm = context.parameterCoordinatesKm;
    sectionInfo(k).parameterKind = context.parameterKind;
    parameterOffset = parameterOffset + localCount;
end

fprintf("Saving east composite Green's function: %s\n", matrixFiles.east);
save(matrixFiles.east, "G_e", "-v7.3");
fprintf("Saving north composite Green's function: %s\n", matrixFiles.north);
save(matrixFiles.north, "G_n", "-v7.3");
fprintf("Saving vertical composite Green's function: %s\n", matrixFiles.vertical);
save(matrixFiles.vertical, "G_u", "-v7.3");

buildInfo = struct();
buildInfo.observationCount = nobs;
buildInfo.parameterCountPerComponent = nparam;
buildInfo.matrixSize = [nobs, 2 * nparam];
buildInfo.columnOrder = "all_faults_strike_slip_then_all_faults_dip_slip";
buildInfo.sectionInfo = sectionInfo;
buildInfo.skipTopDepthLevels = cfg.compositeGreens.skipTopDepthLevels;
buildInfo.pointSourcesPerTriangle = cfg.compositeGreens.pointSourcesPerTriangle;
end

function context = prepare_section(section, cfg, defaultFaultIndex)
grid = section.sourceSurfaceGrid;
required = ["x", "y", "z", "pointX", "pointY", "pointZ"];
for k = 1:numel(required)
    if ~isfield(grid, required(k))
        error("SlipSolve:MissingInput", ...
            "Mesh section %s lacks sourceSurfaceGrid.%s required by the legacy composite method.", ...
            section.name, required(k));
    end
end

[strike, dip] = compute_geometry_pointSource(grid.x, grid.y, grid.z);
strikeInterpolant = scatteredInterpolant( ...
    grid.pointX(:), grid.pointY(:), grid.pointZ(:), strike(:));
dipInterpolant = scatteredInterpolant( ...
    grid.pointX(:), grid.pointY(:), grid.pointZ(:), dip(:));

pointsKm = double(section.nodes);
connectivity = double(section.triangles);
pointsM = pointsKm .* 1e3;
nodeStrike = strikeInterpolant(pointsKm(:, 1), pointsKm(:, 2), pointsKm(:, 3));
nodeDip = dipInterpolant(pointsKm(:, 1), pointsKm(:, 2), pointsKm(:, 3));
faultIndex = get_field(section, "faultIndex", defaultFaultIndex);

slipModel = make_node_geometry(faultIndex, pointsM, connectivity);
slipModel = [slipModel, nodeStrike, nodeDip];
slipModel = alter_nodes_geometry_tri(slipModel);
if size(slipModel, 2) < 27
    difference = 27 - size(slipModel, 2);
    supplement = nan(size(slipModel, 1), difference);
    slipModel = [slipModel(:, 1:8), slipModel(:, 9:16-difference), ...
        supplement, slipModel(:, 16-difference+1:end)];
end

levels = sort(unique(pointsKm(:, 3)), "descend");
[~, levelIndex] = ismember(pointsKm(:, 3), levels);
skipLevels = double(cfg.compositeGreens.skipTopDepthLevels);
deepNodeIndices = find(levelIndex > skipLevels);
topAssignment = assignTrianglesToNodes(slipModel);
topParameterCount = size(topAssignment, 1);

topCoordinatesKm = slipModel(topAssignment(:, 1), 6:8) ./ 1e3;
deepCoordinatesKm = pointsKm(deepNodeIndices, :);

context = struct();
context.name = string(section.name);
context.faultIndex = faultIndex;
context.pointsKm = pointsKm;
context.connectivity = connectivity;
context.strikeInterpolant = strikeInterpolant;
context.dipInterpolant = dipInterpolant;
context.slipModel = slipModel;
context.levelIndex = levelIndex;
context.deepNodeIndices = deepNodeIndices;
context.topParameterCount = topParameterCount;
context.parameterCount = topParameterCount + numel(deepNodeIndices);
context.parameterCoordinatesKm = [topCoordinatesKm; deepCoordinatesKm];
context.parameterKind = [repmat("top_hybrid", topParameterCount, 1); ...
    repmat("layered_node", numel(deepNodeIndices), 1)];
end

function value = get_field(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
