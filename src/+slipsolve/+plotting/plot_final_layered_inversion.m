function figurePaths = plot_final_layered_inversion(cfg, result, fullGridResponse)
%PLOT_FINAL_LAYERED_INVERSION Legacy nodal slip and per-track fit figures.

if nargin < 3
    fullGridResponse = struct();
end

p = slipsolve.project.paths(cfg);
if exist(p.finalInversionFigures, "dir") ~= 7
    mkdir(p.finalInversionFigures);
end
options = cfg.visualization.finalInversion;
figurePaths = struct();
figurePaths.slipFig = fullfile(p.finalInversionFigures, "layered_nodal_slip.fig");
figurePaths.slipPng = fullfile(p.finalInversionFigures, "layered_nodal_slip.png");
figurePaths.triangleFig = fullfile(p.finalInversionFigures, "simple_triangular_slip.fig");
figurePaths.trianglePng = fullfile(p.finalInversionFigures, "simple_triangular_slip.png");
figurePaths.trackFigures = repmat(struct("name", "", "fig", "", "png", ""), ...
    numel(result.tracks), 1);
figurePaths.fullResolutionFigures = struct("name", {}, "fig", {}, "png", {}, ...
    "residualGrid", {});

if ~get_field(options, "enabled", true)
    return
end
visible = "on";
if isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end

hSlip = legacy_nodal_plot(result.completeNodalSlipModel, options, visible, cfg.insar.origin);
save_figure(cfg, hSlip, figurePaths.slipFig, figurePaths.slipPng);
hTriangle = triangle_center_plot(result.triangleModel, options, visible, cfg.insar.origin);
save_figure(cfg, hTriangle, figurePaths.triangleFig, figurePaths.trianglePng);

fitOptions = options.fitPlots;
if get_field(fitOptions, "enabled", true)
    for k = 1:numel(result.tracks)
        track = result.tracks(k);
        name = safe_name(track.name);
        figurePaths.trackFigures(k).name = string(track.name);
        figurePaths.trackFigures(k).fig = fullfile(p.finalInversionFigures, name+"_fit.fig");
        figurePaths.trackFigures(k).png = fullfile(p.finalInversionFigures, name+"_fit.png");
        hTrack = track_fit_plot(track, fitOptions, visible);
        save_figure(cfg, hTrack, figurePaths.trackFigures(k).fig, ...
            figurePaths.trackFigures(k).png);
    end
end
if isfield(options, "fullResolutionFitPlots") && ...
        get_field(options.fullResolutionFitPlots, "enabled", false)
    figurePaths.fullResolutionFigures = plot_full_resolution_fits( ...
        cfg, result, fullGridResponse, options.fullResolutionFitPlots, visible);
end
end

function paths = plot_full_resolution_fits(cfg, result, response, options, visible)
if ~isfield(response, "east")
    error("SlipSolve:MissingFullGridPrediction", ...
        "Full-resolution fit plots require evaluated G_e/G_n/G_u fields.");
end
p = slipsolve.project.paths(cfg);
definitions = options.tracks(:);
paths = repmat(struct("name", "", "fig", "", "png", "", ...
    "residualGrid", ""), numel(definitions), 1);
[faultOld, faultPlot] = read_fault_traces(options.faultTraceFiles, ...
    cfg.insar.origin, options.plotOrigin);
if isempty(faultOld)
    error("SlipSolve:MissingInput", ...
        "Full-resolution legacy plots require at least one fault trace.");
end
barrier = faultOld{1};

for k = 1:numel(definitions)
    definition = definitions(k);
    track = find_result_track(result.tracks, definition.name);
    fprintf("Building full-resolution fit plot %d/%d: %s\n", ...
        k, numel(definitions), definition.name);
    validate_grid_files(definition);
    [modelLon, modelLat, ~] = grdread2(char(definition.modelGridFile));
    [dataLon, dataLat, dataGrid] = grdread2(char(definition.dataGridFile));
    [modelX, modelY] = geographic_grid_to_local(modelLon, modelLat, cfg.insar.origin);
    [dataX, dataY] = geographic_grid_to_local(dataLon, dataLat, cfg.insar.origin);

    projectedComposite = project_full_response(definition, track, response, ...
        modelX, modelY, cfg.insar.origin);
    modelOnPredictionGrid = interpXYZ_barrier( ...
        [response.x, response.y, projectedComposite], ...
        modelX(:), modelY(:), barrier, "Method", "linear", "Extrap", "none");

    sameGrid = isequal(size(modelX), size(dataX)) && ...
        max(abs(modelX(:)-dataX(:)), [], "omitnan") < 1e-6 && ...
        max(abs(modelY(:)-dataY(:)), [], "omitnan") < 1e-6;
    if sameGrid
        modelOnDataGrid = reshape(modelOnPredictionGrid, size(dataGrid));
    else
        threshold = get_field(options, "distanceThresholdKm", []);
        if isempty(threshold), thresholdM = []; else, thresholdM = threshold*1e3; end
        modelOnDataGrid = legacy_near_fault_interpolation( ...
            modelX(:), modelY(:), modelOnPredictionGrid(:), ...
            dataX, dataY, barrier, thresholdM);
    end
    residual = nan(size(dataGrid));
    valid = isfinite(dataGrid);
    residual(valid) = dataGrid(valid)-modelOnDataGrid(valid);
    if lower(string(definition.type)) == "los"
        residualPlot = residual;
    else
        residualPlot = fill_residual_by_fault_side(residual, dataX, dataY, barrier);
    end

    [dataPlotX, dataPlotY] = switch_origin(dataX(:), dataY(:), ...
        cfg.insar.origin.lon, cfg.insar.origin.lat, ...
        options.plotOrigin.lon, options.plotOrigin.lat, options.plotOrigin.refLon);
    dataPlotX = reshape(dataPlotX, size(dataGrid))/1e3;
    dataPlotY = reshape(dataPlotY, size(dataGrid))/1e3;
    h = legacy_full_fit_figure(track, residualPlot, dataPlotX, dataPlotY, ...
        faultPlot, definition, options, visible, cfg.insar.origin, options.plotOrigin);

    name = safe_name(definition.name);
    paths(k).name = string(definition.name);
    paths(k).fig = fullfile(p.finalFullResolutionFigures, name+"_full_fit.fig");
    paths(k).png = fullfile(p.finalFullResolutionFigures, name+"_full_fit.png");
    save_figure(cfg, h, paths(k).fig, paths(k).png);
    if get_field(options, "saveResidualGrids", false)
        residualFile = fullfile(string(options.outputDirectory), name+"_residual.grd");
        grdwrite2(dataLon, dataLat, residual, char(residualFile));
        paths(k).residualGrid = residualFile;
    end
end
end

function projected = project_full_response(definition, track, response, modelX, modelY, origin)
type = lower(string(definition.type));
if type == "azo"
    [sinF, cosF] = heading_from_sample(track.sampleFile);
    projected = sinF*response.east + cosF*response.north;
    return
end
[lookLon, lookLat, lookE] = grdread2(char(definition.lookEFile));
[~, ~, lookN] = grdread2(char(definition.lookNFile));
[~, ~, lookU] = grdread2(char(definition.lookUFile));
[lookX, lookY] = geographic_grid_to_local(lookLon, lookLat, origin);
lookEast = interpXYZ([lookX(:), lookY(:), lookE(:)], ...
    response.x, response.y, "natural", "none");
lookNorth = interpXYZ([lookX(:), lookY(:), lookN(:)], ...
    response.x, response.y, "natural", "none");
lookUp = interpXYZ([lookX(:), lookY(:), lookU(:)], ...
    response.x, response.y, "natural", "none");
projected = response.east.*lookEast + response.north.*lookNorth + response.up.*lookUp;
if all(~isfinite(projected))
    error("SlipSolve:FullGridProjectionFailed", ...
        "LOS look-vector interpolation produced no finite values for %s.", definition.name);
end
% modelX/modelY are intentionally accepted to mirror the legacy projection
% contract and ensure their grid is read before look-vector interpolation.
if isempty(modelX) || isempty(modelY), projected(:) = NaN; end
end

function modelGrid = legacy_near_fault_interpolation( ...
        sourceX, sourceY, sourceModel, queryX, queryY, fault, threshold)
good = isfinite(sourceX) & isfinite(sourceY) & isfinite(sourceModel);
sourceX = sourceX(good);
sourceY = sourceY(good);
sourceModel = sourceModel(good);
[sideSource, ~] = nearest_segment_side(sourceX, sourceY, fault);
[sideQuery, distance] = nearest_segment_side(queryX(:), queryY(:), fault);
globalF = scatteredInterpolant(sourceX, sourceY, sourceModel, "natural", "nearest");
values = globalF(queryX(:), queryY(:));
for side = [false true]
    mask = sideSource == side;
    if nnz(mask) < 3, continue; end
    sideF = scatteredInterpolant(sourceX(mask), sourceY(mask), sourceModel(mask), ...
        "natural", "none");
    if isempty(threshold)
        ux = unique(queryX(:));
        uy = unique(queryY(:));
        thresholdUse = 2*max(median(diff(ux)), median(diff(uy)));
    else
        thresholdUse = threshold;
    end
    queryMask = distance <= thresholdUse & sideQuery == side;
    sideValues = sideF(queryX(queryMask), queryY(queryMask));
    replace = isfinite(sideValues);
    indices = find(queryMask);
    values(indices(replace)) = sideValues(replace);
end
modelGrid = reshape(values, size(queryX));
end

function filled = fill_residual_by_fault_side(residual, x, y, fault)
[side, ~] = nearest_segment_side(x(:), y(:), fault);
raw = residual(:);
filled = nan(size(raw));
for sideValue = [false true]
    source = isfinite(raw) & side == sideValue;
    query = side == sideValue;
    if nnz(source) < 3, continue; end
    interpolant = scatteredInterpolant(x(source), y(source), raw(source), ...
        "natural", "none");
    filled(query) = interpolant(x(query), y(query));
end
filled = reshape(filled, size(residual));
end

function [side, distance] = nearest_segment_side(x, y, fault)
x1 = fault(1:end-1, 1);
y1 = fault(1:end-1, 2);
x2 = fault(2:end, 1);
y2 = fault(2:end, 2);
midpoints = [(x1+x2)/2, (y1+y2)/2];
index = knnsearch(midpoints, [x(:), y(:)]);
dx = x2(index)-x1(index);
dy = y2(index)-y1(index);
crossProduct = dx.*(y(:)-y1(index))-dy.*(x(:)-x1(index));
side = crossProduct >= 0;
distance = abs(crossProduct)./hypot(dx, dy);
end

function h = legacy_full_fit_figure(track, residual, gridX, gridY, ...
        faultPlot, definition, options, visible, oldOrigin, plotOrigin)
[sampleX, sampleY] = switch_origin(track.x, track.y, ...
    oldOrigin.lon, oldOrigin.lat, plotOrigin.lon, plotOrigin.lat, plotOrigin.refLon);
sampleX = sampleX/1e3;
sampleY = sampleY/1e3;
h = figure("Visible", visible, "Color", "w", ...
    "Name", string(definition.name)+" full-resolution layered fit");
layout = tiledlayout(h, 1, 3, "TileSpacing", "compact", "Padding", "compact");
limits = definition.colorbarRange;

ax1 = nexttile;
scatter(ax1, sampleX, sampleY, 30, track.data, "filled");
format_full_map_axis(ax1, limits, definition.axisRange);
plot_faults(ax1, faultPlot);
panel_label(ax1, "a");

ax2 = nexttile;
scatter(ax2, sampleX, sampleY, 30, track.prediction, "filled");
format_full_map_axis(ax2, limits, definition.axisRange);
plot_faults(ax2, faultPlot);
panel_label(ax2, "b");

ax3 = nexttile;
pcolor(ax3, gridX, gridY, residual);
shading(ax3, "flat");
format_full_map_axis(ax3, limits, definition.axisRange);
plot_faults(ax3, faultPlot);
panel_label(ax3, "c");
cb = colorbar(ax3, "eastoutside");
cb.Label.String = "cm";
cb.Label.Rotation = 0;
cb.FontSize = 16;
cb.Label.FontSize = 20;

layout.XLabel.String = "Easting (km)";
layout.YLabel.String = "Northing (km)";
layout.XLabel.FontSize = 20;
layout.YLabel.FontSize = 20;
if get_field(options, "showTitles", true) && strlength(string(definition.title)) > 0
    title(layout, definition.title, "FontSize", 25, "FontWeight", "bold");
end
colormap(h, char(options.colormap));
end

function format_full_map_axis(ax, limits, axisRange)
axis(ax, "equal");
grid(ax, "on");
box(ax, "on");
clim(ax, limits);
if ~isempty(axisRange)
    xlim(ax, axisRange(1:2));
    ylim(ax, axisRange(3:4));
end
set(ax, "FontSize", 20);
hold(ax, "on");
end

function plot_faults(ax, faults)
for k = 1:numel(faults)
    plot(ax, faults{k}(:, 1)/1e3, faults{k}(:, 2)/1e3, ...
        "LineWidth", 2, "Color", "white");
    plot_long_dash_equal(ax, faults{k}(:, 1)/1e3, faults{k}(:, 2)/1e3, ...
        10, 3, "LineWidth", 2, "Color", "black");
end
end

function panel_label(ax, label)
text(ax, 0.1, 1.01, label, "Units", "normalized", ...
    "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
    "FontWeight", "bold", "FontSize", 30);
end

function [faultOld, faultPlot] = read_fault_traces(files, oldOrigin, plotOrigin)
files = string(files(:));
faultOld = cell(numel(files), 1);
faultPlot = cell(numel(files), 1);
for k = 1:numel(files)
    values = readmatrix(files(k), "FileType", "text");
    values = values(all(isfinite(values(:, 1:2)), 2), 1:2);
    [x, y] = ll2xy(values(:, 1), values(:, 2), oldOrigin.refLon);
    [originX, originY] = ll2xy(oldOrigin.lon, oldOrigin.lat, oldOrigin.refLon);
    faultOld{k} = [x-originX, y-originY];
    [plotX, plotY] = ll2xy(values(:, 1), values(:, 2), plotOrigin.refLon);
    [plotOriginX, plotOriginY] = ll2xy(plotOrigin.lon, plotOrigin.lat, plotOrigin.refLon);
    faultPlot{k} = [plotX-plotOriginX, plotY-plotOriginY];
end
end

function [X, Y] = geographic_grid_to_local(longitude, latitude, origin)
[longitudeGrid, latitudeGrid] = meshgrid(longitude, latitude);
[X, Y] = ll2xy(longitudeGrid(:), latitudeGrid(:), origin.refLon);
[originX, originY] = ll2xy(origin.lon, origin.lat, origin.refLon);
X = reshape(X-originX, size(longitudeGrid));
Y = reshape(Y-originY, size(latitudeGrid));
end

function [sinF, cosF] = heading_from_sample(sampleFile)
folder = fileparts(sampleFile);
sinF = readmatrix(fullfile(folder, "sinF.dat"));
cosF = readmatrix(fullfile(folder, "cosF.dat"));
sinF = sinF(1);
cosF = cosF(1);
end

function track = find_result_track(tracks, name)
index = find(string({tracks.name}) == string(name), 1);
if isempty(index)
    error("SlipSolve:MissingInput", ...
        "No final inversion track matches full-grid plot definition %s.", name);
end
track = tracks(index);
end

function validate_grid_files(definition)
files = [string(definition.modelGridFile), string(definition.dataGridFile)];
if lower(string(definition.type)) ~= "azo"
    files = [files, string(definition.lookEFile), ...
        string(definition.lookNFile), string(definition.lookUFile)];
end
for file = files
    if exist(file, "file") ~= 2
        error("SlipSolve:MissingInput", "Full-resolution plot input is missing: %s", file);
    end
end
end

function plot_long_dash_equal(ax, x, y, lengthOn, lengthOff, varargin)
period = lengthOn+lengthOff;
distance = [0; cumsum(hypot(diff(x), diff(y)))];
outputX = [];
outputY = [];
for start = 0:period:distance(end)
    stop = min(start+lengthOn, distance(end));
    if stop <= start, continue; end
    cuts = unique([start; distance(distance > start & distance < stop); stop]);
    outputX = [outputX; interp1(distance, x, cuts); NaN]; %#ok<AGROW>
    outputY = [outputY; interp1(distance, y, cuts); NaN]; %#ok<AGROW>
end
plot(ax, outputX, outputY, "-", varargin{:});
end

function h = legacy_nodal_plot(slipModel, options, visible, sourceOrigin)
data = slipModel;
plotOrigin = get_field(options, "plotOrigin", sourceOrigin);
[x, y] = switch_origin(data(:, 6), data(:, 7), ...
    sourceOrigin.lon, sourceOrigin.lat, ...
    plotOrigin.lon, plotOrigin.lat, plotOrigin.refLon);
vertices = [x, y, data(:, 8)]/1e3;
faces = adjacency_triangles(data(:, 9:15), size(data, 1));
strikeSlip = -data(:, 4)/100;
dipSlip = data(:, 5)/100;

h = figure("Visible", visible, "Color", "w", "Name", "Final layered nodal slip");
ax1 = subplot(2, 1, 1, "Parent", h);
patch(ax1, "Vertices", vertices, "Faces", faces, "FaceColor", "none", ...
    "EdgeColor", "interp", "LineWidth", options.lineWidth, ...
    "FaceVertexCData", strikeSlip);
format_slip_axis(ax1, options);
cb1 = colorbar(ax1, "northoutside", "Orientation", "horizontal");
cb1.Position = get_field(options, "strikeColorbarPosition", ...
    [0.76, 0.60, 0.17, 0.018]);
cb1.Label.String = "Strike-slip (m)";
cb1.Label.FontSize = 16;
cb1.FontSize = 16;
clim(ax1, options.strikeColorbarRange);
if options.showTitles, title(ax1, options.strikeTitle); end

ax2 = subplot(2, 1, 2, "Parent", h);
patch(ax2, "Vertices", vertices, "Faces", faces, "FaceColor", "none", ...
    "EdgeColor", "interp", "LineWidth", options.lineWidth, ...
    "FaceVertexCData", dipSlip);
format_slip_axis(ax2, options);
cb2 = colorbar(ax2, "northoutside", "Orientation", "horizontal");
cb2.Position = get_field(options, "dipColorbarPosition", ...
    [0.76, 0.12, 0.17, 0.018]);
cb2.Label.String = "Dip-slip (m)";
cb2.Label.FontSize = 16;
cb2.FontSize = 16;
clim(ax2, options.dipColorbarRange);
if options.showTitles, title(ax2, options.dipTitle); end
colormap(h, char(options.colormap));
end

function format_slip_axis(ax, options)
axis(ax, "equal");
grid(ax, "on");
view(ax, options.view);
xlabel(ax, "Easting (km)", "FontSize", 18);
ylabel(ax, "Northing (km)", "FontSize", 18);
zlabel(ax, "Depth (km)", "FontSize", 18);
set(ax, "FontSize", 18, "LineWidth", 1);
range = options.axisRange;
if ~isempty(range)
    xlim(ax, range(1:2));
    ylim(ax, range(3:4));
    zlim(ax, range(5:6));
end
end

function faces = adjacency_triangles(connectivity, n)
adjacency = false(n);
for i = 1:n
    connected = connectivity(i, :);
    connected = connected(isfinite(connected) & connected >= 1 & connected <= n);
    adjacency(i, connected) = true;
end
adjacency = adjacency | adjacency.';
faces = zeros(0, 3);
for i = 1:n-2
    neighborsI = find(adjacency(i, i+1:end))+i;
    for j = neighborsI
        neighborsJ = find(adjacency(j, j+1:end))+j;
        common = intersect(neighborsI, neighborsJ);
        for k = common
            faces(end+1, :) = [i j k]; %#ok<AGROW>
        end
    end
end
end

function h = triangle_center_plot(triangles, options, visible, sourceOrigin)
h = figure("Visible", visible, "Color", "w", "Name", "Simple triangular slip");
plotOrigin = get_field(options, "plotOrigin", sourceOrigin);
[x, y] = switch_origin(triangles(:, 1), triangles(:, 2), ...
    sourceOrigin.lon, sourceOrigin.lat, ...
    plotOrigin.lon, plotOrigin.lat, plotOrigin.refLon);
scatter3(x/1e3, y/1e3, triangles(:, 3)/1e3, ...
    36, triangles(:, 4)/100, "filled");
axis equal;
grid on;
xlabel("Easting (km)");
ylabel("Northing (km)");
zlabel("Depth (km)");
if options.showTitles, title(string(options.triangleTitle)); end
view(options.view);
if ~isempty(options.axisRange)
    xlim(options.axisRange(1:2));
    ylim(options.axisRange(3:4));
    zlim(options.axisRange(5:6));
end
colormap(h, char(options.colormap));
colorbar;
end

function h = track_fit_plot(track, options, visible)
h = figure("Visible", visible, "Color", "w", "Name", string(track.name)+" layered fit");
tiledlayout(h, 1, 3, "TileSpacing", "compact", "Padding", "compact");
dataModelLimits = centered_or_configured(options.dataModelColorbarRange, ...
    [track.data; track.prediction]);
residualLimits = centered_or_configured(options.residualColorbarRange, track.residual);
    labels = string(get_field(options, "panelTitles", ["Data" "Model" "Residual"]));
values = {track.data, track.prediction, track.residual};
for k = 1:3
    ax = nexttile;
    scatter(ax, track.x/1e3, track.y/1e3, 18, values{k}, "filled");
    axis(ax, "equal");
    axis(ax, "tight");
    grid(ax, "on");
    xlabel(ax, "Easting (km)");
    ylabel(ax, "Northing (km)");
    if options.showTitles, title(ax, labels(k)); end
    if k < 3, clim(ax, dataModelLimits); else, clim(ax, residualLimits); end
    if ~isempty(options.axisRange)
        xlim(ax, options.axisRange(1:2));
        ylim(ax, options.axisRange(3:4));
    end
    colorbar(ax);
end
colormap(h, char(options.colormap));
end

function limits = centered_or_configured(configured, values)
if ~isempty(configured)
    limits = configured;
    return
end
values = values(isfinite(values));
if isempty(values)
    maximum = 1;
else
    maximum = max(abs(values));
    if maximum == 0, maximum = 1; end
end
limits = [-maximum maximum];
end

function save_figure(cfg, h, figPath, pngPath)
if get_field(cfg.visualization, "saveFig", true)
    savefig(h, figPath);
end
if get_field(cfg.visualization, "savePng", true)
    exportgraphics(h, pngPath, "Resolution", 220);
end
drawnow;
end

function value = get_field(s, name, defaultValue)
if isfield(s, name), value = s.(name); else, value = defaultValue; end
end

function name = safe_name(value)
name = regexprep(string(value), "[^A-Za-z0-9_.-]", "_");
end
