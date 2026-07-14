function figurePaths = plot_fault_geometry(cfg, faultGeometry)
%PLOT_FAULT_GEOMETRY Save interactive fault-geometry diagnostic figures.

p = slipsolve.project.paths(cfg);
if exist(p.faultGeometryFigures, "dir") ~= 7
    mkdir(p.faultGeometryFigures);
end

figurePaths = struct();
figurePaths.geometryFig = fullfile(p.faultGeometryFigures, "fault_geometry_preview.fig");
figurePaths.geometryPng = fullfile(p.faultGeometryFigures, "fault_geometry_preview.png");

visible = "on";
if isfield(cfg, "visualization") && isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end

h = figure("Visible", visible, "Color", "w", "Name", "Fault geometry preview");
hold on;
options = cfg.visualization.faultGeometry;

if isfield(faultGeometry, "faults")
    faults = faultGeometry.faults;
else
    faults = faultGeometry;
end

colors = color_set(options.colormap, max(1, numel(faults)));
maxDepthKm = 0;
for k = 1:numel(faults)
    plot_one_fault(faults(k), colors(k, :));
    maxDepthKm = max(maxDepthKm, faults(k).maxDepthKm);
end

axis equal;
grid on;
xlabel("East (km)");
ylabel("North (km)");
zlabel("Depth (km)");
if get_field(options, "showTitle", true)
    title(string(get_field(options, "title", "Fault Geometry Preview")));
end
view(options.view);
set(gca, "ZDir", "normal");
if ~isempty(options.axisRange)
    xlim(options.axisRange(1:2));
    ylim(options.axisRange(3:4));
    zlim(options.axisRange(5:6));
elseif maxDepthKm > 0
    zlim([-1.1 * maxDepthKm, 2]);
end

savefig(h, figurePaths.geometryFig);
exportgraphics(h, figurePaths.geometryPng, "Resolution", 200);

end

function colors = color_set(name, count)
try
    colors = feval(char(string(name)), count);
catch
    error("SlipSolve:InvalidPlotConfig", "Unknown fault-geometry colormap: %s", name);
end
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end

function plot_one_fault(fault, color)
for k = 1:numel(fault.surfaceTraces)
    tr = fault.surfaceTraces{k};
    plot3(tr.x, tr.y, tr.z, "-", "Color", color, "LineWidth", 3);
end

gridData = fault.surfaceGrid;
mesh(gridData.x, gridData.y, gridData.z, ...
    "LineStyle", "-", "LineWidth", 1, "EdgeColor", 0.25 .* color, ...
    "FaceColor", "none");

if isfield(gridData, "controlPoints")
    scatter3(gridData.controlPoints(:, 1), gridData.controlPoints(:, 2), gridData.controlPoints(:, 3), ...
        10, 0.5 .* color, "filled");
end

for k = 1:numel(fault.depthControls)
    dc = fault.depthControls(k);
    plot3(dc.surfaceX, dc.surfaceY, dc.surfaceZ, "o-", ...
        "Color", color, "LineWidth", 1.5, "MarkerSize", 4);
    plot3(dc.bottomX, dc.bottomY, dc.bottomZ, "o-", ...
        "Color", min(color + 0.35, 1), "LineWidth", 1.5, "MarkerSize", 4);
end
end
