function figurePaths = plot_triangular_mesh(cfg, mesh, faultGeometry)
%PLOT_TRIANGULAR_MESH Save interactive triangular-mesh diagnostic figures.

p = slipsolve.project.paths(cfg);
if exist(p.triangularMeshFigures, "dir") ~= 7
    mkdir(p.triangularMeshFigures);
end

figurePaths = struct();
figurePaths.meshFig = fullfile(p.triangularMeshFigures, "triangular_mesh_preview.fig");
figurePaths.meshPng = fullfile(p.triangularMeshFigures, "triangular_mesh_preview.png");

visible = "on";
if isfield(cfg, "visualization") && isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end

h = figure("Visible", visible, "Color", "w", "Name", "Triangular mesh preview");
hold on;
options = cfg.visualization.triangularMesh;

colorValues = color_values(cfg, mesh);
trisurf(mesh.triangles, mesh.nodes(:, 1), mesh.nodes(:, 2), mesh.nodes(:, 3), colorValues, ...
    "FaceColor", "flat", "FaceAlpha", 0.88, "EdgeColor", [0.18 0.18 0.18], "LineWidth", 0.25);

if should_show_fault_trace(cfg)
    plot_surface_traces(mesh, faultGeometry);
end

axis equal;
grid on;
xlabel("East (km)");
ylabel("North (km)");
zlabel("Depth (km)");
if get_field(options, "showTitle", true)
    titleText = string(get_field(options, "title", ""));
    if strlength(titleText) == 0
        titleText = sprintf("Triangular Mesh: %d Nodes, %d Triangles", ...
            size(mesh.nodes, 1), size(mesh.triangles, 1));
    end
    title(titleText);
end
view(options.view);
set(gca, "ZDir", "normal");
apply_axis_range(gca, options.axisRange);
if ~isempty(options.colorbarRange), clim(options.colorbarRange); end
colormap(h, char(options.colormap));
colorbar;

savefig(h, figurePaths.meshFig);
exportgraphics(h, figurePaths.meshPng, "Resolution", 220);

end

function apply_axis_range(ax, range)
if isempty(range), return; end
xlim(ax, range(1:2));
ylim(ax, range(3:4));
zlim(ax, range(5:6));
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end

function values = color_values(cfg, mesh)
mode = "depth";
if isfield(cfg, "visualization") && isfield(cfg.visualization, "triangularMesh") ...
        && isfield(cfg.visualization.triangularMesh, "colorBy")
    mode = string(cfg.visualization.triangularMesh.colorBy);
end

switch mode
    case "section"
        values = mesh.sectionId;
    case "area"
        values = mesh.areaKm2;
    otherwise
        values = -mesh.centroids(:, 3);
end
end

function tf = should_show_fault_trace(cfg)
tf = true;
if isfield(cfg, "visualization") && isfield(cfg.visualization, "triangularMesh") ...
        && isfield(cfg.visualization.triangularMesh, "showFaultTrace")
    tf = cfg.visualization.triangularMesh.showFaultTrace;
end
end

function plot_surface_traces(mesh, faultGeometry)
if isfield(faultGeometry, "faults")
    faults = faultGeometry.faults;
else
    faults = faultGeometry;
end

colors = lines(max(1, numel(faults)));
for i = 1:numel(mesh.metadata.activeFaults)
    faultIndex = mesh.metadata.activeFaults(i);
    fault = faults(faultIndex);
    for k = 1:numel(fault.surfaceTraces)
        tr = fault.surfaceTraces{k};
        plot3(tr.x, tr.y, tr.z + 0.05, "-", "Color", colors(faultIndex, :), "LineWidth", 2.5);
    end
end
end
