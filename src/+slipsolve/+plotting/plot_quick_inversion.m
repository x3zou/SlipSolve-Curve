function figurePaths = plot_quick_inversion(cfg, quickResult, mesh)
%PLOT_QUICK_INVERSION Save quick-inversion slip diagnostic figures.

p = slipsolve.project.paths(cfg);
if exist(p.quickInversionFigures, "dir") ~= 7
    mkdir(p.quickInversionFigures);
end

figurePaths = struct();
figurePaths.slipFig = fullfile(p.quickInversionFigures, "quick_slip_preview.fig");
figurePaths.slipPng = fullfile(p.quickInversionFigures, "quick_slip_preview.png");

visible = "on";
if isfield(cfg, "visualization") && isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end

slip = hypot(quickResult.slipModel(:, 2), quickResult.slipModel(:, 3));
options = cfg.visualization.quickInversion;
h = figure("Visible", visible, "Color", "w", "Name", "Quick TDE slip preview");
trisurf(mesh.triangles, mesh.nodes(:, 1), mesh.nodes(:, 2), mesh.nodes(:, 3), slip, ...
    "FaceColor", "flat", "FaceAlpha", 0.92, "EdgeColor", [0.18 0.18 0.18], "LineWidth", 0.25);
axis equal;
grid on;
xlabel("East (km)");
ylabel("North (km)");
zlabel("Depth (km)");
if get_field(options, "showTitle", true)
    titleText = string(get_field(options, "title", ""));
    if strlength(titleText) == 0
        titleText = sprintf("Quick TDE Slip: VR %.1f%%", quickResult.varianceReduction);
    end
    title(titleText);
end
view(options.view);
set(gca, "ZDir", "normal");
apply_axis_range(gca, options.axisRange);
if ~isempty(options.colorbarRange), clim(options.colorbarRange); end
colormap(h, char(options.colormap));
colorbar;

savefig(h, figurePaths.slipFig);
exportgraphics(h, figurePaths.slipPng, "Resolution", 220);

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
