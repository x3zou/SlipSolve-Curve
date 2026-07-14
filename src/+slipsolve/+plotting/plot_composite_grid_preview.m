function figureFiles = plot_composite_grid_preview(cfg, X, Y, meshTrace, xRange, yRange, D1, D2, M)
%PLOT_COMPOSITE_GRID_PREVIEW Show the adaptive full-resolution XYZ grid.

figureFiles = strings(0, 1);
if ~isfield(cfg, "visualization") || ~cfg.visualization.enabled
    return
end
options = plot_options(cfg);
if ~options.enabled
    return
end

fig = figure("Name", "Full-resolution composite Green's grid", "Color", "w", ...
    "Visible", figure_visibility(cfg));
ax = axes(fig);
hold(ax, "on");

stride = max(1, round(double(options.pointStride)));
gridX = X(1:stride:end) ./ 1e3;
gridY = Y(1:stride:end) ./ 1e3;
plot(ax, gridX(:), gridY(:), ".", "Color", options.pointColor, "MarkerSize", 4);
plot(ax, meshTrace(:, 1) ./ 1e3, meshTrace(:, 2) ./ 1e3, "-", ...
    "Color", options.faultColor, "LineWidth", 1.8);

axis(ax, "equal");
box(ax, "on");
grid(ax, "on");
xlabel(ax, "X (km)");
ylabel(ax, "Y (km)");
if isempty(options.axisRange)
    xlim(ax, xRange ./ 1e3);
    ylim(ax, yRange ./ 1e3);
else
    xlim(ax, options.axisRange(1:2));
    ylim(ax, options.axisRange(3:4));
end
if options.showTitle
    titleText = string(options.title);
    if strlength(titleText) == 0
        titleText = sprintf("Full XYZ grid: %d points; D1=%.2f km, D2=%.2f km, M=%.2f km", ...
            numel(X), D1 / 1e3, D2 / 1e3, M / 1e3);
    end
    title(ax, titleText);
end

p = slipsolve.project.paths(cfg);
if exist(p.compositeGreensFigures, "dir") ~= 7
    mkdir(p.compositeGreensFigures);
end
base = fullfile(p.compositeGreensFigures, "full_resolution_grid_preview");
savefig(fig, base + ".fig");
exportgraphics(fig, base + ".png", "Resolution", 180);
figureFiles = [base + ".fig"; base + ".png"];
end

function options = plot_options(cfg)
options = struct("enabled", true, "pointStride", 10, ...
    "pointColor", [0.72 0.76 0.80], "faultColor", [0.75 0.10 0.08], ...
    "axisRange", [], "showTitle", true, "title", "");
if isfield(cfg.visualization, "compositeGreens")
    configured = cfg.visualization.compositeGreens;
    fields = fieldnames(configured);
    for k = 1:numel(fields)
        options.(fields{k}) = configured.(fields{k});
    end
end
end

function visible = figure_visibility(cfg)
visible = "on";
if isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end
end
