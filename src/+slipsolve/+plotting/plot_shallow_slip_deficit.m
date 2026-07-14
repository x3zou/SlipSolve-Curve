function files = plot_shallow_slip_deficit(cfg, result)
%PLOT_SHALLOW_SLIP_DEFICIT Plot normalized slip and deficit against depth.

p = slipsolve.project.paths(cfg);
baseName = result.modelType+"_shallow_slip_deficit";
files = struct();
files.fig = fullfile(p.shallowSlipDeficitFigures, baseName+".fig");
files.png = fullfile(p.shallowSlipDeficitFigures, baseName+".png");
if ~isfield(cfg.visualization, "shallowSlipDeficit") || ...
        ~get_field(cfg.visualization.shallowSlipDeficit, "enabled", true)
    return
end
options = cfg.visualization.shallowSlipDeficit;
visible = "on";
if isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end

h = figure("Visible", visible, "Color", "w", ...
    "Name", upper(result.modelType)+" shallow slip deficit", ...
    "Position", [100 100 1200 620]);
layout = tiledlayout(h, 1, 2, "TileSpacing", "compact", "Padding", "compact");
axSlip = nexttile(layout);
axDeficit = nexttile(layout);
hold(axSlip, "on");
hold(axDeficit, "on");
colors = color_set(options.colormap, numel(result.profiles));
handles = gobjects(numel(result.profiles), 1);
for k = 1:numel(result.profiles)
    profile = result.profiles(k);
    handles(k) = plot(axSlip, profile.normalizedSlip, profile.depthKm, "-o", ...
        "Color", colors(k, :), "LineWidth", options.lineWidth, ...
        "MarkerSize", options.markerSize, "MarkerFaceColor", colors(k, :));
    plot(axDeficit, profile.deficitFraction, profile.depthKm, "-o", ...
        "Color", colors(k, :), "LineWidth", options.lineWidth, ...
        "MarkerSize", options.markerSize, "MarkerFaceColor", colors(k, :));
end

format_axis(axSlip, "Normalized slip", options.normalizedSlipRange, options, result);
format_axis(axDeficit, "Shallow slip deficit fraction", options.deficitRange, options, result);
if get_field(options, "showShallowCutoff", true)
    cutoff = result.shallowDepthMaxKm;
    yline(axSlip, cutoff, "--", "Shallow interval", ...
        "Color", [0.35 0.35 0.35], "LineWidth", 1, "HandleVisibility", "off");
    yline(axDeficit, cutoff, "--", "Shallow interval", ...
        "Color", [0.35 0.35 0.35], "LineWidth", 1, "HandleVisibility", "off");
end
legend(axSlip, handles, string({result.profiles.name}), ...
    "Location", "best", "Interpreter", "none");
layout.YLabel.String = "Depth (km)";
layout.YLabel.FontSize = 18;
if get_field(options, "showTitle", true)
    titleText = string(get_field(options, "title", ""));
    if strlength(titleText) == 0
        titleText = upper(result.modelType)+" shallow slip deficit";
    end
    title(layout, titleText, "FontSize", 20, "FontWeight", "bold");
end
if get_field(cfg.visualization, "saveFig", true), savefig(h, files.fig); end
if get_field(cfg.visualization, "savePng", true)
    exportgraphics(h, files.png, "Resolution", 220);
end
drawnow;
end

function format_axis(ax, xLabel, xRange, options, result)
grid(ax, "on");
box(ax, "on");
set(ax, "YDir", "reverse", "FontSize", 14, "LineWidth", 1);
xlabel(ax, xLabel);
xlim(ax, xRange);
depthRange = get_field(options, "depthRangeKm", []);
if isempty(depthRange)
    maximumDepth = max(arrayfun(@(p) max(p.depthKm), result.profiles));
    depthRange = [0, maximumDepth];
end
ylim(ax, depthRange);
end

function colors = color_set(name, count)
try
    palette = feval(char(string(name)), max(count, 3));
catch
    error("SlipSolve:InvalidPlotConfig", "Unknown shallow-deficit colormap: %s", name);
end
colors = palette(1:count, :);
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
