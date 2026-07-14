function files = plot_l_curve(cfg, result)
%PLOT_L_CURVE Plot model roughness against the configured data-misfit metric.

p = slipsolve.project.paths(cfg);
files = struct();
baseName = result.inversionType+"_l_curve";
files.fig = fullfile(p.lCurveFigures, baseName+".fig");
files.png = fullfile(p.lCurveFigures, baseName+".png");
if ~isfield(cfg.visualization, "lCurve") || ...
        ~get_field(cfg.visualization.lCurve, "enabled", true)
    return
end
options = cfg.visualization.lCurve;
visible = "on";
if isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end
h = figure("Visible", visible, "Color", "w", ...
    "Name", upper(result.inversionType)+" smoothness L-curve");
ax = axes(h);
plot(ax, result.modelRoughness, result.plottedMisfit, "-o", ...
    "Color", options.lineColor, "LineWidth", options.lineWidth, ...
    "MarkerSize", options.markerSize, "MarkerFaceColor", options.markerColor, ...
    "MarkerEdgeColor", options.lineColor);
if get_field(options, "showSmoothnessLabels", true)
    xMaximum = max(result.modelRoughness);
    yMaximum = max(result.plottedMisfit);
    for k = 1:numel(result.smoothnessValues)
        horizontal = "left";
        label = compose("  %.4g", result.smoothnessValues(k));
        if result.modelRoughness(k) == xMaximum
            horizontal = "right";
            label = compose("%.4g  ", result.smoothnessValues(k));
        end
        vertical = "bottom";
        if result.plottedMisfit(k) == yMaximum
            vertical = "top";
        end
        text(ax, result.modelRoughness(k), result.plottedMisfit(k), label, ...
            "FontSize", 10, "Color", options.lineColor, ...
            "BackgroundColor", "w", "Margin", 1, ...
            "HorizontalAlignment", horizontal, ...
            "VerticalAlignment", vertical, "Clipping", "on");
    end
end
grid(ax, "on");
box(ax, "on");
xlabel(ax, "Model roughness (RMS)");
ylabel(ax, result.misfitLabel);
set(ax, "XScale", char(options.xScale), "YScale", char(options.yScale), ...
    "FontSize", 16, "LineWidth", 1);
if ~isempty(get_field(options, "xRange", [])), xlim(ax, options.xRange); end
if ~isempty(get_field(options, "yRange", [])), ylim(ax, options.yRange); end
if get_field(options, "showTitle", true)
    titleText = string(get_field(options, "title", ""));
    if strlength(titleText) == 0
        titleText = upper(result.inversionType)+" smoothness L-curve";
    end
    title(ax, titleText);
end
if get_field(cfg.visualization, "saveFig", true), savefig(h, files.fig); end
if get_field(cfg.visualization, "savePng", true)
    exportgraphics(h, files.png, "Resolution", 220);
end
drawnow;
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
