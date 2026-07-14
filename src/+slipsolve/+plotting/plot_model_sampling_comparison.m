function figurePaths = plot_model_sampling_comparison(cfg, trackName, rawData, sampledData)
%PLOT_MODEL_SAMPLING_COMPARISON Compare raw data with model-based samples.

p = slipsolve.project.paths(cfg);
if exist(p.modelSamplingFigures, "dir") ~= 7
    mkdir(p.modelSamplingFigures);
end

safeName = local_safe_name(trackName);
figurePaths = struct();
figurePaths.fig = fullfile(p.modelSamplingFigures, safeName + "_model_sampling_comparison.fig");
figurePaths.png = fullfile(p.modelSamplingFigures, safeName + "_model_sampling_comparison.png");

visible = "on";
if isfield(cfg, "visualization") && isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end
options = cfg.visualization.modelSampling;

h = figure("Visible", visible, "Color", "w", "Name", trackName + " model-based sampling");
tiledlayout(h, 1, 2, "TileSpacing", "compact", "Padding", "compact");
colorLimits = get_field(options, "colorbarRange", []);
if isempty(colorLimits), colorLimits = shared_color_limits(rawData, sampledData); end

axRaw = nexttile;
plot_points(rawData, 5, string(options.beforeTitle), options.showTitles);
axis equal tight;
clim(axRaw, colorLimits);
apply_axis_range(axRaw, options.axisRange);
colorbar;

axSampled = nexttile;
plot_points(sampledData, 18, string(options.afterTitle), options.showTitles);
axis equal tight;
clim(axSampled, colorLimits);
apply_axis_range(axSampled, options.axisRange);
colorbar;

colormap(h, char(options.colormap));
save_requested_figures(cfg, h, figurePaths.fig, figurePaths.png);

end

function plot_points(data, markerSize, titleText, showTitle)
x = get_field(data, "x", []);
y = get_field(data, "y", []);
z = get_field(data, "los", []);

valid = isfinite(x) & isfinite(y) & isfinite(z);
if ~any(valid)
    text(0.5, 0.5, "No finite data", "HorizontalAlignment", "center");
    if showTitle, title(titleText); end
    axis off;
    return
end

scatter(x(valid) ./ 1e3, y(valid) ./ 1e3, markerSize, z(valid), "filled");
xlabel("East (km)");
ylabel("North (km)");
if showTitle, title(titleText); end
grid on;
end

function apply_axis_range(ax, range)
if isempty(range), return; end
xlim(ax, range(1:2));
ylim(ax, range(3:4));
end

function colorLimits = shared_color_limits(rawData, sampledData)
z = [
    get_field(rawData, "los", [])
    get_field(sampledData, "los", [])
];
z = z(isfinite(z));

if isempty(z)
    colorLimits = [0, 1];
    return
end

q = prctile(z, [1, 99]);
if q(1) == q(2)
    colorLimits = [min(z), max(z)];
else
    colorLimits = q;
end
if colorLimits(1) == colorLimits(2)
    colorLimits = colorLimits + [-1, 1];
end
end

function save_requested_figures(cfg, h, figPath, pngPath)
saveFig = true;
savePng = true;

if isfield(cfg, "visualization")
    if isfield(cfg.visualization, "saveFig")
        saveFig = cfg.visualization.saveFig;
    end
    if isfield(cfg.visualization, "savePng")
        savePng = cfg.visualization.savePng;
    end
end

if saveFig
    savefig(h, figPath);
end

if savePng
    exportgraphics(h, pngPath, "Resolution", 200);
end
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name)
    value = s.(name);
else
    value = defaultValue;
end
end

function safeName = local_safe_name(name)
safeName = regexprep(string(name), "[^A-Za-z0-9_\\-]+", "_");
safeName = regexprep(safeName, "_+$", "");
if strlength(safeName) == 0
    safeName = "insar_track";
end
end
