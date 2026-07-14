function figurePaths = plot_insar_subsampling(cfg, track, rawData)
%PLOT_INSAR_SUBSAMPLING Plot raw and quadtree-sampled InSAR observations.

if nargin < 3
    rawData = struct();
end

p = slipsolve.project.paths(cfg);
if exist(p.insarSubsamplingFigures, "dir") ~= 7
    mkdir(p.insarSubsamplingFigures);
end

trackName = get_field(track, "name", "");
if strlength(trackName) == 0
    trackName = get_field(track, "outputFile", "insar_track");
end
trackName = local_safe_name(trackName);

figurePaths = struct();
figurePaths.rawVsSampledFig = fullfile(p.insarSubsamplingFigures, trackName + "_raw_vs_sampled.fig");
figurePaths.rawVsSampledPng = fullfile(p.insarSubsamplingFigures, trackName + "_raw_vs_sampled.png");
figurePaths.rmsHistogramFig = fullfile(p.insarSubsamplingFigures, trackName + "_rms_histogram.fig");
figurePaths.rmsHistogramPng = fullfile(p.insarSubsamplingFigures, trackName + "_rms_histogram.png");

visible = "on";
if isfield(cfg, "visualization") && isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end
options = cfg.visualization.insarSubsampling;

h = figure("Visible", visible, "Color", "w", "Name", trackName + " raw vs sampled");
tiledlayout(h, 1, 2, "TileSpacing", "compact", "Padding", "compact");
colorLimits = get_field(options, "colorbarRange", []);
if isempty(colorLimits), colorLimits = shared_color_limits(rawData, track); end

axRaw = nexttile;
plot_points(rawData, string(options.rawTitle), options.showTitles);
axis equal tight;
clim(axRaw, colorLimits);
apply_axis_range(axRaw, options.axisRange);
colorbar;

axSampled = nexttile;
plot_points(track, string(options.sampledTitle), options.showTitles);
axis equal tight;
clim(axSampled, colorLimits);
apply_axis_range(axSampled, options.axisRange);
colorbar;

colormap(h, char(options.colormap));
save_requested_figures(cfg, h, figurePaths.rawVsSampledFig, figurePaths.rawVsSampledPng);

if isfield(track, "rms") && ~isempty(track.rms)
    hRms = figure("Visible", visible, "Color", "w", "Name", trackName + " RMS histogram");
    histogram(track.rms(isfinite(track.rms)), 40);
    xlabel("Per-cell spread");
    ylabel("Count");
    if options.showTitles
        titleText = string(options.histogramTitle);
        if strlength(titleText) == 0, titleText = trackName + " quadtree cell spread"; end
        title(titleText);
    end
    grid on;
    save_requested_figures(cfg, hRms, figurePaths.rmsHistogramFig, figurePaths.rmsHistogramPng);
else
    figurePaths.rmsHistogramFig = "";
    figurePaths.rmsHistogramPng = "";
end

end

function plot_points(data, titleText, showTitle)
x = get_field(data, "x", []);
y = get_field(data, "y", []);
z = get_field(data, "los", []);

if isempty(x) || isempty(y) || isempty(z)
    text(0.5, 0.5, "No data loaded", "HorizontalAlignment", "center");
    if showTitle, title(titleText); end
    axis off;
    return
end

scatter(x(:) ./ 1e3, y(:) ./ 1e3, 12, z(:), "filled");
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

function colorLimits = shared_color_limits(rawData, track)
z = [
    get_field(rawData, "los", [])
    get_field(track, "los", [])
];
z = z(isfinite(z));

if isempty(z)
    colorLimits = [0, 1];
    return
end

colorLimits = [min(z), max(z)];
if colorLimits(1) == colorLimits(2)
    colorLimits = colorLimits + [-1, 1];
end
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name)
    value = s.(name);
else
    value = defaultValue;
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

function safeName = local_safe_name(name)
safeName = regexprep(string(name), "[^A-Za-z0-9_\\-]+", "_");
safeName = regexprep(safeName, "_+$", "");
if strlength(safeName) == 0
    safeName = "insar_track";
end
end
