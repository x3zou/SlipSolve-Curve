function files = plot_forward_model(cfg, track, fullGrid)
%PLOT_FORWARD_MODEL Plot full-resolution independent data/model/residual.

if nargin < 3
    fullGrid = struct();
end

p = slipsolve.project.paths(cfg);
name = regexprep(string(track.name), "[^A-Za-z0-9_.-]", "_");
files = struct("fig", fullfile(p.forwardModelFigures, name+"_forward_fit.fig"), ...
    "png", fullfile(p.forwardModelFigures, name+"_forward_fit.png"));
options = cfg.visualization.forwardModel;
if ~get_field(options, "enabled", true)
    return
end
visible = "on";
if isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end
h = render_plot(track, fullGrid, options, visible);
if get_field(cfg.visualization, "savePng", true)
    exportgraphics(h, files.png, "Resolution", 220, "BackgroundColor", "white");
end
if get_field(cfg.visualization, "saveFig", true)
    if ~isgraphics(h, "figure")
        h = render_plot(track, fullGrid, options, visible);
    end
    savefig(h, files.fig);
end
drawnow;
end

function h = render_plot(track, fullGrid, options, visible)
if isfield(fullGrid, "data")
    h = full_resolution_plot(track, fullGrid, options, visible);
else
    h = sampled_plot(track, options, visible);
end
end

function h = sampled_plot(track, options, visible)
h = figure("Visible", visible, "Color", "w", ...
    "Name", track.name+" independent forward model", ...
    "Position", [100 100 1450 520]);
layout = tiledlayout(h, 1, 3, "TileSpacing", "compact", "Padding", "compact");
values = {track.data, track.prediction, track.residualDataMinusModel};
labels = string(get_field(options, "panelTitles", ...
    ["Independent data" "Forward model" "Data - model"]));
dataLimits = centered_limits(get_field(options, "dataModelColorbarRange", []), ...
    [track.data; track.prediction]);
residualLimits = centered_limits(get_field(options, "residualColorbarRange", []), ...
    track.residualDataMinusModel);
trace = read_fault_trace(get_field(options, "faultTraceFile", ""));
for k = 1:3
    ax = nexttile(layout);
    scatter(ax, track.x/1e3, track.y/1e3, ...
        get_field(options, "markerSize", 18), values{k}, "filled");
    hold(ax, "on");
    if ~isempty(trace)
        plot(ax, trace(:, 1)/1e3, trace(:, 2)/1e3, "k-", "LineWidth", 1.5);
    end
    axis(ax, "equal");
    axis(ax, "tight");
    grid(ax, "on");
    box(ax, "on");
    if k < 3, clim(ax, dataLimits); else, clim(ax, residualLimits); end
    if ~isempty(options.axisRange)
        xlim(ax, options.axisRange(1:2));
        ylim(ax, options.axisRange(3:4));
    end
    if get_field(options, "showTitles", true), title(ax, labels(k)); end
    colorbar(ax);
    set(ax, "FontSize", 13);
end
colormap(h, char(options.colormap));
xlabel(layout, "Easting (km)", "FontSize", 16);
ylabel(layout, "Northing (km)", "FontSize", 16);
titleText = string(get_field(options, "title", ""));
if strlength(titleText) == 0, titleText = track.name; end
if get_field(options, "showOverallTitle", true)
    title(layout, titleText, "FontSize", 18, "FontWeight", "bold");
end
end

function h = full_resolution_plot(track, raster, options, visible)
h = figure("Visible", visible, "Color", "w", ...
    "Name", track.name+" full-resolution independent forward model", ...
    "Position", [100 100 1550 620]);
layout = tiledlayout(h, 1, 3, "TileSpacing", "compact", "Padding", "compact");
values = {raster.data, raster.model, raster.residual};
labels = string(get_field(options, "panelTitles", ...
    ["Independent data" "Forward model" "Data - model"]));
dataLimits = centered_limits(get_field(options, "dataModelColorbarRange", []), ...
    [raster.data(:); raster.model(:)]);
residualLimits = centered_limits(get_field(options, "residualColorbarRange", []), ...
    raster.residual(:));
for k = 1:3
    ax = nexttile(layout);
    xExtent = [min(raster.x, [], "all"), max(raster.x, [], "all")]/1e3;
    yExtent = [min(raster.y, [], "all"), max(raster.y, [], "all")]/1e3;
    imageHandle = imagesc(ax, xExtent, yExtent, single(values{k}));
    imageHandle.AlphaData = isfinite(values{k});
    set(ax, "YDir", "normal");
    hold(ax, "on");
    plot(ax, raster.fault(:, 1)/1e3, raster.fault(:, 2)/1e3, ...
        "w-", "LineWidth", 2);
    plot(ax, raster.fault(:, 1)/1e3, raster.fault(:, 2)/1e3, ...
        "k--", "LineWidth", 1.5);
    axis(ax, "equal");
    grid(ax, "on");
    box(ax, "on");
    if k < 3, clim(ax, dataLimits); else, clim(ax, residualLimits); end
    if ~isempty(options.axisRange)
        xlim(ax, options.axisRange(1:2));
        ylim(ax, options.axisRange(3:4));
    end
    if get_field(options, "showTitles", true), title(ax, labels(k)); end
    cb = colorbar(ax);
    cb.Label.String = "cm";
    set(ax, "FontSize", 14);
    text(ax, 0.08, 0.98, char('a'+k-1), "Units", "normalized", ...
        "HorizontalAlignment", "left", "VerticalAlignment", "top", ...
        "FontWeight", "bold", "FontSize", 22, "Color", "black");
end
colormap(h, char(options.colormap));
xlabel(layout, "Easting (km)", "FontSize", 17);
ylabel(layout, "Northing (km)", "FontSize", 17);
titleText = string(get_field(options, "title", ""));
if strlength(titleText) == 0, titleText = track.name; end
if get_field(options, "showOverallTitle", true)
    title(layout, titleText, "FontSize", 20, "FontWeight", "bold");
end
end

function trace = read_fault_trace(filePath)
trace = [];
if strlength(string(filePath)) == 0 || exist(filePath, "file") ~= 2
    return
end
trace = readmatrix(filePath);
trace = double(trace(:, 1:2));
trace = trace(all(isfinite(trace), 2), :);
end

function limits = centered_limits(configured, values)
if ~isempty(configured)
    limits = configured;
    return
end
maximum = max(abs(values), [], "omitnan");
if isempty(maximum) || ~isfinite(maximum) || maximum == 0, maximum = 1; end
limits = [-maximum maximum];
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
