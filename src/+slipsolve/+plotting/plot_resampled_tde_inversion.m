function figurePaths = plot_resampled_tde_inversion(cfg, tdeResult, mesh)
%PLOT_RESAMPLED_TDE_INVERSION Plot slip and resampled data fits.

p = slipsolve.project.paths(cfg);
if exist(p.resampledTdeFigures, "dir") ~= 7
    mkdir(p.resampledTdeFigures);
end

figurePaths = struct();
figurePaths.slipFig = fullfile(p.resampledTdeFigures, "resampled_tde_slip.fig");
figurePaths.slipPng = fullfile(p.resampledTdeFigures, "resampled_tde_slip.png");
figurePaths.polarityFig = fullfile(p.resampledTdeFigures, "tde_polarity_sanity.fig");
figurePaths.polarityPng = fullfile(p.resampledTdeFigures, "tde_polarity_sanity.png");
figurePaths.trackFigures = repmat(struct("name", "", "fig", "", "png", ""), numel(tdeResult.tracks), 1);

visible = "on";
if isfield(cfg, "visualization") && isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end
plotOpts = tde_plot_options(cfg);

slip = hypot(tdeResult.slipModel(:, 2), tdeResult.slipModel(:, 3));
h = figure("Visible", visible, "Color", "w", "Name", "Resampled TDE slip");
ax = axes(h);
trisurf(mesh.triangles, mesh.nodes(:, 1), mesh.nodes(:, 2), mesh.nodes(:, 3), slip, ...
    "FaceColor", "flat", "FaceAlpha", 0.92, "EdgeColor", [0.18 0.18 0.18], "LineWidth", 0.25);
axis equal;
grid on;
xlabel("East (km)");
ylabel("North (km)");
zlabel("Depth (km)");
apply_title(ax, plotOpts.slip.showTitle, plotOpts.slip.title, ...
    sprintf("Resampled TDE Slip: VR %.1f%%", tdeResult.varianceReduction));
view(plotOpts.slip.view);
set(gca, "ZDir", "normal");
apply_axis_range(ax, plotOpts.slip.axisRange, true);
apply_colorbar_range(ax, plotOpts.slip.colorbarRange);
colormap(ax, char(plotOpts.slip.colormap));
colorbar;
save_requested_figures(cfg, h, figurePaths.slipFig, figurePaths.slipPng);

if plotOpts.polarity.enabled
    hPolarity = plot_polarity(tdeResult, mesh, visible, plotOpts.polarity);
    save_requested_figures(cfg, hPolarity, figurePaths.polarityFig, figurePaths.polarityPng);
else
    figurePaths.polarityFig = "";
    figurePaths.polarityPng = "";
end

for k = 1:numel(tdeResult.tracks)
    track = tdeResult.tracks(k);
    safeName = safe_file_name(track.name);
    figurePaths.trackFigures(k).name = string(track.name);
    figurePaths.trackFigures(k).fig = fullfile(p.resampledTdeFigures, safeName + "_fit.fig");
    figurePaths.trackFigures(k).png = fullfile(p.resampledTdeFigures, safeName + "_fit.png");
    hTrack = plot_one_track(track, visible, plotOpts.fit);
    save_requested_figures(cfg, hTrack, figurePaths.trackFigures(k).fig, figurePaths.trackFigures(k).png);
end

end

function h = plot_polarity(tdeResult, mesh, visible, opts)
h = figure("Visible", visible, "Color", "w", "Name", "TDE polarity sanity");
tiledlayout(h, 1, 2, "TileSpacing", "compact", "Padding", "compact");

ax = nexttile;
plot_mesh_component(ax, mesh, tdeResult.slipModel(:, 2), opts.axisRange, opts.view);
apply_colorbar_range(ax, opts.strikeColorbarRange);
apply_title(ax, opts.showTitles, opts.strikeTitle, "Strike-slip component (cm)");
colorbar;

ax = nexttile;
plot_mesh_component(ax, mesh, tdeResult.slipModel(:, 3), opts.axisRange, opts.view);
apply_colorbar_range(ax, opts.dipColorbarRange);
apply_title(ax, opts.showTitles, opts.dipTitle, "Dip-slip component (cm)");
colorbar;

colormap(h, char(opts.colormap));
end

function plot_mesh_component(ax, mesh, values, axisRange, viewAngles)
axes(ax);
trisurf(mesh.triangles, mesh.nodes(:, 1), mesh.nodes(:, 2), mesh.nodes(:, 3), values, ...
    "FaceColor", "flat", "FaceAlpha", 0.94, "EdgeColor", [0.18 0.18 0.18], "LineWidth", 0.2);
axis(ax, "equal");
grid(ax, "on");
xlabel(ax, "East (km)");
ylabel(ax, "North (km)");
zlabel(ax, "Depth (km)");
view(ax, viewAngles);
set(ax, "ZDir", "normal");
apply_axis_range(ax, axisRange, true);
clim(ax, centered_color_limits(values));
end

function h = plot_one_track(track, visible, opts)
h = figure("Visible", visible, "Color", "w", "Name", string(track.name) + " resampled TDE fit");
tiledlayout(h, 1, 3, "TileSpacing", "compact", "Padding", "compact");
values = [track.data(:); track.prediction(:)];
limits = configured_or_centered_limits(opts.dataModelColorbarRange, values);
panelTitles = string(opts.panelTitles);

ax = nexttile;
scatter(track.x ./ 1e3, track.y ./ 1e3, 18, track.data, "filled");
axis equal tight;
clim(ax, limits);
apply_axis_range(ax, opts.axisRange, false);
apply_title(ax, opts.showTitles, panelTitles(1), "Data");
xlabel("East (km)");
ylabel("North (km)");
grid on;
colorbar;

ax = nexttile;
scatter(track.x ./ 1e3, track.y ./ 1e3, 18, track.prediction, "filled");
axis equal tight;
clim(ax, limits);
apply_axis_range(ax, opts.axisRange, false);
apply_title(ax, opts.showTitles, panelTitles(2), "Model");
xlabel("East (km)");
ylabel("North (km)");
grid on;
colorbar;

ax = nexttile;
scatter(track.x ./ 1e3, track.y ./ 1e3, 18, track.residual, "filled");
axis equal tight;
clim(ax, configured_or_centered_limits(opts.residualColorbarRange, track.residual));
apply_axis_range(ax, opts.axisRange, false);
apply_title(ax, opts.showTitles, panelTitles(3), "Residual");
xlabel("East (km)");
ylabel("North (km)");
grid on;
colorbar;

colormap(h, char(opts.colormap));
end

function opts = tde_plot_options(cfg)
opts = default_tde_plot_options();
if ~isfield(cfg, "visualization") || ~isfield(cfg.visualization, "resampledTde")
    return
end
user = cfg.visualization.resampledTde;
if isfield(user, "slipPlot")
    opts.slip = merge_struct(opts.slip, user.slipPlot);
end
if isfield(user, "fitPlots")
    opts.fit = merge_struct(opts.fit, user.fitPlots);
end
if isfield(user, "polarityPlot")
    opts.polarity = merge_struct(opts.polarity, user.polarityPlot);
end
end

function opts = default_tde_plot_options()
opts.slip = struct( ...
    "colormap", "parula", ...
    "colorbarRange", [], ...
    "axisRange", [], ...
    "showTitle", true, ...
    "title", "", ...
    "view", [35 25]);
opts.fit = struct( ...
    "colormap", "turbo", ...
    "dataModelColorbarRange", [], ...
    "residualColorbarRange", [], ...
    "axisRange", [], ...
    "showTitles", true, ...
    "panelTitles", ["Data" "Model" "Residual"]);
opts.polarity = struct( ...
    "enabled", true, ...
    "colormap", "turbo", ...
    "strikeColorbarRange", [], ...
    "dipColorbarRange", [], ...
    "axisRange", [], ...
    "showTitles", true, ...
    "strikeTitle", "Strike-slip component (cm)", ...
    "dipTitle", "Dip-slip component (cm)", ...
    "view", [35 25]);
end

function out = merge_struct(out, user)
names = fieldnames(user);
for k = 1:numel(names)
    out.(names{k}) = user.(names{k});
end
end

function limits = configured_or_centered_limits(configuredLimits, values)
if ~isempty(configuredLimits)
    limits = configuredLimits;
else
    limits = centered_color_limits(values);
end
end

function apply_colorbar_range(ax, limits)
if ~isempty(limits)
    clim(ax, limits);
end
end

function apply_axis_range(ax, axisRange, hasZ)
if isempty(axisRange)
    return
end
xlim(ax, axisRange(1:2));
ylim(ax, axisRange(3:4));
if hasZ && numel(axisRange) >= 6
    zlim(ax, axisRange(5:6));
end
end

function apply_title(ax, showTitle, titleText, defaultTitle)
if ~showTitle
    title(ax, "");
    return
end
if strlength(string(titleText)) == 0
    title(ax, defaultTitle);
else
    title(ax, titleText);
end
end

function limits = centered_color_limits(values)
values = values(isfinite(values));
if isempty(values)
    limits = [-1, 1];
    return
end
maxAbs = prctile(abs(values), 99);
if isempty(maxAbs) || maxAbs == 0
    maxAbs = max(abs(values));
end
if isempty(maxAbs) || maxAbs == 0
    maxAbs = 1;
end
limits = [-maxAbs, maxAbs];
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
    exportgraphics(h, pngPath, "Resolution", 220);
end
end

function name = safe_file_name(value)
name = regexprep(string(value), "[^A-Za-z0-9_.-]", "_");
name = regexprep(name, "_+$", "");
if strlength(name) == 0
    name = "track";
end
end
