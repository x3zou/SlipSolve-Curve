function figureFiles = plot_layered_greens(cfg, layeredGreens)
%PLOT_LAYERED_GREENS Plot a compact diagnostic of converted EDGRN databases.

figureFiles = strings(0, 1);
if ~isfield(cfg, "visualization") || ~cfg.visualization.enabled
    return
end

options = plot_options(cfg);
if ~options.enabled
    return
end

p = slipsolve.project.paths(cfg);
fig = figure("Name", "Layered EDGRN Green's functions", "Color", "w", ...
    "Visible", figure_visibility(cfg));
layout = tiledlayout(fig, 1, numel(layeredGreens.products), "TileSpacing", "compact", "Padding", "compact");

for k = 1:numel(layeredGreens.products)
    product = layeredGreens.products(k);
    data = load(product.matFile, "uz", "nr", "nz", "r1", "r2", "z1", "z2");
    ax = nexttile(layout, k);
    values = reshape(data.uz, data.nr, data.nz).';
    if options.log10AbsoluteUz
        values = log10(abs(values) + realmin("double"));
        colorLabel = "log10(abs(uz))";
    else
        colorLabel = "uz";
    end
    rKm = linspace(data.r1, data.r2, data.nr) ./ 1e3;
    zKm = linspace(data.z1, data.z2, data.nz) ./ 1e3;
    imagesc(ax, rKm, zKm, values);
    axis(ax, "xy");
    axis(ax, "tight");
    xlabel(ax, "Radial distance (km)");
    ylabel(ax, "Source depth (km)");
    colormap(ax, char(options.colormap));
    if ~isempty(options.axisRange)
        xlim(ax, options.axisRange(1:2));
        ylim(ax, options.axisRange(3:4));
    end
    if ~isempty(options.colorbarRange), clim(ax, options.colorbarRange); end
    cb = colorbar(ax);
    ylabel(cb, colorLabel);
    if options.showTitles
        panelTitles = string(options.panelTitles);
        if numel(panelTitles) >= k && strlength(panelTitles(k)) > 0
            title(ax, panelTitles(k));
        else
            title(ax, upper(product.type) + " fundamental database");
        end
    end
end
if get_field(options, "showOverallTitle", false)
    title(layout, string(options.title));
end

if exist(p.layeredGreensFigures, "dir") ~= 7
    mkdir(p.layeredGreensFigures);
end
base = fullfile(p.layeredGreensFigures, "layered_edgrn_uz_diagnostic");
savefig(fig, base + ".fig");
exportgraphics(fig, base + ".png", "Resolution", 180);
figureFiles = [base + ".fig"; base + ".png"];
end

function options = plot_options(cfg)
options = struct("enabled", true, "colormap", "turbo", ...
    "log10AbsoluteUz", true, "showTitles", true, ...
    "panelTitles", strings(0, 1), "showOverallTitle", false, ...
    "title", "Layered EDGRN Green's functions", ...
    "axisRange", [], "colorbarRange", []);
if isfield(cfg.visualization, "layeredGreens")
    configured = cfg.visualization.layeredGreens;
    fields = fieldnames(configured);
    for k = 1:numel(fields)
        options.(fields{k}) = configured.(fields{k});
    end
end
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end

function visible = figure_visibility(cfg)
visible = "on";
if isfield(cfg.visualization, "visible") && ~cfg.visualization.visible
    visible = "off";
end
end
