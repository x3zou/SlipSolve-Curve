function compositeGreens = build_composite_greens(cfg, layeredGreens, mesh, observations)
%BUILD_COMPOSITE_GREENS Build full-resolution layered G_e, G_n, and G_u.

if ~isfield(cfg, "compositeGreens") || ~cfg.compositeGreens.enabled
    error("SlipSolve:CompositeGreensDisabled", ...
        "Set cfg.compositeGreens.enabled = true before running composite_greens.");
end
if nargin < 3 || isempty(mesh) || ~isfield(mesh, "sections")
    error("SlipSolve:MissingInput", "A saved triangular mesh is required for composite_greens.");
end
if nargin < 4
    observations = struct();
end

p = slipsolve.project.paths(cfg);
outputDirectory = string(get_field(cfg.compositeGreens, "outputDirectory", p.compositeGreensDirectory));
if strlength(outputDirectory) == 0
    outputDirectory = p.compositeGreensDirectory;
end
if exist(outputDirectory, "dir") ~= 7
    mkdir(outputDirectory);
end

add_legacy_composite_paths(cfg);
validate_layered_products(layeredGreens);
greenName = layered_green_name(layeredGreens);

meshTraceFile = string(cfg.compositeGreens.meshTraceFile);
if exist(meshTraceFile, "file") ~= 2
    error("SlipSolve:MissingInput", "Composite-grid fault trace is missing: %s", meshTraceFile);
end
meshTrace = readmatrix(meshTraceFile);
meshTrace = meshTrace(:, 1:2);
meshTrace = meshTrace(all(isfinite(meshTrace), 2), :);
if size(meshTrace, 1) < 2
    error("SlipSolve:InvalidInput", "Composite-grid fault trace must contain at least two finite [x y] rows.");
end

xRange = validate_range(cfg.compositeGreens.xRange, "xRange");
yRange = validate_range(cfg.compositeGreens.yRange, "yRange");
D1 = validate_positive(cfg.compositeGreens.D1, "D1");
D2 = validate_positive(cfg.compositeGreens.D2, "D2");
M = validate_nonnegative(cfg.compositeGreens.M, "M");
if D2 <= D1
    error("SlipSolve:InvalidCompositeConfig", "cfg.compositeGreens.D2 must be greater than D1.");
end

[X, Y, xVector, yVector] = adaptive_meshgrid(meshTrace, xRange, yRange, D1, D2, M);
dataInsar = [X(:), Y(:)];
gridFile = fullfile(outputDirectory, "full_xyz_grid.mat");
save(gridFile, "X", "Y", "xVector", "yVector", "dataInsar", "meshTrace", ...
    "xRange", "yRange", "D1", "D2", "M", "-v7.3");

figureFiles = slipsolve.plotting.plot_composite_grid_preview( ...
    cfg, X, Y, meshTrace, xRange, yRange, D1, D2, M);
pause_after_preview(cfg);

overwrite = logical(get_field(cfg.compositeGreens, "overwriteExistingFiles", false));
matrixFiles = struct( ...
    "east", fullfile(outputDirectory, "G_e.mat"), ...
    "north", fullfile(outputDirectory, "G_n.mat"), ...
    "vertical", fullfile(outputDirectory, "G_u.mat"));
protect_existing_outputs(matrixFiles, overwrite);

databaseDirectory = string(layeredGreens.outputDirectory);
previousDirectory = pwd;
cleanup = onCleanup(@() cd(previousDirectory));
cd(databaseDirectory);

buildInfo = slipsolve.composite.build_full_xyz( ...
    cfg, mesh, dataInsar, greenName, matrixFiles);

compositeGreens = struct();
compositeGreens.method = "legacy_Comb_Green_Layered_Nodes_xyz";
compositeGreens.status = "complete";
compositeGreens.greenName = greenName;
compositeGreens.layeredDatabaseDirectory = databaseDirectory;
compositeGreens.gridFile = gridFile;
compositeGreens.matrixFiles = matrixFiles;
compositeGreens.grid = struct("xRange", xRange, "yRange", yRange, ...
    "D1", D1, "D2", D2, "M", M, "xCount", numel(xVector), ...
    "yCount", numel(yVector), "observationCount", size(dataInsar, 1));
compositeGreens.buildInfo = buildInfo;
compositeGreens.figureFiles = figureFiles;
compositeGreens.observationMetadata = observations_metadata(observations);
compositeGreens.createdAt = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");
save(p.compositeGreens, "compositeGreens", "-v7.3");

end

function add_legacy_composite_paths(cfg)
legacyRoot = string(cfg.paths.legacyRoot);
paths = [
    legacyRoot
    fullfile(legacyRoot, "TDE_Greens")
    fullfile(legacyRoot, "layered_green")
    fullfile(legacyRoot, "layered_green", "otherfunc")
    string(cfg.paths.legacyOtherFuncRoot)
    fullfile(string(cfg.project.root), "src", "legacy")
];
for k = 1:numel(paths)
    if strlength(paths(k)) > 0 && exist(paths(k), "dir") == 7
        addpath(char(paths(k)));
    end
end
end

function validate_layered_products(layeredGreens)
if ~isfield(layeredGreens, "products") || numel(layeredGreens.products) ~= 3
    error("SlipSolve:MissingInput", ...
        "layered_greens must first create SS, DS, and CL MAT products.");
end
for k = 1:numel(layeredGreens.products)
    if exist(layeredGreens.products(k).matFile, "file") ~= 2
        error("SlipSolve:MissingInput", ...
            "Converted layered Green's MAT file is missing: %s", layeredGreens.products(k).matFile);
    end
end
end

function name = layered_green_name(layeredGreens)
[~, base, extension] = fileparts(string(layeredGreens.edgrnInputFile));
name = base + extension;
end

function protect_existing_outputs(files, overwrite)
names = string(struct2cell(files));
existing = names(isfile(names));
if ~isempty(existing) && ~overwrite
    error("SlipSolve:ExistingOutput", ...
        "Composite Green's output already exists: %s. Set overwriteExistingFiles=true to replace it.", ...
        existing(1));
end
end

function pause_after_preview(cfg)
shouldPause = logical(get_field(cfg.compositeGreens, "pauseAfterGridPreview", false));
visible = ~isfield(cfg.visualization, "visible") || cfg.visualization.visible;
if shouldPause && visible && usejava("desktop")
    input("Inspect the full-resolution grid preview, then press Enter to start the composite Green's calculation: ", "s");
end
end

function value = validate_range(value, name)
value = double(value(:).');
if numel(value) ~= 2 || any(~isfinite(value)) || value(2) <= value(1)
    error("SlipSolve:InvalidCompositeConfig", ...
        "cfg.compositeGreens.%s must be [minimum maximum].", name);
end
end

function value = validate_positive(value, name)
value = double(value);
if ~isscalar(value) || ~isfinite(value) || value <= 0
    error("SlipSolve:InvalidCompositeConfig", ...
        "cfg.compositeGreens.%s must be a positive scalar.", name);
end
end

function value = validate_nonnegative(value, name)
value = double(value);
if ~isscalar(value) || ~isfinite(value) || value < 0
    error("SlipSolve:InvalidCompositeConfig", ...
        "cfg.compositeGreens.%s must be a nonnegative scalar.", name);
end
end

function metadata = observations_metadata(observations)
metadata = struct("provided", ~isempty(fieldnames(observations)), "count", 0);
if isfield(observations, "data")
    metadata.count = numel(observations.data);
end
end

function value = get_field(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
