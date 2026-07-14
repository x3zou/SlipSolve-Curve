function layeredGreens = build_layered_greens(cfg, mesh, observations)
%BUILD_LAYERED_GREENS Convert an externally generated EDGRN database to MAT.
%
% EDGRN execution remains external by design. This preserves the legacy
% workflow where the user owns the Earth model and EDGRN input file, while
% keeping all MATLAB products under the SlipSolve project.

if ~isfield(cfg, "layeredGreens") || ~cfg.layeredGreens.enabled
    error("SlipSolve:LayeredGreensDisabled", ...
        "Set cfg.layeredGreens.enabled = true before running layered_greens.");
end

inputFile = string(cfg.layeredGreens.edgrnInputFile);
if strlength(inputFile) == 0 || exist(inputFile, "file") ~= 2
    error("SlipSolve:MissingInput", ...
        "cfg.layeredGreens.edgrnInputFile does not exist: %s", inputFile);
end

p = slipsolve.project.paths(cfg);
outputDirectory = string(cfg.layeredGreens.outputDirectory);
if strlength(outputDirectory) == 0
    outputDirectory = p.layeredGreensDatabase;
end
overwrite = get_field(cfg.layeredGreens, "overwriteExistingMatFiles", false);

products = slipsolve.edgrn.convert_legacy_getedgrn(inputFile, outputDirectory, overwrite);

layeredGreens = struct();
layeredGreens.method = "legacy_getedgrn_compatible";
layeredGreens.edgrnInputFile = inputFile;
layeredGreens.outputDirectory = outputDirectory;
layeredGreens.edgrnExecutable = string(get_field(cfg.paths, "edgrnExecutable", ""));
layeredGreens.products = products;
layeredGreens.createdAt = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");
layeredGreens.meshPatchCount = mesh_patch_count(mesh);
layeredGreens.observationCount = observation_count(observations);

save(p.layeredGreens, "layeredGreens");
layeredGreens.figureFiles = slipsolve.plotting.plot_layered_greens(cfg, layeredGreens);

end

function value = get_field(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function count = mesh_patch_count(mesh)
count = 0;
if isfield(mesh, "triangles")
    count = size(mesh.triangles, 1);
elseif isfield(mesh, "sections")
    for k = 1:numel(mesh.sections)
        if isfield(mesh.sections(k), "triangles")
            count = count + size(mesh.sections(k).triangles, 1);
        end
    end
end
end

function count = observation_count(observations)
count = 0;
if isfield(observations, "data")
    count = numel(observations.data);
end
end
