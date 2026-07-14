function insarSub = subsample_insar(cfg)
%SUBSAMPLE_INSAR Quadtree subsampling wrapper.

p = slipsolve.project.paths(cfg);
slipsolve.project.ensure_directories(cfg);

if isfield(cfg.insar, "tracks") && ~isempty(cfg.insar.tracks)
    trackDefs = cfg.insar.tracks(:);
else
    trackDefs = sample_list_to_track_defs(cfg);
end

trackCells = cell(numel(trackDefs), 1);
for k = 1:numel(trackDefs)
    trackCells{k} = slipsolve.stages.subsample_insar_track(cfg, trackDefs(k));
end
tracks = vertcat(trackCells{:});

insarSub = assemble_insar_sub(tracks, cfg);
save(p.insarSubsampled, "insarSub", "-v7.3");

end

function trackDefs = sample_list_to_track_defs(cfg)
sampleRows = slipsolve.io.read_sample_list(cfg.paths.sampleListFile);
trackCells = cell(numel(sampleRows), 1);

for k = 1:numel(sampleRows)
    row = sampleRows(k);
    trackDef = struct();
    trackDef.name = sprintf("track_%02d_%s", k, erase(row.outputFile, ".mat"));
    trackDef.dataType = row.type;
    trackDef.dataFile = row.gridFullPath;
    trackDef.boundsLonLat = row.boundsLonLat;
    trackDef.minCellSizeKm = row.minCellSizeKm;
    trackDef.maxCellSizeKm = row.maxCellSizeKm;
    trackDef.rmsThreshold = row.rmsThreshold;
    trackDef.padX = row.padX;
    trackDef.padY = row.padY;
    trackDef.outputFile = row.outputFile;
    trackDef.lookEFile = infer_look_file(row.directoryFullPath, row.type, "e", cfg.insar.inputResolution);
    trackDef.lookNFile = infer_look_file(row.directoryFullPath, row.type, "n", cfg.insar.inputResolution);
    trackDef.lookUFile = infer_look_file(row.directoryFullPath, row.type, "u", cfg.insar.inputResolution);
    trackCells{k} = trackDef;
end
trackDefs = vertcat(trackCells{:});
end

function filePath = infer_look_file(trackDir, dataType, component, inputResolution)
suffix = "";
if inputResolution == "low"
    suffix = "_low";
elseif inputResolution == "high"
    suffix = "_high";
end

if dataType == "azo"
    fileName = "look_" + component + "_azo" + suffix + ".grd";
elseif dataType == "rng"
    fileName = "look_" + component + "_rng" + suffix + ".grd";
else
    fileName = "look_" + component + suffix + ".grd";
end

filePath = fullfile(trackDir, fileName);
end

function insarSub = assemble_insar_sub(tracks, cfg)
insarSub = struct();
insarSub.tracks = tracks;
insarSub.metadata.origin = cfg.insar.origin;
insarSub.metadata.units.xy = "m";
insarSub.metadata.units.displacement = "input_grid_units";
insarSub.metadata.quadtreeFunction = "quadtree_unstructured2";

insarSub.all = struct();
insarSub.all.x = vertcat(tracks.x);
insarSub.all.y = vertcat(tracks.y);
insarSub.all.los = vertcat(tracks.los);
insarSub.all.lookE = vertcat(tracks.lookE);
insarSub.all.lookN = vertcat(tracks.lookN);
insarSub.all.lookU = vertcat(tracks.lookU);
end
