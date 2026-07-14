function shallowSlipDeficitResult = run_shallow_slip_deficit(cfg)
%RUN_SHALLOW_SLIP_DEFICIT Build normalized slip and deficit depth profiles.

if ~isfield(cfg, "shallowSlipDeficit") || ...
        ~get_field(cfg.shallowSlipDeficit, "enabled", false)
    error("SlipSolve:DisabledHelper", ...
        "Set cfg.shallowSlipDeficit.enabled = true before running this helper.");
end

slipsolve.project.ensure_directories(cfg);
p = slipsolve.project.paths(cfg);
settings = cfg.shallowSlipDeficit;
modelTypes = normalize_model_types(settings.modelTypes);
segments = validate_segments(settings.segments);
modelCells = cell(numel(modelTypes), 1);

for k = 1:numel(modelTypes)
    modelType = modelTypes(k);
    fprintf("Computing %s shallow slip deficit.\n", modelType);
    if modelType == "composite"
        [samples, sourceProduct, method] = composite_samples(p, settings);
    else
        [samples, sourceProduct, method] = tde_samples(p, settings);
    end
    profiles = build_profiles(samples, segments, modelType, settings);
    [profileTable, summaryTable] = result_tables(modelType, profiles, settings);
    profileTableFile = fullfile(p.shallowSlipDeficitDirectory, ...
        modelType+"_shallow_slip_profile.csv");
    summaryTableFile = fullfile(p.shallowSlipDeficitDirectory, ...
        modelType+"_shallow_slip_summary.csv");
    writetable(profileTable, profileTableFile);
    writetable(summaryTable, summaryTableFile);

    modelResult = struct();
    modelResult.modelType = modelType;
    modelResult.method = method;
    modelResult.legacyReference = "Shallow_Slip_Deficit_Segmented_New.m";
    modelResult.sourceProduct = sourceProduct;
    modelResult.component = string(settings.component);
    modelResult.referenceDepthRangeKm = settings.referenceDepthRangeKm;
    modelResult.shallowDepthMaxKm = settings.shallowDepthMaxKm;
    modelResult.profiles = profiles;
    modelResult.profileTableFile = profileTableFile;
    modelResult.summaryTableFile = summaryTableFile;
    modelResult.figureFiles = slipsolve.plotting.plot_shallow_slip_deficit( ...
        cfg, modelResult);
    modelCells{k} = modelResult;
end

shallowSlipDeficitResult = struct();
shallowSlipDeficitResult.models = vertcat(modelCells{:});
shallowSlipDeficitResult.definition = ...
    "deficit_fraction = 1 - mean_slip_at_depth/reference_profile_slip";
shallowSlipDeficitResult.createdAt = datetime("now", ...
    "Format", "yyyy-MM-dd HH:mm:ss");
shallowSlipDeficitResult.resultFile = p.shallowSlipDeficitResult;
save(p.shallowSlipDeficitResult, "shallowSlipDeficitResult", "-v7.3");
fprintf("Saved shallow slip deficit result: %s\n", p.shallowSlipDeficitResult);
end

function [samples, sourceProduct, method] = composite_samples(p, settings)
if exist(p.finalResult, "file") ~= 2
    error("SlipSolve:MissingStageProduct", ...
        "Run final_inversion before composite shallow-deficit analysis: %s", ...
        p.finalResult);
end
loaded = load(p.finalResult, "finalResult");
model = loaded.finalResult.completeNodalSlipModel;
if size(model, 2) < 8
    error("SlipSolve:InvalidStageProduct", ...
        "finalResult.completeNodalSlipModel must have at least 8 columns.");
end
samples = struct();
samples.faultId = model(:, 1);
samples.layerId = model(:, 3);
samples.xKm = model(:, 6)/1e3;
samples.yKm = model(:, 7)/1e3;
samples.depthKm = abs(model(:, 8))/1e3;
samples.slipCm = component_amplitude(model(:, 4), model(:, 5), settings.component);
samples.areaM2 = nan(size(model, 1), 1);
sourceProduct = string(p.finalResult);
method = "legacy_nodal_arithmetic_mean_by_native_layer";
end

function [samples, sourceProduct, method] = tde_samples(p, settings)
tdeModel = load_tde_postprocessing_model(p);
model = tdeModel.slipModel;
points = tdeModel.points(:, 1:3);
if size(model, 2) < 6
    error("SlipSolve:InvalidStageProduct", ...
        "The TDE slip model must contain triangle vertex columns 4:6.");
end
vertices = round(model(:, 4:6));
if any(vertices(:) < 1 | vertices(:) > size(points, 1))
    error("SlipSolve:InvalidStageProduct", ...
        "The TDE slip model contains invalid triangle vertex indices.");
end
v1 = points(vertices(:, 1), :);
v2 = points(vertices(:, 2), :);
v3 = points(vertices(:, 3), :);
centers = (v1+v2+v3)/3;
area = 0.5*vecnorm(cross(v2-v1, v3-v1, 2), 2, 2);

samples = struct();
samples.faultId = model(:, 1);
samples.layerId = nan(size(model, 1), 1);
samples.xKm = centers(:, 1)/1e3;
samples.yKm = centers(:, 2)/1e3;
samples.depthKm = abs(centers(:, 3))/1e3;
samples.slipCm = component_amplitude(model(:, 2), model(:, 3), settings.component);
samples.areaM2 = area;
sourceProduct = string(p.resampledTdeResult);
method = "tde_triangle_centroid_depth_profile";
end

function tdeModel = load_tde_postprocessing_model(p)
sourceInfo = dir(p.resampledTdeResult);
if isempty(sourceInfo)
    error("SlipSolve:MissingStageProduct", ...
        "Run resampled_tde_inversion before TDE shallow-deficit analysis: %s", ...
        p.resampledTdeResult);
end
if exist(p.tdePostprocessingModel, "file") == 2
    loaded = load(p.tdePostprocessingModel, "tdeModel");
    if isfield(loaded, "tdeModel") && ...
            get_field(loaded.tdeModel, "sourceBytes", -1) == sourceInfo.bytes && ...
            get_field(loaded.tdeModel, "sourceDatenum", -1) == sourceInfo.datenum
        tdeModel = loaded.tdeModel;
        return
    end
end

fprintf("Creating compact TDE post-processing model from the saved inversion.\n");
loaded = load(p.resampledTdeResult, "tdeResult");
tdeModel = struct();
tdeModel.slipModel = loaded.tdeResult.slipModel;
tdeModel.points = loaded.tdeResult.points;
tdeModel.method = loaded.tdeResult.method;
tdeModel.sourceProduct = string(p.resampledTdeResult);
tdeModel.sourceBytes = sourceInfo.bytes;
tdeModel.sourceDatenum = sourceInfo.datenum;
tdeModel.createdAt = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");
save(p.tdePostprocessingModel, "tdeModel", "-v7.3");
end

function slip = component_amplitude(strikeSlip, dipSlip, component)
switch lower(string(component))
    case {"magnitude", "total", "total_magnitude"}
        slip = hypot(strikeSlip, dipSlip);
    case {"strike", "strike_magnitude", "strike-slip"}
        slip = abs(strikeSlip);
    case {"dip", "dip_magnitude", "dip-slip"}
        slip = abs(dipSlip);
    otherwise
        error("SlipSolve:InvalidShallowDeficitConfig", ...
            "Unknown cfg.shallowSlipDeficit.component: %s", component);
end
end

function profiles = build_profiles(samples, segments, modelType, settings)
profiles = repmat(empty_profile(), numel(segments), 1);
for k = 1:numel(segments)
    segment = segments(k);
    mask = segment_mask(samples, segment);
    if ~any(mask)
        error("SlipSolve:EmptyShallowDeficitSegment", ...
            "Segment '%s' selects no %s model samples.", segment.name, modelType);
    end
    if modelType == "composite"
        profile = composite_profile(samples, mask, settings);
    else
        profile = tde_profile(samples, mask, settings);
    end
    profile.name = string(segment.name);
    profile.faultIds = segment.faultIds;
    profile.xRangeKm = segment.xRangeKm;
    profile.yRangeKm = segment.yRangeKm;
    profile = add_deficit_metrics(profile, settings);
    profiles(k) = profile;
    fprintf("  %s: shallowest deficit %.2f%%; mean deficit to %.2f km %.2f%%\n", ...
        profile.name, 100*profile.shallowestDeficitFraction, ...
        settings.shallowDepthMaxKm, 100*profile.meanShallowDeficitFraction);
end
end

function profile = composite_profile(samples, mask, settings)
averaging = lower(string(settings.composite.averaging));
if averaging ~= "arithmetic"
    error("SlipSolve:InvalidShallowDeficitConfig", ...
        "Composite legacy profiles currently require averaging='arithmetic'.");
end
layerIds = unique(samples.layerId(mask));
n = numel(layerIds);
depth = zeros(n, 1);
meanSlip = zeros(n, 1);
count = zeros(n, 1);
for i = 1:n
    rows = mask & samples.layerId == layerIds(i);
    depth(i) = mean(samples.depthKm(rows));
    meanSlip(i) = mean(samples.slipCm(rows));
    count(i) = nnz(rows);
end
[depth, order] = sort(depth);
profile = empty_profile();
profile.depthKm = depth;
profile.meanSlipCm = meanSlip(order);
profile.sampleCount = count(order);
profile.supportAreaM2 = nan(n, 1);
profile.supportWeight = count(order);
end

function profile = tde_profile(samples, mask, settings)
depth = samples.depthKm(mask);
slip = samples.slipCm(mask);
area = samples.areaM2(mask);
tde = settings.tde;
grouping = lower(string(tde.depthGrouping));
switch grouping
    case {"native_centroid_levels", "native", "centroid_levels"}
        tolerance = double(tde.depthToleranceKm);
        if ~isscalar(tolerance) || ~isfinite(tolerance) || tolerance <= 0
            error("SlipSolve:InvalidShallowDeficitConfig", ...
                "TDE depthToleranceKm must be a positive finite scalar.");
        end
        keys = round(depth/tolerance)*tolerance;
        [~, ~, groups] = unique(keys);
    case {"fixed_bins", "bins"}
        edges = double(get_field(tde, "binEdgesKm", []));
        if isempty(edges)
            width = double(tde.binWidthKm);
            if ~isscalar(width) || ~isfinite(width) || width <= 0
                error("SlipSolve:InvalidShallowDeficitConfig", ...
                    "TDE binWidthKm must be a positive finite scalar.");
            end
            edges = 0:width:(ceil(max(depth)/width)*width+width);
        end
        if numel(edges) < 2 || any(~isfinite(edges)) || any(diff(edges) <= 0)
            error("SlipSolve:InvalidShallowDeficitConfig", ...
                "TDE binEdgesKm must be finite and strictly increasing.");
        end
        groups = discretize(depth, edges);
    otherwise
        error("SlipSolve:InvalidShallowDeficitConfig", ...
            "Unknown TDE depthGrouping: %s", tde.depthGrouping);
end

valid = isfinite(groups);
groups = groups(valid);
depth = depth(valid);
slip = slip(valid);
area = area(valid);
groupIds = unique(groups);
n = numel(groupIds);
profileDepth = zeros(n, 1);
meanSlip = zeros(n, 1);
count = zeros(n, 1);
supportArea = zeros(n, 1);
supportWeight = zeros(n, 1);
averaging = lower(string(tde.averaging));
for i = 1:n
    rows = groups == groupIds(i);
    count(i) = nnz(rows);
    supportArea(i) = sum(area(rows));
    if averaging == "area_weighted"
        weights = area(rows);
    elseif averaging == "arithmetic"
        weights = ones(count(i), 1);
    else
        error("SlipSolve:InvalidShallowDeficitConfig", ...
            "TDE averaging must be 'area_weighted' or 'arithmetic'.");
    end
    supportWeight(i) = sum(weights);
    profileDepth(i) = sum(depth(rows).*weights)/supportWeight(i);
    meanSlip(i) = sum(slip(rows).*weights)/supportWeight(i);
end
[profileDepth, order] = sort(profileDepth);
profile = empty_profile();
profile.depthKm = profileDepth;
profile.meanSlipCm = meanSlip(order);
profile.sampleCount = count(order);
profile.supportAreaM2 = supportArea(order);
profile.supportWeight = supportWeight(order);
end

function profile = add_deficit_metrics(profile, settings)
referenceRange = double(settings.referenceDepthRangeKm(:).');
if numel(referenceRange) ~= 2 || referenceRange(1) > referenceRange(2)
    error("SlipSolve:InvalidShallowDeficitConfig", ...
        "referenceDepthRangeKm must be [minimum maximum].");
end
referenceRows = profile.depthKm >= referenceRange(1) & ...
    profile.depthKm <= referenceRange(2) & isfinite(profile.meanSlipCm);
if ~any(referenceRows)
    error("SlipSolve:InvalidShallowDeficitConfig", ...
        "The reference depth range contains no populated profile levels.");
end
referenceSlip = max(profile.meanSlipCm(referenceRows));
if ~isfinite(referenceSlip) || referenceSlip <= 0
    error("SlipSolve:InvalidShallowDeficitModel", ...
        "Reference profile slip must be positive and finite.");
end
profile.referenceSlipCm = referenceSlip;
profile.normalizedSlip = profile.meanSlipCm/referenceSlip;
profile.deficitFraction = 1-profile.normalizedSlip;
profile.shallowestDepthKm = profile.depthKm(1);
profile.shallowestSlipCm = profile.meanSlipCm(1);
profile.shallowestDeficitFraction = profile.deficitFraction(1);

shallow = profile.depthKm <= settings.shallowDepthMaxKm;
if any(shallow)
    weights = profile.supportWeight(shallow);
    profile.meanShallowSlipCm = ...
        sum(profile.meanSlipCm(shallow).*weights)/sum(weights);
    profile.meanShallowDeficitFraction = ...
        1-profile.meanShallowSlipCm/referenceSlip;
else
    profile.meanShallowSlipCm = NaN;
    profile.meanShallowDeficitFraction = NaN;
end
end

function mask = segment_mask(samples, segment)
mask = samples.xKm > segment.xRangeKm(1) & ...
    samples.xKm <= segment.xRangeKm(2) & ...
    samples.yKm > segment.yRangeKm(1) & ...
    samples.yKm <= segment.yRangeKm(2);
ids = segment.faultIds;
if ~isempty(ids) && ~(isstring(ids) && any(lower(ids) == "all"))
    mask = mask & ismember(samples.faultId, double(ids));
end
mask = mask & isfinite(samples.depthKm) & isfinite(samples.slipCm);
end

function segments = validate_segments(segments)
if isempty(segments) || ~isstruct(segments)
    error("SlipSolve:InvalidShallowDeficitConfig", ...
        "cfg.shallowSlipDeficit.segments must be a nonempty struct array.");
end
required = ["name", "faultIds", "xRangeKm", "yRangeKm"];
for k = 1:numel(segments)
    for field = required
        if ~isfield(segments(k), field)
            error("SlipSolve:InvalidShallowDeficitConfig", ...
                "Shallow-deficit segment %d lacks field %s.", k, field);
        end
    end
    if strlength(string(segments(k).name)) == 0
        error("SlipSolve:InvalidShallowDeficitConfig", ...
            "Every shallow-deficit segment needs a nonempty name.");
    end
    for field = ["xRangeKm", "yRangeKm"]
        range = double(segments(k).(field));
        if numel(range) ~= 2 || any(isnan(range)) || range(1) > range(2)
            error("SlipSolve:InvalidShallowDeficitConfig", ...
                "Segment %s must define ordered two-value %s.", ...
                segments(k).name, field);
        end
        segments(k).(field) = range(:).';
    end
end
segments = segments(:);
end

function [profileTable, summaryTable] = result_tables(modelType, profiles, settings)
modelColumn = strings(0, 1);
componentColumn = strings(0, 1);
segmentColumn = strings(0, 1);
depth = zeros(0, 1);
meanSlip = zeros(0, 1);
normalized = zeros(0, 1);
deficit = zeros(0, 1);
sampleCount = zeros(0, 1);
supportArea = zeros(0, 1);
for k = 1:numel(profiles)
    n = numel(profiles(k).depthKm);
    modelColumn = [modelColumn; repmat(modelType, n, 1)]; %#ok<AGROW>
    componentColumn = [componentColumn; repmat(string(settings.component), n, 1)]; %#ok<AGROW>
    segmentColumn = [segmentColumn; repmat(profiles(k).name, n, 1)]; %#ok<AGROW>
    depth = [depth; profiles(k).depthKm]; %#ok<AGROW>
    meanSlip = [meanSlip; profiles(k).meanSlipCm]; %#ok<AGROW>
    normalized = [normalized; profiles(k).normalizedSlip]; %#ok<AGROW>
    deficit = [deficit; profiles(k).deficitFraction]; %#ok<AGROW>
    sampleCount = [sampleCount; profiles(k).sampleCount]; %#ok<AGROW>
    supportArea = [supportArea; profiles(k).supportAreaM2]; %#ok<AGROW>
end
profileTable = table(modelColumn, componentColumn, segmentColumn, depth, meanSlip, ...
    normalized, deficit, sampleCount, supportArea, ...
    'VariableNames', {'model_type', 'component', 'segment', 'depth_km', 'mean_slip_cm', ...
    'normalized_slip', 'deficit_fraction', 'sample_count', 'support_area_m2'});

n = numel(profiles);
referenceRange = double(settings.referenceDepthRangeKm(:).');
summaryTable = table(repmat(modelType, n, 1), ...
    repmat(string(settings.component), n, 1), reshape([profiles.name], [], 1), ...
    repmat(referenceRange(1), n, 1), repmat(referenceRange(2), n, 1), ...
    repmat(double(settings.shallowDepthMaxKm), n, 1), ...
    reshape([profiles.referenceSlipCm], [], 1), ...
    reshape([profiles.shallowestDepthKm], [], 1), ...
    reshape([profiles.shallowestSlipCm], [], 1), ...
    reshape([profiles.shallowestDeficitFraction], [], 1), ...
    reshape([profiles.meanShallowSlipCm], [], 1), ...
    reshape([profiles.meanShallowDeficitFraction], [], 1), ...
    'VariableNames', {'model_type', 'component', 'segment', ...
    'reference_depth_min_km', 'reference_depth_max_km', ...
    'shallow_depth_max_km', 'reference_slip_cm', ...
    'shallowest_depth_km', 'shallowest_slip_cm', ...
    'shallowest_deficit_fraction', 'mean_shallow_slip_cm', ...
    'mean_shallow_deficit_fraction'});
end

function modelTypes = normalize_model_types(values)
values = lower(string(values(:)));
if any(values == "both")
    modelTypes = ["composite"; "tde"];
else
    aliases = values;
    aliases(ismember(aliases, ["final", "layered", "nodal"])) = "composite";
    aliases(ismember(aliases, ["triangular", "homogeneous_tde"])) = "tde";
    if any(~ismember(aliases, ["composite", "tde"]))
        error("SlipSolve:InvalidShallowDeficitConfig", ...
            "modelTypes must contain 'composite', 'tde', or 'both'.");
    end
    modelTypes = unique(aliases, "stable");
end
end

function profile = empty_profile()
profile = struct("name", "", "faultIds", [], "xRangeKm", [], ...
    "yRangeKm", [], "depthKm", [], "meanSlipCm", [], ...
    "normalizedSlip", [], "deficitFraction", [], "sampleCount", [], ...
    "supportAreaM2", [], "supportWeight", [], "referenceSlipCm", NaN, ...
    "shallowestDepthKm", NaN, "shallowestSlipCm", NaN, ...
    "shallowestDeficitFraction", NaN, "meanShallowSlipCm", NaN, ...
    "meanShallowDeficitFraction", NaN);
end

function value = get_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
