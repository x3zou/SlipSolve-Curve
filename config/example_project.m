function cfg = example_project()
%EXAMPLE_PROJECT User-editable SlipSolve configuration.
%
% This file is the main control panel for the workflow. For the bundled
% precomputed Myanmar example, edit parameters here and then run:
%
%   run examples/myanmar_demo/run_example.m
%
% That example runs ONLY final_inversion. It reuses the project-local mesh,
% seven sampled data sets, full G_e/G_n/G_u matrices, and sampled-Green cache;
% it does not rerun sampling, geometry, TDE, EDGRN, or composite-Green stages.
%
% All MATLAB dependencies used by the workflow are project-local under
% external/legacy_matlab. Generated products stay under cfg.project.outputRoot.

cfg = struct();

%% ------------------------------------------------------------------------
%  1. Project name and output folder
%  ------------------------------------------------------------------------
% cfg.project.name is used in metadata and logs.
cfg.project.name = "myanmar_demo";

% cfg.project.root is the SlipSolve project folder that contains this config.
cfg.project.root = fileparts(fileparts(mfilename("fullpath")));

% cfg.project.outputRoot controls where all generated products are written.
% Keep this inside the SlipSolve project unless you intentionally want a
% separate output directory.
cfg.project.outputRoot = cfg.project.root;

%% ------------------------------------------------------------------------
%  2. Workflow start/stop controls
%  ------------------------------------------------------------------------
% Choose which stage to start from and where to stop. This is usually the
% first thing to change when tuning parameters.
%
% Valid step names:
%   "subsample_insar"
%   "fault_geometry"
%   "triangular_mesh"
%   "quick_inversion"
%   "model_based_sampling"
%   "resampled_tde_inversion"
%   "layered_greens"
%   "composite_greens"
%   "final_inversion"
%
% Bundled-example default: run only the final altered-node layered inversion
% from the copied full-resolution XYZ Greens and seven copied samp3 files.
% Change these values only when intentionally running another workflow stage.
cfg.workflow.startStep = "final_inversion";
cfg.workflow.stopStep = "final_inversion";

% In MATLAB desktop, pause after each completed stage so users can inspect
% figures before continuing.
cfg.workflow.pauseAfterStage = true;

%% ------------------------------------------------------------------------
%  3. Legacy code and input-data paths
%  ------------------------------------------------------------------------
% Internalized legacy numerical-code tree. Keep these relationships intact;
% several original routines locate sibling directories at runtime.
cfg.paths.legacyPackageRoot = fullfile(cfg.project.root, ...
    "external", "legacy_matlab");
cfg.paths.legacyRoot = fullfile(cfg.paths.legacyPackageRoot, "Inversion");
cfg.paths.legacyCurvedGeometryRoot = fullfile(cfg.paths.legacyRoot, ...
    "curved_geometry");
cfg.paths.legacySamplingRoot = fullfile(cfg.paths.legacyPackageRoot, ...
    "geodetic_inversion-master", "sampling");
cfg.paths.legacyOtherFuncRoot = fullfile(cfg.paths.legacyPackageRoot, ...
    "geodetic_inversion-master", "OtherFunc");
cfg.paths.legacyGeodeticFunctionsRoot = fullfile(cfg.paths.legacyPackageRoot, ...
    "geodetic_data");

% Google Drive example-data bundle destination. Users copy the bundle's
% data/ and greens/ folders into the SlipSolve-curve project root.
cfg.paths.exampleGridRoot = fullfile(cfg.project.root, ...
    "data", "raw", "example_grids");
cfg.paths.geodeticRoot = cfg.paths.exampleGridRoot;

% Optional raw InSAR path. Leave empty when tracks are listed explicitly in
% cfg.insar.tracks below.
cfg.paths.rawInsar = "";

% Legacy sample-list path. This is mainly for compatibility with old
% workflows; explicit cfg.insar.tracks entries are preferred for new runs.
cfg.paths.sampleListFile = "";

% Optional single fault-trace alias. The geometry section below is the
% preferred place to define trace files.
cfg.paths.faultTrace = "";

% Fault polyline used by fault-aware quadtree sign cleaning near the trace.
cfg.paths.meshTrace = fullfile(cfg.paths.exampleGridRoot, "mesh_trace.txt");

% EDGRN executable. SlipSolve does not run EDGRN itself: users run it
% externally, then the layered_greens stage converts its completed text output
% into project-local MAT files. This path is retained as a documented record of
% the executable used to create a database.
cfg.paths.edgrnExecutable = ""; % user-installed EDGRN executable, if needed

%% ------------------------------------------------------------------------
%  4. InSAR coordinate convention and fallback quadtree settings
%  ------------------------------------------------------------------------
% Coordinate system and units for sampled observations passed into inversion.
cfg.insar.coordinateSystem = "utm";
cfg.insar.lengthUnit = "m";
cfg.insar.losConvention = "positive_toward_satellite";

% Geographic origin used by ll2xy. Replace all three values for the user's
% own study area. lon/lat define local (0,0), while refLon is the projection
% reference longitude/central meridian expected by the legacy ll2xy routine.
% Every geographic grid, fault trace, saved mesh, and composite Green grid in
% one project must use this same inversion origin and reference longitude.
cfg.insar.origin.lon = 96.05;
cfg.insar.origin.lat = 20.75;
cfg.insar.origin.refLon = 96.05;

% Fallback quadtree parameters used by old configs or legacy sample-list
% tracks that do not define track.initialSampling/modelBasedSampling blocks.
% For new projects, edit each track's initialSampling/modelBasedSampling
% blocks below instead of these fallback values.
cfg.insar.quadtree.minBlockSize = 500;
cfg.insar.quadtree.maxBlockSize = 20000;
cfg.insar.quadtree.varianceThreshold = 0.0025;
cfg.insar.quadtree.maxSamples = 5000;
cfg.insar.quadtree.method = "quadtree";
cfg.insar.quadtree.nanFractionMax = 1;
cfg.insar.quadtree.statistic = "mean";

% Fallback fault-aware options for quadtree_unstructured2.
cfg.insar.quadtree.faultAware = true;
cfg.insar.quadtree.faultTolerance = 10;
cfg.insar.quadtree.minPixelsPerSide = 2;

% "native" means the configured grid files are already the grids to sample.
% Other projects may use "low" or "high" naming conventions.
cfg.insar.inputResolution = "native";

%% ------------------------------------------------------------------------
%  5. InSAR tracks to use
%  ------------------------------------------------------------------------
% Each track entry defines:
%   name          - short label used in output folders and figure names
%   dataType      - "azo", "los", or "rng"
%   dataFile      - displacement/offset grid to sample
%   look*File     - look-vector grids aligned with dataFile
%   initialSampling     - first quadtree pass on the observed data grid
%   modelBasedSampling  - later quadtree pass on the quick model.grd
%
% Sampling parameters are PER TRACK. If different InSAR tracks need
% different density/thresholds, set different values in each track entry.
cfg.insar.tracks = struct( ...
    "name", "SEN_A70_azo_low", ...
    "dataType", "azo", ...
    "dataFile", fullfile(cfg.paths.exampleGridRoot, "SEN", "A70", "azo_ll_low.grd"), ...
    "lookEFile", fullfile(cfg.paths.exampleGridRoot, "SEN", "A70", "look_e_azo_low.grd"), ...
    "lookNFile", fullfile(cfg.paths.exampleGridRoot, "SEN", "A70", "look_n_azo_low.grd"), ...
    "lookUFile", fullfile(cfg.paths.exampleGridRoot, "SEN", "A70", "look_u_azo_low.grd"), ...
    "initialSampling", struct( ...
        "boundsLonLat", [], ...          % [minLon maxLon minLat maxLat], or [] for full grid
        "minCellSizeKm", 0.5, ...        % smallest initial data-based quadtree cell
        "maxCellSizeKm", 20, ...         % largest initial data-based quadtree cell
        "rmsThreshold", 0.0025, ...      % data variability threshold for initial sampling
        "nanFractionMax", 1, ...         % maximum allowed NaN fraction per cell
        "padX", 0, ...                   % legacy NaN padding in grid columns
        "padY", 0, ...                   % legacy NaN padding in grid rows
        "faultToleranceKm", 0.01, ...    % near-fault sign-cleaning tolerance
        "minPixelsPerSide", 2, ...       % minimum pixels on chosen side near fault
        "statistic", "mean"), ...        % "mean" or "median"
    "modelBasedSampling", struct( ...
        "boundsLonLat", [], ...          % usually same area as initial sampling; [] uses full grid
        "minCellSizeKm", 0.5, ...          % smallest model.grd-based quadtree cell
        "maxCellSizeKm", 20, ...        % largest model.grd-based quadtree cell
        "rmsThreshold", 1, ...          % model variability threshold for model-based sampling
        "nanFractionMax", 1, ...         % maximum allowed NaN fraction per cell
        "padX", 0, ...                   % legacy NaN padding in grid columns
        "padY", 0, ...                   % legacy NaN padding in grid rows
        "faultToleranceKm", 0.01, ...    % near-fault sign-cleaning tolerance
        "minPixelsPerSide", 2, ...       % minimum pixels on chosen side near fault
        "statistic", "mean"));           % "mean" or "median"

% To add more tracks, duplicate the first struct entry, for example:
% cfg.insar.tracks(2) = cfg.insar.tracks(1);
% cfg.insar.tracks(2).name = "another_track";
% cfg.insar.tracks(2).dataFile = "/path/to/another_grid.grd";
% cfg.insar.tracks(2).initialSampling.minCellSizeKm = 1;
% cfg.insar.tracks(2).modelBasedSampling.rmsThreshold = 25;

%% ------------------------------------------------------------------------
%  6. Visualization behavior
%  ------------------------------------------------------------------------
% When true, stages create figures and save .fig/.png files.
cfg.visualization.enabled = true;

% Save interactive MATLAB figures for later inspection/tweaking.
cfg.visualization.saveFig = true;

% Save static PNG previews next to the .fig files.
cfg.visualization.savePng = true;

% When true in MATLAB desktop, figures pop up during the run.
cfg.visualization.visible = true;

% Stage-specific figure options.
cfg.visualization.insarSubsampling.showFaultTrace = true;
cfg.visualization.insarSubsampling.showQuadtreeCells = true;
cfg.visualization.insarSubsampling.showRmsHistogram = true;
cfg.visualization.insarSubsampling.colormap = "turbo";
cfg.visualization.insarSubsampling.colorbarRange = [];
cfg.visualization.insarSubsampling.axisRange = []; % [xmin xmax ymin ymax] km
cfg.visualization.insarSubsampling.showTitles = true;
cfg.visualization.insarSubsampling.rawTitle = "Raw cropped data";
cfg.visualization.insarSubsampling.sampledTitle = "Quadtree samples";
cfg.visualization.insarSubsampling.histogramTitle = ""; % empty uses track name

cfg.visualization.faultGeometry.showControlLines = true;
cfg.visualization.faultGeometry.showSurfaceGrid = true;
cfg.visualization.faultGeometry.colormap = "lines";
cfg.visualization.faultGeometry.axisRange = []; % [xmin xmax ymin ymax zmin zmax] km
cfg.visualization.faultGeometry.view = [35 25];
cfg.visualization.faultGeometry.showTitle = true;
cfg.visualization.faultGeometry.title = "Fault Geometry Preview";

cfg.visualization.triangularMesh.colorBy = "depth";
cfg.visualization.triangularMesh.showFaultTrace = true;
cfg.visualization.triangularMesh.colormap = "parula";
cfg.visualization.triangularMesh.colorbarRange = [];
cfg.visualization.triangularMesh.axisRange = []; % [xmin xmax ymin ymax zmin zmax] km
cfg.visualization.triangularMesh.view = [35 25];
cfg.visualization.triangularMesh.showTitle = true;
cfg.visualization.triangularMesh.title = ""; % empty uses node/triangle counts

cfg.visualization.quickInversion.colormap = "parula";
cfg.visualization.quickInversion.colorbarRange = [];
cfg.visualization.quickInversion.axisRange = []; % [xmin xmax ymin ymax zmin zmax] km
cfg.visualization.quickInversion.view = [35 25];
cfg.visualization.quickInversion.showTitle = true;
cfg.visualization.quickInversion.title = ""; % empty uses variance reduction

cfg.visualization.modelSampling.colormap = "turbo";
cfg.visualization.modelSampling.colorbarRange = [];
cfg.visualization.modelSampling.axisRange = []; % [xmin xmax ymin ymax] km
cfg.visualization.modelSampling.showTitles = true;
cfg.visualization.modelSampling.beforeTitle = "Data before sampling";
cfg.visualization.modelSampling.afterTitle = "Data after model-based sampling";

% Resampled TDE inversion plot controls. Empty [] ranges mean automatic.
% Axis ranges are in km:
%   slip/polarity axisRange = [xmin xmax ymin ymax zmin zmax]
%   fit axisRange           = [xmin xmax ymin ymax]
cfg.visualization.resampledTde.slipPlot.colormap = "parula";
cfg.visualization.resampledTde.slipPlot.colorbarRange = [];     % e.g., [0 600]
cfg.visualization.resampledTde.slipPlot.axisRange = [];         % e.g., [-80 80 -320 230 -30 2]
cfg.visualization.resampledTde.slipPlot.showTitle = true;       % false gives no title
cfg.visualization.resampledTde.slipPlot.title = "";             % "" uses automatic title
cfg.visualization.resampledTde.slipPlot.view = [35 25];

cfg.visualization.resampledTde.fitPlots.colormap = "turbo";
cfg.visualization.resampledTde.fitPlots.dataModelColorbarRange = []; % e.g., [-300 300]
cfg.visualization.resampledTde.fitPlots.residualColorbarRange = [];  % e.g., [-80 80]
cfg.visualization.resampledTde.fitPlots.axisRange = [];              % e.g., [-140 140 -320 240]
cfg.visualization.resampledTde.fitPlots.showTitles = true;           % false removes Data/Model/Residual titles
cfg.visualization.resampledTde.fitPlots.panelTitles = ["Data" "Model" "Residual"];

cfg.visualization.resampledTde.polarityPlot.enabled = true;
cfg.visualization.resampledTde.polarityPlot.colormap = "turbo";
cfg.visualization.resampledTde.polarityPlot.strikeColorbarRange = []; % e.g., [-600 600]
cfg.visualization.resampledTde.polarityPlot.dipColorbarRange = [];    % e.g., [-100 100]
cfg.visualization.resampledTde.polarityPlot.axisRange = [];           % e.g., [-80 80 -320 230 -30 2]
cfg.visualization.resampledTde.polarityPlot.showTitles = true;
cfg.visualization.resampledTde.polarityPlot.strikeTitle = "Strike-slip component (cm)";
cfg.visualization.resampledTde.polarityPlot.dipTitle = "Dip-slip component (cm)";
cfg.visualization.resampledTde.polarityPlot.view = [35 25];

%% ------------------------------------------------------------------------
%  7. Fault geometry controls
%  ------------------------------------------------------------------------
% Geometry mode that follows the legacy fit_curved_planeMM.m recipe.
cfg.geometry.mode = "legacy_fit_curved_planeMM";

% Surface fault trace used to construct the curved fault surface.
cfg.geometry.surfaceTraceFiles = fullfile(cfg.paths.exampleGridRoot, "trace.S2");

% Optional explicit depth-control files. Leave empty when using legacy
% segment files plus segment dip angles.
cfg.geometry.depthControlFiles = strings(0, 1);

% Approximate dip used by non-legacy/simple geometry modes.
cfg.geometry.uniformDipDegrees = 60;

% Legacy segment files at depth. Each segment should have a corresponding
% dip value in cfg.geometry.segmentDipDegrees.
cfg.geometry.segmentFiles = [
    fullfile(cfg.project.root, "data", "raw", "example_geometry", "Segment_001.txt")
    fullfile(cfg.project.root, "data", "raw", "example_geometry", "Segment_002.txt")
    fullfile(cfg.project.root, "data", "raw", "example_geometry", "Segment_003.txt")
    fullfile(cfg.project.root, "data", "raw", "example_geometry", "Segment_004.txt")
    fullfile(cfg.project.root, "data", "raw", "example_geometry", "Segment_005.txt")
    fullfile(cfg.project.root, "data", "raw", "example_geometry", "Segment_006.txt")
    fullfile(cfg.project.root, "data", "raw", "example_geometry", "Segment_007.txt")
];
cfg.geometry.segmentDipDegrees = [75 75 70 80 85 90 100];

% Geometry dimensions and smoothing controls. Depth/spacing values are in m.
cfg.geometry.maxDepth = 25000;
cfg.geometry.traceResampleSpacing = 2000;
cfg.geometry.downDipLevels = 21;
cfg.geometry.surfaceFitSmoothness = 0.008;

% Legacy-exact parameters copied from the old curved-geometry scripts.
cfg.geometry.legacyExact.enabled = true;
cfg.geometry.legacyExact.miniKm = 5;
cfg.geometry.legacyExact.depthKm = 25;
cfg.geometry.legacyExact.depthLevels = 30;
cfg.geometry.legacyExact.surfaceInterpolationFactor = 1e5;
cfg.geometry.legacyExact.bottomInterpolationFactor = 1e2;

% Reference products are for QA/debugging only. Keep disabled for normal
% computed geometry runs.
cfg.geometry.legacyExact.useReferenceProducts = false;
cfg.geometry.legacyExact.referenceMeshFile = "";
cfg.geometry.legacyExact.referencePointFile = "";

%% ------------------------------------------------------------------------
%  8. Fault list
%  ------------------------------------------------------------------------
% The workflow supports multiple disconnected faults. Each fault is fit and
% meshed independently, then combined downstream.
cfg.geometry.faults = struct( ...
    "name", "main_fault", ...
    "surfaceTraceFiles", cfg.geometry.surfaceTraceFiles, ...
    "mode", cfg.geometry.mode, ...
    "segmentFiles", cfg.geometry.segmentFiles, ...
    "segmentDipDegrees", cfg.geometry.segmentDipDegrees, ...
    "uniformDipDegrees", cfg.geometry.uniformDipDegrees, ...
    "maxDepth", cfg.geometry.maxDepth, ...
    "traceResampleSpacing", cfg.geometry.traceResampleSpacing, ...
    "downDipLevels", cfg.geometry.downDipLevels, ...
    "surfaceFitSmoothness", cfg.geometry.surfaceFitSmoothness);

% Example for a disconnected second fault:
% cfg.geometry.faults(2) = cfg.geometry.faults(1);
% cfg.geometry.faults(2).name = "second_fault";
% cfg.geometry.faults(2).surfaceTraceFiles = "/path/to/second_trace.txt";
% cfg.geometry.faults(2).segmentFiles = ["/path/to/second_segment_001.txt"];
% cfg.geometry.faults(2).segmentDipDegrees = 70;

%% ------------------------------------------------------------------------
%  9. Triangular mesh controls
%  ------------------------------------------------------------------------
% Mesh method that follows the legacy main_interpolate.m Delaunay routine.
cfg.mesh.method = "legacy_exact_main_interpolate";

% Choose which configured faults to mesh: "all", numeric indices, or names.
cfg.mesh.activeFaults = "all";

% Legacy-style mesh shape controls. biasW > 1 gives depth-dependent spacing.
cfg.mesh.biasL = 1;
cfg.mesh.biasW = 1.15;
cfg.mesh.depthRatio = 1.5;
cfg.mesh.depthKm = cfg.geometry.maxDepth ./ 1e3;
cfg.mesh.interpolationMethod = "cubic";

% General mesh controls kept for future/non-legacy mesh methods.
cfg.mesh.targetEdgeLength = 3000;
cfg.mesh.minPatchArea = 1.0e6;
cfg.mesh.maxPatchArea = 2.5e7;
cfg.mesh.diagonalPattern = "alternating";

% Reference geometry is for QA/debugging only. Keep disabled for normal
% computed mesh runs.
cfg.mesh.legacyExact.useReferenceProducts = false;
cfg.mesh.legacyExact.referenceGeometryFile = "";

%% ------------------------------------------------------------------------
%  10. Quick homogeneous TDE inversion controls
%  ------------------------------------------------------------------------
% Green's-function type and slip parameterization for the quick inversion.
cfg.quickInversion.greenFunction = "homogeneous_tde";
cfg.quickInversion.rakeMode = "fixed";
cfg.quickInversion.rakeDegrees = 90;
cfg.quickInversion.modelType = "okada_curve";

% Regularization and boundary constraints. These follow legacy inversion
% settings and are intentionally exposed for tuning.
cfg.quickInversion.smoothingWeight = 3.3333e-1;
cfg.quickInversion.smoothingDipRatio = 3;
cfg.quickInversion.boundaryConstraintRatio = 5e-4;

% Slip sign constraints. The current default follows the legacy con vector.
cfg.quickInversion.signConstraint = [-1 0 0];
cfg.quickInversion.nonNegativeSlip = false;

% Per-data-type weights in the quick inversion.
cfg.quickInversion.trackWeights.azo = 1;
cfg.quickInversion.trackWeights.los = 1;
cfg.quickInversion.trackWeights.rng = 1;

%% ------------------------------------------------------------------------
%  11. Quick model grid output
%  ------------------------------------------------------------------------
% After the quick inversion, predict displacement back onto the full InSAR
% grid and save model.grd under the SlipSolve project.
cfg.quickInversion.modelGrid.enabled = true;
cfg.quickInversion.modelGrid.source = "insar_tracks";
cfg.quickInversion.modelGrid.chunkSize = 15000;
cfg.quickInversion.modelGrid.fileName = "model.grd";

% Keep false. True would write beside source data, which is not allowed for
% the MyanmarEQ legacy research tree.
cfg.quickInversion.modelGrid.writeToSourceDirectory = false;

%% ------------------------------------------------------------------------
%  12. Model-based sampling stage controls
%  ------------------------------------------------------------------------
% This stage builds a quadtree from model.grd, applies that tree to the real
% data and look vectors, and writes *_samp<iteration>.mat under SlipSolve.
% Numerical quadtree controls for this stage are set per track in:
%   cfg.insar.tracks(k).modelBasedSampling
cfg.quickInversion.modelSampling.enabled = true;
cfg.quickInversion.modelSampling.source = "insar_tracks";

% Output suffix. iterationStep = 99 gives azo_samp99.mat for AZO tracks.
cfg.quickInversion.modelSampling.iterationStep = 99;

% If true, try to use *_high.grd look-vector files when they exist.
cfg.quickInversion.modelSampling.highResolutionLooks = false;

% Keep false for the same read-only reason as modelGrid.writeToSourceDirectory.
cfg.quickInversion.modelSampling.writeToSourceDirectory = false;

%% ------------------------------------------------------------------------
%  13. Resampled-data TDE inversion controls
%  ------------------------------------------------------------------------
% This stage follows the MMInversionTri / make_fault_from_insar_curve setup:
% TDE Greens, smoothingMatrix_test, zero-slip boundary constraints,
% bounds_new polarity constraints, lsqlin.
cfg.resampledTdeInversion.enabled = true;
cfg.resampledTdeInversion.modelType = "okada_curve";

% Input mode for this TDE stage:
%   "legacy_samp3_reference" - use copied legacy samp3 files below; this is
%                              the bit-for-bit reference input set from
%                              make_fault_from_insar_curve.m
%   "model_based_sampling"   - use the *_samp*.mat files produced by the
%                              model_based_sampling stage
cfg.resampledTdeInversion.inputMode = "legacy_samp3_reference";

% Copied reference data root. These files live inside SlipSolve so the legacy
% research tree remains read-only during all new workflow runs.
cfg.resampledTdeInversion.legacySamplesRoot = fullfile(cfg.project.outputRoot, ...
    "data", "raw", "legacy_samp3");

% Exact legacy sample files, in the same order as make_fault_from_insar_curve:
% four Sentinel AZO tracks, then three ALOS2 LOS tracks.
cfg.resampledTdeInversion.legacySamples = [
    struct("name", "SEN_A70_azo_samp3",  "type", "azo", "greenOption", "AZO",   "relativeFile", fullfile("SEN", "A70", "azo_samp3.mat"))
    struct("name", "SEN_A143_azo_samp3", "type", "azo", "greenOption", "AZO",   "relativeFile", fullfile("SEN", "A143", "azo_samp3.mat"))
    struct("name", "SEN_D33_azo_samp3",  "type", "azo", "greenOption", "AZO",   "relativeFile", fullfile("SEN", "D33", "azo_samp3.mat"))
    struct("name", "SEN_D106_azo_samp3", "type", "azo", "greenOption", "AZO",   "relativeFile", fullfile("SEN", "D106", "azo_samp3.mat"))
    struct("name", "ALOS2_A152_los_samp3", "type", "los", "greenOption", "insar", "relativeFile", fullfile("ALOS2", "A152", "los_samp3.mat"))
    struct("name", "ALOS2_D41_los_samp3",  "type", "los", "greenOption", "insar", "relativeFile", fullfile("ALOS2", "D41", "los_samp3.mat"))
    struct("name", "ALOS2_D42_los_samp3",  "type", "los", "greenOption", "insar", "relativeFile", fullfile("ALOS2", "D42", "los_samp3.mat"))
];

% Legacy MMInversionTri smoothness is smoothness = 10e-1.
cfg.resampledTdeInversion.smoothingWeight = 10e-1;
cfg.resampledTdeInversion.smoothingStrikeRatio = 1;
cfg.resampledTdeInversion.smoothingDipRatio = 3;

% Boundary conditions. MMInversionTri constrains bottom, left, and right,
% and allows slip at the surface/top. The legacy zero_slip_boundary_curve
% helper does not define a top/surface option, so keep top.enabled=false for
% legacy-identical runs.
cfg.resampledTdeInversion.boundaryConditions.bottom.enabled = true;
cfg.resampledTdeInversion.boundaryConditions.bottom.ratio = 5e-4;
cfg.resampledTdeInversion.boundaryConditions.left.enabled = true;
cfg.resampledTdeInversion.boundaryConditions.left.ratio = 5e-4;
cfg.resampledTdeInversion.boundaryConditions.right.enabled = true;
cfg.resampledTdeInversion.boundaryConditions.right.ratio = 5e-4;
cfg.resampledTdeInversion.boundaryConditions.top.enabled = false;
cfg.resampledTdeInversion.boundaryConditions.top.ratio = 5e-4;

% Polarity/sign constraints passed to bounds_new. This follows the legacy
% con = [-1 0 0] setting in MMInversionTri.
cfg.resampledTdeInversion.constraints.polarity = [-1 0 0];
cfg.resampledTdeInversion.constraints.nonNegativeSlip = false;

% Optional direct bound overrides. Leave empty to use bounds_new.
cfg.resampledTdeInversion.constraints.customLowerBounds = [];
cfg.resampledTdeInversion.constraints.customUpperBounds = [];

% Relative data weights, matching alpha/beta-style legacy weighting.
cfg.resampledTdeInversion.trackWeights.azo = 1;
cfg.resampledTdeInversion.trackWeights.los = 1;
cfg.resampledTdeInversion.trackWeights.rng = 1;

%% ------------------------------------------------------------------------
%  14. Reserved experimental resampling controls
%  ------------------------------------------------------------------------
% These fields are retained for future experiments but are not currently read
% by model_based_sampling or final_inversion. Per-track modelBasedSampling
% blocks above are the active resampling controls.
cfg.resampling.enabled = true;
cfg.resampling.targetSamples = 12000;
cfg.resampling.modelGradientWeight = 0.7;
cfg.resampling.dataVarianceWeight = 0.3;

%% ------------------------------------------------------------------------
%  15. Layered Green's-function controls (external EDGRN, MATLAB conversion)
%  ------------------------------------------------------------------------
% EDGRN itself is run by the user outside MATLAB. Edit the EDGRN input file
% (receiver depth, distance/depth grids, output directory, and velocity model),
% run the executable, then set edgrnInputFile below and run only
% "layered_greens". See README.md for installation and execution instructions.
%
% The input file's fifth non-comment line names the EDGRN output directory and
% three text files (SS, DS, CL). Those files must exist before this stage runs.
% Paths on that line are interpreted relative to edgrnInputFile, exactly as in
% the legacy getedgrn.m method.
cfg.layeredGreens.enabled = true;
cfg.layeredGreens.edgrnInputFile = fullfile(cfg.project.root, ...
    "data", "raw", "edgrn", "example_model", "edgrnMM");

% Converted MAT files are written here, never next to the external EDGRN input
% or its text database. Use a separate folder for each Earth model/database.
cfg.layeredGreens.outputDirectory = fullfile(cfg.project.outputRoot, "greens", "edgrn", "edgrnMM");
cfg.layeredGreens.overwriteExistingMatFiles = false;

% Reference values only. They are encoded in edgrnInputFile and recorded in
% the converted MAT files; SlipSolve does not use these values to run EDGRN.
cfg.layeredGreens.receiverDepth = 0;
cfg.layeredGreens.sourceDepthSpacing = 1000;
cfg.layeredGreens.radialDistanceSpacing = 1000;
cfg.layeredGreens.maxRadialDistance = 500000;

% Diagnostic plot of log10(abs(uz)) for SS, DS, and CL EDGRN databases.
cfg.visualization.layeredGreens.enabled = true;
cfg.visualization.layeredGreens.colormap = "turbo";
cfg.visualization.layeredGreens.log10AbsoluteUz = true;
cfg.visualization.layeredGreens.showTitles = true;
cfg.visualization.layeredGreens.panelTitles = strings(0, 1); % empty uses SS/DS/CL names
cfg.visualization.layeredGreens.showOverallTitle = false;
cfg.visualization.layeredGreens.title = "Layered EDGRN Green's functions";
cfg.visualization.layeredGreens.axisRange = []; % [rMin rMax depthMin depthMax] km
cfg.visualization.layeredGreens.colorbarRange = [];

%% ------------------------------------------------------------------------
%  16. Composite Green's-function controls
%  ------------------------------------------------------------------------
% Build full-resolution east, north, and vertical layered Green's functions
% using the legacy Comb_Green_Layered_Nodes_xyz.m workflow. This is normally
% the longest-running stage in the workflow.
cfg.compositeGreens.enabled = true;

% Fault trace used to place a fine rectilinear grid near the rupture. The file
% must contain [x y] in meters, in the same local coordinate system as mesh.
cfg.compositeGreens.meshTraceFile = cfg.paths.meshTrace;

% Full-resolution grid controls copied from Comb_Green_Layered_Nodes_xyz.m.
% All five values are in meters. D1 is fine spacing inside the fault bounding
% box expanded by M; D2 is coarser spacing outside it and must exceed D1.
cfg.compositeGreens.xRange = [-240e3 240e3];
cfg.compositeGreens.yRange = [-350e3 290e3];
cfg.compositeGreens.D1 = 1.3e3;
cfg.compositeGreens.D2 = 2.3e3;
cfg.compositeGreens.M = 9e3;

% The legacy node routine treats the top two node-depth levels specially and
% samples each connected triangle with 15 point sources. Keep these defaults
% for reproduction of Comb_Green_Layered_Nodes_xyz.m.
cfg.compositeGreens.skipTopDepthLevels = 2;
cfg.compositeGreens.pointSourcesPerTriangle = 15;

% Full G matrices can be several GB each, so they are saved separately with
% MATLAB -v7.3. Existing files are protected unless overwrite is enabled.
cfg.compositeGreens.outputDirectory = fullfile(cfg.project.outputRoot, ...
    "greens", "composite_full_xyz");
cfg.compositeGreens.overwriteExistingFiles = false;

% Show and optionally pause at the grid preview before starting the expensive
% Green's-function loop. The pause is used only in interactive MATLAB desktop.
cfg.compositeGreens.pauseAfterGridPreview = true;

cfg.visualization.compositeGreens.enabled = true;
cfg.visualization.compositeGreens.pointStride = 10;
cfg.visualization.compositeGreens.pointColor = [0.72 0.76 0.80];
cfg.visualization.compositeGreens.faultColor = [0.75 0.10 0.08];
cfg.visualization.compositeGreens.axisRange = []; % [xmin xmax ymin ymax] km, [] for auto
cfg.visualization.compositeGreens.showTitle = true;
cfg.visualization.compositeGreens.title = "Full-resolution composite Green's grid";

%% ------------------------------------------------------------------------
%  17. Final inversion controls
%  ------------------------------------------------------------------------
% This stage reproduces MMInversion_altered.m, except that track-specific
% Green matrices are interpolated from the shared full-resolution G_e/G_n/G_u
% matrices. To run only this stage, set both workflow steps to
% "final_inversion" at the top of this file.
cfg.finalInversion.enabled = true;

% Sample files used by the final layered inversion. Any number of tracks may
% be supplied. This default reuses the seven copied samp3 definitions above.
% Each MAT file must contain sampled_insar_data = [x y displacement lookE
% lookN lookU ...], with displacement in cm and x/y in meters. AZO tracks use
% sinF.dat and cosF.dat from the sample file's folder instead of look vectors.
cfg.finalInversion.samplesRoot = cfg.resampledTdeInversion.legacySamplesRoot;
cfg.finalInversion.tracks = cfg.resampledTdeInversion.legacySamples;

% Full-resolution component matrices. These project-local files contain
% G_e, G_n, and G_u, with rows ordered as adaptive_meshgrid(X,Y)(:).
cfg.finalInversion.componentFiles.east = fullfile(cfg.compositeGreens.outputDirectory, "G_e.mat");
cfg.finalInversion.componentFiles.north = fullfile(cfg.compositeGreens.outputDirectory, "G_n.mat");
cfg.finalInversion.componentFiles.vertical = fullfile(cfg.compositeGreens.outputDirectory, "G_u.mat");

% Interpolate each component to every sample location, then project it into
% AZO or LOS. "linear" matches the requested interpolation workflow. The
% block size controls memory, not the numerical result. Cached sampled Greens
% make later smoothness/boundary experiments much faster.
cfg.finalInversion.interpolation.method = "linear";
cfg.finalInversion.interpolation.columnBlockSize = 32;
cfg.finalInversion.interpolation.reuseCached = true;
cfg.finalInversion.interpolation.overwriteCached = false;

% MMInversion_altered uses smoothness=7e-1, then multiplies it by 1e-4
% internally before applying smoothingMatrix_laplace2.
cfg.finalInversion.smoothingWeight = 7e-1;
cfg.finalInversion.smoothingInternalScale = 1e-4;
cfg.finalInversion.smoothingStrikeRatio = 1;
cfg.finalInversion.smoothingDipRatio = 3;

% Legacy zero-slip boundary strengths. Surface slip remains unconstrained.
cfg.finalInversion.boundaryConditions.bottom.enabled = true;
cfg.finalInversion.boundaryConditions.bottom.ratio = 6e-4;
cfg.finalInversion.boundaryConditions.left.enabled = true;
cfg.finalInversion.boundaryConditions.left.ratio = 3e-4;
cfg.finalInversion.boundaryConditions.right.enabled = true;
cfg.finalInversion.boundaryConditions.right.ratio = 3e-4;
cfg.finalInversion.boundaryConditions.top.enabled = false;

% Polarity and slip limits in cm. [-1 0 0] constrains strike slip to the
% negative convention used by the legacy right-lateral Myanmar model.
cfg.finalInversion.constraints.polarity = [-1 0 0];
cfg.finalInversion.constraints.maxStrikeSlipCm = 550;
cfg.finalInversion.constraints.maxDipSlipCm = 1000;
cfg.finalInversion.constraints.customLowerBounds = [];
cfg.finalInversion.constraints.customUpperBounds = [];

% Relative track-type weights. With normalizeEachTrack=true, every track's
% Green matrix and data vector are divided by its number of observations,
% exactly as calc_weight_insar_error(ones(N,1)) in the legacy inversion.
cfg.finalInversion.trackWeights.azo = 1;
cfg.finalInversion.trackWeights.los = 1;
cfg.finalInversion.trackWeights.rng = 1;
cfg.finalInversion.normalizeEachTrack = true;

% Velocity model used to assign shear modulus and calculate M0/Mw for the
% simple-triangular export. The copied file is independent of the read-only
% legacy research directory.
cfg.finalInversion.velocityModelFile = fullfile(cfg.project.outputRoot, ...
    "data", "raw", "legacy_reference", "vel2.mat");
cfg.finalInversion.uniformShearModulusPa = []; % [] uses depth-dependent vel2.mat
cfg.finalInversion.simpleTriangularMatFile = fullfile(cfg.project.outputRoot, ...
    "inversions", "simple_triangular_model.mat");
cfg.finalInversion.simpleTriangularTextFile = fullfile(cfg.project.outputRoot, ...
    "inversions", "simple_triangular_model.txt");

% Legacy two-panel nodal-slip visualization. Axis units are km and colorbar
% units are meters. Strike slip is sign-flipped for display so right-lateral
% slip is positive, exactly as plotSlipNodes3D.m.
cfg.visualization.finalInversion.enabled = true;
cfg.visualization.finalInversion.colormap = "parula";
cfg.visualization.finalInversion.axisRange = [-100 100 -420 70 -25 0];
cfg.visualization.finalInversion.view = [45 10];
cfg.visualization.finalInversion.lineWidth = 2;
cfg.visualization.finalInversion.strikeColorbarRange = [0 5];
cfg.visualization.finalInversion.dipColorbarRange = [-0.5 0.5];
% Normalized [left bottom width height] positions. These defaults keep the
% horizontal bars clear of the 3-D y axes; increase left or decrease bottom
% to move a bar farther toward the right-bottom direction.
cfg.visualization.finalInversion.strikeColorbarPosition = [0.76 0.60 0.17 0.018];
cfg.visualization.finalInversion.dipColorbarPosition = [0.76 0.12 0.17 0.018];
% Display origin can differ from the inversion/local-coordinate origin above.
% Set equal to cfg.insar.origin to keep the original local x/y coordinates.
cfg.visualization.finalInversion.plotOrigin = struct( ...
    "lon", 95.936, "lat", 22.011, "refLon", 96.05);
cfg.visualization.finalInversion.showTitles = false;
cfg.visualization.finalInversion.strikeTitle = "Strike-slip";
cfg.visualization.finalInversion.dipTitle = "Dip-slip";
cfg.visualization.finalInversion.triangleTitle = ...
    "Triangle centers colored by slip magnitude (m)";

% Per-track data/model/residual plots. Empty ranges are chosen symmetrically
% about zero; set explicit [min max] values to force a common range.
cfg.visualization.finalInversion.fitPlots.enabled = true;
cfg.visualization.finalInversion.fitPlots.colormap = "turbo";
cfg.visualization.finalInversion.fitPlots.dataModelColorbarRange = [];
cfg.visualization.finalInversion.fitPlots.residualColorbarRange = [];
cfg.visualization.finalInversion.fitPlots.axisRange = [];
cfg.visualization.finalInversion.fitPlots.showTitles = true;
cfg.visualization.finalInversion.fitPlots.panelTitles = ["Data" "Model" "Residual"];

% Full-resolution residual plots following plotGrdModelDifference2.m. Panels
% (a) and (b) show sampled data/model; panel (c) compares the model generated
% from full G_e/G_n/G_u against the full observation grid. Disable this block
% when tuning only inversion parameters and the extra full-grid plots are not
% needed. It is off by default so the bundled example finishes promptly;
% set it true when the seven large-grid comparison figures are wanted.
% Full-resolution source files remain read-only.
cfg.visualization.finalInversion.fullResolutionFitPlots.enabled = false;
cfg.visualization.finalInversion.fullResolutionFitPlots.colormap = "turbo";
cfg.visualization.finalInversion.fullResolutionFitPlots.showTitles = true;
cfg.visualization.finalInversion.fullResolutionFitPlots.faultTraceFiles = ...
    cfg.geometry.surfaceTraceFiles;
cfg.visualization.finalInversion.fullResolutionFitPlots.plotOrigin = ...
    struct("lon", 95.936, "lat", 22.011, "refLon", 96.05);
cfg.visualization.finalInversion.fullResolutionFitPlots.distanceThresholdKm = [];
cfg.visualization.finalInversion.fullResolutionFitPlots.saveResidualGrids = true;
cfg.visualization.finalInversion.fullResolutionFitPlots.outputDirectory = ...
    fullfile(cfg.project.outputRoot, "inversions", "full_resolution_residuals");

% Per-track full-grid definitions copied from the active calls at the end of
% MMInversion_altered.m. modelGridFile is where the full model is evaluated;
% dataGridFile supplies the full observed displacement used for the residual.
% They differ for AZO and are the same file for these LOS tracks.
fullGridRoot = cfg.paths.geodeticRoot;
cfg.visualization.finalInversion.fullResolutionFitPlots.tracks = [
    struct("name", "SEN_A70_azo_samp3", "type", "azo", ...
        "modelGridFile", fullfile(fullGridRoot, "SEN", "A70", "azo_ll_low2.grd"), ...
        "dataGridFile", fullfile(fullGridRoot, "SEN", "A70", "res_deramped_azo_low.grd"), ...
        "lookEFile", "", "lookNFile", "", "lookUFile", "", ...
        "colorbarRange", [-150 150], "axisRange", [-50 110 -410 100], ...
        "title", "Sentinel-1 Ascending Track 70 Azimuth Offset")
    struct("name", "SEN_A143_azo_samp3", "type", "azo", ...
        "modelGridFile", fullfile(fullGridRoot, "SEN", "A143", "azo_ll_low2.grd"), ...
        "dataGridFile", fullfile(fullGridRoot, "SEN", "A143", "res_deramped_azo_low.grd"), ...
        "lookEFile", "", "lookNFile", "", "lookUFile", "", ...
        "colorbarRange", [-150 150], "axisRange", [-100 110 -410 100], ...
        "title", "Sentinel-1 Ascending Track 143 Azimuth Offset")
    struct("name", "SEN_D33_azo_samp3", "type", "azo", ...
        "modelGridFile", fullfile(fullGridRoot, "SEN", "D33", "azo_ll_low2.grd"), ...
        "dataGridFile", fullfile(fullGridRoot, "SEN", "D33", "res_deramped_azo_low.grd"), ...
        "lookEFile", "", "lookNFile", "", "lookUFile", "", ...
        "colorbarRange", [-150 150], "axisRange", [-100 110 -410 100], ...
        "title", "Sentinel-1 Descending Track 33 Azimuth Offset")
    struct("name", "SEN_D106_azo_samp3", "type", "azo", ...
        "modelGridFile", fullfile(fullGridRoot, "SEN", "D106", "azo_ll_low2.grd"), ...
        "dataGridFile", fullfile(fullGridRoot, "SEN", "D106", "res_deramped_azo_low.grd"), ...
        "lookEFile", "", "lookNFile", "", "lookUFile", "", ...
        "colorbarRange", [-150 150], "axisRange", [-100 110 -410 100], ...
        "title", "Sentinel-1 Descending Track 106 Azimuth Offset")
    struct("name", "ALOS2_A152_los_samp3", "type", "los", ...
        "modelGridFile", fullfile(fullGridRoot, "ALOS2", "A152", "los_ll_low.grd"), ...
        "dataGridFile", fullfile(fullGridRoot, "ALOS2", "A152", "los_ll_low.grd"), ...
        "lookEFile", fullfile(fullGridRoot, "ALOS2", "A152", "look_e_low.grd"), ...
        "lookNFile", fullfile(fullGridRoot, "ALOS2", "A152", "look_n_low.grd"), ...
        "lookUFile", fullfile(fullGridRoot, "ALOS2", "A152", "look_u_low.grd"), ...
        "colorbarRange", [-30 30], "axisRange", [-140 140 -410 100], ...
        "title", "ALOS-2 Ascending Track 152 LOS Displacement")
    struct("name", "ALOS2_D41_los_samp3", "type", "los", ...
        "modelGridFile", fullfile(fullGridRoot, "ALOS2", "D41", "los_ll_low.grd"), ...
        "dataGridFile", fullfile(fullGridRoot, "ALOS2", "D41", "los_ll_low.grd"), ...
        "lookEFile", fullfile(fullGridRoot, "ALOS2", "D41", "look_e_low.grd"), ...
        "lookNFile", fullfile(fullGridRoot, "ALOS2", "D41", "look_n_low.grd"), ...
        "lookUFile", fullfile(fullGridRoot, "ALOS2", "D41", "look_u_low.grd"), ...
        "colorbarRange", [-30 30], "axisRange", [-140 140 -410 100], ...
        "title", "ALOS-2 Descending Track 41 LOS Displacement")
    struct("name", "ALOS2_D42_los_samp3", "type", "los", ...
        "modelGridFile", fullfile(fullGridRoot, "ALOS2", "D42", "los_ll_low.grd"), ...
        "dataGridFile", fullfile(fullGridRoot, "ALOS2", "D42", "los_ll_low.grd"), ...
        "lookEFile", fullfile(fullGridRoot, "ALOS2", "D42", "look_e_low.grd"), ...
        "lookNFile", fullfile(fullGridRoot, "ALOS2", "D42", "look_n_low.grd"), ...
        "lookUFile", fullfile(fullGridRoot, "ALOS2", "D42", "look_u_low.grd"), ...
        "colorbarRange", [-30 30], "axisRange", [-140 140 -410 100], ...
        "title", "ALOS-2 Descending Track 42 LOS Displacement")
];

%% ------------------------------------------------------------------------
%  18. Bonus helper: smoothness L-curve
%  ------------------------------------------------------------------------
% Run workflows/run_l_curve.m after the corresponding inversion has been run
% once. This helper does not modify the saved base inversion. It reuses its
% Green matrices, data, geometry, weights, bounds, and boundary conditions.
cfg.lCurve.enabled = true;

% Choose "composite" for the final layered/composite-Green inversion, or
% "tde" for the earlier resampled homogeneous TDE inversion.
cfg.lCurve.inversionType = "composite";

% User-facing smoothness values. Composite trials apply the same internal
% 1e-4 scale as MMInversion_altered; TDE trials use these values directly,
% matching their respective inversion stages.
cfg.lCurve.smoothnessValues = logspace(-1, 1, 12);

% false matches the legacy loop by solving every trial from an empty starting
% model. Set true to reuse each solution as the next lsqlin starting model.
cfg.lCurve.useWarmStart = false;

% The legacy function called this value RMS_misfit, but calculated the raw
% unweighted residual sum of squares: sum((G_raw*m-data_raw).^2). Keep
% "raw_sse" to reproduce that statistic. Other accepted values are
% "raw_rms", "weighted_sse", and "weighted_l2"; all metrics are saved.
cfg.lCurve.misfitMetric = "raw_sse";
cfg.lCurve.saveSolutions = true;

cfg.visualization.lCurve.enabled = true;
cfg.visualization.lCurve.xScale = "log"; % "linear" or "log"
cfg.visualization.lCurve.yScale = "log"; % "linear" or "log"
cfg.visualization.lCurve.lineColor = [0.1 0.1 0.1];
cfg.visualization.lCurve.markerColor = [0.85 0.2 0.12];
cfg.visualization.lCurve.markerSize = 8;
cfg.visualization.lCurve.lineWidth = 1.5;
cfg.visualization.lCurve.showSmoothnessLabels = true;
cfg.visualization.lCurve.showTitle = true;
cfg.visualization.lCurve.title = ""; % empty uses an automatic TDE/composite title
cfg.visualization.lCurve.xRange = [];
cfg.visualization.lCurve.yRange = [];

%% ------------------------------------------------------------------------
%  19. Bonus helper: shallow slip deficit
%  ------------------------------------------------------------------------
% Run workflows/run_shallow_slip_deficit.m after the requested inversion
% result exists. This is post-processing only and does not alter either model.
cfg.shallowSlipDeficit.enabled = true;

% Analyze one or both production models. Accepted values are "composite",
% "tde", "both", or a string array such as ["composite" "tde"].
cfg.shallowSlipDeficit.modelTypes = ["composite" "tde"];

% The attached nodal legacy script uses total slip magnitude. The alternatives
% use absolute component amplitudes so polarity does not cancel spatial means.
% Accepted values: "magnitude", "strike_magnitude", "dip_magnitude".
cfg.shallowSlipDeficit.component = "magnitude";

% Normalize every segment's profile by its maximum mean slip within this
% depth range [minimum maximum] in km. [0 Inf] uses the full model profile.
cfg.shallowSlipDeficit.referenceDepthRangeKm = [0 Inf];

% In addition to the shallowest-bin deficit, report the support-weighted mean
% deficit from the surface through this depth in km.
cfg.shallowSlipDeficit.shallowDepthMaxKm = 5;

% Analyze any number of along-strike regions. Bounds are local model x/y in
% km and use (lower, upper] intervals, matching the non-overlapping legacy
% south/middle/north selections. faultIds=[] means all configured faults.
cfg.shallowSlipDeficit.segments = struct( ...
    "name", "whole_fault", ...
    "faultIds", [], ...
    "xRangeKm", [-Inf Inf], ...
    "yRangeKm", [-Inf Inf]);

% Example matching the commented legacy south/middle/north segmentation:
% cfg.shallowSlipDeficit.segments = [
%     struct("name","south", "faultIds",[], "xRangeKm",[-Inf Inf], "yRangeKm",[-Inf -100])
%     struct("name","middle","faultIds",[], "xRangeKm",[-Inf Inf], "yRangeKm",[-100 50])
%     struct("name","north", "faultIds",[], "xRangeKm",[-Inf Inf], "yRangeKm",[50 Inf])
% ];

% Composite/nodal profiles reproduce Shallow_Slip_Deficit_Segmented_New.m:
% arithmetic mean slip at each native nodal layer.
cfg.shallowSlipDeficit.composite.averaging = "arithmetic";

% TDE profiles default to native triangle-centroid depth levels and an
% area-weighted mean because the triangular patches have unequal areas.
% Set depthGrouping="fixed_bins" to use binEdgesKm or binWidthKm instead.
cfg.shallowSlipDeficit.tde.depthGrouping = "native_centroid_levels";
cfg.shallowSlipDeficit.tde.depthToleranceKm = 1e-4;
cfg.shallowSlipDeficit.tde.averaging = "area_weighted";
cfg.shallowSlipDeficit.tde.binEdgesKm = [];
cfg.shallowSlipDeficit.tde.binWidthKm = 1;

% Shallow-deficit plot controls. Empty depthRangeKm uses the model extent.
cfg.visualization.shallowSlipDeficit.enabled = true;
cfg.visualization.shallowSlipDeficit.colormap = "lines";
cfg.visualization.shallowSlipDeficit.lineWidth = 2;
cfg.visualization.shallowSlipDeficit.markerSize = 7;
cfg.visualization.shallowSlipDeficit.normalizedSlipRange = [0 1.05];
cfg.visualization.shallowSlipDeficit.deficitRange = [0 1];
cfg.visualization.shallowSlipDeficit.depthRangeKm = [];
cfg.visualization.shallowSlipDeficit.showShallowCutoff = true;
cfg.visualization.shallowSlipDeficit.showTitle = true;
cfg.visualization.shallowSlipDeficit.title = "";

%% ------------------------------------------------------------------------
%  20. Bonus helper: forward model on independent data
%  ------------------------------------------------------------------------
% Run workflows/run_forward_model.m after final_inversion. This helper does
% not solve for slip. It projects the saved final composite model onto data
% that were not used in the inversion by interpolating the shared full
% G_e/G_n/G_u matrices at those independent sample coordinates.
cfg.forwardModel.enabled = true;
cfg.forwardModel.modelResultFile = fullfile(cfg.project.outputRoot, ...
    "inversions", "final_result.mat");
cfg.forwardModel.componentFiles = cfg.finalInversion.componentFiles;

% Add any number of independent datasets here. Standard sampled MAT files
% contain sampled_insar_data = [x_m y_m displacement_cm lookE lookN lookU].
% AZO uses sinF/cosF instead of columns 4:6. Direct "east", "north", and
% "up" component types need only the first three columns.
cfg.forwardModel.samplesRoot = fullfile(cfg.project.outputRoot, ...
    "data", "raw", "independent");
cfg.forwardModel.tracks = struct( ...
    "name", "Sentinel2_independent_AZO", ...
    "type", "azo", ...                 % "azo", "los", "rng", "east", "north", or "up"
    "relativeFile", fullfile("SEN2", "azo_samp3.mat"), ...
    "sampleFile", "", ...             % absolute path overrides relativeFile
    "dataGridFile", fullfile("SEN2", "azo_ll.grd"), ... % full observed grid
    "lookEFile", "", ...              % required full grid for LOS/RNG
    "lookNFile", "", ...              % required full grid for LOS/RNG
    "lookUFile", "", ...              % required full grid for LOS/RNG
    "sinF", [], ...                    % [] reads sinF.dat beside the sample file
    "cosF", [], ...                    % [] reads cosF.dat beside the sample file
    "sinFFile", "", ...               % optional explicit heading file
    "cosFFile", "");                   % optional explicit heading file

% Interpolation is numerically the same bilinear full-XYZ sampling used by
% final_inversion. A separate cache is saved for each independent track so
% changing the slip model does not require resampling the Green matrices.
cfg.forwardModel.interpolation.method = "linear";
cfg.forwardModel.interpolation.columnBlockSize = 32;
cfg.forwardModel.interpolation.reuseCached = true;
cfg.forwardModel.interpolation.overwriteCached = false;
cfg.forwardModel.interpolation.saveCached = true;

% Full-resolution comparison following plotForwardModeling2.m. The solved
% model is evaluated on the complete XYZ-Green grid, projected into the
% independent sensor direction, and interpolated through the fault barrier
% onto dataGridFile. The observed, model, and residual panels are all full
% grids. Set enabled=false only when sampled predictions are sufficient.
cfg.forwardModel.fullResolution.enabled = true;
cfg.forwardModel.fullResolution.faultTraceFile = cfg.compositeGreens.meshTraceFile;
cfg.forwardModel.fullResolution.interpolationMethod = "linear";
cfg.forwardModel.fullResolution.extrapolation = "none";
cfg.forwardModel.fullResolution.saveGrids = true;
cfg.forwardModel.fullResolution.reuseComponentCache = true;
cfg.forwardModel.fullResolution.reuseGridCache = true;

% Full-resolution data/model/residual figure controls. Values and colorbars
% are in cm; map coordinates and axisRange are in km. Empty colorbar ranges
% are selected symmetrically around zero.
cfg.visualization.forwardModel.enabled = true;
cfg.visualization.forwardModel.colormap = "turbo";
cfg.visualization.forwardModel.dataModelColorbarRange = [-150 150];
cfg.visualization.forwardModel.residualColorbarRange = [-150 150];
cfg.visualization.forwardModel.axisRange = [-60 60 -320 240];
cfg.visualization.forwardModel.markerSize = 18;
cfg.visualization.forwardModel.showTitles = true;
cfg.visualization.forwardModel.panelTitles = ...
    ["Independent data" "Forward model" "Data - model"];
cfg.visualization.forwardModel.showOverallTitle = true;
cfg.visualization.forwardModel.title = "Sentinel-2 independent forward model";

%% ------------------------------------------------------------------------
%  21. Runtime behavior
%  ------------------------------------------------------------------------
% overwrite=false avoids replacing products unless a stage explicitly does so.
cfg.runtime.overwrite = false;

% saveIntermediate=true keeps stage products so users can restart later.
cfg.runtime.saveIntermediate = true;

% verbose=true prints progress messages during longer calculations.
cfg.runtime.verbose = true;

end
