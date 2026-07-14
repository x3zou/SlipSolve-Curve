function p = paths(cfg)
%PATHS Return canonical project output paths.

root = string(cfg.project.outputRoot);

p = struct();
p.root = root;
p.rawData = fullfile(root, "data", "raw");
p.processedData = fullfile(root, "data", "processed");
p.geometry = fullfile(root, "geometry");
p.greens = fullfile(root, "greens");
p.inversions = fullfile(root, "inversions");
p.figures = fullfile(root, "figures");
p.logs = fullfile(root, "logs");
p.insarSubsamplingFigures = fullfile(p.figures, "insar_subsampling");
p.faultGeometryFigures = fullfile(p.figures, "fault_geometry");
p.triangularMeshFigures = fullfile(p.figures, "triangular_mesh");
p.quickInversionFigures = fullfile(p.figures, "quick_inversion");
p.modelSamplingFigures = fullfile(p.figures, "model_sampling");
p.resampledTdeFigures = fullfile(p.figures, "resampled_tde_inversion");
p.layeredGreensFigures = fullfile(p.figures, "layered_greens");
p.compositeGreensFigures = fullfile(p.figures, "composite_greens");
p.finalInversionFigures = fullfile(p.figures, "final_layered_inversion");
p.finalFullResolutionFigures = fullfile(p.finalInversionFigures, "full_resolution_fits");
p.lCurveFigures = fullfile(p.figures, "l_curve");
p.shallowSlipDeficitFigures = fullfile(p.figures, "shallow_slip_deficit");
p.forwardModelFigures = fullfile(p.figures, "forward_model");

p.insarSubsampled = fullfile(p.processedData, "insar_quadtree.mat");
p.resampledObservations = fullfile(p.processedData, "insar_model_resampled.mat");
p.modelBasedSampling = fullfile(p.processedData, "model_based_sampling");
p.faultGeometry = fullfile(p.geometry, "fault_geometry.mat");
p.triangularMesh = fullfile(p.geometry, "triangular_mesh.mat");
p.quickResult = fullfile(p.inversions, "quick_result.mat");
p.resampledTdeResult = fullfile(p.inversions, "resampled_tde_result.mat");
p.quickModelSampleList = fullfile(p.inversions, "quick_model_sample_list.txt");
p.quickModelGrids = fullfile(p.inversions, "quick_model_grids");
p.layeredGreens = fullfile(p.greens, "layered_greens.mat");
p.layeredGreensDatabase = fullfile(p.greens, "edgrn");
p.compositeGreens = fullfile(p.greens, "composite_greens.mat");
p.compositeGreensDirectory = fullfile(p.greens, "composite_full_xyz");
p.sampledCompositeGreensDirectory = fullfile(p.greens, "composite_sampled_tracks");
p.sampledCompositeGreens = fullfile(p.sampledCompositeGreensDirectory, "sampled_composite_greens.mat");
p.finalResult = fullfile(p.inversions, "final_result.mat");
p.simpleTriangularMat = fullfile(p.inversions, "simple_triangular_model.mat");
p.simpleTriangularText = fullfile(p.inversions, "simple_triangular_model.txt");
p.finalFullResolutionResiduals = fullfile(p.inversions, "full_resolution_residuals");
p.lCurveDirectory = fullfile(p.inversions, "l_curve");
p.shallowSlipDeficitDirectory = fullfile(p.inversions, "shallow_slip_deficit");
p.shallowSlipDeficitResult = fullfile(p.shallowSlipDeficitDirectory, ...
    "shallow_slip_deficit.mat");
p.tdePostprocessingModel = fullfile(p.shallowSlipDeficitDirectory, ...
    "tde_model_for_postprocessing.mat");
p.forwardModelDirectory = fullfile(p.inversions, "forward_model");
p.forwardModelResult = fullfile(p.forwardModelDirectory, "forward_model_result.mat");
p.forwardModelGreensDirectory = fullfile(p.greens, "forward_model_tracks");

end
