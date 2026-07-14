function resampledObs = resample_from_model(cfg, quickResult, insarSub)
%RESAMPLE_FROM_MODEL Resample observations based on the quick model.

contract = [
    "Inputs: cfg.resampling, quickResult, and original/subsampled observations." newline ...
    "Output: resampledObs struct aligned to final Green's function rows." newline ...
    "Legacy target: wrap sample_from_model.m."
];

slipsolve.stages.not_implemented("resample_from_model", contract);
resampledObs = struct();

end

