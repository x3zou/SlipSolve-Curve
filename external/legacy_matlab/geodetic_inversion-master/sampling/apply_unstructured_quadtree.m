function [x_samp, y_samp, data_samp, data_samp_std, nan_frac] = ...
    apply_unstructured_quadtree(x, y, dataVals, data_tree, nan_frac_max)
% APPLY_UNSTRUCTURED_QUADTREE Downsample a gridded dataset using a pre-constructed quadtree.
% Originally coded by Ellis. Translated to MATLAB by Xiaoyu
%   [x_samp, y_samp, data_samp, data_samp_std, nan_frac] = ...
%     apply_unstructured_quadtree(x, y, dataVals, data_tree, nan_frac_max)
%
%   Inputs:
%       x, y       : (n,) gridded x and y coordinates, flattened
%       dataVals   : (n,) gridded data,flattened
%       data_tree  : 1×k cell array where each cell is a vector of indices 
%                    referencing elements in x, y, dataVals
%       nan_frac_max : Threshold for maximum fraction of NaNs allowable 
%                      before returning NaNs for that cell
%
%   Outputs:
%       x_samp, y_samp        : Downsampled coordinates (k×1)
%       data_samp             : Downsampled data values (k×1)
%       data_samp_std         : Standard deviations per cell (k×1)
%       nan_frac              : Fraction of NaNs in each cell (k×1)
%
%   This function mirrors the Python code that iterates over quadtree cells, 
%   checks the fraction of NaNs, and either outputs representative values 
%   or NaNs if the fraction is too high or data are constant.

    % --- Flatten x, y, data if they are 2D (m×n) ---
    % x = x(:);
    % y = y(:);
    % dataVals = dataVals(:);

    % Preallocate output arrays
    nCells = numel(data_tree);
    x_samp        = nan(nCells, 1);
    y_samp        = nan(nCells, 1);
    data_samp     = nan(nCells, 1);
    data_samp_std = nan(nCells, 1);
    nan_frac      = nan(nCells, 1);

    % Loop over cells in the quadtree
    for i = 1:nCells
        % 'idx' is the set of flattened indices for the i-th quadtree cell
        idx = data_tree{i};  

        % Identify NaNs
        i_nan_data = isnan(dataVals(idx));
        n_nan_data = sum(i_nan_data);

        % Compute fraction of NaNs
        if ~isempty(idx)
            nan_frac_data = n_nan_data / numel(idx);
        else
            nan_frac_data = NaN;
        end

        % Decide if we output actual values or NaNs
        if (nan_frac_data > nan_frac_max) || (nanstd(dataVals(idx)) == 0)
            % Too many NaNs or zero standard deviation => fill with NaNs
            x_samp(i)        = NaN;
            y_samp(i)        = NaN;
            data_samp(i)     = NaN;
            data_samp_std(i) = NaN;
            nan_frac(i)      = nan_frac_data;
        else
            % Compute representative values from non-NaN samples
            non_nan_idx      = ~i_nan_data; 
            x_samp(i)        = mean(x(idx(non_nan_idx)));
            y_samp(i)        = mean(y(idx(non_nan_idx)));
            data_samp(i)     = nanmean(dataVals(idx));
            % data_samp_std(i) = nanstd(dataVals(idx));
            data_samp_std(i) = rms(dataVals(idx),"omitnan");
            nan_frac(i)      = nan_frac_data;
        end
    end
end