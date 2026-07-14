function [x_samp, y_samp, data_samp, data_samp_std, nan_frac, data_tree_out] = ...
    apply_unstructured_quadtree2(x, y, dataVals, data_tree, nan_frac_max, fault_xy, varargin)
% APPLY_QUADTREE_UNSTRUCTURED2
% Apply an existing unstructured quadtree (data_tree) to a new dataset,
% with optional sign-cleaning near a provided fault polyline.
%
% Usage:
%   [xs,ys,vs,vs_std,nf,tree2] = apply_quadtree_unstructured2( ...
%       x, y, dataVals, data_tree, nan_frac_max, fault_xy, ...
%       'sign_clean', true, 'fault_tol', 0, 'minPixSide', 6, 'stat', 'mean');
%
% Inputs
%   x, y           : (n,1) or (1,n) vectors of coordinates (flattened)
%   dataVals       : (n,1) data values (flattened; NaN allowed)
%   data_tree      : 1×k cell; each cell is a vector of indices into x,y,dataVals
%   nan_frac_max   : scalar in [0,1]; if NaN fraction in a cell exceeds this,
%                    outputs for that cell are set to NaN (location too, to
%                    preserve original behavior)
%   fault_xy       : (m×2) [x y] polyline for the fault; pass [] if unused
%
% Name-Value options
%   'sign_clean'   : logical (default: true if ~isempty(fault_xy), else false)
%   'fault_tol'    : distance buffer around the fault treated as "on-fault" (default 0)
%   'minPixSide'   : minimum # of points required on chosen side before fallback (default 6)
%   'stat'         : 'mean' (default) or 'median' for within-cell value
%
% Outputs
%   x_samp, y_samp     : (k×1) sampled coordinates (cell-wise means of non-NaN points)
%   data_samp          : (k×1) sampled values (after sign-cleaning if enabled)
%   data_samp_std      : (k×1) RMS of used values in each cell (omit NaNs)
%   nan_frac           : (k×1) NaN fraction in each original cell
%   data_tree_out      : 1×k cell; indices actually used after sign-cleaning
%                        (identical to input data_tree when sign_clean=false)
%
% Notes
%   • Sign-cleaning is sign-agnostic: we classify geometric side relative to the
%     fault polyline; no assumption about +/− of the data. For cells contacting
%     the fault, we prefer the centroid's side; if too few points on that side,
%     we fall back to the majority side within the cell.
%   • We do NOT densify the fault polyline here (per your request).

% -------------------- options --------------------
ip = inputParser;
addParameter(ip,'sign_clean', ~isempty(fault_xy));
addParameter(ip,'fault_tol',  0);
addParameter(ip,'minPixSide', 6);
addParameter(ip,'stat',       'mean');
parse(ip,varargin{:});
sign_clean  = logical(ip.Results.sign_clean);
fault_tol   = ip.Results.fault_tol;
minPixSide  = ip.Results.minPixSide;
stat        = lower(ip.Results.stat);
if ~ismember(stat, {'mean','median'}), stat = 'mean'; end

% -------------------- prealloc -------------------
nCells           = numel(data_tree);
x_samp           = nan(nCells,1);
y_samp           = nan(nCells,1);
data_samp        = nan(nCells,1);
data_samp_std    = nan(nCells,1);
nan_frac         = nan(nCells,1);
data_tree_out    = data_tree;   % will update per cell if sign_clean applies

% -------------------- loop cells -----------------
for i = 1:nCells
    idx = data_tree{i};
    if isempty(idx)
        continue
    end

    % NaN handling (overall)
    i_nan_data    = isnan(dataVals(idx));
    n_nan_data    = sum(i_nan_data);
    nan_frac_data = n_nan_data / numel(idx);
    nan_frac(i)   = nan_frac_data;

    % Decide early exit per original behavior
    if (nan_frac_data > nan_frac_max)
        % too many NaNs -> leave outputs as NaN; keep original indices
        data_tree_out{i} = idx;
        continue
    end

    % Coordinates for location (original behavior: mean over non-NaN data points)
    non_nan_idx  = ~i_nan_data;
    if any(non_nan_idx)
        x_samp(i) = mean(x(idx(non_nan_idx)), 'omitnan');
        y_samp(i) = mean(y(idx(non_nan_idx)), 'omitnan');
    else
        % no valid data; keep NaN location per original behavior
        data_tree_out{i} = idx;
        continue
    end

    % Values (optionally sign-cleaned)
    usedMask = non_nan_idx;  % start with all non-NaN points in this cell

    if sign_clean && ~isempty(fault_xy) && numel(fault_xy)>=2
        % classify geometric side of each non-NaN point
        xg = x(idx(non_nan_idx));
        yg = y(idx(non_nan_idx));
        sP = classifySideMany(xg, yg, fault_xy(:,1), fault_xy(:,2), fault_tol);

        % Does this cell "contact" the fault?
        contactsFault = any(sP==0) || ( any(sP>0) && any(sP<0) );

        if contactsFault
            % Decide side at the REPORTED sample location (x_samp,y_samp)
            sC = pointSideToPolyline(x_samp(i), y_samp(i), fault_xy(:,1), fault_xy(:,2), fault_tol);

            if sC == 0
                % Centroid sits in corridor -> majority geometric side among non-NaN points
                nPos = nnz(sP > 0);  nNeg = nnz(sP < 0);
                useLocal = (nPos >= nNeg) & (sP ~= 0);
            else
                useLocal = (sP == sC);
                if nnz(useLocal) < minPixSide
                    % fallback to majority side (exclude corridor zeros)
                    nPos = nnz(sP > 0);  nNeg = nnz(sP < 0);
                    useLocal = (nPos >= nNeg) & (sP ~= 0);
                end
            end

            % Merge back to cell-wide mask
            tmp = false(size(non_nan_idx));
            tmp(non_nan_idx) = useLocal;
            usedMask = tmp;

            % Update data_tree_out with the *used* indices
            data_tree_out{i} = idx(usedMask);
        else
            % No contact: keep all non-NaN points; keep original indices
            data_tree_out{i} = idx(non_nan_idx);
        end
    else
        % Sign-cleaning disabled: keep non-NaN points only
        data_tree_out{i} = idx(non_nan_idx);
    end

    % If nothing left to use, leave NaNs for this cell
    if ~any(usedMask)
        % keep NaN outputs; tree_out already set
        continue
    end

    % Final value/statistic
    z = dataVals(idx(usedMask));
    switch stat
        case 'median'
            data_samp(i) = median(z, 'omitnan');
        otherwise % 'mean'
            data_samp(i) = mean(z, 'omitnan');
    end

    % Report RMS (match your prior convention)
    try
        data_samp_std(i) = rms(z, "omitnan");
    catch
        % older MATLAB: emulate rms omitnan
        zz = z(isfinite(z));
        data_samp_std(i) = sqrt(mean(zz.^2));
    end

    % Optional: if all used values are identical (std == 0), mimic original behavior
    if numel(z) > 0 && (all(~isfinite(z)) || nanstd(z)==0)
        x_samp(i) = NaN;
        y_samp(i) = NaN;
        data_samp(i) = NaN;
        data_samp_std(i) = NaN;
        % keep nan_frac as computed; keep data_tree_out as set
    end
end

end % main function

% ====================== helpers =======================

function s = pointSideToPolyline(xc,yc,fx,fy,tol)
% +1 or -1 for geometric sides; 0 if within tol of the polyline
if nargin<5, tol=0; end
d2min = inf; s = int8(0);
for k = 1:numel(fx)-1
    x0 = fx(k); y0 = fy(k);
    vx = fx(k+1)-x0; vy = fy(k+1)-y0;
    vv = vx*vx + vy*vy + eps;
    wx = xc-x0; wy = yc-y0;
    t  = max(0, min(1, (wx*vx + wy*vy)/vv));
    px = x0 + t*vx; py = y0 + t*vy;
    dx = xc - px;  dy = yc - py;
    d2 = dx*dx + dy*dy;
    if d2 < d2min
        d2min = d2;
        s = int8(sign(vx*wy - vy*wx)); % cross(v, w)
    end
end
if tol>0 && d2min <= tol^2
    s = int8(0);
end
end

function s = classifySideMany(xp,yp,fx,fy,tol)
% Classify many points by geometric side of a polyline.
if nargin<5, tol=0; end
s = zeros(size(xp),'int8');
for i = 1:numel(xp)
    s(i) = pointSideToPolyline(xp(i), yp(i), fx, fy, tol);
end
end
