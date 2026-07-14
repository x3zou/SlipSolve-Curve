function [X, Y, xvec, yvec] = adaptive_meshgrid(faultXY, xRange, yRange, D1, D2, M)
%ADAPTIVE_MESHGRID Create a rectilinear meshgrid with finer spacing near a fault.
%
% Syntax
%   [X, Y] = adaptive_meshgrid(faultXY, xRange, yRange, D1, D2, M)
%   [X, Y, xvec, yvec] = ...
%
% Inputs
%   faultXY : N×2 array of [x, y] coordinates along the fault polyline.
%   xRange  : [xmin xmax] overall X extent to cover.
%   yRange  : [ymin ymax] overall Y extent to cover.
%   D1      : fine spacing (used near the fault; D1 > 0).
%   D2      : coarse spacing (used elsewhere; D2 > D1).
%   M       : margin distance expanding the fault's bounding box (M ≥ 0).
%
% Outputs
%   X, Y    : meshgrid arrays (nonuniform but rectilinear).
%   xvec    : underlying X coordinates used to form X/Y (optional).
%   yvec    : underlying Y coordinates used to form X/Y (optional).
%
% Notes
%   * Because meshgrid produces rectilinear grids, we approximate the
%     "within M of the fault" region by a bounding box of the fault
%     expanded by ±M in both directions. Inside that box we use D1;
%     outside we use D2.
%   * Endpoints are always included; spacing is adjusted to land on bounds.
%
% Example
%   faultXY = [0 0; 10 5; 20 5];                 % polyline
%   xRange = [-10 40]; yRange = [-10 20];
%   D1 = 0.5; D2 = 2.0; M = 5;
%   [X,Y] = adaptive_meshgrid(faultXY, xRange, yRange, D1, D2, M);
%   plot(faultXY(:,1), faultXY(:,2), 'r-', 'LineWidth', 2); hold on;
%   plot(X, Y, 'k.'); axis equal; box on;
%
% Author: (you)
% -------------------------------------------------------------------------

    % ---- Input checks
    arguments
        faultXY (:,2) double
        xRange (1,2) double {mustBeIncreasing(xRange)}
        yRange (1,2) double {mustBeIncreasing(yRange)}
        D1 (1,1) double {mustBePositive}
        D2 (1,1) double {mustBePositive}
        M  (1,1) double {mustBeNonnegative}
    end
    if ~(D2 > D1)
        error('D2 must be greater than D1.');
    end

    % ---- Expanded fine box around the fault
    fxmin = min(faultXY(:,1)); fxmax = max(faultXY(:,1));
    fymin = min(faultXY(:,2)); fymax = max(faultXY(:,2));

    fineXmin = max(xRange(1), fxmin - M);
    fineXmax = min(xRange(2), fxmax + M);
    fineYmin = max(yRange(1), fymin - M);
    fineYmax = min(yRange(2), fymax + M);

    % If the expanded box collapses (e.g., M=0 and degenerate fault), enforce a tiny width
    if fineXmax < fineXmin, fineXmax = fineXmin; end
    if fineYmax < fineYmin, fineYmax = fineYmin; end

    % ---- Build X vector: coarse outside, fine inside
    x_coarse_left  = stepped_vec(xRange(1), fineXmin, D2, 'left');   % [xmin .. fineXmin]
    x_fine         = stepped_vec(fineXmin, fineXmax, D1, 'middle');  % (fineXmin .. fineXmax)
    x_coarse_right = stepped_vec(fineXmax, xRange(2), D2, 'right');  % [fineXmax .. xmax]

    % Merge and uniquify with tolerance to avoid duplicates
    xvec = merge_with_tol([x_coarse_left, x_fine, x_coarse_right], min(D1, D2));

    % ---- Build Y vector similarly
    y_coarse_bottom = stepped_vec(yRange(1), fineYmin, D2, 'left');
    y_fine          = stepped_vec(fineYmin, fineYmax, D1, 'middle');
    y_coarse_top    = stepped_vec(fineYmax, yRange(2), D2, 'right');

    yvec = merge_with_tol([y_coarse_bottom, y_fine, y_coarse_top], min(D1, D2));

    % ---- Meshgrid (note: meshgrid takes xvec (cols) and yvec (rows))
    [X, Y] = meshgrid(xvec, yvec);
end

% ---------- Helpers ----------

function v = stepped_vec(a, b, step, where)
% Create a vector from a to b using the given step, ensuring endpoints are included
% and avoiding duplicated boundaries when concatenating segments.

    if a == b
        v = a;
        return
    end

    switch where
        case 'left'
            if a < b
                v = a:step:b;
            else
                v = [];
            end
        case 'right'
            if a < b
                v = a:step:b;
            else
                v = [];
            end
        case 'middle'
            if a < b
                % Open interval in theory, but include endpoints for stability.
                v = a:step:b;
            else
                v = [];
            end
        otherwise
            error('Unknown segment location "%s".', where);
    end

    % Ensure exact inclusion of boundaries
    if ~isempty(v)
        v(1)  = a;
        v(end)= b;
    else
        v = [a b];
    end
end

function vout = merge_with_tol(v, baseStep)
% Merge, sort, and uniquify with a numerical tolerance.
    v = sort(v(:).');
    tol = max(1e-9, baseStep*1e-6);
    vout = v(1);
    for k = 2:numel(v)
        if abs(v(k) - vout(end)) > tol
            vout(end+1) = v(k); %#ok<AGROW>
        else
            % Snap to the earlier value to keep exact grid lines aligned
            % (no action needed; we just drop the near-duplicate)
        end
    end
end

function mustBeIncreasing(v)
    if ~(isscalar(v(1)) && isscalar(v(2)) && v(2) > v(1))
        error('Range must be [min max] with max > min.');
    end
end

function mustBeNonnegative(v)
    if ~(isnumeric(v) && isscalar(v) && v >= 0)
        error('Value must be nonnegative.');
    end
end
