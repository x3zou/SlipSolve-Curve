function [x_samp,y_samp,data_samp,data_samp_std,data_tree,data_dims,data_extents,nan_frac] = ...
    quadtree_unstructured2(x,y,data_in,data_index,data_extent,fault_xy,varargin)
%% Unstructured quad-tree sampling with sign-agnostic, fault-aware summarization
% If a finalized cell touches/crosses the fault, we determine the centroid’s
% geometric side relative to the fault polyline and average **only** points
% on that side. No assumption about which side has positive/negative values.
%
% INPUTS (same as your original, plus fault_xy):
%   x, y          : (n,1) vectors of point coordinates
%   data_in       : (n,1) vector of values (NaN allowed)
%   data_index    : (n,1) vector of global indices of the points
%   data_extent   : [xmin xmax ymin ymax] cell extent
%   fault_xy      : (m,2) polyline [x y]; pass [] to ignore fault handling
%
% Name-Value options (kept for compatibility with your original):
%   'rms_min', 'nan_frac_max', 'width_min', 'width_max',
%   'x_samp','y_samp','data_samp','data_samp_std','data_tree',
%   'data_dims','data_extents','nan_frac'
% Additional options:
%   'level'        : recursion depth (internal)
%   'fault_tol'    : half-width buffer around fault treated as “on-fault” (default 0.01)
%   'minPixSide'   : min #points required on chosen side before fallback (default 2)
%   'stat'         : 'mean' (default) or 'median'
%
% OUTPUTS: exactly the same as your original function.

% ---------- defaults (kept gentle to preserve original behavior) ----------
rms_min       = 0;        % user usually overrides
nan_frac_max  = 1;     % typical default
width_min     = 0;        % in same units as x,y
width_max     = inf;
x_samp=[]; y_samp=[]; data_samp=[]; data_samp_std=[];
data_tree={}; data_dims=[]; data_extents=[]; nan_frac=[];
level=0; fault_tol=0.01; minPixSide=2; stat='mean';

% ---------- parse varargin ----------
if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        key = lower(varargin{2*CC-1});
        val = varargin{2*CC};
        switch key
            case 'rms_min',        rms_min = val;
            case 'nan_frac_max',   nan_frac_max = val;
            case 'width_min',      width_min = val;
            case 'width_max',      width_max = val;
            case 'x_samp',         x_samp = val;
            case 'y_samp',         y_samp = val;
            case 'data_samp',      data_samp = val;
            case 'data_samp_std',  data_samp_std = val;
            case 'data_tree',      data_tree = val;
            case 'data_dims',      data_dims = val;
            case 'data_extents',   data_extents = val;
            case 'nan_frac',       nan_frac = val;
            case 'level',          level = val;
            case 'fault_tol',      fault_tol = val;
            case 'minpixside',     minPixSide = val;
            case 'stat'
                stat = lower(val);
                if ~ismember(stat,{'mean','median'}), stat='mean'; end
            otherwise
                error('Unrecognized Keyword: %s', varargin{2*CC-1});
        end
    end
end

% ---------- cell geometry ----------
data_dim = [data_extent(2)-data_extent(1), data_extent(4)-data_extent(3)];
cell_dim = data_dim/2;

% ---------- statistics for split criterion (unchanged intent) ----------
data_mean = nanmean(data_in);
data_rms  = nanstd(data_in - data_mean);

% ---------- split or finalize ----------
if ( (data_rms > rms_min && all(cell_dim > width_min)) || all(cell_dim > width_max) )
    % child extents (TL, TR, BR, BL)
    child_ext = {
        [data_extent(1),               data_extent(1)+cell_dim(1), data_extent(3)+cell_dim(2), data_extent(4)              ]; % TL
        [data_extent(1)+cell_dim(1),   data_extent(2),             data_extent(3)+cell_dim(2), data_extent(4)              ]; % TR
        [data_extent(1)+cell_dim(1),   data_extent(2),             data_extent(3),             data_extent(3)+cell_dim(2)  ]; % BR
        [data_extent(1),               data_extent(1)+cell_dim(1), data_extent(3),             data_extent(3)+cell_dim(2)  ]  % BL
    };

    % masks per child
    child_mask = {
        (x >= data_extent(1))             & (x <= data_extent(1)+cell_dim(1)) & ...
        (y >  data_extent(3)+cell_dim(2)) & (y <= data_extent(4));
        (x >  data_extent(1)+cell_dim(1)) & (x <= data_extent(2))             & ...
        (y >= data_extent(3)+cell_dim(2)) & (y <= data_extent(4));
        (x >= data_extent(1)+cell_dim(1)) & (x <= data_extent(2))             & ...
        (y >= data_extent(3))             & (y <  data_extent(3)+cell_dim(2));
        (x >= data_extent(1))             & (x <  data_extent(1)+cell_dim(1)) & ...
        (y >= data_extent(3))             & (y <= data_extent(3)+cell_dim(2))
    };

    % recurse
    for ii = 1:4
        m = child_mask{ii};
        [x_samp,y_samp,data_samp,data_samp_std,data_tree,data_dims,data_extents,nan_frac] = ...
            quadtree_unstructured2( x(m), y(m), data_in(m), data_index(m), child_ext{ii}, fault_xy, ...
            'level', level+1, 'rms_min', rms_min, 'nan_frac_max', nan_frac_max, ...
            'width_min', width_min, 'width_max', width_max, ...
            'x_samp', x_samp, 'y_samp', y_samp, ...
            'data_samp', data_samp, 'data_samp_std', data_samp_std, ...
            'data_tree', data_tree, 'data_dims', data_dims, ...
            'data_extents', data_extents, 'nan_frac', nan_frac, ...
            'fault_tol', fault_tol, 'minPixSide', minPixSide, 'stat', stat);
    end

    if level==0
        x_samp        = x_samp(:);
        y_samp        = y_samp(:);
        data_samp     = data_samp(:);
        data_samp_std = data_samp_std(:);
    end

else
    % ---------- finalize: compute representative point/value ----------
    i_nan = isnan(data_in);
    nan_frac_data = sum(i_nan)/max(1,numel(data_in));

    if nan_frac_data > nan_frac_max
        x_out = NaN; y_out = NaN; d_out = NaN; dstd_out = NaN;

    else
        good = isfinite(data_in);
        xg = x(good); yg = y(good); zg = data_in(good);

        if isempty(xg)
            x_out = NaN; y_out = NaN; d_out = NaN; dstd_out = NaN;

        else
            useMask = true(size(zg));

            if ~isempty(fault_xy)
                % classify geometric side (+1 / -1) for each point; 0 if within fault_tol
                sP = classifySideMany(xg, yg, fault_xy(:,1), fault_xy(:,2), fault_tol);

                % cell contacts the fault if it has points from both sides or any 0 (corridor)
                contactsFault = any(sP==0) || (any(sP>0) && any(sP<0));

                if contactsFault
                    % centroid (stable, from extent)
                    xc = 0.5*(data_extent(1)+data_extent(2));
                    yc = 0.5*(data_extent(3)+data_extent(4));
                    sC = pointSideToPolyline(xc, yc, fault_xy(:,1), fault_xy(:,2), fault_tol);

                    if sC == 0
                        % centroid on/near trace -> pick majority geometric side
                        nPos = nnz(sP>0); nNeg = nnz(sP<0);
                        if nPos >= nNeg, useMask = (sP>0); else, useMask = (sP<0); end
                    else
                        % use centroid side; if too few points, fall back to majority side
                        useMask = (sP == sC);
                        if nnz(useMask) < minPixSide
                            nPos = nnz(sP>0); nNeg = nnz(sP<0);
                            if nPos >= nNeg, useMask = (sP>0); else, useMask = (sP<0); end
                        end
                    end
                end
            end

            if ~any(useMask)
                x_out = NaN; y_out = NaN; d_out = NaN; dstd_out = NaN;
            else
                xs = xg(useMask); ys = yg(useMask); zs = zg(useMask);

                % representative location (consistent with many implementations)
                x_out = mean(xs);
                y_out = mean(ys);
                % x_out = 0.5*(data_extent(1)+data_extent(2));
                % y_out = 0.5*(data_extent(3)+data_extent(4));


                % statistic
                switch stat
                    case 'median', d_out = median(zs);
                    otherwise,     d_out = mean(zs);
                end
                dstd_out = nanstd(zs);
            end
        end
    end

    % append to outputs (identical structure to original)
    x_samp            = [x_samp; x_out];
    y_samp            = [y_samp; y_out];
    data_samp         = [data_samp; d_out];
    data_samp_std     = [data_samp_std; dstd_out];
    data_tree{end+1}  = data_index;
    data_dims         = [data_dims; data_dim];
    data_extents      = [data_extents; data_extent];
    nan_frac          = [nan_frac; nan_frac_data];
end
end

% ======================= helpers (local functions) ========================

function s = pointSideToPolyline(xc,yc,fx,fy,tol)
% Returns +1 (one geometric side), -1 (the other), or 0 if within tol of the polyline.
% Note: If you reverse the polyline order, all signs flip consistently;
% the selected physical side remains the same.
if nargin<5, tol=0; end
d2min = inf; s = int8(0);
for k=1:numel(fx)-1
    x0=fx(k); y0=fy(k);
    vx=fx(k+1)-x0; vy=fy(k+1)-y0;
    vv=vx*vx + vy*vy + eps;
    wx=xc-x0; wy=yc-y0;
    t = max(0,min(1,(wx*vx + wy*vy)/vv));
    px=x0 + t*vx; py=y0 + t*vy;
    dx=xc-px; dy=yc-py; d2=dx*dx + dy*dy;
    if d2<d2min
        d2min=d2;
        s = int8(sign(vx*wy - vy*wx));  % cross(v, w)
    end
end
if tol>0 && d2min<=tol^2, s=int8(0); end
end

function s = classifySideMany(xp,yp,fx,fy,tol)
% Vectorized over points (loop over points; segments handled inside pointSideToPolyline)
if nargin<5, tol=0; end
s = zeros(size(xp),'int8');
for i=1:numel(xp)
    s(i) = pointSideToPolyline(xp(i),yp(i),fx,fy,tol);
end
end
