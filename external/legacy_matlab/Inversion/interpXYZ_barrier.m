function [z2, meta] = interpXYZ_barrier(XYZ, x2, y2, faultXY, varargin)
%INTERPXYZ_BARRIER Interpolate with a line barrier (e.g., a fault).
%   [z2, meta] = interpXYZ_barrier(XYZ, x2, y2, faultXY, ...)
%
% Inputs
%   XYZ      : N×3 [X,Y,Z] scattered or gridded samples.
%   x2,y2    : query coordinates (same size). Vectors or 2D grids OK.
%   faultXY  : M×2 polyline coordinates of the fault [x y] in map order.
%
% Name-Value options (all optional)
%   'Method'     : 'natural' (default) | 'linear' | 'nearest' | 'v4'
%                  ('v4' = biharmonic, smooth scattered interp via griddata)
%   'Extrap'     : 'none' (default) | 'nearest' | 'linear'
%                  (ignored for 'v4'; outside convex hull returns NaN if 'none')
%   'Buffer'     : 0 (default). Distance (same units as X/Y) around fault within
%                  which input samples are dropped to prevent cross-side bleed.
%   'BlendWidth' : 0 (default). If >0, smoothly blend across a thin strip
%                  (still keeps a sharp transition if small). Set 0 for hard edge.
%   'MinPts'     : 10 (default). Minimum points required per side; if fewer,
%                  that side falls back to 'nearest' or returns NaN.
%   'Fallback'   : 'nan' (default) | 'otherSide' | 'nearest'
%                  Behavior for queries outside a side's hull.
%
% Outputs
%   z2    : interpolated values, same size as x2/y2, without cross-fault smearing.
%   meta  : struct with details (counts per side, interpolants, masks, etc.).
%
% Notes
%   - The fault is treated as an infinitely thin barrier. Samples are split by
%     which side of the *nearest fault segment* they lie on. Set 'Buffer' to drop
%     samples within a small distance of the fault if your data there are noisy.
%   - For sharp discontinuities, use Method='natural' or 'nearest'. 'v4' is smooth.
%   - For multiple faults, call this function once per fault and mask/merge as needed.

    % ---------- Parse inputs
    ip = inputParser; ip.FunctionName = mfilename;
    addParameter(ip,'Method','natural');
    addParameter(ip,'Extrap','none');
    addParameter(ip,'Buffer',0);
    addParameter(ip,'BlendWidth',0);
    addParameter(ip,'MinPts',10);
    addParameter(ip,'Fallback','nan');
    parse(ip,varargin{:});
    opt = ip.Results;

    % ---------- Validate
    validateattributes(XYZ, {'numeric'}, {'2d','ncols',3});
    if ~isequal(size(x2), size(y2)), error('x2 and y2 must have the same size.'); end
    method  = validatestring(lower(opt.Method),  {'natural','linear','nearest','v4'});
    extrap  = validatestring(lower(opt.Extrap),  {'none','nearest','linear'});
    fallback= validatestring(lower(opt.Fallback),{'nan','otherside','nearest'});

    % ---------- Clean data
    X = XYZ(:,1); Y = XYZ(:,2); Z = XYZ(:,3);
    good = isfinite(X) & isfinite(Y) & isfinite(Z);
    X = X(good); Y = Y(good); Z = Z(good);
    if isempty(X), z2 = nan(size(x2)); meta = struct('note','no valid data'); return; end

    % Optional: consolidate exact duplicate (X,Y) by averaging Z
    [~,~,g] = unique([X Y], 'rows', 'stable');
    Zu = accumarray(g, Z, [], @mean);
    XYu = [accumfirst([X Y], g) Zu]; % [X Y Z] unique
    X = XYu(:,1); Y = XYu(:,2); Z = XYu(:,3);

    % ---------- Classify samples and queries by side of the fault
    [sideData, distData] = classifySideToPolyline(X, Y, faultXY);
    [sideQ,    distQ   ] = classifySideToPolyline(x2(:), y2(:), faultXY);

    % Drop samples within buffer distance of the fault (optional)
    if opt.Buffer > 0
        keep = distData > opt.Buffer;
        X = X(keep); Y = Y(keep); Z = Z(keep); sideData = sideData(keep);
    end

    % Split by side
    pos = sideData > 0;
    neg = sideData < 0;

    meta = struct();
    meta.nPos = nnz(pos);
    meta.nNeg = nnz(neg);
    meta.method = method;
    meta.extrap = extrap;

    % ---------- Build per-side interpolants
    buildOK = true;
    Fpos = []; Fneg = [];

    if meta.nPos >= opt.MinPts
        Fpos = makeInterpolant(X(pos), Y(pos), Z(pos), method, extrap);
    else
        buildOK = false;
    end
    if meta.nNeg >= opt.MinPts
        Fneg = makeInterpolant(X(neg), Y(neg), Z(neg), method, extrap);
    else
        buildOK = false;
    end
    meta.Fpos = Fpos; meta.Fneg = Fneg;

    % If a side lacks points, plan fallback behavior
    if ~buildOK
        % Will fallback per 'Fallback' option at evaluation time.
    end

    % ---------- Evaluate
    z = nan(numel(x2),1);

    usePos = sideQ > 0;
    useNeg = sideQ < 0;

    % Evaluate main side(s)
    if ~isempty(Fpos)
        z(usePos) = evalInterpolant(Fpos, x2(usePos), y2(usePos), method, extrap);
    end
    if ~isempty(Fneg)
        z(useNeg) = evalInterpolant(Fneg, x2(useNeg), y2(useNeg), method, extrap);
    end

    % Handle points exactly on the fault: choose nearer side by signed distance 0
    onFault = ~usePos & ~useNeg;
    if any(onFault)
        % pick side based on which interpolant exists; prefer side with more points
        if meta.nPos >= meta.nNeg && ~isempty(Fpos)
            z(onFault) = evalInterpolant(Fpos, x2(onFault), y2(onFault), method, extrap);
        elseif ~isempty(Fneg)
            z(onFault) = evalInterpolant(Fneg, x2(onFault), y2(onFault), method, extrap);
        end
    end

    % Optional soft blend (very thin strip to avoid visual seams if desired)
    if opt.BlendWidth > 0 && ~isempty(Fpos) && ~isempty(Fneg)
        near = distQ <= opt.BlendWidth;
        if any(near)
            zp = evalInterpolant(Fpos, x2(near), y2(near), method, extrap);
            zn = evalInterpolant(Fneg, x2(near), y2(near), method, extrap);
            % Weight toward the side the point is on; center is 50/50 at the fault
            w = 0.5 + 0.5 * tanh( 3 * sideQ(near) .* (1 - distQ(near)/opt.BlendWidth) );
            z(near) = w .* zp + (1-w) .* zn;
        end
    end

    % Fallbacks for NaNs (outside hull or missing side)
    nanQ = isnan(z);
    if any(nanQ)
        switch fallback
            case 'otherside'
                % Try the other side's interpolant
                tryOther = true(size(nanQ));
                tryOther(usePos) = false;  % these used pos; try neg
                tryOther(useNeg) = false;  % these used neg; try pos
                idx = nanQ & (sideQ~=0);
                useOtherNeg = idx & (sideQ>0) & ~isempty(Fneg);
                useOtherPos = idx & (sideQ<0) & ~isempty(Fpos);
                if any(useOtherNeg), z(useOtherNeg) = evalInterpolant(Fneg, x2(useOtherNeg), y2(useOtherNeg), method, extrap); end
                if any(useOtherPos), z(useOtherPos) = evalInterpolant(Fpos, x2(useOtherPos), y2(useOtherPos), method, extrap); end
            case 'nearest'
                % Pure nearest across all samples (ignores barrier only for these)
                if ~isempty(X)
                    idx = find(nanQ);
                    for k = idx(:).'
                        [~,ii] = min( hypot(X - x2(k), Y - y2(k)) );
                        z(k) = Z(ii);
                    end
                end
            case 'nan'
                % leave as NaN
        end
    end

    z2 = reshape(z, size(x2));
    meta.sideQuery = reshape(sign(sideQ), size(x2));
    meta.distQuery = reshape(distQ, size(x2));
end

% ====== Helpers ======

function [side, dist] = classifySideToPolyline(px, py, polyXY)
% Returns side (+1/-1/0) and shortest distance to the polyline for points P.
% Side is computed w.r.t. the orientation of the nearest segment.
    x1 = polyXY(1:end-1,1); y1 = polyXY(1:end-1,2);
    x2 = polyXY(2:end,1);   y2 = polyXY(2:end,2);
    dx = x2 - x1; dy = y2 - y1;
    L2 = dx.^2 + dy.^2;

    nP = numel(px);
    side = zeros(nP,1);
    dist = inf(nP,1);

    % Chunk to limit memory
    blk = max(1, floor(2e5 / max(1, numel(x1)))); % heuristic
    for s = 1:blk:nP
        e = min(nP, s+blk-1);
        Pxs = px(s:e); Pys = py(s:e);

        % Compute projection t for all segments (broadcast over points)
        % t = dot(P-A, v)/|v|^2, clamped to [0,1]
        % Vectorized: points along rows, segments along cols
        PAx = Pxs - x1.';  PAy = Pys - y1.';
        t = (PAx .* dx.' + PAy .* dy.') ./ L2.';
        t = max(0, min(1, t));

        % Nearest point Q on each segment to each point
        Qx = x1.' + t .* dx.';
        Qy = y1.' + t .* dy.';

        % Distances to each segment's nearest point
        d2 = (Pxs - Qx).^2 + (Pys - Qy).^2;

        % Pick best segment
        [d2min, j] = min(d2, [], 2);

        dist(s:e) = sqrt(d2min);

        % Signed side relative to chosen segment j
        jj = sub2ind(size(Qx), (1:numel(Pxs)).', j);
        % Segment vectors at j
        % vx = dx(j).'; vy = dy(j).';
        % ax = x1(j).'; ay = y1(j).';
        % Cross((B-A),(P-A))
        % --- NEW (correct; column vectors matching Pxs/Pys):
        vx = dx(j);
        vy = dy(j);
        ax = x1(j);
        ay = y1(j);

        % Signed side relative to chosen segment j (now column-by-column OK)
        sgn = vx .* (Pys - ay) - vy .* (Pxs - ax);
        side(s:e) = sign(sgn);
    end
end

function F = makeInterpolant(x, y, z, method, extrap)
% Construct per-side interpolant.
    switch method
        case {'natural','linear','nearest'}
            F = scatteredInterpolant(x, y, z, method, 'nearest');
            if ~strcmp(extrap,'none'), F.ExtrapolationMethod = extrap; end
        case 'v4'
            % Marker for v4 path (use griddata at eval time)
            F = struct('x',x,'y',y,'z',z,'method','v4');
    end
end

function z = evalInterpolant(F, xq, yq, method, extrap)
% Evaluate either scatteredInterpolant or the 'v4' path.
    switch method
        case {'natural','linear','nearest'}
            z = F(xq, yq);
            if strcmp(extrap,'none')
                % Mask outside convex hull of this side
                try
                    K = convhull(F.Points(:,1), F.Points(:,2));
                    in = inpolygon(xq, yq, F.Points(K,1), F.Points(K,2));
                    z(~in) = NaN;
                catch
                end
            end
        case 'v4'
            z = griddata(F.x, F.y, F.z, xq, yq, 'v4');
    end
end

function XYu = accumfirst(XY, g)
% Return first occurrence of each group g for matrix XY (cols preserved).
    n = max(g); XYu = zeros(n, size(XY,2));
    first = false(n,1);
    for i = 1:numel(g)
        gi = g(i);
        if ~first(gi), XYu(gi,:) = XY(i,:); first(gi) = true; end
    end
end
