function z2 = interpXYZ(XYZ, x2, y2, method, extrap)
%INTERPXYZ  Interpolate Z at new (x2, y2) from 3-column (X,Y,Z) data.
%
%   z2 = INTERPXYZ(XYZ, x2, y2)
%   z2 = INTERPXYZ(XYZ, x2, y2, method)
%   z2 = INTERPXYZ(XYZ, x2, y2, method, extrap)
%
% Inputs
%   XYZ     : N×3 numeric matrix, columns = [X, Y, Z].
%   x2,y2   : Query coordinates (size must match each other). Can be vectors
%             or 2D grids (e.g., from MESHGRID). z2 will match this shape.
%   method  : (optional) 'linear' (default), 'natural', or 'nearest'.
%   extrap  : (optional) Extrapolation method for points outside convex hull:
%             'none' (default -> returns NaN), 'nearest', or 'linear'.
%
% Output
%   z2      : Interpolated values at (x2, y2), same size as x2/y2.
%
% Notes
%   - Uses scatteredInterpolant (works for scattered or gridded samples).
%   - Duplicate (X,Y) input rows are consolidated by averaging Z.
%
% Examples
%   % Sample scattered data:
%   xy = rand(1000,2)*10; z = peaks(xy(:,1)/3, xy(:,2)/3);
%   XYZ = [xy, z];
%   [Xq,Yq] = meshgrid(0:0.1:10, 0:0.1:10);
%   Zq = interpXYZ(XYZ, Xq, Yq, 'natural', 'none');
%
%   % 1D query vectors (returns column vector):
%   xq = [1; 2; 3]; yq = [4; 5; 6];
%   zq = interpXYZ(XYZ, xq, yq);
%
% Xiaoyu-ready: safe defaults, shape-preserving output, fast.

    % ---- Defaults & validation
    if nargin < 4 || isempty(method), method = 'linear'; end
    if nargin < 5 || isempty(extrap), extrap = 'none';   end

    validateattributes(XYZ, {'numeric'}, {'2d','ncols',3}, mfilename, 'XYZ', 1);
    if ~isequal(size(x2), size(y2))
        error('interpXYZ:SizeMismatch', 'x2 and y2 must have the same size.');
    end
    validMethod  = validatestring(method,  {'linear','natural','nearest'}, mfilename, 'method', 4);
    validExtrap  = validatestring(extrap,  {'none','nearest','linear'},     mfilename, 'extrap', 5);

    % ---- Clean input data: drop NaNs/Inf and consolidate duplicates
    X = XYZ(:,1); Y = XYZ(:,2); Z = XYZ(:,3);
    good = isfinite(X) & isfinite(Y) & isfinite(Z);
    X = X(good); Y = Y(good); Z = Z(good);

    if isempty(X)
        z2 = nan(size(x2)); return;
    end

    % Consolidate duplicate (X,Y) by averaging Z
    XY = [X, Y];
    [~, ~, g] = unique(XY, 'rows', 'stable');
    % Count per unique, sum Z per unique, and average
    nPer  = accumarray(g, 1);
    zSum  = accumarray(g, Z);
    Zu    = zSum ./ nPer;
    XuYu  = XY(accumindexfirst(g), :); % first occurrence per group
    Xu    = XuYu(:,1);
    Yu    = XuYu(:,2);

    % ---- Build interpolant
    F = scatteredInterpolant(double(Xu), double(Yu), double(Zu), validMethod, 'nearest'); % temp extrap
    switch validExtrap
        case 'none'
            % Emulate 'none' by marking outside convex hull as NaN after eval
            doMaskOutside = true;
            F.ExtrapolationMethod = 'nearest'; % compute then mask
        otherwise
            doMaskOutside = false;
            F.ExtrapolationMethod = validExtrap;
    end

    % ---- Evaluate at query points (preserve shape)
    sz = size(x2);
    zvec = F(x2(:), y2(:));

    % If requested 'none', mask points outside the convex hull to NaN
    if doMaskOutside
        % Detect outside via barycentric test using convex hull of input
        try
            K = convhull(Xu, Yu);
            in = inpolygon(x2(:), y2(:), Xu(K), Yu(K));
            zvec(~in) = NaN;
        catch
            % Fallback: trust scatteredInterpolant for most cases
            % (rarely convhull may fail for degenerate data)
        end
    end

    z2 = reshape(zvec, sz);
end

% --- Helper: get index of first occurrence for each group ---
function idxFirst = accumindexfirst(g)
    % Returns the index of the first occurrence for each unique group id in g
    % g must be a column vector of positive integers (from unique(...,'stable')).
    idxFirst = zeros(max(g),1);
    % The first time we see a group, record its index
    seen = false(max(g),1);
    for i = 1:numel(g)
        gi = g(i);
        if ~seen(gi)
            idxFirst(gi) = i;
            seen(gi) = true;
            if all(seen), break; end
        end
    end
end
