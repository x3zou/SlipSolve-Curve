function [xOut, yOut, zOut] = expandGridNaN(x, y, z, padX, padY)
%EXPANDGRIDNAN  Pad a 2-D regular grid with NaNs.
%
%   [xOut,yOut,zOut] = EXPANDGRIDNAN(x,y,z,padX,padY)
%   -------------------------------------------------
%   x      : (n×1) or (1×n) monotonically increasing OR decreasing vector
%   y      : (m×1) or (1×m) monotonically increasing OR decreasing vector
%   z      : (m×n) matrix, z(i,j) ↔ y(i), x(j)
%   padX   : scalar → pad that many columns on BOTH sides
%            1×2 vec [pre post] → pad pre columns left, post columns right
%            0: no padding
%   padY   : 同上，行方向 padding 数
%
%   All added cells in zOut are NaN.  xOut / yOut are extrapolated using
%   the original grid spacing (assumed uniform).

    %% --- orientation & sanity checks ---
    if isrow(x), x = x.'; end
    if isrow(y), y = y.'; end
    [m, n] = size(z);
    assert(numel(x)==n && numel(y)==m, ...
        'Size mismatch: z must be (m × n) with length(x)=n, length(y)=m');

    if isscalar(padX), padX = [padX padX]; end
    if isscalar(padY), padY = [padY padY]; end
    leftX  = padX(1);  rightX = padX(2);
    topY   = padY(1);  botY   = padY(2);

    %% --- extrapolate coordinate vectors ---
    if n > 1
        dx = x(2) - x(1);
    else
        dx = 1;                     % fallback when only one column
    end

    if m > 1
        dy = y(2) - y(1);
    else
        dy = 1;                     % fallback when only one row
    end

    xOut = [x(1) - (leftX:-1:1).' * dx;  x;  x(end) + (1:rightX).' * dx];
    yOut = [y(1) - (topY :-1:1).' * dy;  y;  y(end) + (1:botY ).' * dy];

    %% --- allocate & insert ---
    zOut = NaN(numel(yOut), numel(xOut), class(z));   % preserve numeric type
    rowIdx = topY + (1:m);
    colIdx = leftX + (1:n);
    zOut(rowIdx, colIdx) = z;
end
