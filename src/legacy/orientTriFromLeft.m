function [P1o, P2o, P3o, flipped, sgn] = orientTriFromLeft(P1, P2, P3, varargin)
%ORIENTTRIFROMLEFT Match the triangle orientation helper used by legacy code.

parser = inputParser;
parser.addParameter("desired", "CW", @(s) ischar(s) || isstring(s));
parser.addParameter("tol", 1e-12, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.parse(varargin{:});
desired = upper(string(parser.Results.desired));
tolerance = parser.Results.tol;

P1 = P1(:).';
P2 = P2(:).';
P3 = P3(:).';
normal = [1 0 0];
rawSign = dot(normal, cross(P2 - P1, P3 - P1));
if rawSign > tolerance
    sgn = 1;
elseif rawSign < -tolerance
    sgn = -1;
else
    sgn = 0;
end

flipped = false;
if desired == "CW" && sgn > 0
    [P2, P3] = deal(P3, P2);
    flipped = true;
elseif desired == "CCW" && sgn < 0
    [P2, P3] = deal(P3, P2);
    flipped = true;
end

P1o = P1;
P2o = P2;
P3o = P3;
end
