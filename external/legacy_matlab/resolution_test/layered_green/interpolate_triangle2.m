function [X, Y, Z] = interpolate_triangle2(v1, v2, v3, n)
%INTERPOLATE_TRIANGLE2
%   Near-uniform points strictly INSIDE a 3D triangle using a triangular
%   lattice (barycentric subdivision). Edge/vertex points are excluded.
%
% INPUTS
%   v1, v2, v3 : 1x3 double, triangle vertices [x y z]
%   n          : integer >= 2, subdivision level (edge is cut into n parts)
%                Interior points exist only when n >= 3.
%
% OUTPUTS
%   X, Y, Z : column vectors of 3D coordinates of interior points
%
% NOTES
%   • Spacing along edges is ~ (edge length)/n, producing a uniform
%     triangular lattice inside the triangle.
%   • #points = (n-1)*(n-2)/2  (0 if n < 3)
%   • To INCLUDE edges instead, change the loops to allow zeros for one of
%     (i,j,k) with i+j+k = n (see comment at end).

    % Ensure row vectors
    v1 = v1(:).'; v2 = v2(:).'; v3 = v3(:).';

    % Degenerate triangle check
    if norm(cross(v2 - v1, v3 - v1)) < 1e-14
        X = []; Y = []; Z = [];
        return;
    end

    % If n < 3 there are no strictly interior lattice points
    if n < 3
        X = []; Y = []; Z = [];
        return;
    end

    % Number of interior points in a triangular lattice with i,j,k >= 1 and i+j+k=n
    nPts = (n-1)*(n-2)/2;
    X = zeros(nPts,1);
    Y = zeros(nPts,1);
    Z = zeros(nPts,1);

    % Fill using barycentric integer triples (i,j,k), i+j+k = n, i,j,k >= 1
    t = 0;
    for i = 1:(n-2)
        % For each i, j can go from 1 up to n-1-i (ensures k >= 1)
        for j = 1:(n-1-i)
            k = n - i - j;  % implied, >= 1 here

            % Barycentric weights
            w1 = i / n;  % weight for v1
            w2 = j / n;  % weight for v2
            w3 = k / n;  % weight for v3

            % Interpolate in 3D directly (stays exactly in the plane)
            P = w1 .* v1 + w2 .* v2 + w3 .* v3;

            t = t + 1;
            X(t) = P(1);
            Y(t) = P(2);
            Z(t) = P(3);
        end
    end

end