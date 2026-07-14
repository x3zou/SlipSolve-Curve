function Gnew = resizeG(G,n)
% RESIZEGTONEWROWS   Pad/shift an (a × 2m) matrix G so it matches a
%                    slip matrix that has grown from m → n rows
%                    by prepending t = n-m new triangle rows.
%
% INPUT
%   G   :  a × 2m matrix   (columns 1…m  and m+1…2m map to old rows 1…m)
%   n   :  total rows in the *new* slip matrix  (n ≥ m)
%
% OUTPUT
%   Gnew:  a × 2n matrix
%          • columns 1…t    (and n+1…n+t) are zero   [triangles]
%          • remaining columns keep the original data in order
%
% EXAMPLE
%   Gnew = resizeGtoNewRows(G, nRowsNew);
%
% 2025-06-08  ChatGPT-o3
% -------------------------------------------------------------------------

% --- basic sizes ----------------------------------------------------------
[a,twom] = size(G);
m        = twom/2;
assert(twom == 2*m, 'G must have an even number of columns (2m).');
t        = n - m;                          % number of prepended triangles
assert(t >= 0, 'New row count n must not be smaller than m.');

% --- allocate and copy ----------------------------------------------------
Gnew = zeros(a, 2*n);                      % triangles already zero

% first half
Gnew(:, t+(1:m)) = G(:, 1:m);

% second half
Gnew(:, n+t+(1:m)) = G(:, m+1:2*m);
end
