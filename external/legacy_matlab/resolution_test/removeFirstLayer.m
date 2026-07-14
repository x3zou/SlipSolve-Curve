function Mout = removeFirstLayer(M)
% REMOVEFIRSTLAYER  Drop layer-1 rows and trim to 17 columns.
%
%   Mout = removeFirstLayer(M)
%
% INPUT
%   M      n×27 matrix, columns:
%     1    fault id
%     2    element id
%     3    layer_id (1 = first layer, 2…)
%     4–5  strike-slip, dip-slip
%     6–8  x, y, z
%     9–16 connectivity (row indices, NaN padded)
%     17–18 strike, dip
%     19–27 triangle vertex coords (only layer 1)
%
% OUTPUT
%   Mout   m×17 matrix, columns:
%     1    fault id
%     2    new element id (1…m)
%     3    new layer_id (1…)
%     4–5  strike-slip, dip-slip
%     6–8  x, y, z
%     9–16 connectivity (updated to new row IDs)
%     17   strike
%     18   dip
%
% All layer-1 rows are removed; layer_id decremented; connectivity
% remapped; columns 16 and 19–27 dropped.
%
% 2025-06-09  ChatGPT-o3

% 1. select rows with layer_id ≥ 2
keepMask = M(:,3) > 1;
oldIdxs  = find(keepMask);
nNew     = numel(oldIdxs);

% 2. build map from old row → new row (NaN for dropped)
map = nan(size(M,1),1);
map(oldIdxs) = 1:nNew;

% 3. extract kept rows
K = M(keepMask, :);   % nNew×27

% 4. update element ID (col 2) to new row numbers
K(:,2) = (1:nNew).';

% 5. decrement layer_id (col 3)
K(:,3) = K(:,3) - 1;

% 6. remap connectivity in cols 9–16
for c = 9:16
    conn = K(:,c);
    valid = conn>=1 & conn<=numel(map);
    conn(valid) = map(conn(valid));
    conn(~valid) = NaN;
    K(:,c) = conn;
end

% 7. drop columns 16 and 19–27 to leave 17 columns:
%    keep cols 1–15, then 17–18
colsToKeep = [1:15, 17:18];
Mout = K(:, colsToKeep);

end
