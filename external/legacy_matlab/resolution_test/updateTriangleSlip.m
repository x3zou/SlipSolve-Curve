function M2_new = updateTriangleSlip(M2, M3)
% UPDATETRIANGLESLIP   Set triangle slips in M2 to the slip of their
%                      assigned surface node (from M3).
%
%   M2_new = updateTriangleSlip(M2, M3)
%
% INPUT
%   M2    p×27 mixed model with layer‐1 triangles
%     Col-3     layer_id (1=tri, 2=node)
%     Col-4,5   triangle strike/dip slip (cm)
%   M3    m×17 purely‐nodal model after removeFirstLayer(M2)
%     Col-4,5   strike/dip slip for nodes (cm)
%
% OUTPUT
%   M2_new  = M2 but with Col-4,5 of each triangle row set to the
%            slip values of the node it was assigned to via 
%            assignTrianglesToNodes(M2).
%
% REQUIREMENT: assignTrianglesToNodes(M2) must return A where
%   A(:,1) = nodeRow (global M2 row), 
%   A(:,2:end) = triRows (global M2 rows) assigned to that node.
%
% 2025-06-10  ChatGPT-o3

% 1. Get assignment of triangles to nodes
A = assignTrianglesToNodes(M2);
nodeRows = A(:,1);      % M2 row indices of nodes
triCols  = A(:,2:end);  % triangle row indices or NaN

% 2. Build map from M2 row → M3 row (only for node rows)
oldIdxs = find(M2(:,3)>=2);      % rows kept in M3
map = nan(size(M2,1),1);
map(oldIdxs) = 1:numel(oldIdxs);

% 3. Initialize output
M2_new = M2;

% 4. For each node‐triangle group, update triangle slips
for i = 1:size(A,1)
    nRow = nodeRows(i);
    m3row = map(nRow);
    if isnan(m3row)
        error('Node row %d not found in M3 mapping.', nRow);
    end
    % slip of this node from M3
    slipVals = M3(m3row,4:5);  % [strike, dip]
    % assigned triangles
    tRows = triCols(i, :);
    tRows = tRows(~isnan(tRows));
    % update each triangle’s slip
    for t = tRows
        M2_new(t,4:5) = slipVals;
    end
end
end
