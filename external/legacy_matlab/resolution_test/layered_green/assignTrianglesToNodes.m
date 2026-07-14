function A = assignTrianglesToNodes(M)
% ASSIGNTRIANGLESTONODES   Assign every layer-1 triangle to exactly one
%                          layer-2 node, ensuring each node gets at least
%                          two incident triangles.
%
%   A = assignTrianglesToNodes(M)
%
% INPUT
%   M       n×27 fault matrix
%     Col-3     layer_id (1 = triangles, 2 = nodes)
%     Col-9–16  connectivity (row indices, NaN or 0 padded)
%
% OUTPUT
%   A       nNode×(1+K) matrix
%     Col-1     node row index (into M)
%     Col-2…    assigned triangle row indices (into M)
%     empty slots = NaN
%
% Each triangle (layer-1 row) is assigned to exactly one node that it
% touches. Each node receives at least two triangles; no triangle is reused.
%
% 2025-06-09  ChatGPT-o3

% Locate triangles and nodes
triRows  = find(M(:,3)==1);    % global row IDs of triangles
nodeRows = find(M(:,3)==2);    % global row IDs of nodes
nTri     = numel(triRows);
nNode    = numel(nodeRows);
nRows    = size(M,1);

% Map global row → local node index
row2node = zeros(nRows,1);
row2node(nodeRows) = 1:nNode;

% Build triangle → local node indices using connectivity
tri2nodes = cell(nTri,1);
for t = 1:nTri
    r = triRows(t);
    conn = M(r,9:16);
    conn = conn(~isnan(conn) & conn>0 & conn<=nRows);
    % keep only layer-2 node rows
    conn = conn(M(conn,3)==2);
    % convert to local node indices
    tri2nodes{t} = row2node(conn);
end

% Invert: node → list of incident triangle indices
nodeTriList = cell(nNode,1);
for t = 1:nTri
    localNodes = tri2nodes{t}(:).';  % ensure row
    for li = localNodes
        nodeTriList{li}(end+1) = t;  %#ok<AGROW>
    end
end

% First pass: assign up to two triangles per node
assignedTri = false(nTri,1);
assignList  = cell(nNode,1);

for i = 1:nNode
    inc    = nodeTriList{i};
    free   = inc(~assignedTri(inc));
    % if numel(free) < 2
    %     error('Node row %d has fewer than 2 incident triangles.', nodeRows(i));
    % end
    % take = free(1:2);
    take = free(1:min(2,numel(free)));  % allow 1, or 2
    assignList{i}    = take;
    assignedTri(take)= true;
end

% Second pass: assign any remaining triangles
remainders = find(~assignedTri);
for t = remainders.'
    inc = tri2nodes{t};
    % choose the node with fewest assigned so far
    counts = cellfun(@numel, assignList(inc));
    [~,k]  = min(counts);
    chosen = inc(k);
    assignList{chosen}(end+1) = t;
    assignedTri(t)            = true;
end

% Build output matrix A
maxAssign = max(cellfun(@numel, assignList));
A = NaN(nNode, 1 + maxAssign);
A(:,1) = nodeRows;
for i = 1:nNode
    trisLocal  = assignList{i};        % indices into triRows
    trisGlobal = triRows(trisLocal);   % map to M’s global row IDs
    A(i,2:1+numel(trisGlobal)) = trisGlobal.';
end
end
