function C = getNodeTriangleVertices(M, nodeRow)
% GETNODETRIANGLESCONNECTIVITY
%   C = GETNODETRIANGLESCONNECTIVITY(M, nodeRow)
%
% Given a 27-column fault model M where layer-2 rows are nodes and
% columns 9–16 give each node’s neighbors (row indices, NaN-padded),
% this function finds *all* triangles in the adjacency graph that
% include `nodeRow` as one vertex.  A triangle is any triplet of
% nodes {nodeRow, a, b} such that:
%     • nodeRow connected to a
%     • nodeRow connected to b
%     • a connected to b
%
% INPUT
%   M         n×27 matrix
%     Col-3     layer_id (1=triangles, 2=nodes)
%     Col-9–16  connectivity (row indices, NaN/0 padded)
%   nodeRow   scalar index into M of a node (M(nodeRow,3)>=2)
%
% OUTPUT
%   C         p×2 matrix
%     Each row = [a, b], the two *other* nodeRow indices that,
%     together with nodeRow, form a connectivity triangle.
%     No duplicates; order a<b in each row.
%
% EXAMPLE
%   C = getNodeTrianglesConnectivity(M, 42);
%
% 2025-06-09  ChatGPT-o3

% --- validate input ------------------------------------------------------
n = size(M,1);
if nodeRow<1 || nodeRow>n || M(nodeRow,3)<2
    error('nodeRow must index a node (layer_id >= 2).');
end

% --- get nodeRow’s neighbors --------------------------------------------
nbrs = M(nodeRow,9:16);
nbrs = nbrs(~isnan(nbrs) & nbrs>0 & nbrs<=n);
% keep only layer-2
nbrs = nbrs(M(nbrs,3)>=2);
nbrs = unique(nbrs);

% --- build quick lookup: for each neighbor, its neighbor list in the node set
isNode = M(:,3)>=2;
connBlk = M(:,9:16);

nbrSets = cell(numel(nbrs),1);
for i = 1:numel(nbrs)
    r = nbrs(i);
    c = connBlk(r,:);
    c = c(~isnan(c) & c>0 & c<=n);
    c = c(isNode(c));
    nbrSets{i} = unique(c);
end

% --- find all pairs (i<j) among nbrs that are also connected ------------
C = zeros(0,2);
for i = 1:numel(nbrs)-1
    for j = i+1:numel(nbrs)
        a = nbrs(i);
        b = nbrs(j);
        % check if a and b are neighbors
        % find index of a in nbrs to get its neighbor set
        setA = nbrSets{i};
        if ismember(b, setA)
            C(end+1,:) = sort([a b]); %#ok<AGROW>
        end
    end
end

end
