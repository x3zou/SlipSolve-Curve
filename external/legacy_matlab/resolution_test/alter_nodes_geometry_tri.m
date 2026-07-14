function Mout = alter_nodes_geometry_tri(M)
% TRIANGLELAYERMATRIXFLEX
% -------------------------------------------------------------------------
% Convert a 17-column vertex matrix into a mixed triangle+vertex matrix
% whose width grows automatically to hold *all* connectivity indices.
%
% Key rules
% ---------
% 1.  All original layer-1 vertices are removed.
% 2.  Their triangles become **layer 1** rows and are placed *first*.
%     • Column-16 and -17 of a triangle row are the arithmetic mean of
%       the three source vertices’ column-16/17 values (NaNs ignored).
% 3.  Remaining vertices (layer ≥ 2) follow the triangles, preserving
%     their original layer_id and their own column-16/17 values.
% 4.  The connectivity list (starting at column 9) is expanded so that
%     every row can store its full set of neighbours:
%        – For triangles: all layer-2 vertices *plus* all other triangles
%          sharing at least one vertex.
%        – For vertices:  all still-existing neighbouring vertices *plus*
%          all triangles that contain the vertex.
% 5.  Columns
%       1              constant 1
%       2              running row number (1 … nRows)
%       3              layer_id  (1 = triangles, 2+ = vertices)
%       4–5            zeros
%       6–8            x y z (vertices) / NaN (triangles)
%       9 … 8+K        connectivity indices   (K = max connections ≥ 7)
%       9+K            misc-1  (avg / original)
%       10+K           misc-2  (avg / original)
%       11+K … 19+K    x1 y1 z1 x2 y2 z2 x3 y3 z3  (triangles only)
%
% The final matrix therefore has 19+K columns where K is determined at
% runtime.  All unused cells are filled with NaN.
%
% 2025-06-08

% -------------------------------------------------------------------------
% 0.  Basic bookkeeping
% -------------------------------------------------------------------------
n           = size(M,1);
layerOld    = M(:,3);
toDelete    = (layerOld == 1);          % layer-1 vertices will be removed
toKeep      = ~toDelete;

% -------------------------------------------------------------------------
% 1.  Build an undirected adjacency matrix from columns 9–15
% -------------------------------------------------------------------------
conn        = M(:,9:15);
fromIdx     = repelem((1:n).',1,7);      fromIdx = fromIdx(:);
toIdx       = conn(:);
mask        = ~isnan(toIdx);
A           = sparse(fromIdx(mask), toIdx(mask), true, n, n);
A           = A | A.';                   % make it undirected

% -------------------------------------------------------------------------
% 2.  Enumerate triangles that involve at least one deleted vertex
% -------------------------------------------------------------------------
triSet = [];
delList = find(toDelete);
for a = delList.'
    nbr = find(A(a,:));
    if numel(nbr) < 2,  continue;  end
    pairs = nchoosek(nbr,2);
    keep  = A(sub2ind([n n], pairs(:,1), pairs(:,2)));
    triSet = [triSet;
              sort([repmat(a,sum(keep),1) pairs(keep,:)],2)]; %#ok<AGROW>
end
triSet = unique(triSet,'rows');
nt     = size(triSet,1);               % number of triangle rows

% -------------------------------------------------------------------------
% 3.  Map original row index  ->  new row number
% -------------------------------------------------------------------------
keepList  = find(toKeep);
nKeep     = numel(keepList);
newRow    = zeros(n,1);
newRow(keepList) = nt + (1:nKeep);     % vertices start after triangles

% -------------------------------------------------------------------------
% 4.  Pre-compute full connectivity lists
% -------------------------------------------------------------------------
triConn  = cell(nt,1);
for t = 1:nt
    v   = triSet(t,:);
    v2  = v(layerOld(v)==2);                 % only layer-2 vertices
    triNbr = find(any(ismember(triSet,v),2) & (1:nt).'~=t);
    triConn{t} = [newRow(v2).', triNbr.'];   % vertices first, triangles next
end

vertConn = cell(nKeep,1);
for k = 1:nKeep
    idx        = keepList(k);
    % vertex-to-vertex (exclude deleted vertices)
    v2v        = conn(idx,~isnan(conn(idx,:)));
    v2v        = v2v(~toDelete(v2v));
    v2vNew     = newRow(v2v);  v2vNew = v2vNew(v2vNew>0);
    % vertex-to-triangle
    triTouch   = find(any(triSet==idx,2));
    vertConn{k}= [v2vNew(:).', triTouch(:).'];
end

% -------------------------------------------------------------------------
% 5.  Determine required width
% -------------------------------------------------------------------------
maxConn   = max(7, max([cellfun(@numel,triConn); cellfun(@numel,vertConn)]));
coordCol1 = 9 + maxConn + 2;           % first coordinate column
nCols     = coordCol1 + 9 - 1;         % total columns = 19+maxConn

% -------------------------------------------------------------------------
% 6.  Build triangle rows
% -------------------------------------------------------------------------
triRows              = NaN(nt, nCols);
triRows(:,1)         = 1;              % column-1
triRows(:,2)         = (1:nt).';       % running row numbers
triRows(:,3)         = 1;              % layer_id = 1
triRows(:,4:5)       = 0;              % zeros
triRows(:,6:8)       = NaN;            % triangles have no single xyz

for t = 1:nt
    vIdx = triSet(t,:);
    % connectivity
    triRows(t,9:8+numel(triConn{t})) = triConn{t};
    % misc columns = vertex average
    triRows(t,9+maxConn)  = mean(M(vIdx,16),'omitnan');
    triRows(t,10+maxConn) = mean(M(vIdx,17),'omitnan');
    % triangle vertex coordinates
    triRows(t,coordCol1:coordCol1+8) = reshape(M(vIdx,6:8).',1,9);
end

% -------------------------------------------------------------------------
% 7.  Build vertex rows (layer ≥ 2)
% -------------------------------------------------------------------------
vertRows              = NaN(nKeep, nCols);
origK                 = M(keepList,:);

vertRows(:,1)         = 1;
vertRows(:,2)         = nt + (1:nKeep).';        % row numbers
vertRows(:,3)         = origK(:,3);              % layer_id
vertRows(:,4:5)       = 0;
vertRows(:,6:8)       = origK(:,6:8);            % xyz

for k = 1:nKeep
    vertRows(k,9:8+numel(vertConn{k})) = vertConn{k};
end
vertRows(:,9+maxConn)  = origK(:,16);            % misc-1
vertRows(:,10+maxConn) = origK(:,17);            % misc-2
% coordinate columns stay NaN for vertices (already default)

% -------------------------------------------------------------------------
% 8.  Assemble final matrix and re-index column-2
% -------------------------------------------------------------------------
Mout           = [triRows; vertRows];
Mout(:,2)      = (1:size(Mout,1)).';   % ensure monotone row numbers
end
