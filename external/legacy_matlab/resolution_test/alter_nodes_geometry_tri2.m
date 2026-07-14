function [Mout, patches] = alter_nodes_geometry_tri2(M)
% TRIANGLELAYERMATRIXFLEX
% -------------------------------------------------------------------------
% Convert a 17-column vertex matrix M into a mixed triangle+vertex matrix
% whose width grows automatically to hold all connectivity indices.
% Also outputs a mapping `patches` that links original layer-1 node rows to
% the new triangle row indices in Mout.
%
% USAGE:
%   [Mout, patches] = alter_nodes_geometry_tri(M)
%   patches{i} = vector of Mout-row indices of triangles that originated
%                from original M-row i (layer_id==1).
%
% 2025-06-08  ChatGPT-o3

%% 0. Basic bookkeeping
n           = size(M,1);
layerOld    = M(:,3);
toDelete    = (layerOld == 1);  % layer-1 vertices will be removed
toKeep      = ~toDelete;

%% 1. Build adjacency from cols 9-15
conn        = M(:,9:15);
fromIdx     = repelem((1:n).',1,7); fromIdx = fromIdx(:);
toIdx       = conn(:);
mask        = ~isnan(toIdx);
A           = sparse(fromIdx(mask), toIdx(mask), true, n, n);
A           = A | A.';

%% 2. Enumerate triangles involving deleted vertices
triSet = [];
delList = find(toDelete);
for a = delList.'
    nbr = find(A(a,:)); if numel(nbr)<2, continue; end
    pairs = nchoosek(nbr,2);
    keep  = A(sub2ind([n n],pairs(:,1),pairs(:,2)));
    triSet = [triSet; sort([repmat(a,sum(keep),1) pairs(keep,:)],2)];
end
triSet = unique(triSet,'rows');
nt     = size(triSet,1);

%% 3. Map original row -> new row (0 for dropped)
keepList = find(toKeep);
nKeep    = numel(keepList);
newRow   = zeros(n,1);
newRow(keepList) = nt + (1:nKeep);

%% 4. Precompute full connectivity
triConn  = cell(nt,1);
for t=1:nt
    v   = triSet(t,:);
    % vertices in new numbering
    v2  = v(layerOld(v)==2);
    triNbr = find(any(ismember(triSet,v),2) & (1:nt).'~=t);
    triConn{t} = [newRow(v2).', triNbr.'];
end
vertConn = cell(nKeep,1);
for k=1:nKeep
    idx       = keepList(k);
    v2v       = conn(idx,~isnan(conn(idx,:)));
    v2v       = v2v(~toDelete(v2v));
    v2vNew    = newRow(v2v); v2vNew = v2vNew(v2vNew>0);
    triTouch  = find(any(triSet==idx,2));
    vertConn{k} = [v2vNew(:).', triTouch(:).'];
end

%% 5. Determine width
maxConn   = max(7, max([cellfun(@numel,triConn); cellfun(@numel,vertConn)]));
coordCol1 = 9 + maxConn + 2;
nCols     = coordCol1 + 9 - 1;

%% 6. Build triangle rows & record patches
triRows = NaN(nt,nCols);
triRows(:,1) = 1;
triRows(:,2) = (1:nt).';
triRows(:,3) = 1;
triRows(:,4:5)= 0;
triRows(:,6:8)= NaN;
patches = cell(n,1);
for t=1:nt
    vIdx = triSet(t,:);
    % connectivity
    triRows(t,9:8+numel(triConn{t})) = triConn{t};
    % misc-1, misc-2
    triRows(t,9+maxConn)   = mean(M(vIdx,16),'omitnan');
    triRows(t,10+maxConn)  = mean(M(vIdx,17),'omitnan');
    % coords
    triRows(t,coordCol1:coordCol1+8) = reshape(M(vIdx,6:8).',1,9);
    % record mapping for original deleted vertices
    for a = vIdx(:).'  % original row index
        if toDelete(a)
            patches{a}(end+1) = t;
        end
    end
end

%% 7. Build vertex rows
vertRows = NaN(nKeep,nCols);
vertRows(:,1)=1;
vertRows(:,2)= nt + (1:nKeep).';
vertRows(:,3)= M(keepList,3);
vertRows(:,4:5)=0;
vertRows(:,6:8)= M(keepList,6:8);
for k=1:nKeep
    vertRows(k,9:8+numel(vertConn{k})) = vertConn{k};
end
vertRows(:,9+maxConn)  = M(keepList,16);
vertRows(:,10+maxConn) = M(keepList,17);

%% 8. Assemble and finalize
Mout = [triRows; vertRows];
Mout(:,2) = (1:size(Mout,1)).';
end
