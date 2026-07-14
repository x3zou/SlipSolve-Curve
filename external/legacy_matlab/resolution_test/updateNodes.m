function M1_new = updateNodes(M1, M2, M3, patches)
% UPDATENODES   Update M1’s slip values using M2 triangles and M3 nodal slips.
%
%   M1_new = updateNodes(M1, M2, M3, patches)
%
% INPUT
%   M1       n×17 nodal model before triangle insertion
%            Col-3    layer_id (1=surface, ≥2 deeper)
%            Col-4,5  strike-/dip-slip (cm)
%
%   M2       p×27 mixed model with layer-1 triangles
%            Col-3    layer_id (1=triangle, 2=node)
%            Col-4,5  triangle strike-/dip-slip (cm)
%
%   M3       m×17 nodal model after removeFirstLayer(M2)
%            Rows correspond in order to M1 rows with layer_id>=2
%            Col-4,5  strike-/dip-slip for deeper nodes (cm)
%
%   patches  n×1 cell array from alter_nodes_geometry_tri2
%            patches{i} is a vector of M2-row indices of all layer-1
%            triangles that replaced surface node i.
%
% OUTPUT
%   M1_new   n×17 same as M1, but:
%     • For surface nodes (M1(:,3)==1): 
%         Slip = mean of M2’s updated triangle slips over patches{i}.
%     • For deeper nodes: slip copied from M3.
%
% DEPENDS ON
%   updateTriangleSlip.m
%
% -------------------------------------------------------------------------

% Start with original
M1_new = M1;

% 1. Overwrite triangle slips in M2 so they match their assigned deeper-node slips
M2_new = updateTriangleSlip(M2, M3);

% 2. Surface nodes: use patches mapping to average their triangles’ slips
surf = find(M1(:,3)==1);
for i = surf.'
    triRows = patches{i};
    if isempty(triRows)
        warning('Surface node %d: no assigned triangles.', i);
        continue
    end
    % mean strike and dip over those triangles
    M1_new(i,4) = mean(M2_new(triRows,4));
    M1_new(i,5) = mean(M2_new(triRows,5));
end

% 3. Deeper nodes: copy directly from M3
deepIdx = find(M1(:,3)>=2);
% M3 rows correspond one-to-one to deepIdx in order
M1_new(deepIdx,4:5) = M3(:,4:5);

end
