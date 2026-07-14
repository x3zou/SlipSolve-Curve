function [Ge, Gn, Gp] = comb_green_nodes_xyz(indx, Conn_List, PTs, data_insar, F_strk, F_dip, green_name, varargin)
% Fast: avoid iterative concatenation; use single vertcat at the end.

Vx = PTs(:,1); Vy = PTs(:,2); Vd = PTs(:,3);

% Find all triangles containing the node, and the node's position within each triangle
[tri_rows, node_cols] = find(Conn_List == indx);   % tri_rows: triangle indices; node_cols: 1..3 position
nTri = numel(tri_rows);

% Pre-size a cell to hold each triangle's point-source rows; concat once at end
chunks = cell(nTri,1);

% If your interpolate_triangle2 uses a fixed sample count passed in, keep it here
% (kept at 15 to match your original call); you can expose it via varargin if needed.
npts = 15;

for k = 1:nTri
    tri = Conn_List(tri_rows(k),:);        % three vertex indices
    node_pos = node_cols(k);               % 1, 2, or 3 (position within tri)
    
    % Two "other" vertices without calling setdiff (avoid extra allocs)
    other = tri; 
    other(node_pos) = [];                  % now 1x2

    % Vertex coordinates (x,y,z)
    v_node = [Vx(tri(node_pos)), Vy(tri(node_pos)), Vd(tri(node_pos))];
    v2     = [Vx(other(1)),      Vy(other(1)),      Vd(other(1))];
    v3     = [Vx(other(2)),      Vy(other(2)),      Vd(other(2))];

    % Sample the triangle (returns npts samples because we pass 15)
    [tmpX, tmpY, tmpZ] = interpolate_triangle2(v_node, v2, v3, npts);
    npatch = numel(tmpX);

    % Triangle area & tapered weights
    area    = triangle_area(v_node, v2, v3);
    pmoment = (area / npatch) * 1e6;                % per-sample scalar
    w       = taperTriangle(v_node, v2, v3, [tmpX, tmpY, tmpZ]);

    % Strike & dip at sample points (vectorized if F_* are griddedInterpolant/ scatteredInterpolant)
    tmpStrk = F_strk(tmpX, tmpY, tmpZ);
    tmpDip  = F_dip (tmpX, tmpY, tmpZ);

    % Assemble this triangle's point-source rows (npatch x 11)
    % Columns: [?, ?, ?, x[m], y[m], z[m], 0, 0, strike, dip, moment]
    % (3 leading cols left at 0 to match your original layout)
    T = zeros(npatch, 11);
    T(:,4)  = tmpX * 1e3;
    T(:,5)  = tmpY * 1e3;
    T(:,6)  = tmpZ * 1e3;
    % T(:,7:8) already zero
    T(:,9)  = tmpStrk;
    T(:,10) = tmpDip;
    T(:,11) = pmoment .* w;                 % no extra ones()

    chunks{k} = T;                          % stash without concatenating yet
end

% Single concatenation
TotalModel = vertcat(chunks{:});

% Green's function evaluation
[Ge, Gn, Gp] = calc_green_gps_3d_edcmp_xyz(TotalModel, data_insar, green_name);
end
