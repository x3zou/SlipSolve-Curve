function [H, h1, L, M, W, edgeMode] = smoothingMatrix_laplace2(n1, PTs, varargin)
% smoothingMatrix_laplace_edgeHybrid
% Edge-level hybrid Laplace–Beltrami smoother for triangular vertex meshes.
% Per edge (i,j):
%   - use cotangent weight if BOTH incident triangles to edge (i,j) exist
%     and are ACUTE (opposite angles < 90°);
%   - otherwise use symmetric mean-value weight (average of the two
%     one-sided mean-value weights from i->j and j->i).
%
% Inputs
%   n1   : (nv x K) neighbor list (NaN-padded OK) OR cell{nv} of neighbors
%   PTs  : (nv x 3) vertex coordinates [x y z]
%
% Name-Value options (optional)
%   'Mass'     : 'voronoi' (default) | 'barycentric'
%   'Form'     : 'curvature' (default) | 'dirichlet'
%                'curvature'  -> H = M^{-1/2} * (D - W)   (penalize curvature)
%                'dirichlet'  -> H = D^{-1/2} * (D - W)  (neighbor differences)
%   'ss_ratio' : scalar (default 1)  strike-slip block scale
%   'ds_ratio' : scalar (default 1)  dip-slip block scale
%
% Outputs
%   H        : (2*nv x 2*nv) sparse smoothing matrix (block-diagonal)
%   h1       : number of rows (= 2*nv)
%   L        : Laplacian = D - W (symmetric)
%   M        : (nv x 1) vertex areas (mixed-Voronoi or barycentric)
%   W        : symmetric edge-weight matrix
%   edgeMode : sparse symmetric matrix with 1 on edges where cotan used
%              and 0 where mean-value used
%
% Example:
%   [H,~,L,M,W,edgeMode] = smoothingMatrix_laplace_edgeHybrid(n1,PTs,'Form','curvature');

% ---------- parse ----------
p = inputParser;
p.addParameter('Mass','voronoi',@(s)ischar(s)||isstring(s));
p.addParameter('Form','curvature',@(s)ischar(s)||isstring(s));
p.addParameter('ss_ratio',1,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('ds_ratio',1,@(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
opts = p.Results;

nv = size(PTs,1);

% ---------- neighbors → cell, symmetrize ----------
nbrs = normalizeNeighbors(n1, nv);
nbrs = makeSymmetric(nbrs, nv);

% Precompute ordered 1-ring and local frames for mean-value parts
[ringOrder, P2D, dist2Nbr] = orderRingsAndProject(nbrs, PTs);

% Build W with edge-level hybrid rule
W = spalloc(nv,nv, 16*nv);
edgeMode = spalloc(nv,nv, 16*nv);   % 1=cotan, 0=mean-value

for i=1:nv
    Ni = ringOrder{i}; ki = numel(Ni);
    for t=1:ki
        j = Ni(t);
        if j<=i, continue; end  % undirected edge handling once

        % Incident triangles: common neighbors of i and j
        Sj = intersect(nbrs{i}, nbrs{j});  % third vertices around edge (i,j)
        % Keep only those k where triangles (i,j,k) plausibly exist:
        % require k connected to both i and j (already true by construction)

        useCotan = false;
        w_ij = 0;

        if numel(Sj) == 2
            % Angles opposite edge (i,j) at the two third vertices
            k1 = Sj(1); k2 = Sj(2);
            alpha = angleAtOpp(PTs, i, j, k1);   % at k1 in tri (i,j,k1)
            beta  = angleAtOpp(PTs, i, j, k2);   % at k2 in tri (i,j,k2)

            if isfinite(alpha) && isfinite(beta) && alpha>0 && beta>0 ...
               && alpha < (pi/2) && beta < (pi/2)
                % acute on both sides -> cotan weight
                w_cotan = 0.5*(cotSafe(alpha) + cotSafe(beta));
                if isfinite(w_cotan) && w_cotan > 0
                    w_ij = w_cotan;
                    useCotan = true;
                end
            end
        end

        if ~useCotan
            % Symmetric mean-value fallback
            w_i = meanValueWeightOneSided(i, j, ringOrder, P2D, dist2Nbr);
            w_j = meanValueWeightOneSided(j, i, ringOrder, P2D, dist2Nbr);
            w_ij = 0.5*(w_i + w_j);
            if ~isfinite(w_ij) || w_ij < 0, w_ij = 0; end
        end

        if w_ij > 0
            W(i,j) = w_ij; W(j,i) = w_ij;
            if useCotan
                edgeMode(i,j)=1; edgeMode(j,i)=1;
            end
        end
    end
end

% Degree and Laplacian
d = full(sum(W,2));
L = spdiags(d,0,nv,nv) - W;

% Areas (mass)
switch lower(opts.Mass)
    case 'voronoi'
        M = mixedVoronoiAreas_safe(nbrs, ringOrder, PTs);
    case 'barycentric'
        M = barycentricAreas_safe(nbrs, ringOrder, PTs);
    otherwise
        error('Unknown Mass option.');
end
M = max(M,1e-15);

% Smoothing operator
switch lower(opts.Form)
    case 'curvature'   % penalize curvature: || M^{-1/2} L m ||^2
        H1 = spdiags(1./sqrt(M),0,nv,nv) * L;
    case 'dirichlet'   % penalize neighbor differences: normalized Laplacian
        invsqrtD = spdiags(1./sqrt(max(d,1e-15)),0,nv,nv);
        H1 = invsqrtD * L;
    otherwise
        error('Unknown Form option.');
end

% Block for strike & dip
% Normalization
% sigma = sqrt(eigs(H1'*H1,1,'largestreal'));   % cheap power iteration
% if isfinite(sigma) && sigma>0
%     H1 = H1 / sigma;
% end
[i,j]  = find(triu(W>0,1));
ell    = sqrt(sum((PTs(i,:)-PTs(j,:)).^2,2));
ell0   = median(ell);                     % characteristic length

switch lower(opts.Form)
    case 'dirichlet'
        H1 = (ell0)   * H1;               % H ~ 1/ell → make unitless
    case 'curvature'
        H1 = (ell0^2) * H1;               % H ~ 1/ell^2 → make unitless
end

Z  = sparse(nv,nv);
H  = [opts.ss_ratio*H1, Z; Z, opts.ds_ratio*H1];
h1 = size(H,1);

% ===================== helpers =====================
function nbrs = normalizeNeighbors(n1, nv)
    if iscell(n1)
        nbrs = n1(:);
    else
        nbrs = cell(nv,1);
        for ii=1:nv
            row = n1(ii,:);
            row = row(~isnan(row) & row~=ii);
            nbrs{ii} = unique(row(:))';
        end
    end
end

function nbrs = makeSymmetric(nbrs, nv)
    A = spalloc(nv,nv, 10*nv);
    for ii=1:nv, A(ii, nbrs{ii}) = 1; end
    A = max(A, A.');
    nbrs = cell(nv,1);
    for ii=1:nv, nbrs{ii} = find(A(ii,:)); end
end

function [ringOrder, P2D, dist2Nbr] = orderRingsAndProject(nbrs, PTs)
    % For mean-value weights we need 1-ring neighbors ordered cyclically
    nv = size(PTs,1);
    ringOrder = nbrs;  % will reorder in place
    P2D = cell(nv,1);
    dist2Nbr = cell(nv,1);
    for ii=1:nv
        N = nbrs{ii};
        k = numel(N);
        if k<1
            ringOrder{ii} = N; P2D{ii} = zeros(0,2); dist2Nbr{ii} = zeros(0,1);
            continue;
        end
        V = PTs(N,:) - PTs(ii,:);
        [e1,e2] = localFrame(V);
        P = [V*e1, V*e2];              % 2D projected rays
        ang = atan2(P(:,2), P(:,1));
        [~,ord] = sort(ang);
        ringOrder{ii} = N(ord);
        P2D{ii} = P(ord,:);
        dist2Nbr{ii} = sqrt(sum(V(ord,:).^2,2));
    end
end

function [e1,e2] = localFrame(V)
    C = (V.'*V);
    [Q,D] = eig(C);
    [~,ord] = sort(diag(D),'descend');
    e1 = Q(:,ord(1)); e1 = e1./max(norm(e1),1e-15);
    e2 = Q(:,ord(2)); e2 = e2./max(norm(e2),1e-15);
    if ~all(isfinite(e1)) || ~all(isfinite(e2)), e1=[1;0;0]; e2=[0;1;0]; end
end

function th = angle2D(u,v)
    nu = hypot(u(1),u(2)); nv2 = hypot(v(1),v(2));
    if nu<1e-15 || nv2<1e-15, th = 0; return; end
    c = max(-1,min(1, dot(u,v)/(nu*nv2) ));
    th = acos(c);
end

function th = angle3D(u,v)
    nu = norm(u); nv2 = norm(v);
    if nu<1e-15 || nv2<1e-15, th = 0; return; end
    c = max(-1,min(1, dot(u,v)/(nu*nv2) ));
    th = acos(c);
end

function c = cotSafe(th)
    s = sin(th); c = cos(th)./max(s,1e-15);
end

function a = angleAtOpp(PTs, i, j, k)
    % angle at vertex k opposite edge (i,j) in triangle (i,j,k)
    a = angle3D(PTs(i,:) - PTs(k,:), PTs(j,:) - PTs(k,:));
end

function w = meanValueWeightOneSided(i, j, ringOrder, P2D, dist2Nbr)
    % mean-value contribution computed at vertex i for neighbor j
    N  = ringOrder{i};
    P  = P2D{i};
    r  = dist2Nbr{i};
    k  = numel(N);
    if k==0, w = 0; return; end
    idx = find(N==j, 1, 'first');
    if isempty(idx)
        w = 0; return;
    end
    ip = mod(idx,   k) + 1;        % next of j around i
    im = mod(idx-2, k) + 1;        % prev of j around i
    th_prev = angle2D(P(im,:), P(idx,:));
    th_next = angle2D(P(idx,:),  P(ip,:));
    dij = max(r(idx), 1e-15);
    w = (tan(th_prev/2) + tan(th_next/2)) / dij;
    if ~isfinite(w), w=0; end
end

function M = mixedVoronoiAreas_safe(nbrs, ringOrder, PTs)
    nv = size(PTs,1);
    M = zeros(nv,1);
    for i=1:nv
        N = ringOrder{i}; k = numel(N);
        if k<2, continue; end
        for t=1:k
            j  = N(t);
            jp = N( mod(t, k) + 1 );
            % triangle exists only if edge (j,jp) present
            if ~hasEdge(j,jp,nbrs), continue; end
            vi = PTs(i,:); vj = PTs(j,:); vk = PTs(jp,:);
            [A, alpha, beta, gamma, b, c] = triGeom(vi,vj,vk); %#ok<ASGLU>
            obt = max([alpha,beta,gamma]) > (pi/2 + 1e-14);
            if ~obt
                M(i) = M(i) + 0.125*(cotSafe(gamma)*(b^2) + cotSafe(beta)*(c^2));
            else
                if alpha > pi/2
                    M(i) = M(i) + 0.5*A;
                else
                    M(i) = M(i) + 0.25*A;
                end
            end
        end
    end
    M = max(M,1e-15);
end

function M = barycentricAreas_safe(nbrs, ringOrder, PTs)
    nv = size(PTs,1);
    M = zeros(nv,1);
    for i=1:nv
        N = ringOrder{i}; k = numel(N);
        if k<2, continue; end
        for t=1:k
            j  = N(t);
            jp = N( mod(t, k) + 1 );
            if ~hasEdge(j,jp,nbrs), continue; end
            A = 0.5*norm(cross(PTs(j,:)-PTs(i,:), PTs(jp,:)-PTs(i,:)));
            M(i) = M(i) + A/3;
        end
    end
    M = max(M,1e-15);
end

function tf = hasEdge(u,v,nbrs)
    Nu = nbrs{u};
    tf = any(Nu==v);
end

function [A, alpha, beta, gamma, b, c] = triGeom(vi, vj, vk)
    eij = vj-vi; eik = vk-vi; ejk = vk-vj;
    b = norm(eik);  % |vi - vk|
    c = norm(eij);  % |vi - vj|
    alpha = angle3D(eij, eik);        % at i
    beta  = angle3D(vi-vj, vk-vj);    % at j
    gamma = angle3D(vi-vk, vj-vk);    % at k
    A = 0.5*norm(cross(eij, eik));
end

end
