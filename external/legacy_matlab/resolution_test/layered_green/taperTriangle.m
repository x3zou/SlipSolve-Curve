function w = taperTriangle(A,B,C,P)
% taperTriangle  Compute linear taper inside triangle ABC.
%  A,B,C : 1×3 vectors for triangle vertices (double)
%  P     : n×3 matrix of query points inside the triangle
%  w     : n×1 taper weights, 1 at A, 0 at B and C
%
%  线性权重 w = 1 在顶点 A，沿 A‑B、A‑C 方向线性减到 0

% 边向量 / edge vectors
v0 = B - A;          % A→B
v1 = C - A;          % A→C

% 预计算内积 / pre‑compute dot products
d00 = dot(v0,v0);
d01 = dot(v0,v1);
d11 = dot(v1,v1);
den = d00*d11 - d01*d01;   % 行列式 / denominator

% 向量化处理全部点 / vectorised for all P
v2  = P - A;               % n×3, A→P
d20 = v2*v0.';             % n×1
d21 = v2*v1.';             % n×1

u = (d11.*d20 - d01.*d21) / den;   % λ_B
v = (d00.*d21 - d01.*d20) / den;   % λ_C
w = 1 - u - v;                     % λ_A = taper weight

% 数值微小误差修正 / clamp tiny numerical drift
w = max(min(w,1),0);
end