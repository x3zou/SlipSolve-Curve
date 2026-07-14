function [C, valence] = NodeContact(VID, Points)
%BUILDCONTACTMATRIXAUTO  Contact matrix with auto-sized columns
%   [C, valence] = buildContactMatrixAuto(VID, Points)
%
%   输入 (Inputs)
%   ─────────────────────────────────────────────
%   VID    : (T × 3)   每行存三角形三个顶点的行号 (vertex indices)
%   Points : (N × 3)   顶点坐标，这里只用到 N = size(Points,1)
%
%   输出 (Outputs)
%   ─────────────────────────────────────────────
%   C        : (N × maxVal)  每行列出与该顶点共三角形的邻接顶点；不足列用 NaN
%   valence  : (N × 1)       每个顶点的邻接个数（即度数）
%
%   说明 (Notes)
%   • 函数自动检测最大邻点数 maxVal = max(valence)
%   • 若希望保存“变长”格式，可改用 cell 数组：见下方示例

    % ── Step 0: 基本维度
    N = size(Points,1);

    % ── Step 1: 构造无向邻接稀疏矩阵 A
    %    每个三角形 → 3 条边（AB, BC, AC），再去重
    edges = [VID(:,[1 2]); VID(:,[2 3]); VID(:,[1 3])];
    edges = unique(sort(edges,2),'rows');   % 排序并去重

    %    稀疏矩阵：对称存储
    iIdx = [edges(:,1); edges(:,2)];
    jIdx = [edges(:,2); edges(:,1)];
    A    = sparse(iIdx, jIdx, true, N, N);  % logical sparse

    % ── Step 2: 计算每个顶点的邻接个数（度数）
    valence = full(sum(A,2));               % N×1 向量
    maxVal  = max(valence);                 % 自动检测最大值

    % ── Step 3: 生成固定列数的 NaN 矩阵并填充
    C = NaN(N, maxVal);

    for k = 1:N
        nbrs = find(A(k,:));                % 行 k 的所有邻接顶点索引
        C(k,1:numel(nbrs)) = nbrs;          % 写入；剩余保持 NaN
    end
end