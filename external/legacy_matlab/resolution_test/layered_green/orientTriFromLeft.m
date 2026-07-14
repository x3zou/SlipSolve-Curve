function [P1o,P2o,P3o,flipped,sgn] = orientTriFromLeft(P1,P2,P3,varargin)
% 统一三角形在 y–z 投影中的顺序（从左往右看 = 视线 +x）
% 默认统一为 CW；可用 'desired','CCW' 改为逆时针
% 返回:
%   flipped: 是否交换了2↔3
%   sgn    : dot([1 0 0], cross(P2-P1,P3-P1)) 的符号（>0表示CCW）

ip = inputParser;
ip.addParameter('desired','CW',@(s)ischar(s)||isstring(s));
ip.addParameter('tol',1e-12,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
ip.parse(varargin{:});
want = upper(string(ip.Results.desired));
tol  = ip.Results.tol;

P1 = P1(:).'; P2 = P2(:).'; P3 = P3(:).';
n  = [1 0 0];                         % 从左看
sgnRaw = dot(n, cross(P2-P1, P3-P1)); % >0: CCW in yz when viewed from +x
if     sgnRaw >  tol, sgn =  1;
elseif sgnRaw < -tol, sgn = -1;
else,  sgn =  0;  % 共线/退化
end

flipped = false;
if want=="CW"  && sgn>0      % 现在是CCW，改成CW
    [P2,P3] = deal(P3,P2); flipped = true;
elseif want=="CCW" && sgn<0  % 现在是CW，改成CCW
    [P2,P3] = deal(P3,P2); flipped = true;
end

P1o=P1; P2o=P2; P3o=P3;
end
