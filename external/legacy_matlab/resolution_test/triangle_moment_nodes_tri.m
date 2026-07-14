function [M0_total, Mw, mu] = triangle_moment_nodes_tri(slipMetrics, velProfile,varargin)
% calcMomentFromSlip  Compute seismic moment and moment magnitude from slip model
%
%   [M0_total, Mw] = calcMomentFromSlip(slipMetrics, velProfile)
%
% Inputs:
%   slipMetrics — M×5 array:
%                   [cx, cy, cz, slip_cm, area_m2], where
%                     cx,cy,cz = triangle centroid (m)
%                     slip_cm  = slip magnitude (cm)
%                     area_m2  = triangle area (m^2)
%   velProfile  — N×5 array:
%                   [depth_top_km, depth_bottom_km, Vp_km_s, Vs_km_s, rho_kg_m3]
%                     depth_top_km    = shallow bound of layer (km, ≤ 0)
%                     depth_bottom_km = deep bound of layer (km, ≤ 0)
%                     Vp_km_s, Vs_km_s = P- and S-wave speeds (km/s)
%                     rho_kg_m3       = density (kg/m^3)
%
% Outputs:
%   M0_total — total seismic moment (N·m)
%   Mw       — moment magnitude
%
% Procedure:
%  1. Convert slip from cm to m.
%  2. Convert centroid depth from m to km.
%  3. For each triangle, find its velocity layer (depth_top ≥ depth ≥ depth_bottom).
%  4. Compute shear modulus μ = ρ * (Vs*1000)^2  [Pa].
%  5. Compute triangle moment M0_i = μ * slip_m * area_m2  [N·m].
%  6. Sum M0_i → M0_total.
%  7. Mw = (2/3) * (log10(M0_total) – 9.1).
shear=0;
if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{CC*2-1})
                case 'uniform_shear'
                    shear=varargin{CC*2};
            end
        catch
            error('Unrecognized Keyword');
        end
    end
end

% 1) Slip in meters
slip_m = slipMetrics(:,4) * 1e-2;

% 2) Centroid depth in km
depth_km = slipMetrics(:,3) / 1e3;

% 3) Preallocate shear modulus array
nTri = size(slipMetrics,1);
mu   = zeros(nTri,1);

for i = 1:nTri
    d = depth_km(i);
    % find layer where depth_top >= d >= depth_bottom
    idx = find( velProfile(:,1) >= d & d >= velProfile(:,2), 1, 'first' );
    if isempty(idx)
        error('Depth %.3f km falls outside the velocity profile bounds.', d);
    end
    Vs_km = velProfile(idx,3);
    rho    = velProfile(idx,5)*1e3;
    % 4) shear modulus μ [Pa]
    mu(i) = rho * (Vs_km*1e3)^2;
    if shear~=0
        mu(i)=shear;
    end
end

% 5) triangle moments [N·m]
areas   = slipMetrics(:,5);

M0_tri  = mu .* slip_m .* areas;

% 6) total moment
M0_total = sum(M0_tri);

% 7) moment magnitude
Mw = (2/3) * (log10(M0_total) - 9.1);

end
