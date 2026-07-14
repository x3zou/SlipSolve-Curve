function G=calc_green_insar_curve_AZO(slip_model,PTs,data_insar,varargin)
% Calculate the Green's Function for InSAR AZO observation using Okada's
% solution, using triangular dislocations
%
% Usage: G=calc_green_insar_AZO(data_slip_model,data_insar);
%
%   data_insar=[xe_insar,yn_insar,los_insar,ve_insar,vn_insar,vz_insar];
%
% by Kang Wang on 07/27/2015
% Last Udpated by Kang Wang on 07/27/2015
% Modified by Xiaoyu on 01/09/2025

format long

% Verticies of the triangles
xtemp=PTs(:,1);
ytemp=PTs(:,2);
ztemp=PTs(:,3);
indices=slip_model(:,4:6);

Npatch=length(indices);
Npara=2*Npatch;

xe_insar=data_insar(:,1);
yn_insar=data_insar(:,2);
%zinsar=data_insar(:,3);
zinsar=zeros(length(xe_insar),1);
ve_insar=data_insar(:,4);
vn_insar=data_insar(:,5);
vz_insar=data_insar(:,6);
Nobs=length(xe_insar);

% convert the LOS angle to heading angle
theta_az = -atan2d(vn_insar,ve_insar) - 180;

G=zeros(Nobs,Npara);
nu=0.25;


if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{CC*2-1})
                case 'poisson'
                    nu = varargin{CC*2};
                case 'sinf'
                    cosF=varargin{CC*2};
                case 'cosf'
                    sinF=varargin{CC*2};
            end
        catch
            error('Unrecognized Keyword');
        end
    end
end

for k=1:Npatch
  ind=indices(k,:);
  % x_vert=[xtemp(ind(1)),xtemp(ind(2)),xtemp(ind(3))];
  % y_vert=[ytemp(ind(1)),ytemp(ind(2)),ytemp(ind(3))];
  % z_vert=[ztemp(ind(1)),ztemp(ind(2)),ztemp(ind(3))];

  P1=[xtemp(ind(1)),ytemp(ind(1)),ztemp(ind(1))];
  P2=[xtemp(ind(2)),ytemp(ind(2)),ztemp(ind(2))];
  P3=[xtemp(ind(3)),ytemp(ind(3)),ztemp(ind(3))];

  
  [ue1,un1,~]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,1,0,0,nu);%strike-slip
  [ue2,un2,~]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,0,1,0,nu);%dip-slip
  % ulos1=ue1.*sind(theta_az) + un1.*cosd(theta_az);
  % ulos2=ue2.*sind(theta_az) + un2.*cosd(theta_az);
  ulos1=un1*sinF+ue1*cosF;
  ulos2=un2*sinF+ue2*cosF;

  G(:,k)=ulos1;
  G(:,k+Npatch)=ulos2;

end

end