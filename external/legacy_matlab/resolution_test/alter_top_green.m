function G=alter_top_green(G_raw,slip_model,PTs,data_insar,varargin)
% Alter the top of layered green's function by homogenous green's function
%
% Usage: G=calc_green_insar(data_slip_model,data_insar);
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

%G=zeros(Nobs,Npara);
G=G_raw;
nu=0.25;


if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{CC*2-1})
                case 'poisson'
                    nu = varargin{CC*2};
            end
        catch
            error('Unrecognized Keyword');
        end
    end
end

for k=1:Npatch

  ind=indices(k,:);

  P1=[xtemp(ind(1)),ytemp(ind(1)),ztemp(ind(1))];
  P2=[xtemp(ind(2)),ytemp(ind(2)),ztemp(ind(2))];
  P3=[xtemp(ind(3)),ytemp(ind(3)),ztemp(ind(3))];

  if ztemp(ind(1))==0 || ztemp(ind(2))==0 || ztemp(ind(3))==0

    [ue1,un1,uz1]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,1,0,0,nu);%strike-slip
    [ue2,un2,uz2]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,0,1,0,nu);%dip-slip
    ulos1=ue1.*ve_insar+un1.*vn_insar+uz1.*vz_insar;
    ulos2=ue2.*ve_insar+un2.*vn_insar+uz2.*vz_insar;

    G(:,k)=-ulos1;
    G(:,k+Npatch)=-ulos2;
    
  end

end

end