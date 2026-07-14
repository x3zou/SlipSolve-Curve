function G=calc_green_gps_2d_okada_curve(triangular_model,data_gps,varargin)
%
% Usage:
%     G=calc_green_gps_2d_okada(data_slip_model,data_gps);
%
%     data_gps=[xe_gps,yn_gps,ue,un];

format long

% Verticies of the triangles
load(triangular_model);
points=DT1.Points;
xtemp=points(:,1);
ytemp=Vy1;
ztemp=points(:,2);
indices=DT1.ConnectivityList;

Npatch=length(DT1.ConnectivityList);
Npara=2*Npatch;


xe_gps=data_gps(:,1);
yn_gps=data_gps(:,2);

Nstn=length(xe_gps);
Nobs=2*Nstn;  %only use the horizontal components
G=zeros(Nobs,Npara);

HF=1;
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
    x_vert=[xtemp(ind(1)),xtemp(ind(2)),xtemp(ind(3))];
    y_vert=[ytemp(ind(1)),ytemp(ind(2)),ytemp(ind(3))];
    z_vert=[ztemp(ind(1)),ztemp(ind(2)),ztemp(ind(3))];


    [ue1,un1,~]=TDdispHS(xe_insar,yn_insar,0,x_vert,y_vert,z_vert,1,0,0,nu);%strike-slip
    [ue2,un2,~]=TDdispHS(xe_insar,yn_insar,0,x_vert,y_vert,z_vert,0,1,0,nu);%dip-slip

   U1_green_gps=[ue1;un1];
   U2_green_gps=[ue2;un2];
   G(:,k)=U1_green_gps;
   G(:,k+Npatch)=U2_green_gps;
end