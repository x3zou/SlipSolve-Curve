function G=alter_top_green_nodes_tri_los(G_raw,slip_model,slip_model_unaltered,data_insar,varargin)
% Alter the top of nodes' green's function by Okada's green's function
% input is the altered slip model
%
%
%
%   data_insar=[xe_insar,yn_insar,los_insar,ve_insar,vn_insar,vz_insar];
%
% by Kang Wang on 07/27/2015
% Last Udpated by Kang Wang on 07/27/2015
% Modified by Xiaoyu on 04/27/2025
format long
d2r=pi/180;
HF=1;
% Find the important index of the altered/unaltered slip model
slip_model1=slip_model(find(slip_model(:,1)==1),:);
% slip_model2=slip_model(find(slip_model(:,1)==2),:);



%un-altered number of patch
Npatch_unaltered=length(slip_model_unaltered);

%altered number of patch
N_patch_altered=length(slip_model);


xe_insar=data_insar(:,1);
yn_insar=data_insar(:,2);
%zinsar=data_insar(:,3);
zinsar=zeros(length(xe_insar),1);
ve_insar=data_insar(:,4);
vn_insar=data_insar(:,5);
vz_insar=data_insar(:,6);





%G=zeros(Nobs,Npara);
G=G_raw;
nu=0.25;

% G(:,[1 1+Npatch_unaltered size(slip_model1,1)+1 size(slip_model1,1)+1+Npatch_unaltered])=[];
G=resizeG(G,N_patch_altered);

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

for k=1:N_patch_altered
    if slip_model(k,3)==1
        
        vert1=[slip_model(k,19),slip_model(k,20),slip_model(k,21)];
        vert2=[slip_model(k,22),slip_model(k,23),slip_model(k,24)];
        vert3=[slip_model(k,25),slip_model(k,26),slip_model(k,27)];

        [ue1,un1,uz1]=TDdispHS(xe_insar,yn_insar,zinsar,vert1,vert2,vert3,1,0,0,nu);%strike-slip
        [ue2,un2,uz2]=TDdispHS(xe_insar,yn_insar,zinsar,vert1,vert2,vert3,0,1,0,nu);%dip-slip

        ulos1=ue1.*ve_insar+un1.*vn_insar+uz1.*vz_insar;
        ulos2=ue2.*ve_insar+un2.*vn_insar+uz2.*vz_insar;

        G(:,k)=ulos1;
        G(:,k+N_patch_altered)=ulos2;
    end

end
end