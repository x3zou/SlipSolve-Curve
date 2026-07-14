function G=alter_top_green_nodes(G_raw,slip_model,data_insar,varargin)
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
slip_model2=slip_model(find(slip_model(:,1)==2),:);



%un-altered number of patch
Npatch_unaltered=length(slip_model)+max(slip_model(:,1));

%altered number of patch
N_patch_altered=length(slip_model);


xe_insar=data_insar(:,1);
yn_insar=data_insar(:,2);
%zinsar=data_insar(:,3);
zinsar=zeros(length(xe_insar),1);
ve_insar=data_insar(:,4);
vn_insar=data_insar(:,5);
vz_insar=data_insar(:,6);
Nobs=length(xe_insar);

xp=slip_model(:,6);
yp=slip_model(:,7);
zp=slip_model(:,8);
lp=slip_model(:,18);
wp=slip_model(:,19);

strkp=slip_model(:,16);
dip0=slip_model(:,17);


%G=zeros(Nobs,Npara);
G=G_raw;
nu=0.25;

G(:,[1 1+Npatch_unaltered size(slip_model1,1)+1 size(slip_model1,1)+1+Npatch_unaltered])=[];

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
    if slip_model(k,8)==0
        dxf=lp(k)/2;
        dyf=0;
        dzf=0;
        strike=strkp(k)*d2r;
        theta=(90.0-strkp(k))*d2r;
        [dx,dy]=xy2XY(dxf,dyf,-theta);
        dz=dzf;

        xxo=xp(k)+dx;
        yyo=yp(k)+dy;
        zzo=zp(k)+dz;

        xpt=xe_insar-xxo;
        ypt=yn_insar-yyo;
        delta=dip0(k)*d2r;
        d=-zzo;
        len=lp(k);
        W=wp(k);
        fault_type1=1;
        fault_type2=2;

        U1=1;
        U2=1;
        tp=zeros(size(xe_insar));
        [ue1,un1,uz1]=calc_okada(HF,U1,xpt,ypt,nu,delta,d,len,W,fault_type1,strike,tp);
        [ue2,un2,uz2]=calc_okada(HF,U2,xpt,ypt,nu,delta,d,len,W,fault_type2,strike,tp);
        ulos1=ue1.*ve_insar+un1.*vn_insar+uz1.*vz_insar;
        ulos2=ue2.*ve_insar+un2.*vn_insar+uz2.*vz_insar;

        G(:,k)=ulos1;
        G(:,k+N_patch_altered)=ulos2;

    end

end
end