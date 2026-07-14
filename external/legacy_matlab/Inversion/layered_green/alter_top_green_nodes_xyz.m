%% Alter the top of the Green's function into triangulars, using the new method
function [G1,G2,G3,Gx_ss,Gx_ds,Gy_ss,Gy_ds,Gz_ss,Gz_ds]=alter_top_green_nodes_xyz(model_in,data_insar,green_name,F_strk,F_dip,varargin)
% 1. Load the altered geometry
% 2. For each node on the second layer, assign index for the target
% triangles to compute the homogeneous green's function


if ~isempty(varargin)
    for CC = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{CC*2-1})
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

tri_assignment=assignTrianglesToNodes(model_in); % each row corresponds to a second-layer nodal element
xe_insar=data_insar(:,1);
yn_insar=data_insar(:,2);
zinsar=zeros(length(xe_insar),1);
Nobs=length(xe_insar);
Npara=2*length(tri_assignment);

G1=zeros(Nobs,Npara); % rows: insar observations; columns: slip model (first half: ss/ second half: ds)
G2=zeros(Nobs,Npara); % rows: insar observations; columns: slip model (first half: ss/ second half: ds)
G3=zeros(Nobs,Npara); % rows: insar observations; columns: slip model (first half: ss/ second half: ds)
% Gx_ss=zeros(Nobs,1);
% Gy_ss=zeros(Nobs,1);
% Gz_ss=zeros(Nobs,1);
% Gx_ds=zeros(Nobs,1);
% Gy_ds=zeros(Nobs,1);
% Gz_ds=zeros(Nobs,1);
nu=0.25;
npt=length(tri_assignment);

for i=1:length(tri_assignment)
    element_indx=tri_assignment(i,1);%nodal element index
    v_node=[model_in(element_indx,6),model_in(element_indx,7),model_in(element_indx,8)];

    %% Part 1: combining the homogeneous solutions
    Gx_ss_homo=[];
    Gx_ds_homo=[];

    Gy_ss_homo=[];
    Gy_ds_homo=[];

    Gz_ss_homo=[];
    Gz_ds_homo=[];
    for k=2:size(tri_assignment,2)
        tri_indx=tri_assignment(i,k);
        if ~isnan(tri_indx)
            verts=model_in(tri_indx,19:27);
            P1=verts(:,1:3);
            P2=verts(:,4:6);
            P3=verts(:,7:9);
            [P1,P2,P3,~,~] = orientTriFromLeft(P1,P2,P3);
            % if strcmp(data_type,"azo")
            %     [ue1,un1,~]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,1,0,0,nu);%strike-slip
            %     [ue2,un2,~]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,0,1,0,nu);%dip-slip
            %     ulos1=un1*sinF+ue1*cosF;
            %     ulos2=un2*sinF+ue2*cosF;
            % else
            %     [ue1,un1,uz1]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,1,0,0,nu);%strike-slip
            %     [ue2,un2,uz2]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,0,1,0,nu);%dip-slip
            %     ulos1=ue1.*ve_insar+un1.*vn_insar+uz1.*vz_insar;
            %     ulos2=ue2.*ve_insar+un2.*vn_insar+uz2.*vz_insar;
            % end
            % 
            % if isempty(G_ss_homo)
            %     G_ss_homo=zeros(size(ulos1));
            %     G_ds_homo=zeros(size(ulos2));
            % end
            % 
            % G_ss_homo=G_ss_homo+ulos1;
            % G_ds_homo=G_ds_homo+ulos2;
            [ue_ss,un_ss,uz_ss]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,1,0,0,nu);%strike-slip
            [ue_ds,un_ds,uz_ds]=TDdispHS(xe_insar,yn_insar,zinsar,P1,P2,P3,0,1,0,nu);%strike-slip

            if isempty(Gx_ss_homo)
                Gx_ss_homo=zeros(size(ue_ss));
                Gx_ds_homo=zeros(size(ue_ss));

                Gy_ss_homo=zeros(size(ue_ss));
                Gy_ds_homo=zeros(size(ue_ss));

                Gz_ss_homo=zeros(size(ue_ss));
                Gz_ds_homo=zeros(size(ue_ss));
            end
            
            Gx_ss_homo=Gx_ss_homo+ue_ss;
            Gx_ds_homo=Gx_ds_homo+ue_ds;

            Gy_ss_homo=Gy_ss_homo+un_ss;
            Gy_ds_homo=Gy_ds_homo+un_ds;

            Gz_ss_homo=Gz_ss_homo+uz_ss;
            Gz_ds_homo=Gz_ds_homo+uz_ds;
           
        end
    end



    %% Part 2: combining the layered solutions
    vert_indx=getNodeTriangleVertices(model_in,element_indx);
    TotalModel=[];
    v_node=v_node/1e3;
    for j=1:size(vert_indx,1)
        v2=[model_in(vert_indx(j,1),6),model_in(vert_indx(j,1),7),model_in(vert_indx(j,1),8)];
        v3=[model_in(vert_indx(j,2),6),model_in(vert_indx(j,2),7),model_in(vert_indx(j,2),8)];


        v2=v2/1e3;
        v3=v3/1e3;


        [tmpX, tmpY, tmpZ] = interpolate_triangle2(v_node, v2, v3, 15);%default:110
        area=triangle_area(v_node,v2,v3);
        pmoment=area/length(tmpX)*1e6;
        w=taperTriangle(v_node,v2,v3,[tmpX,tmpY,tmpZ]);



        % interpolate the strike and dip angle
        tmpStrk = F_strk(tmpX, tmpY, tmpZ);
        tmpDip = F_dip(tmpX, tmpY, tmpZ);

        % setup the point source model
        % the length and width of point sources are set 0
        npatch = length(tmpX);
        tmpModel = zeros(npatch, 11);
        tmpModel(:,4) = tmpX .* 1e3;
        tmpModel(:,5) = tmpY .* 1e3;
        tmpModel(:,6) = tmpZ .* 1e3;
        tmpModel(:,7) = 0;
        tmpModel(:,8) = 0;
        tmpModel(:,9) = tmpStrk;
        tmpModel(:,10) = tmpDip;
        tmpModel(:,11) = pmoment * w .* ones(length(tmpX),1);
        TotalModel=[TotalModel;tmpModel];
    end

    % if strcmp(data_type,'insar')
    %     Gtmp = calc_green_gps_3d_edcmp(TotalModel, data_insar, green_name);
    % else
    %     Gtmp = calc_green_azo_edcmp(TotalModel, data_insar, green_name,cosF,sinF);
    % end
    % 
    % npatch=size(Gtmp,2)/2;
    % G_ss_layered = sum(Gtmp(:, 1:npatch), 2);  % strike-slip component
    % G_ds_layered = sum(Gtmp(:, npatch+1:end), 2); % dip-slip component
    % 
    % 
    % G_ss=G_ss_layered+G_ss_homo; % it should be a one-column array
    % G_ds=G_ds_layered+G_ds_homo; % it should be a one-column array
    % 
    % G_raw(:,i) = G_ss;  % strike-slip component
    % G_raw(:,i+npt) = G_ds; % dip-slip component
    [Gx,Gy,Gz] = calc_green_gps_3d_edcmp_xyz(TotalModel, data_insar, green_name);

    npatch=size(Gx,2)/2;

    Gx_ss_layered = sum(Gx(:, 1:npatch), 2);  % strike-slip component
    Gy_ss_layered = sum(Gy(:, 1:npatch), 2);  % strike-slip component
    Gz_ss_layered = sum(Gz(:, 1:npatch), 2);  % strike-slip component


    Gx_ds_layered = sum(Gx(:, npatch+1:end), 2);  % strike-slip component
    Gy_ds_layered = sum(Gy(:, npatch+1:end), 2);  % strike-slip component
    Gz_ds_layered = sum(Gz(:, npatch+1:end), 2);  % strike-slip component

    Gx_ss=Gx_ss_layered+Gx_ss_homo;
    Gx_ds=Gx_ds_layered+Gx_ds_homo;

    Gy_ss=Gy_ss_layered+Gy_ss_homo;
    Gy_ds=Gy_ds_layered+Gy_ds_homo;

    Gz_ss=Gz_ss_layered+Gz_ss_homo;
    Gz_ds=Gz_ds_layered+Gz_ds_homo;

    G1(:,i)=Gx_ss;
    G1(:,i+npt)=Gx_ds;

    G2(:,i)=Gy_ss;
    G2(:,i+npt)=Gy_ds;

    G3(:,i)=Gz_ss;
    G3(:,i+npt)=Gz_ds;


    % Gx_ss=Gx_ss+ue_ss;
    % Gy_ss=Gy_ss+un_ss;
    % Gz_ss=Gz_ss+uz_ss;
    % 
    % Gx_ds=Gx_ds+ue_ds;
    % Gy_ds=Gy_ds+un_ds;
    % Gz_ds=Gz_ds+uz_ds;

end


end