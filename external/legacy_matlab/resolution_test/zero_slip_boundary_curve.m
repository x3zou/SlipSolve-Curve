function [W,d] = zero_slip_boundary_curve(slip_model,PTs,layer_option,ratio)
% add zero-slip constraint at the triangular mesh
% Supplements: add zero-slip constraint at the bottom surface
%           also add zero-slip slip boundary slip at one side of the fault
% options: 'layer_no': add zeros slip for a certain layer (usually top and bottom layer)
%          'left': left side of the fault plane (for patches of the whole layer)
%          'right': right side of the fault plane

    Np = size(slip_model,1);

    %% Flag the left, right and bottom layer for each patch (bottom:-1;left:-2;right:-3)
    points=PTs;
    nflt=max(slip_model(:,1));
    bottom_indx=[];
    inc_indx=0;%adjusting the index by the number of points
    for k=1:nflt
        points_tmp=points(find(points(:,4)==k),:);
        bottom_indx_tmp=find(points_tmp(:,3)==min(points_tmp(:,3)));%Index for bottom points
        if k>1
            inc_indx=length(points(find(points(:,4)<k)));
            bottom_indx_tmp=bottom_indx_tmp+inc_indx;
        end
        bottom_indx=[bottom_indx;bottom_indx_tmp];
    end


    layer_depth=unique(points(:,3));
    layer_depth=flip(layer_depth,1);
    left_indx=[];
    right_indx=[];
    % locate the index of boundary points
    for i = 1:length(layer_depth)
        layerx=points(points(:,3)==layer_depth(i),2);
        temp_left_index=find(points(:,2)==min(layerx) & points(:,3)==layer_depth(i));
        temp_right_index=find(points(:,2)==max(layerx) & points(:,3) == layer_depth(i));
        left_indx=[left_indx;temp_left_index];
        right_indx=[right_indx;temp_right_index];
    end


    conn_list=slip_model(:,4:6);
    flags = zeros(size(conn_list,1),1);

    for i = 1:size(conn_list,1)
        if any(ismember(conn_list(i,:),bottom_indx))
            flags(i)= -1;
        elseif any(ismember(conn_list(i,:),left_indx))
            flags(i)= -2;
        elseif any(ismember(conn_list(i,:),right_indx))
            flags(i)= -3;
        end
    end


    

    % all_layer_id = slip_model(:,3);
    % nL = compute_patch_each_layer(slip_model);
    
    d = zeros(2*Np,1);
    
    % select this fault segment
    % NF = length(segment_ID);
    V = zeros(2*Np,1);
    % for ii = 1:NF
    
        if strcmp(layer_option,'bottom')
            patch_top_this_segment = find(flags==-1);
        elseif strcmp(layer_option,'left')      % start from the first patch of each layer
            patch_top_this_segment = find(flags==-2);
        elseif strcmp(layer_option,'right')
            patch_top_this_segment = find(flags==-3);
        else
            error('There is something wrong with the patch options!');
        end
    
        strike_indx = patch_top_this_segment';
    %     disp(strike_indx);
        dip_indx = strike_indx + Np;
    %     disp(dip_indx);
        zero_slip_indx = [strike_indx,dip_indx];        
        V(zero_slip_indx) = ratio;
    % end
    W = diag(V);      % Use diagonal function to speed it up
    
end