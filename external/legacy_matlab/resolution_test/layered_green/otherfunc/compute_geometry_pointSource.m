function [strike, dip, area,L1,L2,XV,YV,ZV] = compute_geometry_pointSource(xm, ym, zm)
    % [xm, ym, zm] are interpolated mesh coordinates
    % for uniform distributed points, zm are same for each layer
    % each point source is located at the center of grid data
    % that is, size(xm) = size(xp) + 1
    % area is used to plot the slip distribution afterwards
    % Testing: looking for L1 and L2, width W, and top edge vertices XV, YV, ZV
    
    % number of point source along strike and dip
    nL = size(xm, 2) - 1;
    nW = size(xm, 1) - 1;
    L1=[];
    L2=[];
    XV=[];
    YV=[];
    ZV=[];

    
    strike = zeros(nW, nL);
    dip = zeros(nW, nL);
    area = zeros(nW, nL);
    d2r = pi/180;
    
    for i = 1:nW
        % down edge of grid data
        x1 = xm(i,1:nL);
        y1 = ym(i,1:nL);
        
        x2 = xm(i,2:nL+1);
        y2 = ym(i,2:nL+1);
        
        ph1 = -atan2(y1-y2, x1-x2)/d2r + 90;
        
        % top edge of grid data
        x3 = xm(i+1,1:nL);
        y3 = ym(i+1,1:nL);
        z3 = zm(i+1,1:nL);
        
        x4 = xm(i+1,2:nL+1);
        y4 = ym(i+1,2:nL+1);
        z4 = zm(i+1,2:nL+1); 
        
        XV=[XV;x4];
        YV=[YV;y4];
        ZV=[ZV;z4];


        ph2 = -atan2(y3-y4, x3-x4)/d2r + 90;
        
        % using the average strike
        strike(i,:) = 0.5*(ph1 + ph2);       
                   
    end
    
    for i = 1:nW
        topX = movmean(xm(i+1,:),2);  % length nL vector
        topY = movmean(ym(i+1,:),2);
        topZ = movmean(zm(i+1,:),2);
        
        botX = movmean(xm(i,:),2);
        botY = movmean(ym(i,:),2);
        botZ = movmean(zm(i,:),2);
        
        dLayer = mean(topZ - botZ);
        
        for j = 1:nL
%             theta = (90 - strike(i,j)) * d2r;
            theta = -strike(i,j)*d2r;
            [xt, yt] = xy2XY(topX(j), topY(j), theta);
            [xb, yb] = xy2XY(botX(j), botY(j), theta);
            
            relx = xb - xt;
%             rely = yb - yt;  % should be zero
%             disp(['The relative x diff is: ', num2str(relx)]);
%             disp(['The relative y diff is: ', num2str(rely)]);
            dip(i,j) = atan2(dLayer, relx)/d2r;  
            
            botL = sqrt((xm(i,j) - xm(i,j+1))^2 + (ym(i,j) - ym(i,j+1))^2);
            L1=[L1;botL];
            topL = sqrt((xm(i+1,j) - xm(i+1,j+1))^2 + (ym(i+1,j) - ym(i+1,j+1))^2);
            L2=[L2,topL];
            area(i,j) = 0.5*(topL + botL) * dLayer / sind(dip(i,j));
        end
    end
end