function slip_model = make_node_geometry(fault_id,PTs,VID)
%% Geneate the geometry of the slip model for nodes slip element
% 1. faultid 2. node id 3. layer_id, 4.strike-slip 5.dip-slip 6.x 7.y 8.z 
% 9-15. connectivity id of nodes (testing)
% Xiaoyu Zou, 04/24/2025
npt=length(PTs);
z=PTs(:,3);
[~,~,g]=unique(z,'sorted');
layers=max(g)-g+1;
[C,~]=NodeContact(VID,PTs);
slip_model=zeros(npt,15);
if size(C,2)<7 % if less than 7 contact, unify the format by filling with NaN
    diff=7-size(C,2);
    sup=nan(size(slip_model,1),diff);
    C=[C,sup];
end
slip_model(:,9:15)=C;
slip_model(:,1)=fault_id;
slip_model(:,2)=1:npt;
slip_model(:,3)=layers;
slip_model(:,6:8)=PTs(:,1:3);

end



