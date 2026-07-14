function [UX,UY,UZ]=slip2xyz_tri(xpt,ypt,slip_model,PTs)
%function to calculate the surface displacements at given points using slip
%model with triangular solutions
% default input and output values is in cm
%****************INPUT**************
%xpt  ---- vector of  x-coordinates of observations
%ypt  ---- vector of y-coordinates of observations

% Last Updated by Xiaoyu Zou


Npatch=length(slip_model);
Npt=length(xpt);
UX=zeros(Npt,1);
UY=zeros(Npt,1);
UZ=zeros(Npt,1);
zpt=zeros(Npt,1);
nu=0.25;

for i=1:Npatch
    indx1=slip_model(i,4);
    indx2=slip_model(i,5);
    indx3=slip_model(i,6);
    
    vert1=[PTs(indx1,1),PTs(indx1,2),PTs(indx1,3)];
    vert2=[PTs(indx2,1),PTs(indx2,2),PTs(indx2,3)];
    vert3=[PTs(indx3,1),PTs(indx3,2),PTs(indx3,3)];

    [uE,uN,uZ]=TDdispHS(xpt,ypt,zpt,vert1,vert2,vert3,slip_model(i,2),slip_model(i,3),0,nu);%strike-slip
    UX=UX+uE;
    UY=UY+uN;
    UZ=UZ+uZ;
end