%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ke_truss.m - November, 21 2002 %
% author: Jean H. Prevost %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ke,ke,Qe]=Ke_truss(E,A,xn,ien,nen,ndf,nsd)
%------------------------------------------------------------------------
% Purpose:
% compute truss element stiffness
%
% Synopsis:
% [Ke,ke,Qe]=Ke_truss(E,A,xn,ien,nen,ndf,nsd)
%
% Variable Description:
% nen - number of nodes per element
% ndf - number of equations per node
% nsd - number of spatial dimensions
%------------------------------------------------------------------------

node1=ien(1);
node2=ien(2);

% compute the length for each element

v=xn(:,node2)-xn(:,node1);
L=norm(v);

% local stiffness

ke=zeros(2,2);
ke(:,:)=(E*A/L)*[[1 -1]; [-1 1]];

% Rotation matrix to global coordinate system

if ( nsd == 1)        % 1D case
    if (xn(1,node1) < xn(1,node2))Qe=1;
    else Qe=-1;
    end;
else
    Qe=zeros(nen,ndf*nen);
    v=xn(:,node2)-xn(:,node1);
    v=v/norm(v);
    if (nsd ==2)          % 2D case
        Qe(1,1)=v(1);
        Qe(1,2)=v(2);
        Qe(2,3)=v(1);
        Qe(2,4)=v(2);
    elseif (nsd ==3)      % 3D case
        Qe(1,1)=v(1);
        Qe(1,2)=v(2);
        Qe(1,3)=v(3);
        Qe(2,4)=v(1);
        Qe(2,5)=v(2);
        Qe(2,6)=v(3);
    end;
end;

% Global Stiffness Ke(i,j)

Ke=zeros(ndf*nen,ndf*nen);
if (nsd >1)
    Ke(:,:)=Qe(:,:)'*ke(:,:)*Qe(:,:);
else
    Ke(:,:)=Qe*ke(:,:);
end;
