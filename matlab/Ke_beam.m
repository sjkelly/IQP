%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ke,ke,Qe]=Ke_beam(E,A,I,xn,ien,nen,ndf,nsd)
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



% Rotation matrix to global coordinate system

if ( nsd == 1)        % 1D case
    if (xn(1,node1) < xn(1,node2))Qe=1;
    else Qe=-1;
    end;
else
    Qe = zeros(ndf*nen,ndf*nen);
    v = xn(:,node2)-xn(:,node1);
    v = v/norm(v);
    if (nsd == 2)          % 2D case
        ke=zeros(6,6);
        ke(:,:)=[(E*A/L),0          ,0         ,-E*A/L,0          ,0;
            0      ,12*E*I/L^3 ,6*E*I/L^2 ,0     ,-12*E*I/L^3,6*E*I/L^2;
            0      ,6*E*I/L^2  ,4*E*I/L   ,0     ,-6*E*I/L^2 ,2*E*I/L;
            -(E*A/L),0          ,0         , E*A/L,0          ,0;
            0      ,-12*E*I/L^3,-6*E*I/L^2,0     ,12*E*I/L^3 ,-6*E*I/L^2;
            0      ,6*E*I/L^2  ,2*E*I/L   ,0     ,-6*E*I/L^2 ,4*E*I/L];
        
        Qe(1,1) = v(1);
        Qe(1,2) = v(2);
        
        Qe(2,1) = -v(2);
        Qe(2,2) = v(1);
        
        Qe(3,3) = 1;
        
        Qe(4,4) = v(1);
        Qe(4,5) = v(2);
        
        Qe(5,4) = -v(2);
        Qe(5,5) = v(1);
        
        Qe(6,6) = 1;
        
    elseif (nsd ==3)      % 3D case
        ke=zeros(12,12);
        ke(:,:)=[(E*A/L),0          ,0          ,0    ,0          ,0          ,-E*A/L,0           ,0           ,0     ,0          ,0;       
                  0     ,12*E*Iz/L^3,0          ,0    ,0          ,6*E*Iz/L^2 ,0     ,-12*E*Iz/L^3,0           ,0     ,0          ,6*E*Iz/L^2;
                  0     ,0          ,12*E*Iy/L^3,0    ,-6*E*Iy/L^2,0          ,0     ,0           ,-12*E*Iz/L^3,0     ,-6*E*Iy/L^2,0;
                  0     ,0          ,0          ,G*J/L,0          ,0          ,0     ,0           ,0           ,-G*J/L,0          ,0; 
                  0     ,0          ,-6*E*Iy/L^2,0    ,4*E*Iy/L   ,0          ,0     ,0           ,6*E*Iy/L^2  ,0     ,2*E*Iy/L   ,0;
                  0     ,6*E*Iy/L^3 ,0          ,0    ,0          ,4*E*Iz/L   ,0     ,-6*E*Iz/L^2 ,0           ,0     ,0          ,2*E*Iz/L;
                -(E*A/L),0          ,0          ,0    ,0          ,0          ,E*A/L ,0           ,0           ,0     ,0          ,0;       
                  0     ,-12*E*Iz/L^3,0         ,0    ,0          ,-6*E*Iz/L^2,0     ,12*E*Iz/L^3 ,0           ,0     ,0          ,-6*E*Iz/L^2;
                  0     ,0          ,-12*E*Iy/L^3,0   ,6*E*Iy/L^2 ,0          ,0     ,0           ,12*E*Iz/L^3 ,0     ,6*E*Iy/L^2 ,0;
                  0     ,0          ,0          ,-G*J/L,0         ,0          ,0     ,0           ,0           ,G*J/L,0           ,0; 
                  0     ,0          ,-6*E*Iy/L^2,0    ,2*E*Iy/L   ,0          ,0     ,0           ,6*E*Iy/L^2  ,0     ,4*E*Iy/L   ,0;
                  0     ,6*E*Iz/L^3 ,0          ,0    ,0          ,2*E*Iz/L   ,0     ,-6*E*Iz/L^2 ,0           ,0     ,0          ,4*E*Iz/L]
        if v(1) ==0 && v(2) ==0
            if v(3)>0
                r = [0 0 1; 0 1 0; -1 0 0];
            else
                r = [0 0 -1; 0 1 0; 1 0 0];
            end
        else
            CXx = v(1);
            CYx = v(2);
            CZx = v(3);
            D = (CXx*CXx+CYx*CYx)^.5;
            CXy = -CYx/D;
            CYy = CXx/D;
            CZy = 0;
            CXz = -CXx*CZx/D;
            CYz = -CYx*CZx/D;
            CZz = D;
            r = [CXx CYx CZx; CXy CYy CZy; CXz CYz CZz];
        end
        Qe= [r zeros(3) zeros(3) zeros(3);
            zeros(3) r zeros(3) zeros(3);
            zeros(3) zeros(3) r zeros(3);
            zeros(3) zeros(3) zeros(3) r];
        
    end;
end;

% Global Stiffness Ke(i,j)

Ke=zeros(ndf*nen,ndf*nen);
if (nsd >1)
    Ke(:,:) = Qe(:,:)'*ke(:,:)*Qe(:,:);
else
    Ke(:,:) = Qe*ke(:,:);
end;