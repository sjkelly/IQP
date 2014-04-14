%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% truss.m - October 15 2002 %
% author: Jean H. Prevost + David Luet %
% analyses of 1,2 and 3D elastic trusses %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; % removes all variables from the workspace.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
% Material %
%%%%%%%%%%%%

%E=200.; % Young's modulus

%%%%%%%%%%%%%
% Geometric %
%%%%%%%%%%%%%

%A=2; % Area

%%%%%%%%%%%%%%%%%%%%%%%
% Material Properties %
%%%%%%%%%%%%%%%%%%%%%%%
%mat(i,1) = E;  mat(i,1) = A; i = material property number 

mat(1,1) = 200; mat(2,1) = 40e3;
mat(1,2) = 200; mat(2,2) = 25e3;

%%%%%%%%
% Mesh %
%%%%%%%%

nsd=3; % number of space dimension
ndf=nsd; % number of freedom per node
nen=2; % number of element nodes
nel=33; % number of elements/trusses
nnp=18; % number of nodal points

%%%%%%%%%%%%%%%%%%%%%
% Nodal coordinates %
%%%%%%%%%%%%%%%%%%%%%
% xn(i,N):= coordinate i for node N
% N=1,...,nnp
% i=1,...,nsd

xn=zeros(nsd,nnp);

for x = 1:5
    xn(1,x) = (x-1)*20;
    xn(1,x+18) = (x-1)*20;
end

for x = 1:4
    xn(1,x+5) = 10 + (x-1)*20;
    xn(3,x+5) = 8;
end

for x = 1:9
    xn(1,x+9) = (x-1)*10;
    xn(3,x+9) = 16;
end

%%%%%%%%%%%%%%%%
% Connectivity %
%%%%%%%%%%%%%%%%
% ien(a,e)=N
% N: global node number - N=1,...,nnp
% e: element number - e=1,...,nel
% a: local node number - a=1,...,nen
% mp(e) = mat(:,i), i = material property number

ien=zeros(nen,nel);

%% Chord Elements
for x = 1:4
   ien(1,x) = x; ien(2,x) = x+1;  mp(:,x)= mat(:,1);
end

for x = 1:8
   ien(1,x+4) = x+9; ien(2,x+4) = x+10;  mp(:,x+4)= mat(:,1);
end

%%%%%

% Web Elements
for x = 1:5
   ien(1,x+12) = x; ien(2,x+12) = 2*x+8;  mp(:,x+12)= mat(:,2);
end

for x = 1:4
   ien(1,x+17) = x+5; ien(2,x+17) = 2*x+8;  mp(:,x+17)= mat(:,2);
end

for x = 1:4
   ien(1,x+21) = x+5; ien(2,x+21) = 2*x+9;  mp(:,x+21)= mat(:,2);
end

for x = 1:4
   ien(1,x+25) = x+5; ien(2,x+25) = 2*x+10;  mp(:,x+25)= mat(:,2);
end

for x = 1:2
   ien(1,x+29) = x+1; ien(2,x+29) = x+5;  mp(:,x+29)= mat(:,2);
end

for x = 1:2
   ien(1,x+31) = x+2; ien(2,x+31) = x+7;  mp(:,x+31)= mat(:,2);
end

%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%
% prescribed displacement (essential boundary condition)
%
% idb(i,N)=1 if the degree of freedom i of the node N is prescribed
%         =0 otherwise
%
% 1) initialize idb to 0

idb=zeros(ndf,nnp);

% 2) enter the flag for prescribed displacement boundary conditions

idb(1,1) = 1;
idb(2,1) = 1;
idb(2,5) = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prescribed nodal displacement boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g(i,N): prescribed displacement for the dof i of node N
% initialize g

g=zeros(ndf,nnp);

% enter the values

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prescribed nodal forces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(i,N): prescribed force for the dof i of node N
% initialize f

f = zeros(ndf,nnp);

% enter the values

f(2,11) = -100;
f(2,12) = -100;
f(2,13) = -100;

%-----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number the equations; build the id table %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[id,neq] = number_eq(idb,nnp,ndf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the elemental quantities in the elemental coordinate system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for e=1:nel
    E = mp(1,e); A = mp(2,e);
    [Ke(:,:,e),ke(:,:,e),Qe(:,:,e)] = Ke_truss(E,A,xn,ien(:,e),nen,ndf,nsd);
end;

% Contribution of the prescribed displacements to the elemental force vector
% fe=fe-Ke*Ue

fe=zeros(ndf*nen,nel); % fe may be non zero in general
Ue=zeros(ndf*nen,nel);
for e=1:nel
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e) = g(i,ien(n,e));
        end
    end
    fe(:,e)=fe(:,e)-Ke(:,:,e)*Ue(:,e);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly operation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------
% build K and F
%----------------

K=zeros(neq,neq);
F=zeros(neq,1);

% input the prescribed nodal forces in F

for N=1:nnp
    for i=1:ndf
        if (id(i,N) > 0)
            P=id(i,N);
            F(P)=f(i,N);
        end
    end
end

% compute global K and F

if (neq > 0)
    for e=1:nel
        K = addstiff(K,id,Ke(:,:,e),ien(:,e),nen,ndf);
        F = addforce(F,id,fe(:,e),ien(:,e),nen,ndf);
    end
end

% Solve the system

if (neq > 0)
    U=K\F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post processing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% complete U %
%%%%%%%%%%%%%%

Ucomp=zeros(ndf,nnp);
for N=1:nnp
    for i=1:ndf
        if (id(i,N) == 0)
            Ucomp(i,N)=g(i,N);
        else
            P=id(i,N);
            Ucomp(i,N)=U(P);
        end
    end
end

% print results

disp('Nodal Displacements:')
disp(' node d1 d2')
for N=1:nnp
    disp(sprintf('%5d %7g %7g',N,Ucomp(:,N)))
end
disp(' ')

%%%%%%%%%%%%%
% REACTIONS %
%%%%%%%%%%%%%
% build the idb table; overwrite original idb table
% idb(i,N): equation number associated with dof i of node N

ineq=0; % number of equations
for i=1:ndf
    for N=1:nnp
        if (idb(i,N) > 0) % assign an equation number to all prescribed nodes
            ineq=ineq+1;
            idb(i,N)=ineq;
        end;
    end;
end;

% Contribution of the displacement to the elemental force vector
% fe=Ke*Ue

for e=1:nel
    Ue(:,e)=zeros(ndf*nen,1);
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e)=Ucomp(i,ien(n,e));
        end
    end
    fe(:,e)=Ke(:,:,e)*Ue(:,e);
end

% compute reactions R %

R=zeros(ineq,1);
for e=1:nel
    R = addforce(R,idb,fe(:,e),ien(:,e),nen,ndf);
end

% Collect reactions

Rcomp=zeros(ndf,nnp);
for N=1:nnp
    for i=1:ndf
        if (idb(i,N) > 0)
            Rcomp(i,N)=R(idb(i,N));
        end
    end
end

% print results

disp('Nodal Reactions')
disp(' node R1 R2')
for N=1:nnp
    disp(sprintf('%5d %7g %7g',N,Rcomp(:,N)))
end
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%
% AXIAL FORCES/STRESSES %
%%%%%%%%%%%%%%%%%%%%%%%%%

for e=1:nel
    Ue(:,e)=zeros(ndf*nen,1);
    for n=1:nen
        for i=1:ndf
            Ue(i+(n-1)*ndf,e)=Ucomp(i,ien(n,e));
        end
    end
    if (nsd > 1)
        axial(:,e)=ke(:,:,e)*Qe(:,:,e)*Ue(:,e);
    else
        axial(:,e)=ke(:,:,e)*Qe(e)*Ue(:,e);
    end;
    stress(e)=axial(2,e)/A;
    strain(e)=stress(e)/E;
end;

% print results

disp('Element Axial force/stress/strain')
disp('  elem     force      stress      strain')
for e=1:nel
    disp(sprintf('%5d   %7g   %7g   %7g',e,axial(2,e),stress(e),strain(e)))
end
disp(' ')

%%%%%%%%%%%%%%%%%%%%
% plot the results %
%%%%%%%%%%%%%%%%%%%%
plot_results('truss',xn,f,idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp,axial);