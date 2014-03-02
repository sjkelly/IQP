%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_results.m - October,15 2002                         %
% last modified: November 18,2002                          %
% author: David Luet - luet@princeton.edu                  %
% plot results and b.c. for trusses,beams, heat conduction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_results(type,xn,f,Idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       TRUSS AND BEAM                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp(type, 'truss') | strcmp(type,'beam'))
    
    display_factor=0.1;
    eps=1e-4;
    
    % characteristic distances
    xmax=max(max(xn(1,:)));
    xmin=min(min(xn(1,:)));
    
    if (nsd >1)
        ymax=max(max(xn(2,:)));
        ymin=min(min(xn(2,:)));
        
        Lcar=max([xmax-xmin;ymax-ymin]);
    else
        Lcar=xmax-xmin;
    end;
    
    k=0;
    
    while (k~=4)
        k=menu('results','undeformed mesh and BCs', 'undeformed and deformed mesh','reactions','exit');
        
        switch k
            case 1
                figure;
                axis equal;
                title('Undeformed mesh and BCs');
                hold on;
                plot_mesh_underformed(nel,ien,xn,nnp,nsd);
                numbers(nel,ien,xn,nnp,nsd);
                plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xn,nsd);
                plot_bc_force(type,Lcar,f,display_factor,nnp,xn,nsd);
                hold off;
            case 2
                figure;
                axis equal;
                title('Undeformed and deformed mesh');
                hold on;
                plot_mesh_underformed(nel,ien,xn,nnp,nsd);
                numbers(nel,ien,xn,nnp,nsd);
                plot_mesh_deformed(type,xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen);
                hold off;
            case 3
                figure;
                axis equal;
                title('reactions');
                hold on;
                plot_mesh_underformed(nel,ien,xn,nnp,nsd);
                plot_reactions(type,Lcar,Rcomp,Idb,display_factor,xn,nnp,ndf,nsd)
                hold off;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         HEAT                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(type, 'heat'))
    k=0;
    
    while (k~=4)
        k=menu('results','mesh','Temperature with numbers','Temperature w\o numbers','Fluxes','exit');
        
        switch k
            case 1
                figure;
                axis equal;
                title('Mesh');
                hold on;
                if (nen == 3)
                    triplot(ien',xn(1,:),xn(2,:));
                    numbers_tria(nel,ien,xn,nnp,nsd);
                else
                    for e=1:nel,
                        for m=1:nen,
                            node(m)=ien(m,e);
                        end;
                        x=[xn(1,node(1)); xn(1,node(2)); xn(1,node(3)); xn(1,node(4)); xn(1,node(1))];
                        y=[xn(2,node(1)); xn(2,node(2)); xn(2,node(3)); xn(2,node(4)); xn(2,node(1))];
                        plot(x,y,'r-o');
                    end;
                    numbers_quad(nel,ien,xn,nnp,nsd);
                end;
                hold off;
            case 2
                figure;
                axis equal;
                %           title('Temperature Distribution');
                hold on;
                triplot(ien',xn(1,:),xn(2,:));
                plot_temp_tria(Ucomp,nel,ien,xn,nen,nnp,nsd,1);
                hold off;
            case 3
                figure;
                axis equal;
                title('Temperature Distribution');
                hold on;
                %            triplot(ien',xn(1,:),xn(2,:));
                plot_temp_tria(Ucomp,nel,ien,xn,nen,nnp,nsd,0);
                hold off;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       PLOT FUNCTIONS                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% undeformed configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_mesh_underformed(nel,ien,xn,nnp,nsd)
for e=1:nel
    node1=ien(1,e);
    node2=ien(2,e);
    if (nsd == 3)
        
    end
    x0=[xn(1,node1);xn(1,node2)];
    if (nsd > 1)
        y0=[xn(2,node1);xn(2,node2)];
        if nsd ==3
            z0=[xn(3,node1);xn(3,node2)];
        end
    else
        y0=[0;0];
    end;
    if nsd < 3
        plot(x0,y0,'b-o');
    else
        line(x0,y0,z0);
    end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbers the nodes and elements - truss and beam %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numbers(nel,ien,xn,nnp,nsd)
% number the nodes in the undeformed configuration
for n=1:nnp
    if (nsd > 1)
        text(xn(1,n),xn(2,n), num2str(n),'FontSize',14);
    else
        text(xn(1,n),0, num2str(n),'FontSize',14);
    end;
end
% number the elements in the undeformed configuration
for e=1:nel
    node1=ien(1,e);
    node2=ien(2,e);
    %compute the coordinates of the middle point
    xg=0.5*(xn(:,node1)+xn(:,node2));
    s=sprintf('(%d)', e);
    if (nsd > 1)
        text(xg(1),xg(2), s,'FontSize',14);
    else
        text(xg(1),0, s,'FontSize',14);
    end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbers the nodes and elements - triangular mesh %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numbers_tria(nel,ien,xn,nnp,nsd)
% number the nodes
for n=1:nnp
    if (nsd > 1)
        text(xn(1,n),xn(2,n), num2str(n),'FontSize',12);
    else
        text(xn(1,n),0, num2str(n),'FontSize',12);
    end;
end
% number the elements
for e=1:nel
    node1=ien(1,e);
    node2=ien(2,e);
    node3=ien(3,e);
    %compute the coordinates of the center of mass
    xg=(1/3)*(xn(:,node1)+xn(:,node2)+xn(:,node3));
    s=sprintf('(%d)', e);
    text(xg(1),xg(2), s,'FontSize',12);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbers the nodes and elements - quad mesh       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numbers_quad(nel,ien,xn,nnp,nsd)
% number the nodes
for n=1:nnp
    if (nsd > 1)
        text(xn(1,n),xn(2,n), num2str(n),'FontSize',12);
    else
        text(xn(1,n),0, num2str(n),'FontSize',12);
    end;
end
% number the elements
for e=1:nel
    node1=ien(1,e);
    node2=ien(2,e);
    node3=ien(3,e);
    node4=ien(4,e);
    %compute the coordinates of the center of mass
    xg=(1/4)*(xn(:,node1)+xn(:,node2)+xn(:,node3)+xn(:,node4));
    s=sprintf('(%d)', e);
    text(xg(1),xg(2), s,'FontSize',12);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% deformed configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_mesh_deformed(type,xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen);

switch type
    case 'truss'
        scale=display_factor*Lcar/max(max(abs(Ucomp))); % scale factor for the displacements
        for e=1:nel
            node1=ien(1,e);
            node2=ien(2,e);
            xt(:,1)=xn(:,node1);
            xt(:,2)=xn(:,node2);
            
            if (nsd == 1)
                xt(2,1)=0;
                xt(2,2)=0;
            end;
            
            for i=1:ndf
                if Idb(i,node1)~=0
                    xt(i,1)=xt(i,1)+scale*Ucomp(i,node1);
                else
                    xt(i,1)=xt(i,1)+scale*Ucomp(i,node1);
                end;
                if Idb(i,node2)~=0
                    xt(i,2)=xt(i,2)+scale*Ucomp(i,node2);
                else
                    xt(i,2)=xt(i,2)+scale*Ucomp(i,node2);
                end;
            end;
            plot(xt(1,:),xt(2,:)','r-o');
            if nsd == 3
                line(xt(1,:),xt(2,:),xt(3,:));
            end
        end
    case 'beam'
        xi=-1:0.1:1;
        nxi=size(xi,2);
        
        Umax=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                   SCALE FACTOR                                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e=1:nel
            node1=ien(1,e);
            node2=ien(2,e);
            if nsd == 2 
            % displacement + rotation
            Un(1,1)=Ucomp(1,node1);
            Un(1,2)=Ucomp(2,node1);
            Un(1,3)=Ucomp(3,node1);
            Un(2,1)=Ucomp(1,node2);
            Un(2,2)=Ucomp(2,node2);
            Un(2,3)=Ucomp(3,node2);
            elseif nsd ==3
            Un(1,1)=Ucomp(1,node1);
            Un(1,2)=Ucomp(2,node1);
            Un(1,3)=Ucomp(3,node1);
            Un(1,4)=Ucomp(1,node1);
            Un(1,5)=Ucomp(2,node1);
            Un(1,6)=Ucomp(3,node1);
            Un(2,1)=Ucomp(1,node2);
            Un(2,2)=Ucomp(2,node2);
            Un(2,3)=Ucomp(3,node2);
            Un(2,4)=Ucomp(1,node2);
            Un(2,5)=Ucomp(2,node2);
            Un(2,6)=Ucomp(3,node2);
            end
            
            % intial postion
            x0=zeros(nxi,nsd);
            for n=1:nxi
                x0(n,:)=((1-xi(n))/2)*xn(:,node1)'+((1+xi(n))/2)*xn(:,node2)';
            end
            
            Le=sqrt((xn(1,node1)-xn(1,node2))^2+(xn(2,node1)-xn(2,node2))^2);
            
            % g(i,j) - vector i coordinate j
            g=zeros(nsd,nsd);
            g(1,:)=(1/Le)*(xn(:,node2)-xn(:,node1))';
            
            g(2,1)=-g(1,2);
            g(2,2)=g(1,1);
            
            % to make g_3 = z - so that the rotation are the same in both systems of coordinates
            cross_prod=g(1,1)*g(2,2)-g(1,2)*g(2,1);
            if cross_prod < 0
                g(2,:)=-g(2,:);
            end;
            
            %transformation matrices
            beta=g;
            betap=inv(beta);
            
            % transform displacement coordinates in the global system
            un=zeros(nen,nsd);
            
            for n=1:nen
                for i=1:nsd
                    for j=1:nsd
                        un(n,i)=un(n,i)+Un(n,j)*betap(j,i);
                    end;
                end;
            end;
            
            for n=1:nen
                un(n,3)=Un(n,3);
            end;
            
            % transverse displacement in the local system of coordinates
            for i=1:nxi
                N1(i)=(1/4)*(1-xi(i))^2*(2+xi(i));
                N2(i)=(Le/8)*(1-xi(i)^2)*(1-xi(i));
                N3(i)=(1/4)*(1+xi(i))^2*(2-xi(i));
                N4(i)=(Le/8)*(-1+xi(i)^2)*(1+xi(i));
                u(i,2)=N1(i)*un(1,2)+N2(i)*un(1,3)+N3(i)*un(2,2)+N4(i)*un(2,3);
            end;
            
            % horizontal displacement in the local system of coordinates
            for i=1:nxi
                u(i,1)=(1/2)*((1-xi(i))*un(1,1)+(1+xi(i))*un(2,1));
            end;
            
            % transform from local to global system
            U=zeros(nxi,nsd);
            for i=1:nsd
                for n=1:nxi
                    for j=1:nsd
                        U(n,i)=U(n,i)+u(n,j)*beta(j,i);
                    end;
                end;
            end;
            
            % scale factor
            if (max(max(abs(U)))>Umax)
                Umax=max(max(abs(U)));
            end;
        end;
        
        scale=display_factor*Lcar/Umax;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                       PLOT                                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e=1:nel
            node1=ien(1,e);
            node2=ien(2,e);
            
            % displacement + rotation
            Un(1,1)=Ucomp(1,node1);
            Un(1,2)=Ucomp(2,node1);
            Un(1,3)=Ucomp(3,node1);
            Un(2,1)=Ucomp(1,node2);
            Un(2,2)=Ucomp(2,node2);
            Un(2,3)=Ucomp(3,node2);
            
            % intial postion
            x0=zeros(nxi,nsd);
            for n=1:nxi
                x0(n,:)=((1-xi(n))/2)*xn(:,node1)'+((1+xi(n))/2)*xn(:,node2)';
            end
            
            Le=sqrt((xn(1,node1)-xn(1,node2))^2+(xn(2,node1)-xn(2,node2))^2);
            
            % g(i,j) - vector i coordinate j
            g=zeros(nsd,nsd);
            g(1,:)=(1/Le)*(xn(:,node2)-xn(:,node1))';
            
            g(2,1)=-g(1,2);
            g(2,2)=g(1,1);
            
            % to make g_3 = z - so that the rotation are the same in both systems of coordinates
            cross_prod=g(1,1)*g(2,2)-g(1,2)*g(2,1);
            if cross_prod < 0
                g(2,:)=-g(2,:);
            end;
            
            %transformation matrices
            beta=g;
            betap=inv(beta);
            
            % transform displacement coordinates in the global system
            un=zeros(nen,nsd);
            
            for n=1:nen
                for i=1:nsd
                    for j=1:nsd
                        un(n,i)=un(n,i)+Un(n,j)*betap(j,i);
                    end;
                end;
            end;
            
            for n=1:nen
                un(n,3)=Un(n,3);
            end;
            
            % transverse displacement in the local system of coordinates
            for i=1:nxi
                N1(i)=(1/4)*(1-xi(i))^2*(2+xi(i));
                N2(i)=(Le/8)*(1-xi(i)^2)*(1-xi(i));
                N3(i)=(1/4)*(1+xi(i))^2*(2-xi(i));
                N4(i)=(Le/8)*(-1+xi(i)^2)*(1+xi(i));
                u(i,2)=N1(i)*un(1,2)+N2(i)*un(1,3)+N3(i)*un(2,2)+N4(i)*un(2,3);
            end;
            
            % horizontal displacement in the local system of coordinates
            for i=1:nxi
                u(i,1)=(1/2)*((1-xi(i))*un(1,1)+(1+xi(i))*un(2,1));
            end;
            
            % transform from local to global system
            U=zeros(nxi,nsd);
            for i=1:nsd
                for n=1:nxi
                    for j=1:nsd
                        U(n,i)=U(n,i)+u(n,j)*beta(j,i);
                    end;
                end;
            end;
            
            % deformed configuration in the global system
            xt=x0+scale*U;
            
            plot(xt(:,1),xt(:,2),'b-');
        end;
        
    otherwise
        disp('type not supported by plot_results');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary conditions on displacements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xn,nsd)
alpha=display_factor*Lcar; % scale factor  for bc symbols
switch type
    case 'truss'
        for P=1:nnp
            if (nsd >1)
                if ((Idb(1,P) ~= 0) & (Idb(2,P) ~= 0))
                    bc_symbols(xn(:,P),alpha,3);
                end;
                if ((Idb(1,P) ~= 0) & (Idb(2,P) == 0))
                    bc_symbols(xn(:,P),alpha,2);
                end;
                if ((Idb(1,P) == 0) & (Idb(2,P) ~= 0))
                    bc_symbols(xn(:,P),alpha,1);
                end;
            else
                if (Idb(1,P) ~= 0) bc_symbols([xn(1,P),0],alpha,2);
                end;
            end;
        end;
    case 'beam'
        for P=1:nnp
            if ((Idb(1,P) ~= 0) & (Idb(2,P) ~= 0) & (Idb(3,P) == 0))
                bc_symbols(xn(:,P),alpha,3);
            end;
            if ((Idb(1,P) ~= 0) & (Idb(2,P) == 0) & (Idb(3,P) == 0))
                bc_symbols(xn(:,P),alpha,2);
            end;
            if ((Idb(1,P) == 0) & (Idb(2,P) ~= 0) & (Idb(3,P) == 0))
                bc_symbols(xn(:,P),alpha,1);
            end;
            if ((Idb(1,P) == 0) & (Idb(2,P) ~= 0) & (Idb(3,P) ~= 0))
                bc_symbols(xn(:,P),alpha,4);
            end;
            if ((Idb(1,P) ~= 0) & (Idb(2,P) == 0) & (Idb(3,P) ~= 0))
                bc_symbols(xn(:,P),alpha,5);
            end;
            if ((Idb(1,P) ~= 0) & (Idb(2,P) ~= 0) & (Idb(3,P) ~= 0))
                bc_symbols(xn(:,P),alpha,6);
            end;
        end;
    otherwise
        disp('type not supported by plot_results');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary conditions on force %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_bc_force(type,Lcar,f,display_factor,nnp,xn,nsd)
delta=display_factor*Lcar/max(max(abs(f))); % scale factor for force b.c.
alpha=2*display_factor*Lcar; % scale factor for moment
for N=1:nnp
    if ( nsd >1)
        if ( (f(1,N) ~=0 ) | (f(2,N) ~= 0))
            quiver(xn(1,N),xn(2,N),f(1,N),f(2,N),delta,'k');
        end;
    else
        if (f(1,N) ~=0) quiver(xn(1,N),0,f(1,N),0,delta,'k');
        end;
    end;
end;

if (strcmp(type, 'beam'))
    for N=1:nnp
        if (f(3,N) > eps)
            plot_moments([xn(1,N);xn(2,N)],alpha,0);
        end;
        if (f(3,N) < -eps)
            plot_moments([xn(1,N);xn(2,N)],alpha,1);
        end;
    end;
end;


%%%%%%%%%%%%%
% reactions %
%%%%%%%%%%%%%
function plot_reactions(type,Lcar,Rcomp,Idb,display_factor,xn,nnp,ndf,nsd)
beta=display_factor*Lcar/max(max(abs(Rcomp))); % scale factor for reactions
alpha=2*display_factor*Lcar; % scale factor for moment
for N=1:nnp
    RN=zeros(ndf);
    if (nsd >1)
        if ((Idb(1,N) ~= 0) | (Idb(2,N) ~= 0))
            quiver(xn(1,N),xn(2,N),Rcomp(1,N),Rcomp(2,N),beta,'k');
        end;
        if (strcmp(type, 'beam') & (Idb(3,N) ~= 0))
            if (Rcomp(3,N) > eps)
                plot_moments([xn(1,N);xn(2,N)],alpha,0);
            end;
            if (Rcomp(3,N) < -eps)
                plot_moments([xn(1,N);xn(2,N)],alpha,1);
            end;
        end;
    else
        if (Idb(1,N) ~= 0)
            quiver(xn(1,N),0,Rcomp(1,N),0,beta,'k');
        end;
    end;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              BOUNDARY CONDITIONS SYMBOLS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bc_symbols(xp,alpha,symbol);

switch symbol
    case 1
        % v fixed
        x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
        y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
        
        line(x,y,'Color','k','LineWidth',1.2);
        
        for i=0:3,
            circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(7/8)*alpha],alpha/8);
        end;
        
        
    case 2
        % u fixed
        x=[xp(1);xp(1)-(3/4)*alpha;xp(1)-(3/4)*alpha;xp(1)];
        y=[xp(2);xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha;xp(2)];
        
        line(x,y,'Color','k','LineWidth',1.2);
        
        for i=0:3,
            circle([xp(1)-(7/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
        end;
        
    case 3
        % u and v fixed
        x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
        y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];
        
        line(x,y,'Color','k','LineWidth',1.2);
        
        for i=0:3,
            line([xp(1)-(alpha/4)+i*alpha/4;xp(1)-(alpha/2)+i*(alpha/4)], ...
                [xp(2)-(3/4)*alpha;xp(2)-alpha],'Color','k','LineWidth',1.2);
        end;
        
    case 4
        % v and theta fixed
        x=[xp(1)-alpha/2;xp(1)+alpha/2];
        y=[xp(2);xp(2)];
        
        line(x,y,'Color','k','LineWidth',1.2);
        
        for i=0:3,
            circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(1/8)*alpha],alpha/8);
        end;
        
    case 5
        % u and theta fixed
        x=[xp(1);xp(1)];
        y=[xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha];
        
        line(x,y,'Color','k','LineWidth',1.2);
        
        for i=0:3,
            circle([xp(1)-(1/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
        end;
        
    case 6
        % u, v and theta fixed
        line([xp(1)-alpha/2;xp(1)+alpha/2],[xp(2),xp(2)],'Color','k','LineWidth',1.2);
        for i=0:3,
            line([xp(1)-alpha/2+(i+1)*alpha/4, xp(1)-alpha/2+i*alpha/4],[xp(2),xp(2)-alpha/4]...
                ,'Color','k','LineWidth',1.2);
        end;
        
end;



function circle(x0,r);
theta=0:0.1:2*pi;
x=r*cos(theta)+x0(1);
y=r*sin(theta)+x0(2);

plot(x,y,'k','LineWidth',1.2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       MOMENTS                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_moments(xp,alpha,sign);
switch sign
    case 0  %positive moment
        r=alpha/4;
        theta=0:0.1:3*pi/2;
        x=r*cos(theta)+xp(1);
        y=r*sin(theta)+xp(2);
        
        plot(x,y,'k','LineWidth',1.2);
        line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],...
            'Color','k','LineWidth',1.2);
        line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],...
            'Color','k','LineWidth',1.2);
        
    case 1  % negative moment
        r=alpha/4;
        theta=pi:-0.1:-pi/2;
        x=r*cos(theta)+xp(1);
        y=r*sin(theta)+xp(2);
        
        plot(x,y,'k','LineWidth',1.2);
        line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],...
            'Color','k','LineWidth',1.2);
        line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],...
            'Color','k','LineWidth',1.2);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             HEAT TEMPERATURE DISTRIBUTION                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_temp_tria(Ucomp,nel,ien,xn,nen,nnp,nsd,flag_number)
%%%%%%%%%%%%%
% tria mesh %
%%%%%%%%%%%%%
if (nen == 3)
    for e=1:nel,
        for k=1:nen,
            node(k)=ien(k,e);
        end
        X(e,:,:)=[[xn(1,node(1))   xn(1,node(3))]; [xn(1,node(2))  xn(1,node(3))]];
        Y(e,:,:)=[[xn(2,node(1))   xn(2,node(3))]; [xn(2,node(2))  xn(2,node(3))]];
        Z(e,:,:)=[[Ucomp(node(1))  Ucomp(node(3))];[Ucomp(node(2)) Ucomp(node(3))]];
        
    end
    for e=1:nel,
        for i=1:2,
            for j=1:2,
                XX(i,j)=X(e,i,j);
                YY(i,j)=Y(e,i,j);
                ZZ(i,j)=Z(e,i,j);
            end
        end
        contourf(XX,YY,ZZ,5);
        %         if (flag_number == 1)
        %            [C,h]=contourf(XX,YY,ZZ,5);
        %             clabel(C,h);
        %         end
    end
    colorbar;
    
else
    %%%%%%%%%%%%%
    % quad mesh %
    %%%%%%%%%%%%%
    for e=1:nel,
        for k=1:nen,
            node(k)=ien(k,e);
        end
        X(e,:,:)=[[xn(1,node(1))  xn(1,node(4))]; [xn(1,node(2))  xn(1,node(3))]];
        Y(e,:,:)=[[xn(2,node(1))  xn(2,node(4))]; [xn(2,node(2))  xn(2,node(3))]];
        Z(e,:,:)=[[Ucomp(node(1)) Ucomp(node(4))];[Ucomp(node(2)) Ucomp(node(3))]];
        
    end
    for e=1:nel,
        for i=1:2,
            for j=1:2,
                XX(i,j)=X(e,i,j);
                YY(i,j)=Y(e,i,j);
                ZZ(i,j)=Z(e,i,j);
            end
        end
        contourf(XX,YY,ZZ,5);
        %         if (flag_number == 1)
        %             [C,h]=contour(XX,YY,ZZ,5);
        %             clabel(C,h);
        %        end
    end
    colorbar;
    
end

if (flag_number == 1)
    for n=1:nnp
        if (nsd > 1)
            text(xn(1,n),xn(2,n), num2str(Ucomp(n)),'FontSize',12);
        else
            text(xn(1,n),0, num2str(Ucomp(n)),'FontSize',12);
        end
    end
end

function N = shape_tria(n,xi,eta)
switch n
    case 1
        N=xi;
    case 2
        N=eta;
    case 3
        N=1-xi-eta;
end




%     n=2; % number of points on the parent element
%     h=1/(n-1);
%
%     % points on the parent element
%     for i=1:n,
%         xi(1,i)=1-(h*(i-1));
%     end;
%
%     for j=1:n
%         eta(2,j)=h*(j-1);
%     end;
%
%     for k=1:(2*n-1),
%         xi(3,k)=1-(h/2)*(k-1);
%         eta(3,k)=1-xi(3,k);
%     end;
%
%     m=0;
%     alpha=((2*n-1)-rem(2*n-1,2))/2
%     for p=1:alpha,
%         for q=1:n,
%             t=(q-1)*h;
%             m=m+1;
%             xi(4,m)=(1-t)*xi(1,p)+t*xi(3,p);
%             eta(4,m)=t*eta(3,p);
%         end;
%     end;
%
%     for p=alpha+1:(2*n-1),
%         for q=1:n,
%             m=m+1;
%             t=(q-1)*h;
%             xi(4,m)=t*xi(3,p);
%             eta(4,m)=(1-t)*eta(2,p-alpha)+t*eta(3,p);
%         end;
%     end;
%
%     m=0;
%     for i=1:2*n-1,
%         for j=1:n,
%             m=m+1;
%             XI(i,j)=xi(4,m);
%             ETA(i,j)=eta(4,m);
%         end;
%     end;
%
%     % actual plot
%
%     for e=1:nel,
%         X=zeros(2*n-1,n);
%         Y=zeros(2*n-1,n);
%         Z=zeros(2*n-1,n);
%
%         for k=1:nen,
%             node(k)=ien(k,e);
%         end
%
%         for i=1:2*n-1,
%             for j=1:n,
%                 for k=1:nen,
%                     X(i,j)=X(i,j)+xn(1,node(k))*shape_tria(k,XI(i,j),ETA(i,j));
%                     Y(i,j)=Y(i,j)+xn(2,node(k))*shape_tria(k,XI(i,j),ETA(i,j));
%                     Z(i,j)=Z(i,j)+Ucomp(node(k))*shape_tria(k,XI(i,j),ETA(i,j));
%                 end;
%             end;
%         end;
%
%
%         contourf(X,Y,Z,5);
%         if (flag_number == 1)
%             [C,h]=contour(X,Y,Z,5);
%             clabel(C,h);
%         end
%         colorbar;
%     end 
