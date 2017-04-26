%--------------------------------------------------------------------------
%Program: Splitting2D
%Description: Fluid flow solver based on splitting method, in 2D.
%Version: 2.0
%Author: Ilyas Ramdzan
%Last Updated: 29/11/16
%Uses external modules for derivatives of functions
%--------------------------------------------------------------------------
%Version 2 Notes:
%Try to apply code with immersed boundary method before moving on with FSI

%% initialize program

% Clear workspace and command window
clear; clc; tic
disp([datestr(now),' START'])
disp(' Runtime   Simtime    P-iter    U-iter       res')
dispformat='%8.2f  %8.4f  %8d  %8d  %8.5f \n';
format short

% input variables
Re=100; %for 1 unit length / 1 unit speed
nu=1/Re;
dt=0.01;
sor_w=0.5;

% fluid domain
dx=0.0125;
dy=0.0125;
x=0:dx:1;
y=0:dy:1;
M=length(x);
N=length(y);
[x,y]=ndgrid(x,y);

% preallocating solution variables
u=zeros(M,N);
v=u;
P=u;
xN=u;
yN=u;

%% Start time iterations
time=0;
for itime=1:46
    %% Output contour plot while running
%     if(rem(itime,10)==0)
%         quiver(x,y,u,v,0.6,'k-')
%         axis equal
%         axis([0 1 0 1])
%         hold on 
%         pcolor(x,y,P);
%         colormap(jet)
%         colorbar
%         shading interp
%         axis equal
%         axis([0 1 0 1])
%         title({['2D Cavity Flow with Re = ',num2str(Re)];['time(\itt) = ',num2str(time)]})
%         xlabel('Spatial co-ordinate (x) \rightarrow')
%         ylabel('Spatial co-ordinate (y) \rightarrow')
%         drawnow;
%         hold off
%     end
    %% Fluid Solver code
    time=time+dt;
    uold=u;
    vold=v;
    
    % Velocity BC here
    u(:,end)=1;
    
    
    % Calculate Non Linear terms, xN, yN
    % Solve explicitly using second order Adam Bashforth formula

    xN1=xN; %store previous iteration values
    yN1=yN;

    uuddx=ddx(u.*u,dx);
    vuddy=ddy(v.*u,dy);
    xN=-0.5*(uuddx+vuddy);
    uddx=ddx(u,dx);
    uddy=ddy(u,dy);
    xN=xN-0.5*(u.*uddx+v.*uddy);

    uvddx=ddx(u.*v,dx);
    vvddy=ddy(v.*v,dy);
    yN=-0.5*(uvddx+vvddy);
    vddx=ddx(v,dx);
    vddy=ddy(v,dy);
    yN=yN-0.5*(u.*vddx+v.*vddy);

    % Calculate U hat
    uhat=(dt/2)*(3*xN-xN1)+u;
    vhat=(dt/2)*(3*yN-yN1)+v;

    % U hat BC


    % Pressure Poisson Eqn Gauss Seidel with SOR

    divUhat=(ddx(uhat,dx)+ddy(vhat,dy)); % divergence of U hat

    M=size(P,1);
    N=size(P,2);
    isorP=0;
    sor_res=1;
    while sor_res>1e-7
        isorP=isorP+1;
        Po=P; % store previous iteration values
        for i=2:M-1
            for j=2:N-1
                P(i,j)=(-1/4)*(dx*dy/dt*divUhat(i,j)...
                    -P(i,j+1)-P(i,j-1)-P(i+1,j)-P(i-1,j));% need to correct for dx=/=dy
                P(i,j)=Po(i,j)+sor_w*(P(i,j)-Po(i,j));
            end
        end
        % Pressure BC here
        sor_res=max(max(abs(Po-P)));
    end
    P(1,:)=P(2,:);
    P(end,:)=P(end-1,:);
    P(:,1)=P(:,2);
    P(:,end)=P(:,end-1);

    % Calculate U double hat
    uhat2=-ddx(P,dx)*dt+uhat;
    vhat2=-ddy(P,dy)*dt+vhat;

    % U double hat BC

    % Velocity marching
    % Now comes the linear terms in the N-S equations
    % Solve for Un+1 implicitly which results in the Helmholtz equation
    % Using Gauss Seidel with SOR
    % Differencing scheme used for time is the 2nd order Adam Moulton formula
    % Similar to Crank Nicholson

    % First, calculate the laplacian of U (del2U) for the diffusive term of 
    % the nth time index
    del2u=ddx2(u,dx)+ddy2(u,dy);
    del2v=ddx2(v,dx)+ddy2(v,dy);

    M=size(u,1);
    N=size(u,2);
    sor_res=1;
    isorU=0;
    while sor_res>1e-7
        isorU=isorU+1;
        uo=u; % store previous iteration values
        vo=v;
        for i=2:M-1
            for j=2:N-1
                u(i,j)=((((u(i-1,j)+u(i+1,j))*dy^2+(u(i,j-1)+u(i,j+1))*dx^2)...
                    /(dy^2*dx^2)+del2u(i,j))*0.5*nu*dt+uhat2(i,j))...
                    /(1+2*(dx^2+dy^2)/(dx^2*dy^2)*nu*dt*0.5);
                u(i,j)=uo(i,j)+sor_w*(u(i,j)-uo(i,j));

                v(i,j)=((((v(i-1,j)+v(i+1,j))*dy^2+(v(i,j-1)+v(i,j+1))*dx^2)...
                    /(dy^2*dx^2)+del2v(i,j))*0.5*nu*dt+vhat2(i,j))...
                    /(1+2*(dx^2+dy^2)/(dx^2*dy^2)*nu*dt*0.5);
                v(i,j)=vo(i,j)+sor_w*(v(i,j)-vo(i,j));
            end
        end
        sor_res=max(max([abs(uo-u) abs(vo-v)]));
    end

    % calculate residual from previous time step
    res=max(max([abs(uold-u) abs(vold-v)]))...
        /max(max([abs(u) abs(v)]));
    
    % Display iteration information
    fprintf(dispformat,toc,time,isorP,isorU,res)
    
end

disp([datestr(now),' FINISH ',num2str(toc,3),' seconds'])

%% Output final contour plot

figure
hold
contourf(x,y,P)
quiver(x,y,u,v,10)
xlim([0,1])
ylim([0,1])
% figure
% plot(u(ceil(end/2),:),y(ceil(end/2),:))
% figure
% plot(x(:,ceil(end/2)),v(:,ceil(end/2)))

%seedx=0.1:dx*10:0.9;
%seedy=0.1:dy*10:0.9;

%[seedy,seedx]=meshgrid(seedy,seedx);

%figure
%streamline(x,y,u,v,seedx,seedy)
%streamline(x,y,u,v,0.5,0.5)
