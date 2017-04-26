% nodal coordinates and element connectivity matrix are taken from 
% nodes.txt and elem.txt 
% solve 8 node hexahedral elements only
% stress and strain calculation not included yet

clear
clc
tic
disp('time')

solve_for='static'; %solver flag, either 'static' or 'propagation'

nodes=load('nodes.txt');%file containing node coordinates
elem=load('elem.txt');  %file containing element connectivity information
[num_nodes,dimension]=size(nodes);
[num_elem,order_elem]=size(elem);

disp([num2str(toc,3),'  files loaded, number of nodes ',...
    num2str(num_nodes),'  number of elements ',num2str(num_elem)])



coord=zeros(8,3); %initialize local coord matrix for storing coordinate of nodes for a hexahedron
Km=zeros(num_nodes*3,num_nodes*3); % initialize global stiffness matrix
Mm=zeros(num_nodes*3,num_nodes*3); % initialize global mass matrix

%% Generating the Stiffness Matrix [Km] and Mass Matrix [Mm]
disp([num2str(toc,3),'  compute stiffness matrix and mass matrix'])
for elemnum=1:num_elem;
    for nodenum=1:8; % nodenum is the local node number
        coord(nodenum,:)=nodes(elem(elemnum,nodenum),:); % adding the coordinate information
    end
    
    ng=2; %number of gauss point for integration through the legendre quadrature
    xg=zeros(3,3);
    wg=zeros(3,3);

    % gauss points and weights, the first index denotes the current gauss 
    % point while the second index is the total number of integration 
    % points
    xg(1,1)=0; 
    wg(1,1)=2;
    
    xg(1,2)=-0.5773502691896;
    xg(2,2)=+0.5773502691896;
    wg(1,2)=1;
    wg(2,2)=1;
    
    xg(1,3)=-0.7745966692415;
    xg(2,3)=0;
    xg(3,3)=+0.7745966692415;
    wg(1,3)=0.5555555555555;
    wg(2,3)=0.8888888888888;
    wg(3,3)=0.5555555555555;
    
    % loop for every gauss point i,j,k represent the indices for the gauss
    % points
    km=zeros(24,24);
    for i=1:ng; 
        for j=1:ng;
            for k=1:ng;
                % a matrix of derivatives of shape function w.r.t local axis
                % taken at the gauss points 
                der=0.125*[...
                    
              -(1-xg(j,ng))*(1-xg(k,ng)), (1-xg(j,ng))*(1-xg(k,ng)),...
               (1+xg(j,ng))*(1-xg(k,ng)),-(1+xg(j,ng))*(1-xg(k,ng)),...
              -(1-xg(j,ng))*(1+xg(k,ng)), (1-xg(j,ng))*(1+xg(k,ng)),...
               (1+xg(j,ng))*(1+xg(k,ng)),-(1+xg(j,ng))*(1+xg(k,ng));...

              -(1-xg(i,ng))*(1-xg(k,ng)),-(1+xg(i,ng))*(1-xg(k,ng)),...
               (1+xg(i,ng))*(1-xg(k,ng)), (1-xg(i,ng))*(1-xg(k,ng)),...
              -(1-xg(i,ng))*(1+xg(k,ng)),-(1+xg(i,ng))*(1+xg(k,ng)),...
               (1+xg(i,ng))*(1+xg(k,ng)), (1-xg(i,ng))*(1+xg(k,ng));...

              -(1-xg(i,ng))*(1-xg(j,ng)),-(1+xg(i,ng))*(1-xg(j,ng)),...
              -(1+xg(i,ng))*(1+xg(j,ng)),-(1-xg(i,ng))*(1+xg(j,ng)),...
               (1-xg(i,ng))*(1-xg(j,ng)), (1+xg(i,ng))*(1-xg(j,ng)),...
               (1+xg(i,ng))*(1+xg(j,ng)), (1-xg(i,ng))*(1+xg(j,ng))];
  
                jac= der*coord;  % jacobian of the transformation matrix
                detjac=det(jac); % determinant of the jacobian
                deriv=jac\der;% a matrix of derivatives of the shape function w.r.t to global coordinates

                %% [B] matrix
                % computing the [B] matrix by assembling the derivatives 
                B=zeros(6,24);
                for m=1:8;
                    n=3*m-2;
                    B(1,n)=deriv(1,m);
                    B(4,n)=deriv(2,m);
                    B(6,n)=deriv(3,m);
                    B(2,n+1)=deriv(2,m);
                    B(4,n+1)=deriv(1,m);
                    B(5,n+1)=deriv(3,m);
                    B(3,n+2)=deriv(3,m);
                    B(5,n+2)=deriv(2,m);
                    B(6,n+2)=deriv(1,m);
                end
    
                %% Material Matrix [D]
                % compute the stress-strain relationship matrix, otherwise known as the
                % material matrix. 
                E=7e7; % young's modulus
                nu=0.3; % poisson's ratio
                rho=7800; %density

                c1=(E*(1-nu))/((1+nu)*(1-2*nu));
                c2=nu/(1-nu);
                c3=(1-2*nu)/(2*(1-nu));

                D=c1.*[...
                    1,  c2, c2, 0,  0,  0;...
                    c2, 1,  c2, 0,  0,  0;...
                    c2, c2, 1,  0,  0,  0;...
                    0,  0,  0,  c3, 0,  0;...
                    0,  0,  0,  0,  c3, 0;...
                    0,  0,  0,  0,  0,  c3];

                %% [B]'[D][B]dxdydz 
                % compute the sum of the integrands at all the gauss points
                BtDB=B'*D*B;
                km=km+detjac*wg(i,ng)*wg(j,ng)*wg(k,ng)*BtDB;
            end
        end
    end
    %km(:,:,elemnum)=km; %not sure why this line of code is here
    %% assemble the global stiffness matrix [Km]
    for i=1:8
        localnum1=i*3;
        globalnum1=elem(elemnum,i)*3;
        for j=1:8
            localnum2=j*3;
            globalnum2=elem(elemnum,j)*3;
            Km(globalnum1-2:globalnum1,globalnum2-2:globalnum2)=...
                Km(globalnum1-2:globalnum1,globalnum2-2:globalnum2)+...
                km(localnum1-2:localnum1,localnum2-2:localnum2);
        end
    end
    
    %% Generating the element Mass Matrix [M]
    mm=zeros(24,24);
    for i=1:ng;
        for j=1:ng;
            for k=1:ng;
                % N is the shape function in the form f a matrix for the
                % three variables
                N=0.125*[...
                (1-xg(i,ng))*(1-xg(j,ng))*(1-xg(k,ng)),0,0,...
                (1+xg(i,ng))*(1-xg(j,ng))*(1-xg(k,ng)),0,0,...
                (1+xg(i,ng))*(1+xg(j,ng))*(1-xg(k,ng)),0,0,...
                (1-xg(i,ng))*(1+xg(j,ng))*(1-xg(k,ng)),0,0,...
                (1-xg(i,ng))*(1-xg(j,ng))*(1+xg(k,ng)),0,0,...
                (1+xg(i,ng))*(1-xg(j,ng))*(1+xg(k,ng)),0,0,...
                (1+xg(i,ng))*(1+xg(j,ng))*(1+xg(k,ng)),0,0,...
                (1-xg(i,ng))*(1+xg(j,ng))*(1+xg(k,ng)),0,0;...

                0,(1-xg(i,ng))*(1-xg(j,ng))*(1-xg(k,ng)),0,...
                0,(1+xg(i,ng))*(1-xg(j,ng))*(1-xg(k,ng)),0,...
                0,(1+xg(i,ng))*(1+xg(j,ng))*(1-xg(k,ng)),0,...
                0,(1-xg(i,ng))*(1+xg(j,ng))*(1-xg(k,ng)),0,...
                0,(1-xg(i,ng))*(1-xg(j,ng))*(1+xg(k,ng)),0,...
                0,(1+xg(i,ng))*(1-xg(j,ng))*(1+xg(k,ng)),0,...
                0,(1+xg(i,ng))*(1+xg(j,ng))*(1+xg(k,ng)),0,...
                0,(1-xg(i,ng))*(1+xg(j,ng))*(1+xg(k,ng)),0;...

                0,0,(1-xg(i,ng))*(1-xg(j,ng))*(1-xg(k,ng)),...
                0,0,(1+xg(i,ng))*(1-xg(j,ng))*(1-xg(k,ng)),...
                0,0,(1+xg(i,ng))*(1+xg(j,ng))*(1-xg(k,ng)),...
                0,0,(1-xg(i,ng))*(1+xg(j,ng))*(1-xg(k,ng)),...
                0,0,(1-xg(i,ng))*(1-xg(j,ng))*(1+xg(k,ng)),...
                0,0,(1+xg(i,ng))*(1-xg(j,ng))*(1+xg(k,ng)),...
                0,0,(1+xg(i,ng))*(1+xg(j,ng))*(1+xg(k,ng)),...
                0,0,(1-xg(i,ng))*(1+xg(j,ng))*(1+xg(k,ng))];

                NtN=N'*N;
                mm=mm+rho*detjac*wg(i,ng)*wg(j,ng)*wg(k,ng)*NtN;
            end
        end
    end
        
    %% Assemble the global Mass Matrix [Mm]
    for i=1:8
        localnum1=i*3;
        globalnum1=elem(elemnum,i)*3;
        for j=1:8
            localnum2=j*3;
            globalnum2=elem(elemnum,j)*3;
            Mm(globalnum1-2:globalnum1,globalnum2-2:globalnum2)=...
                Mm(globalnum1-2:globalnum1,globalnum2-2:globalnum2)+...
                mm(localnum1-2:localnum1,localnum2-2:localnum2);
        end
    end    
        
end
%gKm=Km;


%% Loads
% constructing the force matrix
disp([num2str(toc,3),'  applying loads'])

loaded=load('faceyz2.txt');
F=zeros(num_nodes*3,1);
Fx=zeros(num_nodes,1);
Fy=zeros(num_nodes,1);
Fz=zeros(num_nodes,1);

Fy(loaded)=-100/numel(loaded);

for i=loaded;
    F(i*3-2)=Fx(i);  
    F(i*3-1)=Fy(i);
    F(i*3)  =Fz(i);
end

%% Boundary Conditions
% apply displacement constraint=0 for bottom nodes
disp([num2str(toc,3),'  applying constraints'])

fixed=load('faceyz.txt'); % add the node number of fixed nodes here

trimvector=zeros((num_nodes-numel(fixed))*3,1);
trimcounter=0;
for i=1:num_nodes
    if any(i==fixed)
        continue
    end
    trimcounter=trimcounter+1;
    trimvector(trimcounter*3-2:trimcounter*3)=(i*3-2):i*3;
end
Km=Km(trimvector,trimvector);
F=F(trimvector);
Mm=Mm(trimvector,trimvector);

%{
% uses matlab variable matrix coder. not efficient.
sort(fixed);
for i=fixed(end:-1:1);
    Km(i*3-2:i*3,:)=[];
    Km(:,i*3-2:i*3)=[];
    F(i*3-2:i*3)=[];
    Mm(i*3-2:i*3,:)=[];
    Mm(:,i*3-2:i*3)=[];
end
%}
disp([num2str(toc,3),'  Number of equations = ',num2str(length(trimvector))])

%% Initial Conditions
% for dis and velocity
disp([num2str(toc,3),'  applying initial conditions'])

Uo=zeros(length(F),1);
dUdto=zeros(length(F),1);



%% solve [K]{U}={F} for static case
% currently using full inverse matrix computation. other methods should be 
% considered if given enough time
if strcmp(solve_for,'static');
    disp([num2str(toc,3),'  Solving [K]{u}={F}'])
    U=Km\F;
end

%% Time Marching using Newmark's B=1/4 (theta=1/2)
if strcmp(solve_for,'propagation');
    disp([num2str(toc,3),'  start time marching'])

    dt=0.05;
    endtime=1;
    time=0;

    MplusK=Mm/(0.5*dt)+Km*0.5*dt;
    Frames(endtime/dt) = struct('cdata',[],'colormap',[]);
    step=0;
    while time<=endtime;
        step=step+1;
        time=time+dt;
        disp([num2str(toc,3),'  time is currently ',num2str(time)])
        RHS=dt*F+(2/dt)*Mm*Uo+2*Mm*dUdto-(0.5*dt)*Km*Uo;
        U=MplusK\RHS;
        dUdt=(2/dt)*(U-Uo)-dUdto;
        % d2Udt2 is not computed as it is not used in the iterations
        Uo=U;
        dUdto=dUdt;
        
        %% Plot Animation
        % the code down here is a copy of the plot script further down with
        % added getframe to create an animation
        dis=zeros(num_nodes,dimension);
        count=0;
        for i=1:num_nodes;
            if any(i==fixed) % skip the boundary condition nodes
                count=count+1;
                continue
            end
            dis(i,:)=U(3*(i-count)-2:3*(i-count))';
        end
        dis_mag=sqrt(dis(:,1).^2+dis(:,2).^2+dis(:,3).^2);
        vis_scale=100; % visualization scale
        dis_scaled=dis*vis_scale;
        displaced=dis_scaled+nodes;
        
        face(1,:)=[1,2,3,4];
        face(2,:)=[1,2,6,5];
        face(3,:)=[1,4,8,5];
        face(4,:)=[2,6,7,3];
        face(5,:)=[3,4,8,7];
        face(6,:)=[5,6,7,8];

        figure
        hold on
        for elemnum=1:num_elem;
            for facenum=1:6;
                x=displaced(elem(elemnum,face(facenum,:)),1);
                y=displaced(elem(elemnum,face(facenum,:)),2);
                z=displaced(elem(elemnum,face(facenum,:)),3);
                c=dis_mag(elem(elemnum,face(facenum,:)),1);
                patch(x,y,z,c)
            end    
        end
        daspect([1,1,1])
        view(3)
        xlim([min(nodes(:,1))-0.3 max(nodes(:,1))+0.3])
        ylim([min(nodes(:,2))-0.3 max(nodes(:,2))+0.3])
        zlim([min(nodes(:,3)) max(nodes(:,3))])
        quiver3(displaced(loaded,1),displaced(loaded,2),displaced(loaded,3),...
            Fx(loaded),Fy(loaded),Fz(loaded))
        scatter3(displaced(fixed,1),displaced(fixed,2),displaced(fixed,3),'filled')

        drawnow
        Frames(step)=getframe;
        close all
    end
end
    
%% sorting the computed displacements into a matrix of x y z coordinates
dis=zeros(num_nodes,dimension);
count=0;
for i=1:num_nodes;
    if any(i==fixed) % skip the boundary condition nodes
        count=count+1;
        continue
    end
    dis(i,:)=U(3*(i-count)-2:3*(i-count))';
end
dis_mag=sqrt(dis(:,1).^2+dis(:,2).^2+dis(:,3).^2);
vis_scale=0.3*max(max(nodes)); % visualization scale
if max(max(abs(dis)))==0;
    dis_scaled=dis;
else
    dis_scaled=dis*vis_scale/max(max(abs(dis)));
end
displaced=dis_scaled+nodes;


%% visualization
disp([num2str(toc,3),'  processing plot'])

face(1,:)=[1,2,3,4];
face(2,:)=[1,2,6,5];
face(3,:)=[1,4,8,5];
face(4,:)=[2,6,7,3];
face(5,:)=[3,4,8,7];
face(6,:)=[5,6,7,8];

figure
hold on
for elemnum=1:num_elem;
    for facenum=1:6;
        x=displaced(elem(elemnum,face(facenum,:)),1);
        y=displaced(elem(elemnum,face(facenum,:)),2);
        z=displaced(elem(elemnum,face(facenum,:)),3);
        c=dis_mag(elem(elemnum,face(facenum,:)),1);
        patch(x,y,z,c)
    end    
end
daspect([1,1,1])
% view(3)
xlim([min(displaced(:,1)) max(displaced(:,1))])
ylim([min(displaced(:,2)) max(displaced(:,2))])
zlim([min(displaced(:,3)) max(displaced(:,3))])
quiver3(displaced(loaded,1),displaced(loaded,2),displaced(loaded,3),...
    Fx(loaded),Fy(loaded),Fz(loaded))
scatter3(displaced(fixed,1),displaced(fixed,2),displaced(fixed,3),'filled')

if strcmp(solve_for,'propagation');
    movie(Frames,2)
end

disp([num2str(toc,3),'  end of run'])

I=max(nodes(:,2))^3*max(nodes(:,3))/12;
L=max(nodes(:,1));
w=numel(loaded)*Fy((loaded(1)))/L;

format long
calc=abs(100*L^3/(3*E*I));
code=max(abs(dis(:,2)));
pcterr=(calc-code)*100/calc;
fprintf('handcalc       code   pcterror \n')
fprintf('%8.6f   %8.6f   %8.6f \n',calc,code,pcterr)



