% solves the problem of a 2x3x2 element beam that is held fixed at y=0 and 
% applied nodal forces at x=0,y=3, using Finite Element Method based on
% Galerkin's method of weighted residual, with linear shape functions and
% one gauss lagrange integration point. nodal coordinates and element 
% connectivity matrix are taken from nodes.txt and elem.txt respectively.

clear
clc

nodes=load('nodes.txt');%file containing node coordinates
elem=load('elem.txt');  %file containing element connectivity information

coord=zeros(8,3); %initialize coord matrix for storing coordinate of nodes for a hexahedron
Km=zeros(length(nodes)*3,length(nodes)*3); % initialize global stiffness matrix
Mm=zeros(length(nodes)*3,length(nodes)*3);

%% Generating the Stiffness Matrix [Km] and Mass Matrix [Mm]
for elemnum=1%:length(elem);
    for nodenum=1:8; % nodenum is the local node number
        coord(nodenum,:)=nodes(elem(elemnum,nodenum),:); % adding the coordinate information
    end

    %{
    currently the code solves for 1 gauss point integration
    development for multiple gauss points will be done later

    ng=1; %number of gauss point for gauss lagandre quadrature
    xg=zeros(ng,ng);
    wg=zeros(ng,ng);

    xg(1,1)=0; % gauss points and weight
    wg(1,1)=2;

    for i=1:ng; % loop for every gauss point
        for j=1:ng;
            for k=1:ng;
                der=0.125*[...
                    -(1-xg(ng,j)(1-xg(ng,k),-(1-xg(ng,j)(1+xg(ng,k),...
                    (1-xg(ng,j)(1+xg(ng,k),(1-xg(ng,j)(1-xg(ng,k),...
                    -(1+xg(ng,j)(1-xg(ng,k),-(1+xg(ng,j)(1+xg(ng,k)...
                    (1+xg(ng,j)(1+xg(ng,k),(1+xg(ng,j)(1-xg(ng,k);...
    %}                
    % gauss point 0 and weight 2 is used throughout the code from here on
    wg=2;
    der=0.125*[-1,-1, 1, 1,-1,-1, 1, 1;... % a matrix of derivatives of the shape function wrt local coord using gauss point 0,0,0 as integrand
               -1,-1,-1,-1, 1, 1, 1, 1;...
               -1, 1, 1,-1,-1, 1, 1,-1];

    jac= der*coord;  % jacobian of the transformation matrix
    detjac=det(jac); % determinant of the jacobian
    invjac=inv(jac); % inverse of the jacobian
    deriv=invjac*der;% a matrix of derivatives of the shape function wrt to global coordinates at the gauss point 0,0,0

    %% [B] matrix
    % computing the [B] matrix by assembling the derivatives 
    B=zeros(6,24);
    for j=1:8;
        i=3*j-2;
        B(1,i)=deriv(1,j);
        B(4,i)=deriv(2,j);
        B(6,i)=deriv(3,j);
        B(2,i+1)=deriv(2,j);
        B(4,i+1)=deriv(1,j);
        B(5,i+1)=deriv(3,j);
        B(3,i+2)=deriv(3,j);
        B(5,i+2)=deriv(2,j);
        B(6,i+2)=deriv(1,j);
    end
    
    %% Material Matrix [D]
    % compute the stress-strain relationship matrix, otherwise known as the
    % material matrix
    E=7e7; % young's modulus
    v=0.3; % poisson's ratio
    rho=7800; %density

    const1=(E*(1-v))/((1+v)*(1-2*v));
    const2=v/(1-v);
    const3=(1-2*v)/(2*(1-v));
    
    D=const1.*[1,const2,const2,0,0,0;
        const2,1,const2,0,0,0;
        const2,const2,1,0,0,0;
        0,0,0,const3,0,0;
        0,0,0,0,const3,0;
        0,0,0,0,0,const3];


    BtDB=B'*D*B;

    km=detjac*wg^3*BtDB;

    %% assemble the global stiffness matrix [Km]
    for i=1:8
        localnum1=i*3;
        globalnum1=elem(elemnum,i)*3;
        for j=1:8
            localnum2=j*3;
            globalnum2=elem(elemnum,j)*3;
            Km(globalnum1-2:globalnum1,globalnum2-2:globalnum2)=...
                km(localnum1-2:localnum1,localnum2-2:localnum2);
        end
    end
    
    %% Generating the element Mass Matrix [M]
    
    N=0.125*[1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0;...
        0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,;...
        0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1];
    NtN=N'*N;
    mm=rho*detjac*wg^3*NtN;
    
    %% Assemble the global Mass Matrix [Mm]
    for i=1:8
        localnum1=i*3;
        globalnum1=elem(elemnum,i)*3;
        for j=1:8
            localnum2=j*3;
            globalnum2=elem(elemnum,j)*3;
            Mm(globalnum1-2:globalnum1,globalnum2-2:globalnum2)=...
                mm(localnum1-2:localnum1,localnum2-2:localnum2);
        end
    end    
        
end



%% Loads
% constructing the force matrix. 10 unit load is applied to the three nodes
% at x=0 and y=3
F=zeros(length(nodes)*3,1);
for i=[10,22,34];
    F(i*3-2)=1000000;
end

%% Boundary Conditions
% apply displacement constraint=0 for bottom nodes 1,2,3,13,14,15,25,26,27
% note that this can only be used with Matlab coder that supports variable matrix
for i=[27,26,25,15,14,13,3,2,1];
    Km(i*3-2:i*3,:)=[];
    Km(:,i*3-2:i*3)=[];
    F(i*3-2:i*3)=[];
end

%% solve Ax=b
% currently using full inverse matrix computation. other methods should be 
% considered if given enough time
u=Km\F;

%% sorting the computed displacements into a matrix of x y z coordinates
displacement=zeros(size(nodes));
count=0;
for i=1:36;
    if any(i==[1,2,3,13,14,15,25,26,27])
        count=count+1;
        continue
    end
    displacement(i,:)=u(3*(i-count)-2:3*(i-count))';
end

threed_to_vtk ( nodes, elem, displacement, 'results.vtk', 'displacement' )

