%--------------------------------------------------------------------------
%Program: CreateMesh
%Description: A function to generate a simple grid of mesh for 
%             Structure2D solver
%Version: 2.0
%Author: Ilyas Ramdzan
%Last Updated: 30/11/16
%Notes:
% the resulting data is a nodal coordinate file nodes.txt, an element
% connectivity matrix file elem.txt, a node list for the loads 
% faceloaded.txt, and a node list for the fixed constraint.
% the parameters that make up the mesh are lengths x y z and 
% number of elements nx ny nz
%--------------------------------------------------------------------------
 
clear
clc
tic

disp('time')

x=1;
y=10;
nx=5;
ny=50;
dx=x/nx;
dy=y/ny;

fid=fopen('nodes.txt','w');
for j=0:ny;
    for i=0:nx;
        fprintf ( fid, '  %f  %f \n', [i*dx,j*dy]);
    end
end
fclose(fid);

disp([num2str(toc,3),'   nodes.txt saved'])
fid=fopen('elem.txt','w');
nxy=(nx+1)*(ny+1);
for j=1:ny;
    for i=1:nx;
        num=i+(nx+1)*(j-1);
        fprintf ( fid, '  %d  %d  %d  %d\n', ...
            [num, num+(nx+1), num+1+(nx+1), num+1]);
    end
end
fclose(fid);
disp([num2str(toc,3),'   elem.txt saved'])

fid=fopen('faceloaded.txt','w');
num=(nx+1)*(ny)+1;
for i=0:nx;
    fprintf ( fid, '  %d', num+i);
end
fclose(fid);
disp([num2str(toc,3),'   faceloaded.txt saved'])

fid=fopen('facefixed.txt','w');
for i=1:nx+1;
    fprintf ( fid, '  %d', i);
end
fclose(fid);
disp([num2str(toc,3),'   facefixed.txt saved'])
