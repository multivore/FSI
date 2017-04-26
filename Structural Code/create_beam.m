% a simple function to generate a grid of mesh for fea of a beam
% the resulting data is a nodal coordinate file nodes, an element
% connectivity matrix file elem, a node list for the various faces
% the parameters that make up the beam are lengths x y z and 
% number of elements nx ny nz
clear
clc
tic


disp('time')

x=10;
y=1;
z=1;
nx=50;
ny=5;
nz=5;
dx=x/nx;
dy=y/ny;
dz=z/nz;

fid=fopen('nodes.txt','w');
for k=0:nz;
    for j=0:ny;
        for i=0:nx;
            fprintf ( fid, '  %f  %f  %f\n', [i*dx,j*dy,k*dz]);
        end
    end
end
fclose(fid);
disp([num2str(toc,3),'   nodes.txt saved'])

fid=fopen('elem.txt','w');
nxy=(nx+1)*(ny+1);
for k=1:nz;
    for j=1:ny;
        for i=1:nx;
            num=i+(nx+1)*(j-1)+nxy*(k-1);
            fprintf ( fid, '  %d  %d  %d  %d  %d  %d  %d  %d\n', [num, num+1, num+1+(nx+1), num+(nx+1),...
                num+nxy,num+nxy+1,num+nxy+1+(nx+1),num+nxy+(nx+1)]);
        end
    end
end
fclose(fid);
disp([num2str(toc,3),'   elem.txt saved'])

fid=fopen('facexy.txt','w');
count=0;
for j=0:ny;
    for i=0:nx;
        count=count+1;
        fprintf ( fid, '  %d', count);
    end
end
fclose(fid);
disp([num2str(toc,3),'   facexy.txt saved'])

fid=fopen('faceyz.txt','w');
for k=0:nz;
    for j=0:ny;
        nodenum=1+j*(nx+1)+k*nxy;
        fprintf ( fid, '  %d', nodenum);
    end
end
fclose(fid);
disp([num2str(toc,3),'   faceyz.txt saved'])

fid=fopen('faceyz2.txt','w');
for k=0:nz;
    for j=0:ny;
        nodenum=(j+1)*(nx+1)+k*nxy;
        fprintf ( fid, '  %d', nodenum);
    end
end
fclose(fid);
disp([num2str(toc,3),'   faceyz2.txt saved'])

fid=fopen('facexz.txt','w');
for k=0:nz;
    for i=0:nx;
        nodenum=1+i+ny*(nx+1)+k*nxy;
        fprintf ( fid, '  %d', nodenum);
    end
end
fclose(fid);
disp([num2str(toc,3),'   facexz.txt saved'])
