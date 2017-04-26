function [Uddx2]=ddx2(U,dx)

% computes the second partial derivative of the variable U wrt x
% x is the first index of the array
% dx is constant, uniform grid
% for 2D
% using higher order central difference formula
% need to improve on the x boundary

Uddx2=U; %to copy the array size

Uddx2(3:end-2,:)=(-U(5:end,:)+16*U(4:end-1,:)-30*U(3:end-2,:)...
    +16*U(2:end-3,:)-U(1:end-4,:))/(12*dx^2);
Uddx2(end-1:end,:)=(-U(end-4:end-3,:)+4*U(end-3:end-2,:)...
    -5*U(end-2:end-1,:)+2*U(end-1:end,:))/(dx^2);
Uddx2(1:2,:)=(2*U(1:2,:)-5*U(2:3,:)+4*U(3:4,:)-U(4:5,:))/(dx^2);
