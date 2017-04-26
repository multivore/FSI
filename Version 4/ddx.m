function [Uddx]=ddx(U,dx)

% computes the partial derivative of the variable U wrt x
% x is the first index of the array
% dx is constant, uniform grid
% for 2D
% using higher order central difference formula
% need to improve on the x boundary

Uddx=U; %to copy the array size

Uddx(3:end-2,:)=(-U(5:end,:)+8*U(4:end-1,:)-8*U(2:end-3,:)+U(1:end-4,:))...
    /(12*dx);
Uddx(end-1:end,:)=(U(end-3:end-2,:)-4*U(end-2:end-1,:)+3*U(end-1:end,:))...
    /(2*dx);
Uddx(1:2,:)=(-3*U(1:2,:)+4*U(2:3,:)-U(3:4,:))/(2*dx);

