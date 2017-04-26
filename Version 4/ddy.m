function [Uddy]=ddy(U,dy)

% computes the partial derivative of the variable U wrt x
% x is the first index of the array
% dx is constant, uniform grid
% for 2D
% using higher order central difference formula
% need to improve on the x boundary

Uddy=U; %to copy the array size

Uddy(:,3:end-2)=(-U(:,5:end)+8*U(:,4:end-1)-8*U(:,2:end-3)+U(:,1:end-4))...
    /(12*dy);
Uddy(:,end-1:end)=(U(:,end-3:end-2)-4*U(:,end-2:end-1)+3*U(:,end-1:end))...
    /(2*dy);
Uddy(:,1:2)=(-3*U(:,1:2)+4*U(:,2:3)-U(:,3:4))/(2*dy);
