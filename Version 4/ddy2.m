function [Uddy2]=ddy2(U,dy)

% computes the second partial derivative of the variable U wrt y
% y is the second index of the array
% dy is constant, uniform grid
% for 2D
% using higher order central difference formula
% need to improve on the y boundary

Uddy2=U; %to copy the array size

Uddy2(:,3:end-2)=(-U(:,5:end)+16*U(:,4:end-1)-30*U(:,3:end-2)...
    +16*U(:,2:end-3)-U(:,1:end-4))/(12*dy^2);

Uddy2(:,end-1:end)=(-U(:,end-4:end-3)+4*U(:,end-3:end-2)...
    -5*U(:,end-2:end-1)+2*U(:,end-1:end))/(dy^2);
Uddy2(:,1:2)=(2*U(:,1:2)-5*U(:,2:3)+4*U(:,3:4)-U(:,4:5))/(dy^2);