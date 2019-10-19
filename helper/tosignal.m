function [ xc ] = tosignal( prob,x,y )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
[xSize,ySize] = size(prob.signal);
xc = zeros(size(prob.signal));
for k = 1:prob.os
    if prob.ndim == 1
        xc = xc + sqrt(xSize*ySize)*conj(prob.mask(:,k)).*ifft(y(:,k));
    else
        xc = xc + sqrt(xSize*ySize)*conj(prob.mask(:,:,k)).*ifft2(y(:,:,k));
    end
end
xc = (x+xc)/(1+prob.os);
end

