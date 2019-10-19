function [ y ] = toimage( prob,x )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
y = zeros(size(prob.data));
[xSize,ySize] = size(prob.signal);
for k = 1:prob.os
    if prob.ndim == 1
        y(:,k) = fft(prob.mask(:,k).*x)/sqrt(xSize*ySize);
    else
    y(:,:,k) = fft2(prob.mask(:,:,k).*x)/sqrt(xSize*ySize);
    end
end

end

