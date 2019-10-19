function [ x,y ] = graph_calfourier(c,d,prob)
%graph_cal graph projection with the graph Ax=y
%   Detailed explanation goes here
% w = L\(A*c+A*A'*d);
% y = transpose(L)\w;
% x = c+A'*(d-y);

[xSize,ySize] = size(prob.tsignal);
x = 0.5*sqrt(xSize*ySize)*ifft2(fft2(c)/sqrt(xSize*ySize)+d);
y = fft2(x)/sqrt(xSize*ySize);
end

