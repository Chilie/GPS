function [ X ] = octanary( xsize )
%OCTANARY construct octanary pattern mask
% the helper function in phase retrieval.
%   xsize : the signal size
%   X: output mask generated by rand 'd=b_1b_2'

B = zeros(xsize);
B1 = rand(xsize);
mask = B1<=0.25;
B(mask)=1;
mask = B1>0.25&B1<=0.5;
B(mask)=-1;
mask = B1>0.5&B1<=0.75;
B(mask)=-1i;
mask = B1>0.75;
B(mask) = 1i;

% C = zeros(xsize);
% C1 = rand(xsize);
% mask = C1<=0.8;
% C(mask) = sqrt(2)/2;
% mask = C1>0.8;
% C(mask) = sqrt(3);
% 
% X = B.*C;
X = B;
end

