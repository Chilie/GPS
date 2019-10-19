function [ X ] = zerosadd( A, mode, ratio )
% Zerosadd before do fft, zeros-padding for matrix or cell A(element is matrix)
% Input:
%       A - data need to zeros-padding or zeros-cancled, A can be a
%       matrix/vector and cell, whose element is matrix/vector
%       and 
%       mode - 1 , to implement zeros-padding
%              2,  to implement zeros-cancling
%       ratio  , oversampling ratio for fft
%   Kee Lee <liji597760593@126.com> from Peking Univ.
if nargin < 3
    ratio = 2;
end
if isa(A,'numeric')
    switch mode
        case 1
            X = zeros(ratio*size(A));
            X(1:size(A,1),1:size(A,2)) = A;
        case 2
            X = zeros(size(A)/ratio);
            X = A(1:size(X,1),1:size(X,2));
        otherwise
            error('mode == 1 or mode == 2 error,check the input.');
    end       
elseif isa(A,'cell')
        switch mode
        case 1
             X = cellfun(@(x) zerosadd(x, 1, ratio), A, 'UniformOutput',false);
        case 2
                X = cellfun(@(x) zerosadd(x, 2, ratio), A, 'UniformOutput',false);
        otherwise
            error('mode == 1 or mode == 2 error,check the input.');
    end
else
    error('the data structure is not supported with the function.');
end
end

