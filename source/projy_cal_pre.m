function [ x ] = projy_cal_pre(L,A,c)
%graph_cal graph projection with the graph Ax=y
%   Detailed explanation goes here
% w = L\(A*c+A*A'*d);
% y = ctranspose(L)\w;
% x = c+A'*(d-y);

[m,n] = size(A);
if m < n
w = forwardSubstitution(L,c);
y = backSubstitution(ctranspose(L),w);
x = A'*y;

else
% w = L\(c+A'*d);
% x = ctranspose(L)\w;

w = forwardSubstitution(L,A'*c);
x = backSubstitution(ctranspose(L),w);
% 
% % x = (eye(size(A,2)) + A'*A)\(c+A'*d);
% % fprintf('',);
end

% fprintf('inner rel:Ax-b %.3e\n',norm(A*x-y));


end

