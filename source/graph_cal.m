function [ x,y ] = graph_cal(L,A,c,d)
%graph_cal graph projection with the graph Ax=y
%   Detailed explanation goes here
% w = L\(A*c+A*A'*d);
% y = ctranspose(L)\w;
% x = c+A'*(d-y);

[m,n] = size(A);
if m < n
w = forwardSubstitution(L,(A*(c+A'*d)));
y = backSubstitution(ctranspose(L),w);
x = c+A'*(d-y);

else
% w = L\(c+A'*d);
% x = ctranspose(L)\w;

w = forwardSubstitution(L,(c+A'*d));
x = backSubstitution(ctranspose(L),w);
% 
% % x = (eye(size(A,2)) + A'*A)\(c+A'*d);
% % fprintf('',);
y = A*x;
end

% fprintf('inner rel:Ax-b %.3e\n',norm(A*x-y));


end

