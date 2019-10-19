function x=forwardSubstitution(U,b)
% Solving an upper triangular system by back-substitution
% Input matrix L is an n by n lower triangular matrix
% Input vector b is n by 1
% Input scalar n specifies the dimensions of the arrays
% Output vector x is the solution to the linear system
% U x = b
% K. Ming Leung, 01/26/03
n = size(U,1);
x=zeros(n,1);
for j=1:1:n
    if (U(j,j)==0) error('Matrix is singular!'); end;
    x(j)=b(j)/U(j,j);
    b(j+1:n)=b(j+1:n)-U(j+1:n,j)*x(j);
end