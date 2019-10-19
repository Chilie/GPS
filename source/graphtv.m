function [ x ] = graphtv(TV,TVt,d,x)
%GRAPHTV solve linear equation with conjuate gradient method
%   Detailed explanation goes here

n1 = size(d,1);
n2 = size(d,2);
A = @(x) 2*x+ TVt(TV(x));
r = d - A(x);
p = r;
rsold = norm(r,'fro')^2;

for i = 1:numel(d)
%     disp(['iter ' i]);
    Ap = A(p);
    alpha = rsold/sum(sum(p.*Ap));
    x = x+alpha*p;
    r = r-alpha*Ap;
    rsnew = norm(r,'fro')^2;
%     fprintf('residual %.2e\n',sqrt(rsnew));
    if sqrt(rsnew) < 1e-2
        
        break;
    end
    p = r + (rsnew/rsold)*p;
    rsold = rsnew;
end

end

