function [ x, output ] = GPSPR( prob,option )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% the default parameter
opt.max_iter = 5000;
opt.verbosity = 1;
opt.log = 1; % to trace the progress
opt.tolx = 1e-6;
opt.toldiff = 1e-6;
opt.balance_para = 1;

fields =fieldnames(option);
for k = 1:size(fields,1)
   opt = setfield(opt,fields{k},getfield(option,fields{k}));
end

type = prob.type;
n = prob.n;
m = prob.m;
A = prob.A;
data = prob.data;
x0 = prob.x0;
y0 = opt.y0;
verbosity = opt.verbosity;
disp_skip = opt.disp_skip;
max_iter = opt.max_iter;
pri = prob.pri;
% initialization
xk = y0;
yk = A*xk;
lambda = zeros(size(A,2),1);
nu = zeros(size(data));

% balance-para
balance_para = opt.balance_para; % 0< <1;
g = @(x) applypriori(x,pri);
mode = 'nonoise';
f = @(y) prox_ampl(y,data,mode,type,1); % support noise mode: 'nonoise','gaussian','outlier'
% Cholsky decomposition, according to magnitude of m and n
if m>= n
    L = chol(eye(n)+A'*A,'lower');
else
    L = chol(eye(m)+A*A','lower');
end
graph_proj = @ (x,y) graph_cal(L,A,x,y);

% relerr = @(x,y) min(norm(x+y),norm(x-y));
relerr = @(x) computerelerror(x,x0); 
if verbosity
   fprintf(' Iter || relerror_x  || relerror_y ||   DRS  || graph_veri\n'); 
end
if nargout > 1
   output = struct();
end
flag = 0;
if opt.log
   relerror = [];
   reldr = [];
end
for k = 1:max_iter
    % step 1 graph_proj
    tempx = xk-lambda;
    tempy = yk-nu;
    [xkb,ykb] = graph_proj(tempx,tempy);
    if balance_para <= 1
       xkb = (1-balance_para)*tempx + balance_para*xkb;
       ykb = (1-balance_para)*tempy + balance_para*ykb;
    end
%     fprintf('rel:Ax-b %.3e\n',norm(A*xkb-ykb));
% step 2 proximity operator for priopri information and maginitude
    
    xkt = g(xkb+lambda);
    ykt = f(ykb+nu);
    if k > 2 & norm(xkt-xk,'fro')/norm(xk,'fro') < opt.toldiff
        flag = 1;
    end

     xk = xkt;
     yk = ykt;

%     fprintf('rel:Ax-b %.3e\n',norm(abs(yk)-data));
   
% step 3 dual update
    lambda = lambda +xkb-xk;
    nu = nu+ykb-yk; 
    
    rel_x = relerr(xk);
    if opt.log
        
    DRS = computerelerror([xk+lambda;yk+nu],[x0;A*x0]);
    end
    if verbosity & mod(k,disp_skip)==0
        rel_y = norm(abs(yk)-data);
        
        graph_ver = norm(A*xkb-ykb);
     fprintf('  %5d  || %.2e  || %.2e  || %.2e  || %.2e\n',k,rel_x,rel_y,DRS,graph_ver);
    end
    if rel_x < opt.tolx || flag == 1
        break;
    end
        
    if opt.log
       relerror = [relerror rel_x]; 
       reldr = [reldr DRS];
    end
end
x = xk;
if nargout >1
    if k< max_iter
        output.msg = 'break due to epsilon';
    else
        output.msg = 'attained maximum iteration';
    end
    output.trace = relerror;
    output.tracedr = reldr;
end

