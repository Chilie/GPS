function [ x, output ] = DRPR( prob,option )
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

% balance-para
balance_para = opt.balance_para; % 0< <1;
mode = 'nonoise';
f = @(y) prox_ampl(y,data,mode,type,1); % support noise mode: 'nonoise','gaussian','outlier'
% Cholsky decomposition, according to magnitude of m and n
% Cholsky decomposition, according to magnitude of m and n
if m>= n
    L = chol(A'*A,'lower');
else
    L = chol(A*A','lower');
end
proj_Mpart = @ (x) projy_cal_pre(L,A,x);

% relerr = @(x,y) min(norm(x+y),norm(x-y));
relerr = @(x) computerelerror(x,x0); 
if verbosity
   fprintf(' Iter || relerror_x  || relerror_y \n'); 
end
if nargout > 1
   output = struct();
end
flag = 0;
if opt.log
   relerror = [];
end
for k = 1:max_iter
    % step 1 graph_proj
        ykb = f(yk);
    medy0 = 2*ykb-yk;
    tempmedy = proj_Mpart(medy0);
    medy = A*tempmedy;
    if balance_para < 1
        medy = (1-balance_para)*(medy0) + balance_para*medy;
    end
    xk = proj_Mpart(yk);
    yk = yk + medy -ykb;
    xkt = proj_Mpart(yk);
%     fprintf('toldiff: %.2e\n',norm(ykt-yk,'fro')/norm(yk,'fro'))
    if k > 2 & norm(xkt-xk,'fro')/norm(xk,'fro') < opt.toldiff
        flag = 1;
    end
%      yk = ykt;
     rel_x = relerr(xkt);
    if verbosity & mod(k,disp_skip)==0
        
        rel_y = norm(abs(yk)-data);
     fprintf('  %5d  || %.2e  || %.2e \n',k,rel_x,rel_y);
    end
    
    if rel_x < opt.tolx || flag == 1
        break;
    end
        
    if opt.log
       relerror = [relerror rel_x]; 
    end
end
x = xkt;
if nargout >1
    if k< max_iter
        output.msg = 'break due to epsilon';
    else
        output.msg = 'attained maximum iteration';
    end
    output.trace = relerror;
end

