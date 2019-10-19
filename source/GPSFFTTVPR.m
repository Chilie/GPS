function [ x, output ] = GPSFFTTVPR( prob,option )
%DRSFFTTVPR routine for Fourier phase retrieval with TV minimization term
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
[xSize,ySize] = size(prob.tsignal);

A = @(x) fft2(x)/sqrt(xSize*ySize);
At = @(y) sqrt(xSize*ySize)*ifft2(y);
data = prob.data;
x0 = prob.x0;
y0 = opt.y0;
verbosity = opt.verbosity;
disp_skip = opt.disp_skip;
max_iter = opt.max_iter;
pri = prob.pri;
% initialization
xk = y0;
yk = A(xk);
lambda = zeros([xSize,ySize]);
nu = zeros(size(data));
nu2 = zeros([xSize,ySize]);

% balance-para
balance_para = opt.balance_para; % 0< <1;
tvp = opt.tvp; % parameter to control the tv.
g = @(x) applypriori(x,pri);
mode = 'nonoise';
fintprox = @(y) prox_ampl(y,data,mode,type,1); % support noise mode: 'nonoise','gaussian','outlier'
ftvm = @(y) exp(1i*angle(y)).*max(abs(y)-tvp,0);
%%% Define the TV forward and conjuate operator in matrix form
n1 = xSize;
n2 = ySize;
mat = @(x) reshape(x,n1,n2);
action = 'variation'; % 'variation'
if strcmp(action,'matrix')
    e = ones(max(n1,n2),1);
    e2 = e;
    e2(n2:end) = 0;
    J = spdiags([-e2,e], 0:1,n2,n2);
    I = eye(n1);
    Dh = kron(J,I);  % horizontal differences, sparse matrix
    % see also blkdiag
    
    e2 = e;
    e2(n1:end) = 0;
    J = spdiags([-e2,e], 0:1,n1,n1);
    I = eye(n2);
    Dv = kron(I,J);  % vertical differences, sparse matrix
    
    op = Dh + 1i*Dv;
elseif strcmp(action,'variation')
    Dh     = @(X) [diff(X,1,2),  zeros(n1,1)];
    diff_h = @(X) [zeros(n1,1),X(:,1:end-1)] - [X(:,1:end-1),zeros(n1,1) ];
    Dv     = @(X) [diff(X,1,1); zeros(1,n2)];
    %         diff_v = @(X) [zeros(1,n2);X(1:end-1,:)] - [X(1:end-1,:);zeros(1,n2) ];
    % sometimes diff_v is much slower than diff_h
    % We can exploit data locality by working with transposes
    diff_v_t = @(Xt) ([zeros(n2,1),Xt(:,1:end-1)] - [Xt(:,1:end-1),zeros(n2,1)])';
    
    Dh_transpose = @(X)      diff_h(mat(X))  ;
    %     Dv_transpose = @(X)      diff_v(mat(X))  ;
    Dv_transpose = @(X)      diff_v_t(mat(X)')  ; % faster
    
    TV  = @(x) ( Dh(mat(x)) + 1i*Dv(mat(x)) );     % real to complex
    TVt = @(z) ( Dh_transpose(real(z)) + Dv_transpose(imag(z)) );
    
end

zk = TV(xk);
% if verbosity
%    fprintf(' Iter || relerror_x  || relerror_y ||   DRS  || graph_veri\n');
% end
if nargout > 1
    output = struct();
end
flag = 0;
if opt.log
    iterX = struct();
end
for k = 1:max_iter
    % step 1 graph_proj
    tempx = xk-lambda;
    tempy = yk-nu;
    tempz = zk-nu2;
    d = tempx + TVt(tempz) + real(At(tempy));
    xkb= graphtv(TV,TVt,d,tempx);
    ykb = A(xkb);
    zkb = TV(xkb);
    if balance_para <= 1
        xkb = (1-balance_para)*tempx + balance_para*xkb;
        ykb = (1-balance_para)*tempy + balance_para*ykb;
        zkb = (1-balance_para)*tempz + balance_para*zkb;
    end
    %     fprintf('rel:Ax-b %.3e\n',norm(A*xkb-ykb));
    % step 2 proximity operator for priopri information and maginitude
    
    xkt = g(xkb+lambda);
    ykt = fintprox(ykb+nu);
    zkt = ftvm(zkb+nu2);
    if k > 2 & norm(xkt-xk,'fro')/norm(xk,'fro') < opt.toldiff
        flag = 1;
    end
    
    xk = xkt;
    yk = ykt;
    zk = zkt;
    
    %     fprintf('rel:Ax-b %.3e\n',norm(abs(yk)-data));
    
    % step 3 dual update
    lambda = lambda +xkb-xk;
    nu = nu+ykb-yk;
    nu2 = nu2 +zkb-zk;
    %     rel_x = relerr(xk);
    fprintf('Iteration  %5d\n',k);
    
    if opt.log
        %        relerror = [relerror rel_x];
        iterX(k).X = zerosadd(xk,2,prob.os);
    end
end
x = xk;
if nargout >1
    if k< max_iter
        output.msg = 'break due to epsilon';
    else
        output.msg = 'attained maximum iteration';
    end
    %     output.trace = relerror;
    output.x = iterX;
end

