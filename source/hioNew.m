function [x,output] = hioNew(prob,option)
% HIO Hybird Input-Output method used in PR,
% here it is used to generate the initial image.
% Input:
%       data: structure of data info
% Output:
%       x: return by the HIO algorithm.

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

A = @(x) fft2(x);
data = prob.data;
x0 = prob.x0;
yk = opt.y0;
verbosity = opt.verbosity;
disp_skip = opt.disp_skip;
max_iter = opt.max_iter;
pri = prob.pri;
% initialization
% xk = y0;
% yk = A(xk);
lambda = zeros([xSize,ySize]);


% balance-para
balance_para = opt.balance_para; % 0< <1;
g = @(x) applypriori(x,pri);
mode = 'nonoise';
f = @(y) sqrt(xSize*ySize)*ifft2(prox_ampl(y,data,mode,type,1));

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
    temp = yk-lambda;
    xk = g(temp);
    
    yk = f(A(xk+lambda));
    lambda = lambda + balance_para*(xk-yk);
    fprintf('Iteration %5d\n',k);
    
    if opt.log
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
    output.x = iterX;
end


% rel = norm(x-z,'fro')/norm(x,'fro');
% h = figure;
% imagesc(z);
% colormap(gray(256));
% axis off;
% tightfig(h);
% print(h, '-depsc2',['mapflip',num2str(kk), '.eps']);
end

