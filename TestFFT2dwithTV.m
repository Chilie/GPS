% This is a generic test file for calling the algorithm
% for the low-rank matrix completion problem.
%
% Algorithm: Alternating Projection methods via random projection method to
% approximate the low rank matrix

% To accerlate the computation, the observation constraints are implemented
% by mex file updateM.

% Author : Kee Lee, Beijing Computational Science and Research Center
% 8/8/2017


clear all; clc; close all;

% Set randstream based on clock
% s = RandStream.create('mt19937ar','seed',sum(100*clock));
% RandStream.setDefaultStream(s);


%% Problem specifications
prob.ndim = 2;
prob.type = 'comp';
prob.addnoise = 0; % whether containing noise
prob.addnoisetype = 'gaussian'; % 'poisson'
% input parameter at the initial stage
prioris = {'real-valuedness','nonnegativity','support'};%e.g. {'sparsity',40,'real-valuedness','nonnegativity','support'};
maxiter = 200; % +INF; % Maximum iterations in the inner subproblem
% epsilon = 100; % 'ls' 1 ,'map' 1e-1% the samllest eigenvalue gap to zero tolerance.
%% Generating data

if strcmp(prob.type,'real')
        prob.signal = load('img_data/camera128.mat','-ascii');
        [s1,s2] = size(prob.signal);
        prob.d1 = s1;
        prob.d2 = s2;
        prob.os = 2;
 elseif strcmp(prob.type,'comp')
        prob.signal =  load('camera128.mat','-ascii');%randi(10,prob.d1,prob.d2);
load('s.mat');
prob.signal =  s;
        prob.signal = load('img_data/camera128.mat','-ascii');
        [s1,s2] = size(prob.signal);
        prob.d1 = s1;
        prob.d2 = s2;
        prob.os = 2;
end

if prob.d2 == 1
    prob.data = zeros(prob.os*prob.d1,1);
    prob.tsignal = zerosadd(prob.signal,1,prob.os);
    [xSize,ySize] = size(prob.tsignal);
    prob.data = abs(fft(prob.tsignal))/sqrt(xSize*ySize);
else
    prob.data = zeros(prob.os*prob.d1,prob.os*prob.d2);
    prob.tsignal = zerosadd(prob.signal,1,prob.os);
    [xSize,ySize] = size(prob.tsignal);
    prob.data = abs(fft2(prob.tsignal))/sqrt(xSize*ySize);
end


if prob.addnoise
    SNR = prob.addnoise;
   if strcmp(prob.addnoisetype, 'gaussian')
       prob.data= max(awgn(prob.data,SNR,'measured'),realmin);
   else
        alpha = (norm(prob.data,'fro')^2)/(sum(sum((prob.data)))*(10^(SNR/10)));
        prob.data= max(alpha*poissrnd(prob.data/alpha),realmin);
   end
end

%% Initialization for algorithms
% from the rank one beginning

Npeat = 1;
Recov = struct();

for repeat = 1:Npeat
    
fprintf('**** First stage: Random initialization\n');
fprintf('The Iteration: %d\n------------',repeat);
% the initialization stage value of x
if strcmp(prob.type,'real')
    x0 = zerosadd(randi(100,size(prob.signal)),1,prob.os);%randi(10,prob.d1,prob.d2);
else
    x0 = zerosadd(randi(100,size(prob.signal)),1,prob.os);%randi(100,prob.d1,prob.d2)+1i*randi(100,prob.d1,prob.d2);
end

prob.x0 = prob.signal;
opt.y0 = x0;%prob.x0+0.01*prob.x0;%randn(prob.n,1);

opt.max_iter = 200;
opt.verbosity = 1;
opt.disp_skip = 20;
opt.balance_para = 1;
opt.tvp = 10;
opt.log = 1; % to trace the progress
opt.tolx = 1e-10;
opt.toldiff = -Inf;
prob.pri = prioris;
fprintf('++++++++++++++++++DRSPR+++++++++++++++++++++++++++++++++++++++++++++\n');
%             fprintf('**** First stage: Random initialization\n');
% [x1,output1] = hioNew(prob,opt);
[x1,output1] = hioNew(prob,opt);
% [x1,output1] = DRSFFTPR(prob,opt);
opt.y0 = x1;
opt.balance_para = 0.9;
opt.max_iter = 300;
[x2,output2] = DRSFFTTVPR(prob,opt);

Recov(repeat).X1 = output1.x(end).X;
Recov(repeat).X2 = output2.x(end).X;
% pause;
% x1 = rot90(output1.x(end).X,2);
% x2 = rot90(output2.x(end).X,2);
% x1 = output1.x(end).X;
% x2 = output2.x(end).X;
% relhio = norm(x1-prob.signal,'fro')/norm(prob.signal,'fro')
% reltv = norm(x2-prob.signal,'fro')/norm(prob.signal,'fro')

% [z]=hioNew(prob,x0);        
% 
% RecoveryTace1(repeat).second(:,:,end) = reshape(y,prob.d1,prob.d2);
% 
% [z]=hioNew(prob.signal,x0);
% RecoveryTace2(repeat).final = reshape(z,prob.d1,prob.d2);
end
