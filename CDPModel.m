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



%% Problem specifications
prob.ndim = 2;
% prob.type = 'real';
prob.type = 'comp';
% prob.model = 'ls'; % 'ls','map'
prob.d1 = 128;
prob.d2 = 128;
prob.os = 2; % oversampling ratio, 2 is the limit case for real, and 4 for complex case.
prob.addnoise = 0; % whether containing noise
prob.addnoisetype = 'gaussian'; % 'poisson'
% input parameter at the initial stage
prioris = {'real-valuedness','nonnegativity'};%e.g. {'sparsity',40,'real-valuedness','nonnegativity','support'};
% lambda0 = 0; % the initial point of regualrization path 
maxiter = 1200; % +INF; % Maximum iterations in the inner subproblem
% epsilon = 100; % 'ls' 1 ,'map' 1e-1% the samllest eigenvalue gap to zero tolerance.
%% Generating data

fprintf('generating data... the dim: %d, the length: %d \n',prob.ndim,prob.d1);

verbosity = 0; % Must be an integer; Display information

fprintf('the measurement is sampling L: %d',prob.os);
if strcmp(prob.type,'real')
%         prob.signal = load('camera128.mat','-ascii');%randn(prob.d1,prob.d2);
        prob.signal = randn(prob.d1,prob.d2);
 elseif strcmp(prob.type,'comp')
        prob.signal = randn(prob.d1,prob.d2) + 1i*randn(prob.d1,prob.d2);
        prob.signal = load('camera128.mat','-ascii');%randn(prob.d1,prob.d2);
end

if prob.d2 == 1
    prob.data = zeros(prob.d1,prob.os);
    prob.mask = zeros(prob.d1,prob.os);
    [xSize,ySize] = size(prob.signal);
    for k = 1:prob.os
        if k == 1
            prob.mask(:,k) = ones([xSize,ySize]);
        else
            prob.mask(:,k) = octanary([xSize,ySize]);
        end
        prob.data(:,k) = abs(fft(prob.mask(:,k).*prob.signal))/sqrt(xSize*ySize);
    end
else
    prob.data = zeros(prob.d1,prob.d2,prob.os);
    prob.mask = zeros(prob.d1,prob.d2,prob.os);
    [xSize,ySize] = size(prob.signal);
    for k = 1:prob.os
        if k == 1
            prob.mask(:,:,k) = ones([xSize,ySize]); % octanary
        else
            prob.mask(:,:,k) = octanary([xSize,ySize]);
        end
        prob.data(:,:,k) = abs(fft2(prob.mask(:,:,k).*prob.signal))/sqrt(xSize*ySize);
    end
end


if prob.addnoise
    SNR = prob.addnoise;
    for k = 1:prob.os
        if prob.d2 ==1
            if strcmp(prob.addnoisetype, 'gaussian')
                prob.data(:,k) = max(awgn(prob.data(:,k),SNR,'measured'),realmin);
            else
                alpha = (norm(prob.data(:,k),'fro')^2)/(sum(prob.data(:,k))*(10^(SNR/10)));
                prob.data(:,k) = max(alpha*poissrnd(prob.data(:,k)/alpha),realmin);
            end
        else
            if strcmp(prob.addnoisetype, 'gaussian')
                prob.data(:,:,k) = max(awgn(prob.data(:,:,k),SNR,'measured'),realmin);
            else
                alpha = (norm(prob.data(:,:,k),'fro')^2)/(sum(sum((prob.data(:,:,k))))*(10^(SNR/10)));
                prob.data(:,:,k) = max(alpha*poissrnd(prob.data(:,:,k)/alpha),realmin);
            end
        end
    end
end

%% Initialization for algorithms
% from the rank one beginning


    
fprintf('**** First stage: Random initialization\n');
% the initialization stage value of x
if strcmp(prob.type,'real')
      x0 = randi(100,prob.d1,prob.d2);% randn(prob.d1,prob.d2);
%     x = randn(prob.d1,prob.d2);
else
    x0 = randn(prob.d1,prob.d2)+1i*randn(prob.d1,prob.d2);
%     x0 = (1+0.01)*prob.signal;
end

prob.x0 = prob.signal;
opt.y0 = x0;%prob.x0+0.01*randn(prob.n,1);

opt.max_iter = 400;
opt.verbosity = 1;
opt.disp_skip = 20;
opt.balance_para = 1;
opt.log = 1; % to trace the progress
opt.tolx = 1e-10;
opt.toldiff = -Inf;
prob.pri = prioris;
fprintf('++++++++++++++++++DRSPR+++++++++++++++++++++++++++++++++++++++++++++\n');
%             fprintf('**** First stage: Random initialization\n');
[x1,output1] = GPSCDPPR(prob,opt);
