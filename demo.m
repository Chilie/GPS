% this is a demo for 1-d Gaussian phase retrieval via GPS
% Author Li Ji <matliji@nus.edu.sg>
clc;clear;
% parameter setting
prob.ndim = 1; % 1-d PR or 2-d PR
prob.type = 'real'; % real Gaussian or complex Gaussian
prob.d1 = 400; % the length of signal
prob.n = prob.d1; 
prob.os = 1.7; % oversampling ratio
prob.m = round(prob.os*prob.d1);
prob.addnoise = 0; % noise SNR, 0 for the case of noiseless
prob.addnoisetype = 'gaussian'; % noise type
% a-priori information of signal
prioris = {};%e.g. {'sparsity',40,'real-valuedness','nonnegativity','support'};
% generate gaussian random signal and data w/ noise
if prob.ndim == 1
    fprintf('the measurment m is :%d\n', prob.m);
    if strcmp(prob.type,'real')
        prob.x0 = randn(prob.d1,1);
        prob.A = randn(prob.m, prob.d1);
    else
        prob.x0 = randn(prob.d1,1)+1i*randn(prob.d1,1);
        prob.A = sqrt(2)/2*(randn(prob.m, prob.d1)+1i*randn(prob.m,prob.d1));
    end
    prob.data = abs(prob.A*prob.x0);
end
if prob.addnoise
    SNR = prob.addnoise;
    if strcmp(prob.addnoisetype, 'gaussian')
        prob.data = max(awgn(prob.data,SNR,'measured'),realmin);
    else
        alpha = (norm(prob.data,'fro')^2)/(sum(prob.data(:))*(10^(SNR/10)));
        prob.data = max(alpha*poissrnd(prob.data/alpha),realmin);
    end
end
% solving parameter
opt.max_iter = 5000;
opt.verbosity = 1;
opt.disp_skip = 20;
opt.balance_para = 1; % the parameter corresponding to (1-t) in paper, 1 for GPS, <1 for RGPS.
opt.log = 1; % to trace the progress
opt.tolx = 1e-10;
opt.toldiff = -Inf;
prob.pri = prioris;

if strcmp(prob.type,'real')
    x0 = randn(prob.d1,1);
else
    x0 = randn(prob.d1,1)+1i*randn(prob.d1,1);
end

opt.y0 = x0;%prob.x0+0.01*randn(prob.n,1);
fprintf('++++++++++++++++++DRSPR+++++++++++++++++++++++++++++++++++++++++++++\n');
[x2,output2] = GPSPR(prob,opt);
