% This is a generic test file for calling the algorithm
% for the low-rank matrix completion problem.
%
% Algorithm: increamental Phase Retrieval Algorithm

% Author : Kee Lee, Beijing Computational Science and Research Center
% 8/8/2017

addpath(genpath('.'));
warning off
clear all; clc; close all;

% Set randstream based on clock
% s = RandStream.create('mt19937ar','seed',sum(100*clock));
% RandStream.setDefaultStream(s);

TestAlg = {'GPSPR-1','GPSPR-0-9','DRPR-1','DRPR-0-9','RAF','TAF','TWF','WirtFlow'}; % to compare true/false.
TestNalg = length(TestAlg);
TestOS = [3 4];%[1:0.1:5];
Ntestos = length(TestOS);

Nprob = 1;
Nsln = 1;
Recov = zeros(Nprob,TestNalg,Ntestos);

Recovredius = zeros(Nprob,Nsln,TestNalg,Ntestos);
% frel = zeros(Nprob,Nsln,1+TestNalg,Ntestos);
Iterhist = zeros(Nprob,Nsln,4,Ntestos);
for nt = 1:Ntestos
    %% Problem specifications
    prob.ndim = 1;
    % prob.type = 'real';
    prob.type = 'comp';
    prob.d1 = 100;
    prob.n = prob.d1;
    prob.d2 = 1;
    
    prob.os = TestOS(nt); % oversampling ratio, 2 is the limit case for real, and 4 for complex case.
    prob.m = round(prob.os*prob.d1);
    prob.addnoise = 0; % whether containing noise
    prob.addnoisetype = 'gaussian'; % 'poisson'
    % input parameter at the initial stage
    
    
    %% Generating data
    
    fprintf('generating data... the dim: %d, the length: %d \n',prob.ndim,prob.d1);
    
    for nk = 1:Nprob
        disp('=================================================================');
        fprintf('==============the prob # %d============================\n',nk);
        if prob.ndim == 1
            fprintf('the measurment m is :%d\n', prob.m);
            if strcmp(prob.type,'real')
                prob.signal = randn(prob.d1,1);
                prob.A = randn(prob.m, prob.d1);
            else
                prob.signal = randn(prob.d1,1)+1i*randn(prob.d1,1);
                prob.A = sqrt(2)/2*(randn(prob.m, prob.d1)+1i*randn(prob.m,prob.d1));
            end
            prob.data = abs(prob.A*prob.signal);
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
        
        %% Initialization for algorithms
        % from the rank one beginning
        
        
        for isln = 1:Nsln
            fprintf('===========================the prob #%d, solution #%d=========\n',nk,isln);
            fprintf('**** First stage: Random initialization\n');
            % the initialization stage value of x
            %             if strcmp(prob.type,'real')
            %                 x0 = randn(prob.d1,1);
            %             else
            %                 x0 = randn(prob.d1,1)+1i*randn(prob.d1,1);
            %             end
            opts.verbose = 0;
            opts.initMethod = 'weighted';%'truncated','amplitude','weighted';
            opts.isComplex = true;
            A = @(x) prob.A*x;
            At = @(x) prob.A'*x;
            b0 = prob.data;
            n = size(prob.signal,1);
            
            % other methods---------------------------------------------------------
            x0 = initX(A, At, b0, n, opts);
            if strcmp(prob.type,'real')
                x0 = real(x0);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Our Method
            prob.x0 = prob.signal;
            opt.y0 = x0;%prob.x0+0.01*randn(prob.n,1);
            
            opt.max_iter = 5000;
            opt.verbosity = 0;
            opt.disp_skip = 20;
            opt.balance_para = 1;
            opt.log = 1; % to trace the progress
            opt.tolx = 1e-10;
            opt.toldiff = -Inf;
            prob.pri = 'none';
            fprintf('++++++++++++++++++GPSPR+++++++++++++++++++++++++++++++++++++++++++++\n');
            %             fprintf('**** First stage: Random initialization\n');
            [x1,output1] = GPSPR(prob,opt);
            
            relredius = output1.trace(end);
            Iterhist(nk,isln,1,nt) = length(output1.trace);
            Recovredius(nk,isln,1,nt) = relredius;
            if relredius < 1e-3
                Recov(nk,1,nt) = Recov(nk,1,nt) + 1;
            end
            opt.balance_para = 0.9;
            fprintf('++++++++++++++++++GPSPR-0-9+++++++++++++++++++++++++++++++++++++++++++++\n');
            [x1,output1] = GPSPR(prob,opt);
            
            relredius = output1.trace(end);
            Iterhist(nk,isln,2,nt) = length(output1.trace);
            Recovredius(nk,isln,2,nt) = relredius;
            if relredius < 1e-3
                Recov(nk,2,nt) = Recov(nk,2,nt) + 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Our Method
            fprintf('++++++++++++++++++DRPR+++++++++++++++++++++++++++++++++++++++++++++\n');
            %             fprintf('**** First stage: Random initialization\n');
            opt.balance_para = 1;
            [x2,output2] = DRPR(prob,opt);
            relredius = output2.trace(end);
            Iterhist(nk,isln,3,nt) = length(output2.trace);
            Recovredius(nk,isln,3,nt) = relredius;
            if relredius < 1e-3
                Recov(nk,3,nt) = Recov(nk,3,nt) + 1;
            end
            fprintf('++++++++++++++++++DRPR-0-9+++++++++++++++++++++++++++++++++++++++++++++\n');
            opt.balance_para = 0.9;
            [x2,output2] = DRPR(prob,opt);
            relredius = output2.trace(end);
            Iterhist(nk,isln,4,nt) = length(output2.trace);
            Recovredius(nk,isln,4,nt) = relredius;
            if relredius < 1e-3
                Recov(nk,4,nt) = Recov(nk,4,nt) + 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for k = 5:TestNalg
                opts.algorithm = TestAlg{k};
                opts.tol = 1e-10;
                %     opts.recordReconErrors = true;
                opts.maxIters = opt.max_iter;
                opts = manageOptions(opts);
                opts.xt = prob.signal;
                [sol, outs] = solveX(A, At, b0, x0, opts);
                innertemp = dot(prob.signal,sol);
                c = innertemp/abs(innertemp);
                relredius = norm(c*prob.signal-sol)/norm(prob.signal);
                Recovredius(nk,isln,k,nt) = relredius;
                if relredius < 1e-3
                    Recov(nk,k,nt) = Recov(nk,k,nt) + 1;
                end
            end
        end
    end
end
% save('Mat.mat','Recov','Recovredius','Iterhist');
