function [ y ] = prox_ampl(x,data,mode,type,regu_para)
%prox_amplitude calculate the proximity of amplitude arised in phase
%retrieval.
% mode varies with the noise model, and regu_para is the regularization
% parameter.
% data is the measurement amplitude.

% flag = 1; % 1 is for real

if strcmp(type,'real')
    switch mode
        case 'nonoise'
            y = real(data.*exp(1i*angle(x)));
%             y = A*x0;
        case 'gaussian'
            y = real((data+regu_para*abs(x)).*exp(1i*angle(x)))/(1+regu_para);
        case 'outlier'
            y = real((data+max(data-abs(x)-1/regu_para,0)).*exp(1i*angle(x)));
    end
elseif strcmp(type,'comp')
    switch mode
        case 'nonoise'
            y = data.*exp(1i*angle(x));
        case 'gaussian'
            y = (data+regu_para*abs(x)).*exp(1i*angle(x))/(1+regu_para);
        case 'outlier'
            y = (data+max(data-abs(x)-1/regu_para,0).*exp(1i*angle(x)));
    end
end
end

