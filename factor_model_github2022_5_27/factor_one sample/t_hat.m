function [t_fdr_hat,t_fval,index] = t_hat(pai1,T_k,sigma_hat)
% function to get critical value

% Input： 
%   pai1：the estimation of the proportion of non-nulls
%   T_k: the mean level after removing factor
%   sigma_hat: the estimation of covariance of T_k

% output：
%   t_fdr_hat：critical value
%   t_fval: estimated FDP given t_fdr_hat
%   index: transition variables in the calculation of t_fdr_hat

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

global gama
%% solve critical value
aaaaaa = 0.01:0.001:5;
f_t_hat = @(t_fdr_hat) 2*(1-pai1)*cdf('Normal',-t_fdr_hat,0,1) - ...
    gama.*mean(bsxfun(@ge,abs(T_k./sqrt(diag(sigma_hat))),t_fdr_hat),1); 

f_t_hatt = f_t_hat(aaaaaa);
[~,index] = find(f_t_hatt<=0);

if isempty(index) == 1
    t_fdr_hat = Inf;
    t_fval = 0;
else
    index = index(1);
    t_fdr_hat = aaaaaa(index);
    t_fval = f_t_hat(t_fdr_hat);
end