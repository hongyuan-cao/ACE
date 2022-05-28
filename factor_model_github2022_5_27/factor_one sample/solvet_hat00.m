function [t_fdr_hat, T_k, pai1, sigma_hat]...
    = solvet_hat00(Y,W0_hat,T0_hat)
% function to get FAME statistics, critical value

% Input： 
%   Y: generated data followed factor model by users
%   W0_hat：the estimation of factor
%   T0_hat: sqrt sample number times sample mean of Y

% output：
%   t_fdr_hat：critical value
%   T_k: the mean level after removing factor
%   sigma_hat: the estimation of covariance T_k
%   pai1：the estimation of the proportion of non-nulls

%   Call description:
%   solvet_hat00(Y,W0_hat,T0_hat): call other functions, such as pai.m, t_hat.m,
%   and get the results together

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

global n h_hat
%% get mu_hat, B_hat, sigma_hat,
W_piao = [ones(1,n);W0_hat];
Px = (W_piao')*((W_piao*(W_piao'))\W_piao);%inv(A)*b = A\b
muB_hat = Y*(W_piao')/(W_piao*(W_piao'));
B_hat = muB_hat(:,2:end);
sigma_hat = Y*(eye(n)-Px)*Y'/(n-h_hat-1);

%% calculate statistics
T_k = T0_hat - B_hat*sum(W0_hat,2)/sqrt(n);

%% for small p, we need add a constant s0
if length(T_k) < 1000
    pai1 = pai(T_k, sigma_hat);
    zhong = sqrt(diag(sigma_hat));
    co = 0.1*mean(prctile(zhong,[0:0.001:100*(1-pai1)]));
    sigma_hat = diag((zhong + co).^2);
end

%% solve \pi_1
pai1 = pai(T_k, sigma_hat);
%% solve critical value
[t_fdr_hat, ~, ~] = t_hat(pai1,T_k,sigma_hat);