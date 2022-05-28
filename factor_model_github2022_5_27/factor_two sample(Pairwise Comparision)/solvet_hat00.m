function [t_fdr_hat, T_k, pai1, sigma_hat]...
    = solvet_hat00(Y,X,T0_hat,W0_hat1,W0_hat2,gama)
% function to get FAME statistics, critical value

% Input： 
%   Y, X: generated data followed factor model by users
%   W0_hat1, W0_hat2：the estimation of factor
%   T0_hat: sqrt sample number times sample mean of Y
%   gama: FDR control level

% output：
%   t_fdr_hat：critical value
%   T_k: the mean level after removing factor
%   sigma_hat: the estimation of covariance T_k
%   pai1：the estimation of the proportion of non-nulls

%   Call description:
%   solvet_hat00(Y,X,T0_hat,W0_hat1,W0_hat2,gama): call other functions, 
%   such as pai.m, t_hat.m, and get the results together

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

global n1 n2 h_hat1 h_hat2
%% get mu_hat, B_hat, sigma_hat,
W_piao1 = [ones(1,n1);W0_hat1];
Px1 = (W_piao1')*((W_piao1*(W_piao1'))\W_piao1);%inv(A)*b = A\b
muB_hat1 = Y*(W_piao1')/(W_piao1*(W_piao1'));
B_hat1 = muB_hat1(:,2:end);
sigma_hat1 = Y*(eye(n1)-Px1)*Y'/(n1-h_hat1-1);

W_piao2 = [ones(1,n2);W0_hat2];
Px2 = (W_piao2')*((W_piao2*(W_piao2'))\W_piao2);%inv(A)*b = A\b
muB_hat2 = X*(W_piao2')/(W_piao2*(W_piao2'));
B_hat2 = muB_hat2(:,2:end);
sigma_hat2 = X*(eye(n2)-Px2)*X'/(n2-h_hat2-1);

sigma_hat = sigma_hat1/n1 + sigma_hat2/n2;


%% calculate statistics
T_k = T0_hat - [B_hat1,-B_hat2]*[mean(W0_hat1,2);mean(W0_hat2,2)];

%% for pairwise comparision, we always need add a constant s0
pai1 = pai(T_k, sigma_hat);
zhong = sqrt(diag(sigma_hat));
co = 0.1*mean(prctile(zhong,[0:0.001:100*(1-pai1)]));
sigma_hat = diag((zhong + co).^2);

%% solve pai1 estimate, critic value
pai1 = pai(T_k, sigma_hat);
[t_fdr_hat, ~, ~] = t_hat(pai1,T_k,gama,sigma_hat);