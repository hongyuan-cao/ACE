function [t_fdr_hat, T_k, pai1, sigma_hat]...
    = solvet_hat00(Z,W0_hat,T0_hat,gama)
global n h_hat
%% get mu_hat, B_hat, sigma_hat,
%mean(W0_hat,2)
W_piao = [ones(1,n);W0_hat];
Px = (W_piao')*((W_piao*(W_piao'))\W_piao);%inv(A)*b = A\b
muB_hat = Z*(W_piao')/(W_piao*(W_piao'));
B_hat = muB_hat(:,2:end);
sigma_hat = Z*(eye(n)-Px)*Z'/(n-h_hat-1);

%% calculate statistics
T_k = T0_hat - B_hat*sum(W0_hat,2)/sqrt(n);

%% for small p, we need add a constant s0
if length(T_k) < 1000
    pai1 = pai(T_k, sigma_hat);
    zhong = sqrt(diag(sigma_hat));
    co = 0.1*mean(prctile(zhong,[0:0.001:100*(1-pai1)]));
    sigma_hat = diag((zhong + co).^2);
end

%% solve pai1
pai1 = pai(T_k, sigma_hat);
%% solve critical value
[t_fdr_hat, ~, ~] = t_hat(pai1,T_k,gama,sigma_hat);