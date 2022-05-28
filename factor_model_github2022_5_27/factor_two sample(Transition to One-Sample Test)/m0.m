function [pai1,FDP_hat,NDR,true_FDP,R] = ...
    m0(h, h2,gama,berlii)
% main function to get all important variables we need 

% Input： 
%   h, h2：number of factors
%   gama: FDR control level
%   berlii: A vector only containing 0 and 1. 0 represents null，1 represents
%   non-null

% output：
%   pai1：the estimation of the proportion of non-nulls
%   FDP_hat: the estimation of the FDP
%   true_FDP: true FDP
%   S: the number of truly rejections
%   R: the number of rejections
%   NDR: truly rejections over total non-nulls, i.e. power

%   Call description:
%   m0(h,h2,Y,berlii): call other functions, such as solvet_hat00.m, solveW.m,
%   and get the results together

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

global p n1 n2 Y n
n = n1;
[Y,Z,X] = generate_Y(p, n1, n2 , h , h2 ,berlii);
Yba = mean(Y,2);
T = (sqrt(n)* Yba);

%% get factor estimation
deltaaZ = cov(Z'); deltaaX = cov(X');
deltaa = deltaaZ + n1*deltaaX/n2;
% deltaa = cov(Y');%%%%%%%%%%%%%%%
T0_hat = T;
W0_hat = solveW(Y,deltaa);
%% get CK's pai1 estimation, Statistics, sigma_hat
[t_fdr_hat, T_k, pai1, sigma_hat] = solvet_hat00(Y,W0_hat,T0_hat,gama);

%% get FDP_hat, number of rejection and power
R = sum((abs(T_k./sqrt(diag(sigma_hat)))>=t_fdr_hat));

FDP_hat = 2*(1-pai1)*cdf('Normal',-t_fdr_hat,0,1)/(R/p);

index_0 = berlii==0;
bbb = sqrt(diag(sigma_hat));

[~, false_reject]=find(abs(T_k(index_0)./...
    bbb(index_0)) >= t_fdr_hat);

if R == 0
    S = 0; true_FDP = 0
else
    S = R - length(false_reject);  true_FDP = length(false_reject)/R
end
NDR = S/sum(berlii==1);
end