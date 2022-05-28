% An example for our FAME method in one-sample test

% Input£º 
%   p: number of variables
%   n: number of samples
%   h£ºnumber of factors
%   gama: FDR control level
%   P: the proportion of non-nulls

% output£º
%   Y£ºgenerate data followed factor model
%   pai1£ºthe estimation of the proportion of non-nulls
%   true_FDR: true FDR
%   FDR_hat: the estimation of the FDR
%   R: the number of rejections
%   power: truly rejections over total non-nulls, i.e. power

%   Call description:
%   p n true_pai1 h gama: need to be set as global variables in MATLAB

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

clear;clc;close all
global p n true_pai1 h gama

tic
p = 2000; n = 200; h = 3; gama = 0.05; P = 0.2;
for jj = 1:10
    %% generate data followed factor model
    Y = zeros(p,n); mu = zeros(p,1);
    berlii = binornd(1,P,[p,1]);%0 represents H0£¬1 represents H1
    true_pai1 = sum(berlii==1)/p;
    index_pai1 = find(berlii==1);
    mu(index_pai1) = unifrnd(0.2, 0.5, length(index_pai1), 1);    
    B = 2*rand(p,h) - 1;
    for j = 1:n
        w = randn(h,1);
        e = randn(p,1);
        Y(:,j) = mu + B*w + e;%
    end
    %% main function, get results, such as FDR, the number of rejections and power
    [pai1,FDP_hat,NDR,true_FDP,R,S] = m0(Y,berlii);
    
    ppai1(jj) = pai1; FFDP_hat(jj) = FDP_hat; NNDR(jj) = NDR;
    ttrue_FDP(jj) = true_FDP; RR(jj) = R; SS(jj) = S;
end
pai1 = mean(ppai1); FDR_hat = mean(FFDP_hat); power = mean(NNDR);
true_FDR = mean(ttrue_FDP); R = mean(RR); S = mean(SS);
clc
fprintf('pai1 = %f\n FDR_hat = %f\n true_FDR = %f\n R = %f\n S = %f\n power = %f\n',...
    pai1,FDR_hat,true_FDR,R,S,power);
toc