% An example for our FAME method in two-sample test(Pairwise Comparision) 
% which is a less powerful way, we do not suggest to use it. 

% Input£º 
%   p: number of variables
%   n1, n2: number of samples
%   h, h2£ºnumber of factors
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
%   p n1 n2 true_pai1: need to be set as global variables in MATLAB

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

clear;clc;close all
tic
global p n1 n2 true_pai1                                                         
p = 2000; n1 = 200; n2 = 200; h = 2; h2 = 2; gama = 0.05;
P = 0.2;%%%%% true pai1
for jj = 1:100
    berlii = binornd(1,P,[p,1]);%0 represents H0£¬1 represents H1
    true_pai1 = sum(berlii==1)/p;  
    [pai1,FDP_hat,NDR,true_FDP,R] = m0(h, h2,gama,berlii);
    
    ppai1(jj) = pai1; FFDP_hat(jj) = FDP_hat; NNDR(jj) = NDR;
    ttrue_FDP(jj) = true_FDP; RR(jj) = R;
end
pai1 = mean(ppai1); FDR_hat = mean(FFDP_hat); power = mean(NNDR);
true_FDR = mean(ttrue_FDP); R = mean(RR);
clc
fprintf('pai1 = %f\n FDR_hat = %f\n true_FDR = %f\n R = %f\n power = %f\n',...
    pai1,FDR_hat,true_FDR,R,power);
toc