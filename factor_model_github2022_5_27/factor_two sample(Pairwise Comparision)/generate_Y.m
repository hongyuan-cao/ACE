function [Y,X] = generate_Y(p, n1 ,n2, h, h2 ,berlii)
% function of generating data followed factor model

% Input： 
%   p: number of variables
%   n1, n2: number of samples
%   h, h2：number of factors
%   berlii: A vector only containing 0 and 1. 0 represents null，1 represents
%   non-null

% output：
%   Y, X：generate data followed factor model

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

%% generate data followed factor model
Y = zeros(p,n1);
X = zeros(p,n2);
mu = zeros(p,1);
mu2 = zeros(p,1);
index_pai1 = find(berlii==1);
mu(index_pai1) = unifrnd(0.3, 0.6, length(index_pai1), 1);%%%%%%%%%%

B = 2*rand(p,h) - 1;
B2 = 2*rand(p,h2) - 1;
%% independent
for jy = 1:n1
    w = randn(h,1);
    e = randn(p,1);
    Y(:,jy) = mu + B*w + e;%
end
for j = 1:n2
    w2 = randn(h2,1);
    e2 = randn(p,1);
    X(:,j) = mu2 + B2*w2 + e2;%
end

%% AR
% mu_e = zeros(1, p);
% load('sigma_e.mat');
% % sigma_e = sigma_e(1:200,1:200);%%% for p = 200
% ee = mvnrnd(mu_e,sigma_e,n1); ee2 = mvnrnd(mu_e,sigma_e,n2);
% ee = ee'; ee2 = ee2';
% for j = 1:n1  
%     w = randn(h,1);  
%     e = ee(1:end,j);   
%     Y(:,j) = mu + B*w + e;%
% end
% for j = 1:n2  
%     w2 = randn(h,1);  
%     e2 = ee2(1:end,j);   
%     X(:,j) = mu2 + B2*w2 + e2;%
% end