function [W0_hat1,W0_hat2]  = solveW(Y, X ,deltaa1,deltaa2)
% estimate factor model

% Input： 
%   Y, X: generated data followed factor model by users
%   deltaa1, deltaa2: estimation of covariance Y and X

% output：
%   W0_hat1, W0_hat2：the estimation of factor

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

%% estimate factor loading matrix B_hat, factor W0_hat
global n1 n2 h_hat1 h_hat2
%%
[gama_norm1, lambdaa1] = svd(deltaa1);
lam_sort1 = diag(lambdaa1);
[gama_norm2, lambdaa2] = svd(deltaa2);
lam_sort2 = diag(lambdaa2);
%% get factor number estimation h_hat
if ~(isreal(lam_sort1)) || (~(isreal(lam_sort2)))
    error('特征根是复数')
end
% bizhi = lam_sort(1:(end-1))./lam_sort(2:end);
h_max = 10;
bizhi1 = lam_sort1(1:(h_max-1))./lam_sort1(2:h_max);
[~,h_hat1] = max(bizhi1);
bizhi2 = lam_sort2(1:(h_max-1))./lam_sort2(2:h_max);
[~,h_hat2] = max(bizhi2);
%% estimate factor loading matrix B_hat1,2
B_hat1 = gama_norm1(:,1:h_hat1)*diag(sqrt(lam_sort1(1:h_hat1)));
B_hat2 = gama_norm2(:,1:h_hat2)*diag(sqrt(lam_sort2(1:h_hat2)));
%% get W0_hat，
x01 = zeros(h_hat1,1); lb1 = repmat(-3,h_hat1,1); ub1 = repmat(3,h_hat1,1);
x02 = zeros(h_hat2,1); lb2 = repmat(-3,h_hat2,1); ub2 = repmat(3,h_hat2,1);

Y_mean = mean(Y,2); X_mean = mean(X,2);
eee1 = (Y_mean + X_mean)/2;
for j1 = 1:n1
    objfun1 =@(W) sum(abs(Y(:,j1) - eee1 - B_hat1*W));
    problem1 = createOptimProblem('fmincon','objective',objfun1,'x0',x01,'lb',lb1,...
        'ub',ub1,'options',optimset('Algorithm','SQP','Disp','none'));
    gs1 = GlobalSearch;
    W0_hat1(:,j1) = run(gs1,problem1);
    
    objfun2 =@(W) sum(abs(X(:,j1) - eee1 - B_hat2*W));
    problem2 = createOptimProblem('fmincon','objective',objfun2,'x0',x02,'lb',lb2,...
        'ub',ub2,'options',optimset('Algorithm','SQP','Disp','none'));
    gs2 = GlobalSearch;
    W0_hat2(:,j1) = run(gs2,problem2);
end