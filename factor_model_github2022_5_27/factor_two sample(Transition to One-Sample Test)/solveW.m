function W0_hat = solveW(Y,deltaa)
% estimate factor model

% Input： 
%   Y: construct new statistics followed factor model by FAME method
%   deltaa: estimation of covariance Y

% output：
%   W0_hat：the estimation of factor

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

%% estimate factor loading matrix B_hat, factor W0_hat
global n1 h_hat

%%
[gama_norm, lambdaa] = svd(deltaa);
lam_sort = diag(lambdaa);

%% get factor number estimation h_hat
if ~(isreal(lam_sort))
    error('特征根是复数')
end
% bizhi = lam_sort(1:(end-1))./lam_sort(2:end);
h_max = 10;
bizhi = lam_sort(1:(h_max-1))./lam_sort(2:h_max);
[~,h_hat] = max(bizhi);
%% estimate factor loading matrix B_hat
B_hat = gama_norm(:,1:h_hat)*diag(sqrt(lam_sort(1:h_hat)));

%% get factor W0_hat
x0 = zeros(h_hat,1); lb = repmat(-3,h_hat,1); ub = repmat(3,h_hat,1);
for jjj = 1:n1
    objfun =@(W) sum(abs(Y(:,jjj) - B_hat*W));
    problem = createOptimProblem('fmincon','objective',objfun,'x0',x0,'lb',lb,...
        'ub',ub,'options',optimset('Algorithm','SQP','Disp','none'));
    gs = GlobalSearch;
    W0_hat(:,jjj) = run(gs,problem);
end