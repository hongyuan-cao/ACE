function W0_hat = solveW(Z,deltaa)
%% estimate factor loading matrix B_hat, factor W0_hat
global n h_hat

%% （2）求bigsigma的特征值特征向量，正交化，排序
[gama_norm, lambdaa] = svd(deltaa); %对于对称矩阵有 deltaa = UU*VV*UU'
lam_sort = diag(lambdaa);

%% get factor number estimation h_hat
if ~(isreal(lam_sort))
    error('特征根是复数')
end
bizhi = lam_sort(1:10)./lam_sort(2:11);
[~,h_hat] = max(bizhi);
%% estimate factor loading matrix B_hat
B_hat = gama_norm(:,1:h_hat)*diag(sqrt(lam_sort(1:h_hat)));

%% get factor W0_hat，
x0 = zeros(h_hat,1); lb = repmat(-3,h_hat,1); ub = repmat(3,h_hat,1);
W0_hat = zeros(h_hat,n);
for jjj = 1:n
    objfun =@(W) sum(abs(Z(:,jjj) - B_hat*W));
    problem = createOptimProblem('fmincon','objective',objfun,'x0',x0,'lb',lb,...
        'ub',ub,'options',optimset('Algorithm','SQP','Disp','none'));
    gs = GlobalSearch;
    W0_hat(:,jjj) = run(gs,problem);
end