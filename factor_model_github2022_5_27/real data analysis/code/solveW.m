function W0_hat = solveW(Z,deltaa)
%% estimate factor loading matrix B_hat, factor W0_hat
global n h_hat
[gama_norm, lambdaa] = svd(deltaa); % deltaa = UU*VV*UU'
lam_sort = diag(lambdaa);

%% get factor number estimation h_hat
% bizhi = lam_sort(1:(end-1))./lam_sort(2:end);
h_max = 10;
bizhi = lam_sort(1:(h_max-1))./lam_sort(2:h_max);
[~,h_hat] = max(bizhi);
%% estimate factor loading matrix B_hat
B_hat = gama_norm(:,1:h_hat)*diag(sqrt(lam_sort(1:h_hat)));

%% get factor W0_hat，
W0_hat = zeros(h_hat,n);
if h_hat == 1
    for jjj = 1:n
        ab = quantreg(B_hat, Z(:,jjj), 0.5, 1);
        W0_hat(:,jjj) = ab(1);
    end
else
    for jjj = 1:n
        W0_hat(:,jjj) = quantreg(B_hat, Z(:,jjj), 0.5);
    end    
end
