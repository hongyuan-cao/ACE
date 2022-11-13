clear;clc;tic
global p n aaaaaa
aaaaaa = 0.01:0.001:10;
% p = 3523; n = 96;
gama = 0.01;
xe = readtable('D:\matlab\bin\SAVE\real_data\code\Ebench.csv');%% change with your pathway
xe = cell2mat(table2cell(xe));
xv = readtable('D:\matlab\bin\SAVE\real_data\code\Vbench.csv');%% change with your pathway
xv = cell2mat(table2cell(xv));

vv = find(mean(xv,2)<6);
ee = find(mean(xe,2)<6);
ve = intersect(vv, ee);
xv(ve,:) = [];
xe(ve,:) = [];

% ve = table(ve);
% writetable(ve,'D:\R\SAVE\real data with 6 threshold\file\index_remove.csv');

Z = xv - xe;
[p, n] = size(Z);
Zba = mean(Z,2);
T0_hat = (sqrt(n)*Zba);

deltaa = cov(Z');

for jj = 1:length(gama)
%%
W0_hat = solveW(Z,deltaa);

%% Step 2, 3, 4
[t_fdr_hat, T_k, pai1, sigma_hat, mu_hat] = solvet_hat00(Z,W0_hat,T0_hat,gama(jj));

%% Step 5 solution FDP_hat
R = sum((abs(T_k./sqrt(diag(sigma_hat)))>=t_fdr_hat));

FDP_hat = 2*(1 - pai1)*cdf('Normal',-t_fdr_hat,0,1)/(R/p);
ppai1(jj) = pai1;
RR(jj) = R;
FFDP_hat(jj) = FDP_hat;
end
clc
fprintf('pai1 = %f\n FDP_hat = %f\n R = %f\n',ppai1,FFDP_hat,RR);
toc

%% record results
% dat = importdata('D:\matlab\bin\SAVE\real_data\code\data.benchmark.csv');
% miRNA = dat.textdata(2:end,2);
% miRNA(ve) = [];
% 
% % index reject and rejected miRNA among 507 miRNA
% Index_reject = find((abs(T_k./sqrt(diag(sigma_hat)))>=t_fdr_hat));
% miRNA_reject = miRNA(Index_reject);
% 
% % adj pvalues
% adj_pvalues = 2*cdf('Normal',-abs(T_k./sqrt(diag(sigma_hat))),0,1);
% adj_pvalues = adj_pvalues(Index_reject);
% 
% % adjust mean difference
% mu_hat = mu_hat(Index_reject);
% 
% % estimated standard deviation
% sd = sqrt(diag(sigma_hat/n));
% sd = sd(Index_reject);
% 
% % statistic
% statistic = mu_hat./sd;
% 
% % up or down regulate
% regulate = repmat('Up',length(Index_reject),1);
% regulate = string(regulate);
% regulate(statistic < 0) = 'Down';
% 
% 
% Index_reject = table(Index_reject); miRNA_reject = table(miRNA_reject);
% adj_pvalues = table(adj_pvalues); mu_hat = table(mu_hat); sd = table(sd);
% statistic = table(statistic); regulate = table(regulate);
% 
% writetable([Index_reject,miRNA_reject,adj_pvalues,mu_hat,sd,statistic,regulate],...
%     'D:\R\SAVE\real data with 6 threshold\file\gene_reject_by_ACE.csv');
