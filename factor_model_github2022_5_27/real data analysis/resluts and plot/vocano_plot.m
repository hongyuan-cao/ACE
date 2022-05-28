clear;clc
global n aaaaaaa aaaaaa
p = 3523; n = 96; gamma = 0.05; aaaaaa = 0.01:0.001:5; aaaaaaa = 0.3:0.001:4;
xe = readtable('E:\RR\SAVE\graduate\pfa\R\Ebench.csv');%% change with your pathway
Y = cell2mat(table2cell(xe));
xv = readtable('E:\RR\SAVE\graduate\pfa\R\Vbench.csv');%% change with your pathway
X = cell2mat(table2cell(xv));

%% Marginal
[~,b,~] = ttest2(X',Y');
pvalues_M = b';
mean_diff_M = mean(X,2) - mean(Y,2);
% T = table(-log10(pvalues));
% writetable(T,'E:\RR\SAVE\Fan\real_data_analysis\log10_pvalue.csv')

%% PP remove factor
Z = X - Y;
Zba = mean(Z,2);
T0_hat = (sqrt(n)* Zba);
p = length(T0_hat);

deltaa = cov(Z');

W0_hat = solveW(Z,deltaa);

[t_fdr_hat, T_k, pai1, sigma_hat] = solvet_hat00(Z,W0_hat,T0_hat,gamma);
T = T_k;
% mean_diff = table(T);
% writetable(mean_diff,'E:\RR\SAVE\Fan\real_data_analysis\mean_diff_PP.csv')

pvalues = 2*cdf('Normal',-abs(T_k./sqrt(diag(sigma_hat))),0,1);
% Pvalue = table(-log10(pvalues));
% writetable(Pvalue,'E:\RR\SAVE\Fan\real_data_analysis\log10_pvalue_PP.csv')

scatter(T/sqrt(96),-log10(pvalues))
hold on
scatter(mean_diff_M,-log10(pvalues_M))