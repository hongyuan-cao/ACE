clear;clc
tic
global p n aaaaaa aaaaaaa Z
aaaaaa = 0.01:0.001:5; aaaaaaa = 0.3:0.001:4;
p = 3523; n = 96; %gama = [0.001,0.005,0.01:0.01:0.04];
gama = 0.05;
xe = readtable('E:\RR\SAVE\graduate\pfa\R\Ebench.csv');%% change with your pathway
xe = cell2mat(table2cell(xe));
xv = readtable('E:\RR\SAVE\graduate\pfa\R\Vbench.csv');%% change with your pathway
xv = cell2mat(table2cell(xv));

% remove poorly expressed markers
ve = readtable('E:\RR\SAVE\factor_model_real_data2022\real_data_analysis2022_5_7\index_remove.csv');
% change with your pathway
ve = cell2mat(table2cell(ve));
xv(ve,:) = [];
xe(ve,:) = [];

Z = xv - xe;
Zba = mean(Z,2);
T0_hat = (sqrt(n)* Zba);
p = length(T0_hat);
deltaa = cov(Z');

for jj = 1:length(gama)
%%
W0_hat = solveW(Z,deltaa);

%% Step 2, 3, 4
[t_fdr_hat, T_k, pai1, sigma_hat] = solvet_hat00(Z,W0_hat,T0_hat,gama(jj));

%% Step 5 solution FDP_hat
R = sum((abs(T_k./sqrt(diag(sigma_hat)))>=t_fdr_hat));
[~,index_max] = sort(abs(T_k./sqrt(diag(sigma_hat))),'descend');
index_max = index_max(1:20);
Index_reject = find((abs(T_k./sqrt(diag(sigma_hat)))>=t_fdr_hat));

FDP_hat = 2*(1-pai1)*cdf('Normal',-t_fdr_hat,0,1)/(R/p);
ppai1(jj) = pai1;
RR(jj) = R;
FFDP_hat(jj) = FDP_hat;
end
clc
fprintf(' pai1 = %f\n FDP_hat = %f\n R = %f\n',ppai1,FFDP_hat,RR);
toc
