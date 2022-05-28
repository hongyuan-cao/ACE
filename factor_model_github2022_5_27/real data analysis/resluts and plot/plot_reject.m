clear;clc;close all
%% FDR vs Number of Rejection after removing the poorly-expressed markers
gamma = [0.001 0.005 0.01 0.02 0.03 0.04 0.05];
PP_b = [185  241  270  297  322  338 349];

Fan_b = [57  99  113  141  162  179  183];

CK_b = [52  89  112  147  178  188  204];

ST_b = [41  80  109  138  163  180  194];

% subplot(3,1,1)
b = plot(gamma,PP_b,'-*' ...
    ,gamma,Fan_b,'--o',gamma,CK_b,'k:^',gamma,ST_b,'-.+');
for jj =1:4
    b(jj).LineWidth = 1;
end
title('\bf Benchmark')

xlabel('FDR control level \gamma')
ylabel('Number of rejection')
legend('PP','Fan','CK','Storey');
xlabel('FDR control level \gamma')
ylabel('Number of rejection')

