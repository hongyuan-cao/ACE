function [t_fdr_hat,t_fval,index] = t_hat(pai1,T_k,gama,sigma_hat)
%% solve critical value
global aaaaaa
f_t_hat = @(t_fdr_hat) 2*(1-pai1)*cdf('Normal',-t_fdr_hat,0,1) - ...
    gama.*mean(bsxfun(@ge,abs(T_k./sqrt(diag(sigma_hat))),t_fdr_hat),1); 

f_t_hatt = f_t_hat(aaaaaa);
[~,index] = find(f_t_hatt<=0);

if isempty(index) == 1
    t_fdr_hat = Inf;
    t_fval = 0;
else
    index = index(1);
    t_fdr_hat = aaaaaa(index);
    t_fval = f_t_hat(t_fdr_hat);
end