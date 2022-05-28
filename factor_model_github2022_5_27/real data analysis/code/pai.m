function pai1 = pai(T_k, sigma_hat)
global aaaaaaa
%% CK's pai1 estimation
x = T_k./sqrt(diag(sigma_hat));
pai =@(c) (mean(bsxfun(@min,abs(x),c),1)./c - 2*(1-exp(-c.^2/2))./(c*sqrt(2*pi)) -...
    2 * cdf('Normal',-c,0,1))./(1 -2*(1-exp(-c.^2/2))./(c*sqrt(2*pi)) - ...
    2 * cdf('Normal',-c,0,1));

paii = pai(aaaaaaa);
pai1 = max(paii);
if pai1>1
    pai1 = 1;
else if pai1<0
        pai1 = 0;
    end
end
end