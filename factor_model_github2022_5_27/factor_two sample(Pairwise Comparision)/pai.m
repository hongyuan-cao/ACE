function pai1 = pai(T_k, sigma_hat)
% function to get the estimation of the proportion of non-nulls

% Input： 
%   T_k: the mean level after removing factor
%   sigma_hat: the estimation of covariance of T_k

% output：
%   pai1：the estimation of the proportion of non-nulls

%   Call description:
%   This method is based on Cao, H., Kosorok, M. (2011). Simultaneous 
%   critical values for t-tests in very high dimensions. Bernoulli, 17, 347–394.

%   VersionV1.0, the code was written in 2022, May 27, revised in 2022,
%   May, 28, author: Peng Wang

%% CK's pai1 estimation
aaaaaaa = 0.3:0.001:4;
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