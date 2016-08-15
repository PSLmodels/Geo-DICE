function [ prc ] = lognorm_trunc_inv( dist, low, up, mu, sig  )

%Truncated Lognorm Distribution
%   Calculates the Lognorm distribution for values between dt_l and dt_u 

cdf2 = cdf('lognorm', up , mu, sig);
cdf1 = cdf('lognorm', low , mu, sig);
pp = cdf('lognorm', dist,  mu, sig);
prc = (pp - cdf1) / (cdf2 - cdf1);

end

