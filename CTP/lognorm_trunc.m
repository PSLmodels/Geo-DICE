function [ dt ] = lognorm_trunc( p, low, up, mu, sig  )

%Truncated lognormal Distribution
%   Calculates the lognormal distribution for values between dt_l and dt_u

cdf2 = cdf('lognorm', up , mu, sig);
cdf1 = cdf('lognorm', low , mu, sig);
pc = p * (cdf2 - cdf1) + cdf1;
dt = icdf('lognorm', pc,  mu, sig);

end