function Out = NE_cdf(x,param)
% CDF of the normal-exponential distribution defined by
% exponential distribution with parameter alpha and
% N(mu,sigma), where param=[mu,sigma,alpha]
mu=param(1);
sigma=param(2);
alpha=param(3);
if (sigma<0)||(alpha<=0);
    Out = NaN*x;
elseif alpha == Inf; % normal distribution
    Out = normcdf(x,mu,sigma);
elseif sigma==0; % shifted exponential distribution
    Out = (1-exp(-alpha*abs(x-mu))).*(x>=mu);
else % normal-exponential distribution
    temp1=exp(alpha^2*sigma^2/2-alpha*(x-mu))...
        .*normcdf((x-mu)/sigma-alpha*sigma);
    Out=normcdf((x-mu)/sigma)-temp1;
end
end

