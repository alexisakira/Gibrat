function Out = normal_MLbinObj(params,counts,cutoff)
% log likelihood of normal-Laplace from binned data

if cutoff(1) == -Inf   %THIS PRODUCES THE SAME RESULTS AS IN THE OLD CODE SO LONG AS THE FIRST CUTOFF IS 0 (log(0)=-INF)
    cutoff=cutoff(2:end);
    
    mu = params(1);
    sigma =params(2);
    
    N = sum(counts); % number of observations
    K = length(cutoff);
    
    P = 0*counts;
    P(1) = normcdf(cutoff(1),mu,sigma);
    
    for k=2:K-1
        P(k) = normcdf(cutoff(k),mu,sigma) - normcdf(cutoff(k-1),mu,sigma);
    end
    
    P(K) = 1-normcdf(cutoff(K-1),mu,sigma);
    P = max(P,eps); % avoid very small numbers
    Out = sum(counts.*log(P))/N;

elseif cutoff(1) ~= -Inf    %THIS ACCOMMADATES FIRST CUTOFF THAT IS NOT ZERO
    mu = params(1);
    sigma =params(2);
    
    N = sum(counts); % number of observations
    K = length(cutoff);
    
    P = 0*counts;
    
    for k=1:K-2 %-2, because the first element is not removed
        P(k) = normcdf(cutoff(k+1),mu,sigma) - normcdf(cutoff(k),mu,sigma);
    end
    
    P(K-1) = 1-normcdf(cutoff(K-1),mu,sigma);
    P = max(P,eps); % avoid very small numbers
    Out = sum(counts.*log(P))/N;
end

