function [g,dg] = NEbinObj(m,sigma,a,counts,cutoff)
% log likelihood of normal-Laplace from binned data

if size(counts,1) > size(counts,2)
    counts = counts';
end
if size(cutoff,1) > size(cutoff,2)
    cutoff = cutoff';
end

N = sum(counts); % number of observations
K = length(counts); % number of categories

if length(cutoff) ~= K+1
    error('length of cutoff must be +1 of counts')
end

[F,dF] = NE_help(cutoff,m,sigma,a);

P = F(2:K+1)-F(1:K);
g = -sum(counts.*log(P))/N;

temp = dF(:,2:K+1)-dF(:,1:K);
dg = -sum(bsxfun(@times,counts./P,temp),2)/N;


end

