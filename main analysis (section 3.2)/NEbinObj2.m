function [g,dg] = NEbinObj2(m,sigma,a,counts,cutoff)
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

P = (F(2:K+1)-F(1:K))/(1-F(1)); % probability conditional on truncation
g = -sum(counts.*log(P))/N;

temp1 = dF(:,2:K+1)-dF(:,1:K);
temp2 = sum(bsxfun(@times,counts./P,temp1),2)/N;
if cutoff(1) == -Inf
    dg = -temp2;
else
    dg = temp2 - dF(:,1)/(1-F(1)); % last term adjusts for conditioning on truncation
end

end

