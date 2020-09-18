function [param,se,logL] = NEbinML(counts,cutoff)

N = sum(counts); % number of observations
K = length(counts); % number of categories
if length(cutoff) ~= K+1
    error('length of cutoff must be +1 of counts')
end

[param_lN,logL0] = normal_MLbin(counts,cutoff);  % CHANGE. OLD VERSION: log(cutoff(2:end))
mu0 = param_lN(1);
sigma0 = param_lN(2);

obj = @(theta)(NEbinObj(theta(1),theta(2),theta(3),counts,cutoff));
theta0 = [mu0,sigma0,2];

% do maximum likelihood

options = optimset('TolFun',1e-8,'TolX',1e-8,'Display','off','GradObj','on');
[theta,fval] = fmincon(obj,theta0,[],[],[],[],[-Inf,0,0],[],[],options);
logL = -fval*sum(counts);

m = theta(1);
sigma = theta(2);
a = theta(3);
alpha = a/sigma;
param = [m,sigma,alpha];

if logL < logL0
    param = [param_lN,Inf];
    logL = logL0;
end

% compute standard errors

[F,dF] = NE_help(cutoff,m,sigma,a);

P = F(2:K+1)-F(1:K);
temp = dF(:,2:K+1)-dF(:,1:K);
S1 = temp(1,:);
S2 = temp(2,:) + alpha*temp(3,:);
S3 = sigma*temp(3,:);
S = [S1; S2; S3];

M = zeros(3);
for k=1:K
    s = S(:,k);
    M = M + counts(k)/P(k)^2*(s*s');
end
V = inv(M);
se = sqrt(diag(V))';

end

