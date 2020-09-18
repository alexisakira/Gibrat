close all

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',14);

% get data
N = length(data); % number of data sets

for n=1:N
    X = data{n};
    counts = X(:,3);
    lcutoff = X(:,1);
    %lcutoff(1) = 0;
    cutoff = [lcutoff;Inf];
    
    F = flipud(cumsum(flipud(counts)))/sum(counts); % empirical complementary CDF
    z = norminv(1-F);
    z = z(2:end);
    x = log(lcutoff(2:end));
    K = length(z);
    
    X = [ones(K,1),z];
    betaOLS = X\x; % OLS estimator

    W = diag(counts(2:end))/sum(counts(2:end));
    betaWLS = (X'*W*X)\(X'*W*x); % WLS estimator
    
    [param_lN,logL] = normal_MLbin(counts,log(cutoff(2:end))); % ML estimator
    betaML = param_lN';
    
    xHatOLS = X*betaOLS;
    xHatWLS = X*betaWLS;
    xHatML = X*betaML;
    
    %plot results
    figure
    plot(z,x,'ok'); hold on
    plot(z,xHatOLS,z,xHatWLS); hold off
    xlabel('Quantile of standard normal')
    ylabel('Log cutoff')
    legend('Data','OLS','WLS','Location','Best')
    title(figuretitles{n})
end

