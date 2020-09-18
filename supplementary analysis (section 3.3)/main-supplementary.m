close all

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',14);

%Specify the directory
cd('/Users/')

filenames = {'068 establishments.xls','071 income (England).xls','072 income (England).xls',...
    '106a income (Saxony).xls','106b income (Saxony).xls','123 income (Prussia).xls',...
    '124 income (Oldenburg).xls','128 income (France).xls','133 income (Japan).xls',...
    '139a wealth (Basel).xls','139b wealth (Basel).xls','147 rental income.xls',...
    '151a estates (France).xls','151b estates (England).xls','151c estates (Italy).xls',...
    '152 inheritance.xls','165 rent.xls','166 property.xls','167 bank deposits.xls',...
    '168 real estate.xls','185 corporate profit.xls','251 communities.xls',...
    '252 cities.xls','253 children.xls'};

figuretitles = {'p.68 Establishments (France, 1921)','p.71 Income (England, 1893)','p.72 Income (England, 1911)',...
    'p.106 Income (Saxony, 1879)','p.106 Income (Saxony, 1890)','p.123 Income (Prussia, 1911)',...
    'p.124 Income (Oldenburg, 1890)','p.128 Income (France, 1906)','p.133 Income (Japan, 1904)',...
    'p.139 Wealth (Basel, 1454)','p.139 Wealth (Basel, 1887)','p.147 Rental Income (France, 1927)',...
    'p.151 Estates (France, 1926)','p.151 Estates (England, 1913)','p.151 Estates (Italy, 1910)',...
    'p.152 Inheritance (France, 1911)','p.165 Rent (Paris, 1878)','p.166 Property Values (France, 1899)',...
    'p.167 Bank Deposits (France, 1896)','p.168 Real Estate (France, 1914)','p.185 Corporate Profit (France, 1927)',...
    'p.251 Communities (France, 1926)','p.252 Cities (Europe, 1927)','p.253 Children (France, 1926)'};

n = length(filenames);

data = cell(n,1);
sums = zeros(n,1);
for idx = 1:n 
    data{idx} = xlsread(filenames{idx}); 
    sums(idx) = sum(data{idx}(:));
end

%getData % import data sets
N = length(data);

alpha_PL = zeros(1,N);
bmin_PL = zeros(1,N);
probbmin_PL = zeros(1,N);
LogL = zeros(1,N);
stder_alpha_PL = zeros(1,N);
stder_bmin_PL = zeros(1,N);
pval_PL = zeros(1,N);
blim_store= zeros(1,N);

mu_LN = zeros(1,N);
sigma_LN = zeros(1,N);

lambda_E = zeros(1,N);

lambda_SE = zeros(1,N);
beta_SE = zeros(1,N);

LR_PL_LN = zeros(1,N);
p_PL_LN = zeros(1,N);

LR_PL_E = zeros(1,N);
p_PL_E = zeros(1,N);

LR_PL_SE = zeros(1,N);
p_PL_SE = zeros(1,N);

tic %Takes 8307.179699 seconds (2.5 hours) to complete the loop
for n=1:N
    X = data{n};
    
    if n==15      
        X(1,1)=1;   % set the first cutoff to 1 (instead of 0) for data set 15 (Rental Income - France)
    end
    
    counts = round(X(:,3));
    lcutoff = X(:,1);
    cutoff = [lcutoff;Inf];
    cCDF = flipud(cumsum(flipud(counts))/sum(counts)); % complementary CDF
    bLim = max(cutoff(cCDF >= 0.1)); % restrict search for bmin to contain 10% of tail observations
    blim_store(1,n)=bLim;
    
    %Fitting a binned power-law distribution
    [alpha,bmin,L] = bplfit(counts,cutoff,'limit',bLim);
    
    if ((sum(X(X(:,1)>=bmin,3))/sum(X(:,3)))*100)==100  % if bmin includes 100% of the sample, set bmin = second smallest cutoff
        bLim=X(2,1);
        [alpha,bmin,L] = bplfit(counts,cutoff,'bmin',bLim);
    end
            
    alpha_PL(1,n) = alpha-1; % subtract 1 because definition in Virkar & Clauset is different
    bmin_PL(1,n) = bmin;
    probbmin_PL(1,n) = (sum(X(X(:,1)>=bmin,3))/sum(X(:,3)))*100;
    LogL(1,n) = L;
    
    %Visualizing the plotting function
    %bplplot2(counts, cutoff, bmin, alpha)
    %title(figuretitles{n},'FontSize', 22) %changed: 14
    
    %Calculating p-value for the fitted power-law mode (KS statistic)
    %[p, d] = bplpva(counts,cutoff,bmin,'limit',bLim);
    %pval_PL(1,n) = p;

    %Estimating uncertainty in the fitted parameters
    %[alpha, bmin] = bplvar(counts,cutoff,'limit',bLim);
    %stder_alpha_PL(1,n) = alpha;
    %stder_bmin_PL(1,n) = bmin;
    
    %Estimate alternative distirbutions
    
    %1.Estimate lognormal
    [mu_est, sigma_est, L] = blgnormfit(counts, cutoff, [], []);
    mu_LN(1,n) = mu_est;
    sigma_LN(1,n) = sigma_est;
    
    %2.Estimate exponential
    [lambda_est, L] = bexpnfit(counts, cutoff);
    lambda_E(1,n) = lambda_est;
    
    %3.Estimate stretched exponential
    [lambda_est2, beta_est2, L] = bstexpfit(counts, cutoff);
    lambda_SE(1,n) = lambda_est2;
    beta_SE(1,n) = beta_est2;
    
    %Compare power law to alternative distirbutions
    
    %1. Power law (PL) vs Lognormal (LN)
    p1 = getPDF(cutoff, bmin, 'pl', alpha);
    p2 = getPDF(cutoff, bmin, 'lgnorm', mu_est, sigma_est);
    [normR, p] = blrtest(p1, p2, counts, cutoff, bmin, 0);
    LR_PL_LN(1,n) = normR;
    p_PL_LN(1,n) = p;
    
    %2. Power law (PL) vs Exponential (E)
    p1 = getPDF(cutoff, bmin, 'pl', alpha);
    p2 = getPDF(cutoff, bmin, 'expn', lambda_est);
    [normR, p] = blrtest(p1, p2, counts, cutoff, bmin, 0);
    LR_PL_E(1,n) = normR;
    p_PL_E(1,n) = p;
    
    %3. Power law (PL) vs Stretched Exponential (SE)
    p1 = getPDF(cutoff, bmin, 'pl', alpha);
    p2 = getPDF(cutoff, bmin, 'stexp', lambda_est2, beta_est2);
    [normR, p] = blrtest(p1, p2, counts, cutoff, bmin, 0);
    LR_PL_SE(1,n) = normR;
    p_PL_SE(1,n) = p;
end
toc
