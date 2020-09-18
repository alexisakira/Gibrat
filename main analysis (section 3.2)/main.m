close all

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',14);

%Specify the directory
cd('/Users/')

%Get Data
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

N = length(data); % number of data sets

%Create storage for estimation results
Param_PlN = zeros(3,N);
Param_lN = zeros(2,N);
LogL = zeros(2,N);
stderr = zeros(3,N);

for n=1:N
    X = data{n};
    counts = X(:,3);
    lcutoff = X(:,1);
    cutoff = [lcutoff;Inf];
    
    tic
    [param_PlN,se,logL] = NEbinML2(counts,log(cutoff));
    logL = logL*sum(counts);
    toc
    Param_PlN(:,n) = param_PlN';
    stderr(:,n) = se';
    LogL(1,n) = logL;
    
    [param_lN,logL] = normal_MLbin(counts,log(cutoff));
    logL = logL*sum(counts);
    Param_lN(:,n) = param_lN';
    LogL(2,n) = logL;
    temp1 = 10^ceil(log10(lcutoff(end)));
    if lcutoff(1)>0
        temp2 = 10^floor(log10(lcutoff(1)));
    else
        temp2 = 10^floor(log10(lcutoff(2)));
    end
    grid = linspace(log(temp2),log(temp1),101); % construct grid for PlN
    
    PlNCDF = NE_cdf(grid,param_PlN); % CDF of PlN evaluated at grid
    lNCDF = normcdf(grid,param_lN(1),param_lN(2)); % lognormal
    cumNcounts = flipud(cumsum(flipud(counts)));
    cumNpred_PlN = sum(counts)*(1-PlNCDF); % predicted number of firms above cutoffs using PlN
    cumNpred_lN = sum(counts)*(1-lNCDF); % lognormal
    
     figure
     loglog(lcutoff,cumNcounts,'ok','markersize', 8); hold on
     loglog(exp(grid),cumNpred_PlN,'-');  
     loglog(exp(grid),cumNpred_lN,'--'); hold off
     xlabel('Size cutoff','FontSize', 22)
     ylabel('Number of larger units','FontSize', 22)
     legend('Data','PlN','lognormal','Location','SW')
     lgd.FontSize = 22;                           
     xlim([temp2,temp1]);
     ylim([10^floor(min(log10(cumNcounts))),10^ceil(log10(cumNcounts(2)))])
     title(figuretitles{n},'FontSize', 22)       
end

chi2stat = 2*(LogL(1,:)-LogL(2,:));
format long g
DeltaL = LogL(1,:)-LogL(2,:);
Pval = 1-chi2cdf(chi2stat,1);  
Alpha = Param_PlN(3,:);
Alphaerr = stderr(3,:);

