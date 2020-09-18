function [param,logL] = normal_MLbin(counts,cutoff)

%K = length(cutoff);
%if length(counts) ~= K
%    error('counts and cutoff must have same length')
%end

if cutoff(1) == -Inf   %THIS PRODUCES THE SAME RESULTS AS IN THE OLD CODE SO LONG AS THE FIRST CUTOFF IS 0 (log(0)=-INF)
    obj=@(theta)(-normal_MLbinObj(theta,counts,cutoff));
    
    % compute mean and standard deviation from binned data
    temp = counts(1:end-1); % delete last observation because it's infinity
    mu0 = sum(temp.*cutoff(2:end-1))/sum(temp);
    sigma0 = sqrt(sum(temp.*cutoff(2:end-1).^2)/sum(temp)-mu0^2);
    
    options = optimset('TolFun',1e-8,'TolX',1e-6,'Display','off');
    [param,fval] = fmincon(obj,[mu0,sigma0],[],[],[],[],[-Inf,0],[],[],options);
    logL = -fval*sum(counts);
    
elseif cutoff(1) ~= -Inf  %THIS ACCOMMADATES FIRST CUTOFF THAT IS NOT ZERO
    obj=@(theta)(-normal_MLbinObj(theta,counts,cutoff));
    
    % compute mean and standard deviation from binned data
    temp = counts; 
    mu0 = sum(temp.*cutoff(1:end-1))/sum(temp);
    sigma0 = sqrt(sum(temp.*cutoff(1:end-1).^2)/sum(temp)-mu0^2);
    
    options = optimset('TolFun',1e-8,'TolX',1e-6,'Display','off');
    [param,fval] = fmincon(obj,[mu0,sigma0],[],[],[],[],[-Inf,0],[],[],options);
    logL = -fval*sum(counts);
end
end

