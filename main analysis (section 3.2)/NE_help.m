function [f,df] = NE_help(x,m,sigma,a)

if size(x,1) > size(x,2)
    x = x'; % convert to row vector
end

z = (x-m)/sigma;
eaz = exp(a^2/2-a*z);
Paz = 1-normcdf(a-z);
paz = normpdf(a-z);

f = normcdf(z) - eaz.*Paz; 
    
df1 = (1/sigma)*(-normpdf(z)-eaz.*(a*Paz-paz)); %df/dm
    
df2 = z.*df1;                                   %df/dsigma
    
df3 = eaz.*(-(a-z).*Paz+paz);                   %df/da
    
%%% Changed the above code (commented out) to accommodate truncation (first cutoff not zero)
%f=zeros(size(x));
%df1=zeros(size(x));
%df2=zeros(size(x));
%df3=zeros(size(x));

%for i=1:size(x)
    %if x(1) == 0      
%        f(i) = normcdf(z(i)) - eaz(i)*Paz(i);
    
%        df1(i) = (1/sigma)*(-normpdf(z(i))-eaz(i)*(a*Paz(i)-paz(i))); %df/dm
    
%        df2(i) = z(i)*df1(i);                                   %df/dsigma
    
%        df3(i) = eaz(i)*(-(a-z(i))*Paz(i)+paz(i));                   %df/da

    %elseif x(1) ~= 0
    %    f(i) = (normcdf(z(i)) - eaz(i)*Paz(i))/(1-(normcdf(z(1)) - eaz(1)*Paz(1)));
    
    %    df1(i) = (1/sigma)*(-normpdf(z(i))-eaz(i)*(a*Paz(i)-paz(i))); %df/dm
    
    %    df2(i) = z(i)*df1(i);                                   %df/dsigma
    
    %    df3(i) = eaz(i)*(-(a-z(i))*Paz(i)+paz(i));                   %df/da
    %end
%end

%%%%%%%

df = [df1; df2; df3];
    
f(x==Inf) = 1;
f(x==-Inf) = 0;
df(:,x==Inf) = 0;
df(:,x==-Inf) = 0;

end

