function rbf = sliding_slope_fit(P, cfs, df, lin)
%rbf = sliding_slope_fit(P, cfs, df, lin)
% performs linear (regular or robust) fit over a PSD (P) over all the
% windows, as defined by cfs (window edges)
%
% P: PSD
% cfs: freq window edges
% df: freq resolution of PSD
% (optional) lin: lin = 1 - regular fit, lin = 0 - robust fit (default)

if nargin==3
    lin = 0;
end

rbf = zeros(1,length(cfs));
for i=1:length(cfs)
    frange = cfs(i,1):df:cfs(i,2);
    if lin
        %linear fit
        R = linfitPSD(P,frange,df);
    else
        %robust fit
        R = robfitPSD(P,frange,df);
    end
    rbf(i) = R(2);    
end