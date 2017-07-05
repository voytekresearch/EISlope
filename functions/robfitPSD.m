function [fitparam, fiterr] = robfitPSD(PSD,freqs,dF)
%[fitparam, fiterr]= robfitPSD(PSD,freqs,dF)
% performs robust fit on PSD, over specifed frequency regions
%
%   PSD: can be 1,2 or 3D, must start from 0Hz
%       if 1D, will automatically rotate to have long axis as frequency
%       if 2D, needs to be []
%       if 3D, needs to be [X,Y,freq], where X & Y can be trial or time
%   freq: the frequency range to fit, should be at dF resolution
%   dF: is the frequency resolution of PSD [default 1]

if nargin == 2
    dF=1;
end
if length(size(PSD))==2
    %2D matrix, check if really 1D
    if any(size(PSD)==1)
       %really 1D, long axis is freq 
       if size(PSD,2)>size(PSD,1)
           PSD = PSD';
       end
    end
    PSD = permute(PSD,[2 3 1]); %expand into 3D array
end
[LX LY ~] = size(PSD);
fitparam = zeros(LX, LY, 2);
fiterr = zeros(LX, LY);
%fitInds = (round(freqs(1)/dF):round(freqs(end)/dF))+1;
fitInds = round(freqs/dF)+1; % find indices of PSD to fit over, corresponding to freqs
for xx = 1:LX
    for yy = 1:LY
        % perform fit and retrieve error
        [fitparam(xx,yy,:), temp] = robustfit(log10((fitInds-1)*dF),log10(squeeze(PSD(xx,yy,fitInds))));
        fiterr(xx,yy) = temp.ols_s;
    end
end
fitparam = squeeze(fitparam);
fiterr = squeeze(fiterr);