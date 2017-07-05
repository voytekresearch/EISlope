function [fitparam, fiterr] = linfitPSD(PSD,freqs,dF)
%[fitparam, fiterr] = linfitPSD(PSD,freqs,dF)
% performs simple linear fit on PSD, over specifed frequency regions
%
%   PSD: power spectrum, can be 1,2 or 3D, must start from 0Hz
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
    else
        %2D, long axis is time
        if size(PSD,1)>size(PSD,2)
            %PSD = PSD';
        end
    end
    PSD = permute(PSD,[2 3 1]); %expand into 3D array
end
[LX LY ~] = size(PSD);
fitparam = zeros(LX, LY, 2);
fiterr = zeros(LX, LY);
%fitInds = (round(freqs(1)/dF):round(freqs(end)/dF))+1;
fitInds = round(freqs/dF)+1;
for xx = 1:LX    
    for yy = 1:LY        
        [temp_param, temp_err] = polyfit(log10((fitInds-1)*dF)',log10(squeeze(PSD(xx,yy,fitInds))),1);        
        fitparam(xx,yy,:) = temp_param(end:-1:1); %fit parameters [offset, slope]
        fiterr(xx,yy) = temp_err.normr; %error
    end
end
fitparam = squeeze(fitparam);
fiterr = squeeze(fiterr);