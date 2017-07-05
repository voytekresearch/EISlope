function [PSDcat, PSDcatv, PSDall] = phase_split_PSD(data, fs, ph_freq, ph_bins)
%[PSDcat, PSDcatv, PSDall] = phase_split_PSD(data, fs, ph_freq, ph_bins)
% computes very short-window PSD, based on the phase of a specified
% frequency. It first computes the Hilbert phase of the low-frequency
% oscillation, then bin data for a given phase range, and computes the PSD
% (or FFT, really) for each small segment of data. The phase-bin aggregate
% PSD is the average of all the FFTs of data that falls within the same
% phase bin (peak /trough)
% 
% data: time series
% fs: sampling frequency of data
% ph_freq: [low high] passband frequency for oscillation
% ph_bin: vector of left bin edges of oscillation phase
%    -pi/2 to pi/2 is peak, i.e. cosine phase


disp('Filtering..')
osc = eegfilt(data,fs,ph_freq(1),ph_freq(2));
disp('Scanning..')

%butterworth filt
%osc = butterpass(data,fs,ph_freq,2);

%get phase
ph = angle(hilbert(osc-mean(osc)));

%initialize all phase labels as N in case there's a dump bin
phase_cat = ones(size(ph))*length(ph_bins);
for i = 1:length(ph_bins)-1
    phase_cat(find(ph>=ph_bins(i) & ph<ph_bins(i+1))) = i;
end

%find indices of segment beginnings
cat_inds = [1 find(diff(phase_cat))+1];
cat_inds = [cat_inds; [cat_inds(2:end)-1 length(data)]];
num_segs = length(cat_inds);
labels = phase_cat(cat_inds(1,:));
lfp_segs=cell(num_segs,1);

% gather lfp segments by phase
for i=1:num_segs
    lfp_segs{i} = data(cat_inds(1,i):cat_inds(2,i));
end

% take psd
Nfft = fs;
categories = unique(labels);
PSDall = cell(length(categories),1);
PSDcat = zeros(length(categories),Nfft/2+1);
PSDcatv = zeros(length(categories),Nfft/2+1);
lfp_split = cell(length(categories),1);
for i=1:length(categories)
    lfp_split{i} = lfp_segs(labels==categories(i));
    PSDcur = ones(length(lfp_split{i}), Nfft/2+1);
    for j=1:length(lfp_split{i})
        L = length(lfp_split{i}{j});
        if L>1
            curP = abs(fft([((lfp_split{i}{j}-mean(lfp_split{i}{j})).*hamming(L)') zeros(1,Nfft-L)])).^2;
            PSDcur(j,:) = curP(1:Nfft/2+1);
        end
    end    
    PSDall{i} = PSDcur;
    PSDcat(i,:) = median(PSDall{i});
    PSDcatv(i,:) = std(PSDall{i});
end