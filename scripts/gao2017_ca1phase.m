%% analysis pipeline and figure plotting for CA1 phase-slope analysis (Figure 3)
%%    necessary variables are computed from cell 1 of gao2017_ca1depth.m

%% theta-phase single cycle fits (30-50Hz fit)
no_beta = 1;
dF=1;
fit_range = 30:dF:50;
files = dir('phasecat*');
for file = 1:length(files)            
    load(files(file).name)
    disp(files(file).name)
    
    if no_beta
        slpfit = slpfit_nb;
        Pcat = Pcat_nb;
        PSD = PSDnb;
    end
    
    % ttest and median from single cycle fits
    for chan = 1:length(slpfit)
       [h(chan), pv(chan)] = ttest2(slpfit{chan}{1},slpfit{chan}{2});
       mslope(:,chan) = cellfun(@median,slpfit{chan});
    end
    % fit over averaged PSD
    for i=1:2
        temp = robfitPSD(squeeze(Pcat(i,:,:)),fit_range,dF);
        slope(i,:) = temp(:,2);
    end
    H{file} = h;
    PV{file} = pv;
    MSLOPE{file} = mslope;
    SLOPE{file} = slope;  
    %high frequency activity (high gamma)
    HFA{file} = squeeze(mean(log10(PSD(100:300,:)),1));
    HFA_ph{file} = squeeze(mean(log10(Pcat(:,100:300,:)),2))';
    LN{file} = PSD(61,:)-PSD(55,:); %excessive line noise
    
    clear h pv mslope slope
end    
%% correlation and t-test analysis for single theta cycle PSD slope
clc
S = cell2mat(SLOPE); %slope fit to average PSD
MS = cell2mat(MSLOPE); %median slope of single trial (cycle) PSD
dG = cell2mat(HFA_ph')';
diag(corr(S',MS'));
disp(sprintf('Slope diff: %i out of %i channels significant.', sum(cell2mat(H)),length(cell2mat(H))))
[h pv] = ttest(S(1,:)',S(2,:)');
disp(sprintf('Mean Aggregate Fit ---\nTrough: %.3f, Peak: %.3f, p=%f', mean(S(2,:)),mean(S(1,:)), pv))
disp(sprintf('Mean Median of Individual Fit ---\nTrough: %.3f, Peak: %.3f, p=%f', mean(MS(2,:)),mean(MS(1,:)), pv))
[h pv] = ttest(dG(1,:),dG(2,:));
disp('-----')
disp(sprintf('HFA diff: %.4f, p=%.4f', mean(diff(dG)), pv))


%% FIGURES
figure_folder = '/Users/rgao/Dropbox/Research/Reports/PSPPSD/newdraft/figures/';
%% time series schematic
load('~/Documents/Data/CA1_Buzsaki/Raw/lfp_spiketimes_ec013.527.mat')
data_win = 5500+(1:srate);
t= (0:srate-1)/srate;

%get phase bin categories (peak or valley)
ph = angle(hilbert(eegfilt(data(1,:),srate,5,10)));
ph = ph(data_win);
ph_bins = [-pi/2, pi/2];
phase_cat = ones(size(ph))*length(ph_bins);
for i = 1:length(ph_bins)-1
    phase_cat(find(ph>=ph_bins(i) & ph<ph_bins(i+1))) = i;
end
indtr = find(phase_cat==1);
indpk = find(phase_cat==2);

peaks = phase_cat; peaks(peaks==2)=nan;
troughs = phase_cat; troughs(troughs==1)=nan;

figure
colors = get(gca,'colororder');
plot(t,data(1,data_win), 'k', 'linewidth', 1);
hold on
plot(t,500*ph, '.', 'color', [1 1 1]*0.6);
plot(t,troughs*-1000,'color', colors(1,:))
plot(t,peaks*2000,'color', colors(2,:))
hold off
ylim([-2500 2500])
legend({'LFP' 'Phase' 'Troughs' 'Peaks'})
set(gca,'YTick',[])
set(gca,'XTick',0:0.5:1)
xlabel('Time (s)')
box off
nice_figure(gcf, [figure_folder 'CAphase_TS'],[6 3])

%% PSDs
load('~/Documents/Data/CA1_Buzsaki/results_phase/nb/phasecat_nb_lfp_spiketimes_ec013.527.mat');
if no_beta
    Pcat = Pcat_nb;
    PSD = PSDnb;
end
figure
loglog(0:625, squeeze(Pcat(2,:,1)))
hold on
loglog(0:625, squeeze(Pcat(1,:,1)))
loglog(0:500, PSD(:,1)*4e4, 'k')
hold off
xlim([2 200])
xlabel('Frequency')
ylabel({'Power'})
set(gca,'YTick',[])
legend({'Trough' 'Peak', 'All (raw)'})
box off
nice_figure(gcf, [figure_folder 'CAphase_PSD'],[4 4])

%% slope change histograms
slope_chan = S;
bins = -5:0.1:0;
figure
hold on
stairs(bins,hist(slope_chan(2,:),bins), 'linewidth',2)
stairs(bins,hist(slope_chan(1,:),bins),'linewidth',2)
hold off
legend({'Trough (Exc)' 'Peak (Inh)'}, 'location', 'northwest')
xlabel('Slope')
ylabel('Count')
set(gca,'ytick',[])
box off
nice_figure(gcf, [figure_folder 'CAphase_phslope'],[4 4])

figure
bins = -0.5:0.1:1.5;
stairs(bins,hist(diff(slope_chan), bins), 'linewidth', 2, 'color', 'k')
xlabel('Slope Change (Trough - Peak)')
ylabel('Count')
set(gca,'ytick',[])
%box off
nice_figure(gcf, [figure_folder 'CAphase_slopediff'],[4 4])

%% gamma change histograms
hfa_chan = dG;
bins = 4.5:0.05:6;
figure
hold on
stairs(bins,hist(hfa_chan(2,:),bins), 'linewidth',2)
stairs(bins,hist(hfa_chan(1,:),bins),'linewidth',2)
hold off
xlabel('High Freq. Power (100-300Hz)')
ylabel('Count')
legend({'Trough (Exc)' 'Peak (Inh)'}, 'location', 'northwest')
set(gca,'ytick',[])
box off
nice_figure(gcf, [figure_folder 'CAphase_phhfa'],[4 4])

figure
bins = -0.05:0.01:0.2;
stairs(bins,hist(diff(hfa_chan), bins), 'linewidth', 2, 'color', 'k')
xlabel('HFA Change (Trough - Peak)')
ylabel('Count')
set(gca,'ytick',[])
box off
nice_figure(gcf, [figure_folder 'CAphase_hfadiff'],[4 4])

%% connected lines (changeogram? idk)
slope_chan = S;
figure
colors = get(gca,'colororder');
plot([2 1], slope_chan, 'color', [0 0 0 0.1], 'linewidth', 1)
hold on
scatter(1*ones(1,length(slope_chan)), slope_chan(2,:), 20, colors(1,:), 'filled')
scatter(2*ones(1,length(slope_chan)), slope_chan(1,:), 20, colors(2,:), 'filled')

xlim([0.5 2.5])
ylabel('Slope')
set(gca,'xtick',[1 2])
set(gca,'ytick',-4:2:2)
set(gca,'xticklabel', {'Trough (Exc)' 'Peak (Inh)'})
hold off
box off
nice_figure(gcf, [figure_folder 'CAphase_slpline'],[3 4])