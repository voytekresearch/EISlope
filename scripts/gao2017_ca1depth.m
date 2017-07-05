%% analysis pipeline and figure plotting for CA1 depth-slope analysis (Figure 2)
%% -----------------------------------------------
% process all CA1 LFP data for PSD & phase-split PSD
files = dir('lfp*');
sv_root =  './results/CA1_depth/';
%parameters
warning('off'); %turn off robustfit warning
for file = 1:length(files)
    load(files(file).name)
    disp(files(file).name)    
    tic
    fs = srate;
    winLen = fs*2;
    stepLen = fs/2;
    endfreq = 400;
    osc_band = [5 10];
    dF = 0.5;    
    numchan = size(data,1);

    for chan = 1:numchan
        disp(chan)        
        % calculate PSD & HFA from STFT
        [F, Ft, Fa] = stft([], data(chan,:)', fs, winLen, stepLen, endfreq);
        PSD(:,chan) = median(abs(F).^2,3);
        %HFA calculated from 140 to 230Hz, avoiding 180Hz power harmonics
        HFA(:,chan) = squeeze(sum(log10(abs(F([140:175 185:230]*2,:,:)).^2),1))'; %140-230HzHz

        %phase-split PSD
        %raw data
        disp('raw')
        [Pcat(:,:,chan),Pcatv(:,:,chan), Ps{chan}]=phase_split_PSD(data(chan,:),srate,[5 10],[-pi/2, pi/2]);                
        %stopband around beta peak
        disp('beta filtered')
        nobeta = data(chan,:)-eegfilt(data(chan,:),srate,16,21);
        [Pcat_nb(:,:,chan),Pcatv_nb(:,:,chan),Ps_nb{chan}]=phase_split_PSD(nobeta,srate,[5 10],[-pi/2, pi/2]);        
        
        %single trial fitting, PSDs not saved due to memory
        disp('fitting...')
        for ph = 1:2
            tmp =robfitPSD(Ps{chan}{ph}',30:50,1);
            slpfit{chan}{ph} = tmp(:,2);
            tmp =robfitPSD(Ps_nb{chan}{ph}',30:50,1);
            slpfit_nb{chan}{ph} = tmp(:,2);
        end
        
        disp('computing PSDs...')
        PSD(:,chan) = mPSD(data(chan,:)',srate,srate, srate/2, 500);
        PSDnb(:,chan) = mPSD(nobeta',srate,srate, srate/2, 500);
    end
    
    %save_file = [sv_root files(file).name];    
    %save(save_file, 'PSD', 'HFA');
    %clear data PSD HFA
    
    %save(save_file, 'PSD', 'PSDnb', 'Pcat', 'Pcat_nb', 'Pcatv', 'Pcatv_nb', 'Ps', 'Ps_nb', 'slpfit', 'slpfit_nb');
    %clear data PSD PSDnb Pcat Pcat_nb Pcatv Pcatv_nb Ps Ps_nb slpfit slpfit_nb
    toc    
end
%% analysis for depth-slope change
dF = 0.5;
fit_range = [30:dF:50];
files = dir('*lfp*');
all_corr = {};
all_pv = {};
for file = 1:length(files)
    load(files(file).name)
    disp(files(file).name)        
    
    numchan = size(PSD,2);
    %fitting slope and align shanks by depth
    SWR = mean(HFA,1)';    
    temp = robfitPSD(PSD, fit_range, dF);
    EItr = temp(:,2);    
    
    %shank configuration in hc2 doc on CRCNS database
    %pad missing data channels with NaN to make 8-by-{4,7,8} grid
    switch numchan
        case 31
            SWR = [SWR; nan];
            EItr = [EItr; nan];
            num_shank = 4;
        case 64
            num_shank = 8;
        case 55
            SWR = [SWR(1:15); nan; SWR(16:22); nan; SWR(23:29); nan; SWR(30:51); nan; nan; SWR(52:end); nan*zeros(4,1)];
            EItr = [EItr(1:15); nan; EItr(16:22); nan; EItr(23:29); nan; EItr(30:51); nan; nan; EItr(52:end); nan*zeros(4,1)];
            num_shank = 8;
        case 56
            SWR = [SWR(1:23); nan; SWR(24:30); nan; SWR(31:52); nan; nan; SWR(53:end); nan*zeros(4,1)];
            EItr = [EItr(1:23); nan; EItr(24:30); nan; EItr(31:52); nan; nan; EItr(53:end); nan*zeros(4,1)];
            num_shank = 8;
    end
    
    EI_exp = nan(15,num_shank);
    EIz = nan(15,num_shank);
    HFA_exp = nan(15,num_shank);
    HFA_grid = reshape(SWR,8,[]);
    [~, offset] = max(HFA_grid);           
    EI_grid = reshape(EItr,8,num_shank);           
    for shank=1:num_shank
        %absolute slope
        EI_exp((8-offset(shank))+(1:length(EI_grid(:,shank))),shank) = EI_grid(:,shank);
        %normalized slope
        EIz((8-offset(shank))+(1:length(EI_grid(:,shank))),shank) = (EI_grid(:,shank)-nanmean(EI_grid(:,shank)))/nanstd(EI_grid(:,shank));
        %high frequency activity
        HFA_exp((8-offset(shank))+(1:length(EI_grid(:,shank))),shank) = HFA_grid(:,shank);
    end    
    Ir{file} = EI_exp;
    Iz{file} = EIz;
    Hr{file} = HFA_exp;
end

%% slope-inhibition correlation & linear model (slope averaged across 124 shank)
%cell synapse density from megias et al 2001
% densities (E,I): 
%   rad: (0.03, 1.7)/um 
%   pyr: (0, 20)/100 sq um --> sqrt(20/100)/um 
%   ori: (0.64, 0.61)/um & (3.4, 0.1)/um (proximal/distal)
segs = [5,3,4,3]; %rad, pyr, ori.p, ori,dis
E_val = [0.03 0.00001     0.64 3.4];
I_val = [1.7 sqrt(20/100) 0.61 0.1];
E_density = [];
I_density = [];
for i=1:length(segs)
    E_density = [E_density E_val(i)*ones(1,segs(i))];
    I_density = [I_density I_val(i)*ones(1,segs(i))];
end
mask = gausswin(5);
mask = mask./sum(mask);
E = conv(E_density, mask, 'same')';
I = conv(I_density, mask, 'same')';
EI = log10(E./I);

%grand average correlation
slope = nanmean(cell2mat(Ir),2);
synv = ([E I EI]);
clc
[rh_p pv_p] = corr(slope,synv);
[rh_s pv_s] = corr(slope,synv, 'type','spearman');
disp('Grand Average Correlations (E, I, EI)')
disp('---Pearson---')
disp(sprintf('Rho: %f, P-val: %f \n',[rh_p' pv_p']'));
disp('---Spearman---')
disp(sprintf('Rho: %f, P-val: %f \n',[rh_s' pv_s']'));

% individual shank correlation (124 shanks)
slope = reshape(cell2mat(Ir),1,[]);
Ev = repmat(E, length(slope)/length(E),1);
Iv = repmat(I, length(slope)/length(E),1);
EIv = log10(Ev./Iv);
inds = ~isnan(slope); %every shank is shifted, so has NaN at some depths
slope = slope(inds);
synv = [Ev(inds) Iv(inds) EIv(inds)];
[rh_p_all, pv_p_all] = corr(slope',synv, 'type', 'pearson');
[rh_s_all, pv_s_all] = corr(slope',synv, 'type', 'spearman');
disp('Individual shank correlations (E, I, EI)')
disp('---Pearson---')
disp(sprintf('Rho: %f, P-val: %f \n',[rh_p_all' pv_p_all']'));
disp('---Spearman---')
disp(sprintf('Rho: %f, P-val: %f \n',[rh_s_all' pv_s_all']'));

%% multivariate linear model
disp('Multivariate Linear Models')
comp_label = {'E', 'I', 'E:I'};
mod_inds = {1,2,3,1:2,[1 3],2:3,1:3};
model_p = ones(length(mod_inds),1);
coeff = NaN*ones(length(mod_inds),4);
coeff_p = NaN*ones(length(mod_inds),4);
Rs = ones(length(mod_inds),2);
for i=1:length(mod_inds)
    disp(comp_label(mod_inds{i}));
    mdl = LinearModel.fit(synv(:,mod_inds{i}), slope');   
    Rs(i,1) = mdl.Rsquared.Ordinary; %regular R^2
    Rs(i,2) = mdl.Rsquared.Adjusted; %adjusted R^2
    %coefficient estimates
    tmp = mdl.Coefficients.Estimate;
    coeff(i,1) = tmp(1); %constant coefficient
    coeff(i,mod_inds{i}+1) = tmp(2:end); %factor coefficients
    
    %coefficient p-values
    tmp = mdl.Coefficients.pValue;
    coeff_p(i,1) = tmp(1); %constant coefficient
    coeff_p(i,mod_inds{i}+1) = tmp(2:end); %factor coefficients
    
    model_p(i) = coefTest(mdl); %model p-value against null (constant) model
end

%% -------SUPPLEMENTARY ANALYSIS-----------
%% per rat correlation
clc
rh_s=[];
pv_s=[];
rats = {1:11, 12:15, 16:17 18:21};
disp('Per rat correlation (E, I, EI)')
for i=1:4    
    slope = reshape(cell2mat(Ir(rats{i})),1,[]);
    Ev = repmat(E, length(slope)/length(E),1);
    Iv = repmat(I, length(slope)/length(E),1);
    EIv = log10(Ev./Iv);
    inds = ~isnan(slope);
    synv = [Ev(inds) Iv(inds) EIv(inds)];
    [rh_s(i,:), pv_s(i,:)] = corr(slope(inds)',synv, 'type', 'spearman');
    disp(sprintf('Rat %i ----------', i))
    disp(sprintf('Rho: %f, P-val: %f \n',[rh_s(i,:)' pv_s(i,:)']'));
end

%% rat averaged shank correlation
clc
slope=[];
for i=1:4    
    slope(i,:) = nanmedian(cell2mat(Ir(rats{i})),2);    
end
slope = reshape(slope',1,[]);
Ev = repmat(E, length(slope)/length(E),1);
Iv = repmat(I, length(slope)/length(E),1);
EIv = log10(Ev./Iv);
inds = ~isnan(slope);
synv = [Ev(inds) Iv(inds) EIv(inds)];
disp('Rat averaged correlations (E, I, EI)')
disp('---Spearman---')
[rh_s, pv_s] = corr(slope(inds)',synv, 'type', 'spearman');
disp(sprintf('Rho: %f, P-val: %f \n',[rh_s' pv_s']'));



%% FIGURES -----------------------
%save figure directory
figure_folder = '/Users/rgao/Dropbox/Research/Reports/PSPPSD/newdraft/figures/';
%% single-shank PSD
load('~/Documents/Data/CA1_Buzsaki/results_depth/lfp_spiketimes_ec013.527.mat')
for i=1:8
    loglog(0:dF:400,PSD(:,8+i), 'color', [1 1 1]*(8-i)/10, 'linewidth', 1)
    [1 1 1]*(8-i)/10
    hold on
end
hold off
xlim([1 400])
ylim([0.5 4e4])
xlabel('Frequency (Hz)')
set(gca,'YTickLabel',[])
ylabel('Power')
C = linspace(0.7,0,100);
colormap([C; C; C]')
hc = colorbar('location', 'eastoutside');
set(hc,'YTick',[])
set(hc,'YTickLabel',[])
ylabel(hc,'CA1 Depth');
hc.Box = 'off';
box off
nice_figure(gcf, [figure_folder 'CAdepth_PSD'],[5 4])

%% depth plot of slope
figure
colors = get(gca,'colororder');
M = nanmean(cell2mat(Ir),2);
S = nanstd(cell2mat(Ir)')';
line([M+S M-S]', (8-[1:15;1:15])*20, 'color', 'k', 'linewidth',3)
hold on
scatter(M, (7:-1:-7)*20, 100, [1 1 1]*0.5, 'filled');
hold off
xlabel('Slope')
ylabel('Depth (um)')
set(gca,'YTick',-100:100:100)
%nice_figure(gcf, [figure_folder 'CAdepth_slope'],[4 4])

%% synapse density plot
figure
colors = get(gca,'colororder');
barh((7:-1:-7)*20,E, 'facecolor', colors(1,:)*1.1, 'linewidth',1)
hold on
barh((7:-1:-7)*20,-I, 'facecolor', colors(2,:)*1.1, 'linewidth',1);
hold off
xlim([-3.5 3.5])
legend({'Excitatory' 'Inhibitory'})
set(gca,'XTick',-2:2:2)
set(gca,'XTickLabel',[2 0 2])
set(gca,'YTick',-120:40:120)
xlabel('Synapse Density (#/um)')
ylabel('Depth (um)')
box off
%nice_figure(gcf, [figure_folder 'CAdepth_syncount'],[4 4])

%% slope-EI correlation
M = nanmean(cell2mat(Ir),2);
S = nanstd(cell2mat(Ir)')';
figure
line([EI EI]',[M-S M+S]', 'color', 'k')
hold on
scatter(EI,M,80,[1 1 1]*0.5,'filled')
hold off
xlim([-2.5 2])
set(gca,'XTick',-2:2:2)
set(gca,'YTick',-3:1:-1)
xlabel('log_{10}EI ratio')
ylabel('Slope')
title([sprintf('Spearman''s Rho = %.3f, p = %.3f', rh_s_all(3), pv_s_all(3))])
box off
nice_figure(gcf, [figure_folder 'CAdepth_slopeEI'],[4 4])

%% slope-I correlation
M = nanmean(cell2mat(Ir),2);
S = nanstd(cell2mat(Ir)')';
figure
line([I I]',[M-S M+S]', 'color', 'k')
hold on
scatter(I,M,80,[1 1 1]*0.5,'filled')
hold off
xlim([0 2])
set(gca,'XTick',0:0.5:2)
set(gca,'YTick',-3:1:-1)
xlabel('Inhibitory Density (#/um)')
ylabel('Slope')
title([sprintf('Spearman''s Rho = %.3f, p = %.5f', rh_s_all(2), pv_s_all(2))])
box off
nice_figure(gcf, [figure_folder 'CAdepth_slopeI'],[4 4])

