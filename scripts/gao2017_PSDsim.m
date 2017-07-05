%% --- simulation and analysis of resulting data for Figure 1
%% set up simulation
dt = 0.001; %simulation time step
fs = 1/dt; %sampling rate
tk = 0:dt:1; %PSC kernel time vector
t = 0:dt:60*2; %simulation time vector

%spike train parameters
FR_E = 2;
FR_I = 5;
N_E = 8000;
N_I = 2000;

%ampa/gaba PSC kernels
Vr = -65;
Ee = 0;
Ei = -80;
AMPA_tau = [0.1 2]/1000;
GABA_tau = [0.5 10]/1000;
kA = syn_kernel(tk,AMPA_tau);
kG = syn_kernel(tk,GABA_tau);

%oscillation mask
o_filt = ones(fs/2+1,1);
o_filt(9:13) = 1+4*gausswin(5); %alpha
o_filt(18:25) = 1+1.5*gausswin(8); %beta

%save figure directory
figure_folder = '/Users/rgao/Dropbox/Research/Reports/PSPPSD/newdraft/figures/';

%% simulate with different EI ratios
EI_ratio = 2:0.2:6;
num_trs = 5;
for i=1:length(EI_ratio)
    boost = EI_ratio(i)./((N_I*FR_I*sum(kG))/(N_E*FR_E*sum(kA)));
    disp(EI_ratio(i))   
    for tr = 1:num_trs        
        spk_E = pois_spikes(t(end)+tk(end), dt,N_E,FR_E);
        spk_I = pois_spikes(t(end)+tk(end), dt,N_I,FR_I);

        GE(:,i,tr) = conv(spk_E,kA, 'valid');
        GI(:,i,tr) = conv(spk_I,kG, 'valid')*boost;
        LFP_E(:,i,tr) = detrend(GE(:,i,tr),'constant')*(Ee-Vr);
        LFP_I(:,i,tr) = detrend(GI(:,i,tr),'constant')*(Ei-Vr);
        LFP(:,i,tr) = LFP_E(:,i,tr) + LFP_I(:,i,tr); 
    end    
end
%% compute PSDs & fit
o_mask = repmat(o_filt,[1, length(EI_ratio)]); %oscillation mask
for i = 1:size(LFP,3)
    %include oscillation
    %PSD(:,:,i) = mPSD(LFP(:,:,i),fs,fs,fs/2,fs/2).*o_mask;
         
    %no oscillation
    PSD(:,:,i) = mPSD(LFP(:,:,i),fs,fs,fs/2,fs/2);
    
    %loglog linear fit 
%    lif = linfitPSD(PSD(:,:,i),30:50,1); %regular linear fit
%    slope(:,i) = lif(:,2);
    
    rbf = robfitPSD(PSD(:,:,i),30:50,1); %robust fit    
    slope(:,i) = rbf(:,2);
end
%% sliding freq fit
rr = [];
df=1;
win_len = 20;
fit_win = [20:5:160];
fit_win = [fit_win-win_len/2; fit_win+win_len/2]';
for ei = 1:size(PSD,2)
    for tr = 1:size(PSD,3)
        rr(:,ei,tr) = sliding_slope_fit(squeeze(PSD(:,ei,tr)),fit_win,df);
    end    
end

%% --------------------- FIGURES
%% PSC kernel plot
%PSC kernels
tk_p = 0:0.0001:0.2;
kA_p = syn_kernel(tk_p,AMPA_tau);
kG_p = syn_kernel(tk_p,GABA_tau);
figure
plot(tk_p*1e3,kA_p)
hold on
plot(tk_p*1e3,-kG_p)
hold off
line([-5 20], [0 0], 'linewidth', 1, 'color', 'k')
xlim([-5 20])
ylim([-1.1 1.1])
set(gca,'xTick', 0:10:20)
xlabel('Time (ms)')
set(gca,'yTick', 0)
set(gca,'yTickLabel', [])
legend({'AMPA (Exc.)' 'GABA (Inh.)'})
box off
nice_figure(gcf, [figure_folder 'sim_psc'],[3 3])

%% current & LFP time series
figure
subplot('position', [0.05 0.55 0.9 0.40])
plot(t, LFP_E(:,11,1), 'linewidth',1);
hold on
plot(t, LFP_I(:,11,1), 'linewidth',1);
hold off
xlim([0 0.2])
ylabel('Current')
legend({'Excitatory' 'Inhibitory'}, 'Location', 'southeast')
set(gca,'xTick', [])
set(gca,'yTick', 0)
set(gca,'yTickLabel', [])
box off

subplot('position', [0.05 0.1 0.9 0.40])
plot(t,LFP(:,11,1),'k','linewidth',1);
xlabel('Time (s)')
xlim([0 0.2])
set(gca,'xTick', 0:0.1:0.2)
set(gca,'yTick', 0)
set(gca,'yTickLabel', [])
ylabel('Voltage')
legend({'LFP'}, 'Location', 'southeast')
box off
%nice_figure(gcf, [figure_folder 'sim_time'],[5 4])

%% current & LFP PSDs
figure
PE = mPSD(LFP_E(:,11,1),fs,fs,fs/2,fs/2);
PI = mPSD(LFP_I(:,11,1),fs,fs,fs/2,fs/2);
loglog(0:500,PSD(:,11,1), 'k','linewidth',1);
hold on
loglog(0:500,PE,'linewidth',1);
loglog(0:500,PI,'linewidth',1);
hold off
xlim([0 200])
ylim([3 1e3])
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'yTickLabel', [])
legend({'LFP' 'Excitatory' 'Inhibitory'}, 'Location', 'southwest')
box off
%nice_figure(gcf, [figure_folder 'sim_PSD'],[4 4])
%% PSD with different EI ratios
figure
loglog(0:500, PSD(:,1,1)./sum(PSD(:,1,1)), 'k', 'linewidth',1)
hold on
loglog(0:500, PSD(:,end,1)./sum(PSD(:,end,1)), '-.', 'color', [1 1 1]*0.5,'linewidth',1)
xlim([0 200])
hold off
xlim([0 200])
ylim([10e-4 3e-2])
xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'yTickLabel', [])
legend({'g_E:g_I = 1:2','g_E:g_I = 1:6'}, 'location','southwest')
box off
%nice_figure(gcf, [figure_folder 'sim_PSD_ei'],[4 4])
%% 30-50Hz fit vs. EI ratio
figure
x_ax = 20:-1:0;
plot_filled(x_ax, mean(slope,2), std(slope')', 'k')
hold on;
scatter(x_ax, mean(slope,2), 20, 'ok', 'filled')
hold off
set(gca,'xtick', [0:10:20])
set(gca,'xticklabel', {'1:6' '1:4' '1:2'})
set(gca,'ytick', [-1:0.5:0.5])
xlabel('g_E : g_I ratio')
ylabel('Slope (30-50Hz)')
box off
%nice_figure(gcf, [figure_folder 'sim_slope_ei'],[4 4])

%% sliding freq fit
corrEI = zeros(length(fit_win),size(PSD,3));
for ff = 1:size(rr,1)
    corrEI(ff,:) = corr(1./EI_ratio',squeeze(rr(ff,:,:)), 'type','spearman');
end
figure
plot_filled(fit_win(:,1)+win_len/2,mean(corrEI,2),std(corrEI')',[1 1 1]*0)
hold on
scatter(fit_win(:,1)+win_len/2,mean(corrEI,2),20, 'ok', 'filled')
hold off
xlabel('Fit Center Frequency (20Hz window)')
ylabel('Spearman Correlation')
set(gca,'yTick', -1:0.5:1)
box off
%nice_figure(gcf, [figure_folder 'sim_freqslidecorr'],[4 4])
