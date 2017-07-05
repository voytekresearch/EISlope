%% analysis pipeline and figure plotting for monkey ECoG analysis (Figure 4)
%% batch process all sub-session, get PSD and time-resolved slopes
fs = 1000;
df = 1;
frange = [40:df:48 52:df:60]; %skip 50 & 60Hz powerline (Japan)
hfa_fr = [125:df:145 155:df:175];
PSD = {};
SLP = {};
OFF = {};
ERR = {};
HFA = {};
CONDT = {};
CONDIND = {};
FT = {};
CONDLABEL = {};
for sesh = 1:3
    tic
    disp(sprintf('Session%i/',sesh))
    cd(sprintf('Session%i/',sesh))
    load Condition.mat
    load('ECoGTime.mat')  
    num_chan = length(dir('ECoG_*'));       
    [F, Ft, Fa] = stft([],ECoGTime',fs,fs,fs/4,400);
    %initialize matrices
    SLP{sesh} = zeros(length(Ft),num_chan);
    OFF{sesh} = zeros(length(Ft),num_chan);
    HFA{sesh} = zeros(length(Ft),num_chan);
    ERR{sesh} = zeros(length(Ft),num_chan);
    for chan = 1:num_chan     
        fname = sprintf('ECoG_ch%i.mat',chan);
        fdata = sprintf('ECoGData_ch%i',chan);
        data = load(fname, fdata);
        disp(fname);

        %compute STFT over file
        [F, Ft, Fa] = stft([],data.(fdata)',fs,fs,fs/4,400);        
        P = squeeze(abs(F).^2);        
        
        %fit slope
        
        %robust fit: ~2hr run time per monkey, TOO SLOW, not much gain
        [tmp err]= robfitPSD(P,frange,1);
        SLP{sesh}(:,chan) = tmp(:,2);
        OFF{sesh}(:,chan) = tmp(:,1);
        ERR{sesh}(:,chan) = err;

        %linear fit: ~30min run time        
%         slp = zeros(length(Ft),2);
%         for i=1:length(slp)
%             [slp(i,:) err] = polyfit(log10(frange),log10(P(frange+1,i))',1);
%         end        
%         SLP{sesh}(:,chan) = slp(:,1);
%         OFF{sesh}(:,chan) = slp(:,2);
%         HFA{sesh}(:,chan) = mean(log10(P(hfa_fr+1,:)),1);
%         ERR{sesh}(:,chan) = err.normr;
        
        %get condition-averaged PSD        
        for conds = 1:length(ConditionTime)-1
            %find STFT times within condition 
            cur_inds = [find(Ft>ConditionTime(conds),1, 'first')...
                find(Ft<ConditionTime(conds+1),1,'last')];
            PSD{sesh}(:,chan,conds) = median(P(:,cur_inds(1):cur_inds(2)),2);
        end
    end
    CONDLABEL{sesh} = ConditionLabel;
    if sesh~=1
        CONDT{sesh} = ConditionTime + FT{sesh-1}(end);
        CONDIND{sesh} = ConditionIndex + CONDIND{sesh-1}(end);
        FT{sesh} = Ft+FT{sesh-1}(end);
    else
        CONDT{sesh} = ConditionTime;
        CONDIND{sesh} = ConditionIndex;
        FT{sesh} = Ft;
    end
    cd ..
    toc
end
save results4060_rbf.mat PSD SLP HFA OFF ERR CONDT CONDIND FT CONDLABEL frange fs

%% fit slope over block averaged PSDs
warning('off'); %turn off robustfit warning
df = 1;
frange = [40:df:48 52:df:60]; %skip 50Hz powerline
off_avg = {};
slp_avg = {};
mse = {};
num_chan = size(PSD{1},2);
for sesh = 1:3
    for cond = 1:size(PSD{sesh},3)       
        %robustfit
        [tmp err] = robfitPSD(PSD{sesh}(:,:,cond),frange,df);
        off_avg{sesh}(:,cond) = tmp(:,1);
        slp_avg{sesh}(:,cond) = tmp(:,2);
        mse{sesh}(:,cond) = err;
        
%         %linear fit
%         for chan = 1:num_chan
%             [tmp err] = polyfit(log10(frange), log10(PSD{sesh}(frange+1,chan,cond))',1);         
%             slp_avg{sesh}(chan,cond) = tmp(1);
%             off_avg{sesh}(chan,cond) = tmp(2);
%             normR{sesh}(chan,cond) = err.normr;
%             mse{sesh}(chan,cond) = sum(((tmp(1)*log10(frange)+tmp(2))-log10(PSD{sesh}(frange+1,chan,cond)')).^2)/length(frange);            
%         end
    end
end
%t-test between awake and anesthesized conditions
% {1}(1): eyes open; 
% {1}(3): eyes closed; 
% {2}(1): injection; 
% {2}(2): anes
S = [slp_avg{1}(:,3),slp_avg{2}(:,2)];
[h, p] = ttest2(S(:,1), S(:,2));
%clc
disp(sprintf('Mean Slopes ---\nAwake: %.3f, Anes: %.3f. p=%.4f',...
    mean(S(:,1)), mean(S(:,2)), p))
%%





%% FIGURES
figure_folder = '/Users/rgao/Dropbox/Research/Reports/PSPPSD/newdraft/figures/';
%% PSD (mean over all 128 channels)
figure
colors = get(gca,'colororder');
chan = 110;
loglog(0:400, squeeze(mean(PSD{1}(:,:,1),2)), 'k')
hold on
loglog(0:400, squeeze(mean(PSD{2}(:,:,2),2)), 'color',colors(2,:))
hold off
xlim([2 100])
set(gca,'ytick', [])
ylabel('Power')
xlabel('Frequency (Hz)')
legend({'Awake (eyes closed)' 'Anesthesia'})
box off
nice_figure(gcf, [figure_folder 'tycho_PSD'],[4 4])

%% Slope comparison between conditions - connected lines (changeogram? idk)
slope_chan = S;
figure
colors = get(gca,'colororder');
plot([1 2], S, 'color', [0 0 0 0.1], 'linewidth', 1)
hold on
scatter(ones(1,length(S)), S(:,1), 20, 'k', 'filled')
scatter(2*ones(1,length(S)), S(:,2), 20, colors(2,:), 'filled')
hold off
hold off
ylabel('Slope')
xlim([0.5 2.5])
set(gca,'xtick', [1 2])
set(gca,'ytick', [-5:1:-2])
set(gca,'xticklabel', {'Awake' 'Anes'})
box off
nice_figure(gcf, [figure_folder 'tycho_slopecomp'],[2.5 4])


%% plot distributions of slope
bins = -4.5:0.1:-2;
H = hist(S,bins);
figure
stairs(bins,H(:,1), 'linewidth',2,'color', 'k');
hold on
stairs(bins,H(:,2), 'linewidth',2,'color',colors(2,:));
hold off
xlabel('Slope')
ylabel('Count')
legend({'Awake' 'Anes'})
box off
nice_figure(gcf, [figure_folder 'tycho_slopehist'],[4 4])

%% plot slope timecourse
smo_len = 4*15;
S_combined = cell2mat(SLP');
%S_combined = cell2mat(OFF');
T_combined = cell2mat(FT');
M = mean(S_combined,2);
ST = std(S_combined')';
CT_combined = cell2mat(CONDT);
CT_combined = CT_combined([4 5 6 7 8 9]);
figure
plot(T_combined, M, 'color', 0.7*[1 1 1 0.5], 'linewidth',1)
hold on
plot(T_combined, smooth(M,smo_len), 'linewidth',1)
plot((CT_combined'*[1 1])', (ones(size(CT_combined))'*[-8 2])', 'k--', 'linewidth',1)
hold off
xlim([0 6000])
ylim([-5 -1])
set(gca,'xtick', [0:600:6000])
set(gca,'xtickLabel', [0:600:6000]/60)
set(gca,'ytick', [-5:1:-1])
xlabel('Time (min)')
ylabel('Slope')
box off
nice_figure(gcf, [figure_folder 'tycho_slopetime'],[8 3])

%% plot spatial topology
load('~/Documents/Data/NeuroTycho/Propofol/GridLocations/20110621KTMD_Anesthesia+and+Sleep_Chibi_Toru+Yanagawa_mat_2Dimg/ChibiMAP.mat');
%load('~/Documents/Data/NeuroTycho/Propofol/GridLocations/20110112KTMD_Anesthesia+and+Sleep_George_Toru+Yanagawa_mat_2Dimg/GeorgeMAP.mat');
S = [slp_avg{1}(:,3),slp_avg{2}(:,2)];
%S = [off_avg{1}(:,3),off_avg{2}(:,2)];
ds = S(:,2)-S(:,1);
figure
N_step = 200;
cmp = redblue(N_step);
image(I);axis equal
hold on
for i=1:128
    cur_color = cmp(round(99*ds(i)./max(abs(ds)))+100,:);
    scatter(X(i),Y(i),40, cur_color, 'filled', 'markeredgecolor','k')
end
hold off
colormap(cmp)
h=colorbar;
set(h,'ticks',[0 100 200]);
set(h,'ticklabels',[-2:2:2])
ylabel(h,'Slope Difference (Anes - Awake)')
axis off
nice_figure(gcf, [figure_folder 'tycho_slopebrain'],[5 4])

%% slope - offset plot & zoomed time course
smo_len = 4*15;
CT_combined = cell2mat(CONDT);
S_combined = smooth(mean(cell2mat(SLP'),2),smo_len);
O_combined = smooth(mean(cell2mat(OFF'),2),smo_len);
Twin = round(4*(CT_combined(4)-60*3)):round(4*(CT_combined(6)+60*4));
cmp = redblue(length(Twin));
cmp = repmat([1 1 1],length(cmp),1)-cmp;

figure
%plot slope time course
subplot(2,1,1)
scatter(Twin/4, S_combined(Twin),25,cmp,'.')
hold on
plot((CT_combined(3:6)'*[1 1])', (ones(size(CT_combined(3:6)))'*[-5 -2])', 'k--')
hold off
xlim([Twin(1) Twin(end)]/4)
set(gca,'ytick',-5:1:-2)
set(gca,'xtick', [0:300:6000])
set(gca,'xtickLabel', [0:300:6000]/60)
ylabel('Slope')
xlabel('Time (min)')

%plot in 2D slope and offset
subplot(2,1,2)
scatter(S_combined(Twin),O_combined(Twin),10,cmp, 'filled')
xlabel('Slope')
ylabel('Offset')
set(gca,'xtick',-5:1:-2)
nice_figure(gcf, [figure_folder 'tycho_slopeoffset'],[5 6])


%% SUPPLEMENTARY
%% slope & offset timecourse around anesthesia
figure
colors = get(gca,'colororder');
F = dir('*PF_*');
smo_len = 4*15; %4 samples per sec, X secs
for f=1:length(F)
    load([F(f).name '/results4060.mat'])   
    CT_combined = cell2mat(CONDT);
    S_combined = smooth(mean(cell2mat(SLP'),2),smo_len);
    O_combined = smooth(mean(cell2mat(OFF'),2),smo_len);
    Twin = round(4*(CT_combined(4)-60*3)):round(4*(CT_combined(6)+60*4));
    cmp = redblue(length(Twin));
    cmp = repmat([1 1 1],length(cmp),1)-cmp;
    
    %plot in 2D slope and offset
    subplot(4,2,f*2-1)
    scatter(S_combined(Twin),O_combined(Twin,:),10,cmp, 'filled')
    
    %plot slope time course
    subplot(4,2,f*2)
    scatter(Twin/4, S_combined(Twin),20,cmp,'.')
    hold on
    plot((CT_combined(3:6)'*[1 1])', (ones(size(CT_combined(3:6)))'*[-5 -2])', 'k--')
    hold off
    xlim([Twin(1) Twin(end)]/4)
end
ylabel('Slope')
xlabel('Time (s)')
subplot(4,2,7)
xlabel('Slope')
ylabel('Offset')

%% worst fit channels 30-50 and 40-60
sesh = 1; cond = 1;
for chan = 1:num_chan
    frange = 30:48;
    [tmp err] = polyfit(log10(frange), log10(PSD{sesh}(frange+1,chan,cond))',1);
    err3050(chan) = err.normr;
    
    frange = [40:48 52:60];
    [tmp err] = polyfit(log10(frange), log10(PSD{sesh}(frange+1,chan,cond))',1);    
    err4060(chan) = err.normr;
end
[val ind30] = sort(err3050);
[val ind40] = sort(err4060);
figure
subplot(121)
loglog(0:400,PSD{sesh}(:,ind30(end-(0:19)),cond), 'color', [1 0 0 0.6], 'linewidth', 1);
hold on
loglog(0:400,PSD{sesh}(:,ind30(1:20),cond), 'color', [0 0 0 0.6], 'linewidth', 1);
hold off
xlim([0 400])
title('Best and worst fits: 30-50Hz')
subplot(122)
loglog(0:400,PSD{sesh}(:,ind40(end-(0:19)),cond), 'color', [1 0 0 0.6], 'linewidth', 1);
hold on
loglog(0:400,PSD{sesh}(:,ind40(1:20),cond), 'color', [0 0 0 0.6], 'linewidth', 1);
hold off
xlim([0 400])
xlabel('Frequency (Hz)')
title('Best and worst fits: 40-60Hz')
figure
plot(err3050,err4060,'ok')
hold on
plot([0 0.3], [0 0.3], 'k--')
hold off
ylim([0 0.15])
xlabel('30-50Hz error')
ylabel('40-60Hz error')