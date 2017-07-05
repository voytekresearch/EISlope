function discretized = pois_spikes(sim_t, dt, N_neu, FR)
%discretized = pois_spikes(sim_t, dt, N_neu, FR)
%   simulate population spiking of N neurons firing at FR each, return a
%   single spike train that is the total spiking

%mu parameter for exponential distribution
MU = 1./(N_neu*FR);

%draw ISI from exp RV 
ISI = exprnd(MU, [((sim_t+2)/MU) 1]); %pad 2s of sim time
spk_times = cumsum(ISI);
spk_times(spk_times>sim_t)=[];

%discretize
bins = (0:dt:sim_t)+dt/2; %make discretizing bins
discretized = hist(spk_times,bins);