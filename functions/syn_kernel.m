function kernel = syn_kernel(t, tau, type)
%function kernel = syn_kernel(t, tau, type)
% given a specific synaptic kernel type and time constant, this returns a
% time series of the kernel that spans the time defined (t) in seconds
%
% t: time vector in seconds (e.g. t=0:0.001:5)
% tau: t_decay or [t_rise t_decay] in seconds
% type: single- or double-exponential, or alpha synapse

if length(tau)==2
    type='2exp';
end
switch type
    case 'exp'
        %single decay exponential (e^-t)
        kernel = exp(-t./tau(1));
        
    case '2exp'
        %double exponential -(e^t/t_rise) + e^t/t_decay
        if length(tau)==2
            % compute the normalization factor
            tpeak = tau(2)*tau(1)/(tau(2)-tau(1))*log(tau(2)/tau(1));
            normf = 1/(-exp(-tpeak./tau(1))+exp(-tpeak./tau(2)));            
            kernel = (-exp(-t./tau(1))+exp(-t./tau(2)))*normf;
        else
            kernel = [];
            disp('Need two time constants for double exponential.')
        end
                
    case 'alpha'
        %alpha kernel (t*e^t)
        kernel = (t/tau).*exp(-t./tau);        
end
end

