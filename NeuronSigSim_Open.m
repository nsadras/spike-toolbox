function [rec_spike,rec_exp_RV] = NeuronSigSim_Open(rates,delta)

% This program is a simulation of spike neuron signal under open loop framework. Only one trial and one neuron.
% The first spike value is always zero (consider it as a determinant initial state)
% 
% Input:
% rates: instantaneous firing rate over time
% delta: the spike time interval in second.
%        
% Output:
% rec_spike:  Records the neuron spike signal along the time so its length must be the same as "time length". It is a row vector.
% rec_exp_RV: exponential R.V. record. I use -1 as the initial value for distinguishing the ones been assigned with the ones not.

% First, I set up some initial parameters.

n_samples = size(rates,2);
rec_spike = zeros(1,n_samples); % I assume the spike value at time=0 (rec_spike(1)) is always 0.
rec_exp_RV = repmat(-1,[1,n_samples]);

spike_index = 1;

% Now, I use a loop to find the next spike index continuously till it over the number of state.  

two_points_integral = [0,0]; % set up the vector for recording two boundary points in the trapezoid integral.
two_points_integral(1) = rates(spike_index);
while(spike_index < n_samples)

    exp_RV = exprnd(1);
    rec_exp_RV(spike_index+1) = exp_RV;
    partial_sum = 0;
    
    while (partial_sum < exp_RV && spike_index < n_samples)
        spike_index = spike_index+1;
        two_points_integral(2) = rates(spike_index); % I update the second term in the integral one step.
        partial_sum = partial_sum+trapz(two_points_integral)*delta;
        two_points_integral(1) = two_points_integral(2);
    end
    
    if (partial_sum >= exp_RV)
        rec_spike(spike_index)=1;
    end
    
end


end