function [x_hat, W_hat] = ppf_strf(spikes, srf, A, Q, trf, trf_basis_fns, trf_start, stim_indic, reset_times, delta)
%PPF_STRF Summary of this function goes here
%   Detailed explanation goes here


n_neurons = size(spikes, 1);
state_size = size(srf, 1) - 1;
x0_hat = 0*ones(1, state_size);   % Initial state estimate
W0 = .0005*eye(state_size);
W_prev = W0;                    % Initial state estimate covariance
%x_hat = x0_hat.';               % Matrix of all estimated states over time
x_hat = [];
x_prev = x0_hat.';
%W_hat = zeros([state_size, state_size, size(spikes, 2)]);
W_hat = nan;

lambda_est = zeros(n_neurons, size(spikes, 2));
stim_times = find(stim_indic);


Rc = zeros(n_neurons, size(spikes, 2));
for neuron_idx = 1:n_neurons
    cur_trf = trf(neuron_idx, :); % 1 x n_weights
    r = cur_trf*trf_basis_fns;
    R = conv(stim_indic, r, 'same');
    
    R = conv(stim_indic, r, 'full');
    R = R(1:size(spikes,2));
    R = circshift(R, floor(trf_start/delta));
    Rc(neuron_idx, :) = R;
end

for spike_idx = 1:size(spikes,2)       
    % reset state if necessary
        % reset state if necessary
    if ismember(spike_idx, reset_times)
       fprintf('[PPF]: reset state\n')
       x_prev = x0_hat';
    end

 
    x_predicted = A*x_prev;         % Random walk prior   
    W_predicted = A*W_prev*A.' + Q;
    %W_predicted = A*W_prev*A.' + Q.*diag(Rc); % j a n k
   
    W_update = zeros(state_size);

    for neuron_idx = 1:n_neurons
        Rt = Rc(neuron_idx, spike_idx);
        lambda = exp([1 x_predicted.']*[srf(1,neuron_idx) ; Rt*srf(2:state_size+1,neuron_idx)]);
        W_update = W_update + (Rt^2)*srf(2:state_size+1,neuron_idx)*srf(2:state_size+1,neuron_idx).'*lambda*delta;
    end
    
    W_updated_inv = inv(W_predicted) + W_update;
    W_updated = inv(W_updated_inv);
    
    %W_hat(:,:,spike_idx) = W_updated;
        
    x_update = zeros(state_size,1);
    
    for neuron_idx = 1:n_neurons
        Rt = Rc(neuron_idx, spike_idx);
        lambda = exp([1 x_predicted.']*[srf(1,neuron_idx) ; Rt*srf(2:state_size + 1,neuron_idx)]);
        x_update = x_update + Rt*srf(2:state_size+1,neuron_idx)*(spikes(neuron_idx,spike_idx) - lambda*delta); 
        lambda_est(neuron_idx, spike_idx) = lambda;
    end
    
    x_update = W_updated_inv\x_update;   
    x_updated = x_predicted + x_update;
   
    if any(spike_idx == (stim_times+500))
        x_hat = [x_hat x_updated];
    end
    
    W_prev = W_updated;
    x_prev = x_updated;
       
end

