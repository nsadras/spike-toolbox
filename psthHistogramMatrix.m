function psthHistogramMatrix(spikes_matrix, grid_h, grid_w, N_bins, delta)
% PSTHHISTOGRAMMATRIX - Generates histogram for each neuron in the spikes_matrix
% spikes_matrix -   matrix of spiking activity, size is n_trials x n_neurons
%                   x n_samples
% grid_h:           height of plot grid
% grid_w:           width of plot grid
% N_bins:           number of bins in
% delta:            sampling period
    
    [n_trials, n_neurons, n_samples] = size(spikes_matrix);
    window_size = n_samples/N_bins;
    
    figure()
    title('Firing rate histograms')
    
    neuron_indices = reshape(1:n_neurons, [grid_h grid_w]).'; % assume 5x5 neuron arrangement
    neuron_indices = neuron_indices(:);
    for plot_idx = 1 : n_neurons
        neuron_idx = neuron_indices(plot_idx);
        samples_matrix = squeeze(spikes_matrix(:,neuron_idx,:));
        
        bins = zeros(1, N_bins);
        for bin_idx = 0 : N_bins-1
            j = bin_idx*window_size + 1;
            while j <= window_size*(bin_idx+1)
                bins(bin_idx+1) = bins(bin_idx+1) + sum(samples_matrix(:,j));
                j = j + 1;
            end
        end

        bins = bins/(n_trials*window_size*delta);
        
        subplot(grid_h, grid_w, plot_idx);
        bar(bins)
        set(gca, 'XTick', [0 50 100], 'XTickLabel', {'0','500', '1000'});
        ylim([0 80])
        title(strcat('Neuron ', num2str(neuron_idx))) 
        xlabel('Time (ms)')
        ylabel('Hz')    
    end
end

