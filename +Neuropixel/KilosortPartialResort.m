classdef KilosortPartialResort < handle

    properties(Transient)
        ks  % Neuropixel.KilosortDataset
    end
    
    properties
        sample_mask(:, 1) logical % [nSamples] sparse logical matrix indicating the timepoints for which resorted data is available
        
        batch_mask(:, 1) logical % [nBatches] sparse matrix indicating which batches were resorted
        
        % amplitudes.npy - [nSpikes, ] double vector with the amplitude scaling factor that was applied to the template when extracting that spike
        amplitudes(:,1) double;

        % pc_features.npy - [nSpikes, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike.
        % The channels that those features came from are specified in pc_features_ind.npy. E.g. the value at pc_features[123, 1, 5]
        % is the projection of the 123rd spike onto the 1st PC on the channel given by pc_feature_ind[5].
        pc_features(:, :, :) single

        % spike_templates.npy - [nSpikes, ] uint32 vector specifying the identity of the template that was used to extract each spike
        spike_templates(:, 1) uint32
        
        % spike_templates.npy - [nSpikes, ] uint32 vector specifying the identity of the template that was originally used to extract each spike, before splitAllClusters
        spike_templates_preSplit(:, 1) uint32

        % spike_times.npy - [nSpikes, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        spike_times(:, 1) uint64

        % template_features.npy - [nSpikes, nTemplateRank] single matrix giving the magnitude of the projection of each spike onto nTemplateRank other features.
        % Which other features is specified in template_feature_ind.npy
        template_features(:, :) single

        % spike_clusters.npy - [nSpikes, ] uint32 vector giving the cluster identity of each spike. This file is optional and
        % if not provided will be automatically created the first time you run the template gui, taking the same values as
        % spike_templates.npy until you do any merging or splitting.
        spike_clusters(:, 1) uint32

        % [nSpikesCutoff, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        cutoff_spike_times(:, 1) uint64
        cutoff_amplitudes(:, 1) double
        cutoff_spike_templates(:, 1) uint32
        cutoff_spike_templates_preSplit(:, 1) uint32
        cutoff_spike_clusters(:, 1) uint32
        
        % [nSpikesCutoff, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike (from .cProjPC_cutoff_invalid)
        cutoff_pc_features(:, :, :) uint32
    end
    
    methods
        function ksr = KilosortPartialResort()
            
        end
        
        
    end
    
    methods(Static)
        
    end
end