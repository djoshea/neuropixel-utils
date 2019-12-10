classdef KilosortPartialResort < handle

    properties(Transient)
        ks  % Neuropixel.KilosortDataset
    end
    
    properties
        fsAP % stored from ks
        sort_windows (:, 2) uint64 % nWindows x 2 list of start, stop samples which were actually resorted
        
        % spike_times.npy - [nSpikes, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        spike_times(:, 1) uint64
        
        % amplitudes.npy - [nSpikes, ] double vector with the amplitude scaling factor that was applied to the template when extracting that spike
        amplitudes(:,1) single;

        % spike_templates.npy - [nSpikes, ] uint32 vector specifying the identity of the template that was used to extract each spike
        spike_templates(:, 1) uint32
        
        % spike_templates.npy - [nSpikes, ] uint32 vector specifying the identity of the template that was originally used to extract each spike, before splitAllClusters
        spike_templates_preSplit(:, 1) uint32

        
        % spike_clusters.npy - [nSpikes, ] uint32 vector giving the cluster identity of each spike. This file is optional and
        % if not provided will be automatically created the first time you run the template gui, taking the same values as
        % spike_templates.npy until you do any merging or splitting.
        spike_clusters(:, 1) uint32

        % [nSpikesCutoff, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        cutoff_spike_times(:, 1) uint64
        cutoff_amplitudes(:, 1) single
        cutoff_spike_templates(:, 1) uint32
        cutoff_spike_templates_preSplit(:, 1) uint32
        cutoff_spike_clusters(:, 1) uint32
    end
    
    properties(Transient)
        % pc_features.npy - [nSpikes, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike.
        % The channels that those features came from are specified in pc_features_ind.npy. E.g. the value at pc_features[123, 1, 5]
        % is the projection of the 123rd spike onto the 1st PC on the channel given by pc_feature_ind[5].
        pc_features(:, :, :) single
        
        % template_features.npy - [nSpikes, nTemplateRank] single matrix giving the magnitude of the projection of each spike onto nTemplateRank other features.
        % Which other features is specified in template_feature_ind.npy
        template_features(:, :) single

        cutoff_template_features(:, :) single % [nSpikesCutoff, nTemplateRank]
        
        % [nSpikesCutoff, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike (from .cProjPC_cutoff_invalid)
        cutoff_pc_features(:, :, :) uint32
    end
    
    properties(Dependent)
        nSpikes
        nSpikesCutoff
        nSortWindows
        nSortSamples
    end
    
    methods
        function kspr = KilosortPartialResort()
            
        end
        
        function n = get.nSpikes(kspr)
            n = numel(kspr.spike_times);
        end
        
        function n = get.nSpikesCutoff(kspr)
            n = numel(kspr.cutoff_spike_times);
        end
        
        function n = get.nSortWindows(kspr)
            n = size(kspr.sort_windows, 1);
        end
        
        function n = get.nSortSamples(kspr)
            durations = kspr.sort_windows(:, 2) - kspr.sort_windows(:, 1) + uint64(1);
            n = sum(durations);
        end
        
        function splice_into_ks(kspr, ks)
            assert(isa(ks, 'Neuropixel.KilosortDataset'));
            
            % compute mask of spikes to keep
            mask_keep = true(ks.nSpikes, 1);
            mask_keep_cutoff = true(ks.nSpikesCutoff, 1);
            for iW = 1:kspr.nSortWindows
                % keep if not in this window
                mask_keep = mask_keep & ~(ks.spike_times >= kspr.sort_windows(iW, 1) & ks.spike_times <= kspr.sort_windows(iW, 2));
                mask_keep_cutoff = mask_keep_cutoff & ~(ks.cutoff_spike_times >= kspr.sort_windows(iW, 1) & ks.cutoff_spike_times <= kspr.sort_windows(iW, 2));
            end
           
            debug('Removing %d spikes / %d cutoff, inserting %d / %d in %d sort windows\n', nnz(~mask_keep), nnz(~mask_keep_cutoff), kspr.nSpikes, kspr.nSpikesCutoff, kspr.nSortWindows);
            ks.mask_spikes(mask_keep, mask_keep_cutoff);
            ks.append_spikes(kspr);
            ks.sort_spikes();
        end
        
        function [spike_idx_segmented, cutoff_spike_idx_segmented] = segment_into_windows_clusters(kspr, cluster_ids)
            % used mostly for evaluating the response, quickly segments spike times into sort_windows and cluster_ids
            % spike_idx_segmented is a nSortWindows x nClusters 
            
            windows = kspr.sort_windows;
            nWindows = size(windows, 1);
            nClusters = numel(cluster_ids);
            
            spike_idx_segmented = do_segment(kspr.spike_times, kspr.spike_clusters);
            cutoff_spike_idx_segmented = do_segment(kspr.cutoff_spike_times, kspr.cutoff_spike_clusters);
            
            function idx_segmented = do_segment(times, clusters)
                idx = (1:numel(times))';
                [mask_in_cluster, cluster_ind] = ismember(clusters, cluster_ids);
                window_ind = Neuropixel.Utils.discretize_windows(times, windows);
                mask_in_window = ~isnan(window_ind);
                
                mask = mask_in_cluster & mask_in_window;
                subs = [window_ind(mask), cluster_ind(mask)];
                idx_segmented = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(idx(mask), 1, [nWindows, nClusters], subs);
            end
        end
        
        function [spike_counts_segmented, cutoff_spike_counts_segmented] = count_by_window_cluster(kspr, cluster_ids)
            windows = kspr.sort_windows;
            nWindows = size(windows, 1);
            nClusters = numel(cluster_ids);
            
            spike_counts_segmented = do_count(kspr.spike_times, kspr.spike_clusters);
            cutoff_spike_counts_segmented = do_count(kspr.cutoff_spike_times, kspr.cutoff_spike_clusters);
            
            function counts_segmented = do_count(times, clusters)
                [mask_in_cluster, cluster_ind] = ismember(clusters, cluster_ids);
                window_ind = Neuropixel.Utils.discretize_windows(times, windows);
                mask_in_window = ~isnan(window_ind);
                
                mask = mask_in_cluster & mask_in_window;
                subs = [window_ind(mask), cluster_ind(mask)];
                counts_segmented = accumarray(subs, 1, [nWindows, nClusters]);
            end
        end
    end
    
    methods(Static)
        function kspr = construct_reextractSpikesWithFixedTemplates(ks, varargin)
            p = inputParser();
            % these are used to replace specific time windows in the raw data with 
            p.addParameter('data_replace_windows', zeros(0, 2), @(x) ismatrix(x) && size(x, 2) == 2);
            p.addParameter('data_replace', {}, @iscell);
            p.addParameter('debug_ignore_replace_data', false, @islogical);
            p.addParameter('pad_spike_extract_windows', [0 0], @isvector); % pad spike-sorted windows by this amount [pre post] relative to the data_replace_windows
            p.parse(varargin{:});
            
            assert(isa(ks, 'Neuropixel.KilosortDataset'));
            kspr = Neuropixel.KilosortPartialResort();
            kspr.ks = ks;
            kspr.fsAP = ks.fsAP;
            
            data_replace = p.Results.data_replace;
            data_replace_windows = p.Results.data_replace_windows;
            pad_spike_extract_windows = p.Results.pad_spike_extract_windows;
            spike_extract_windows = [data_replace_windows(:, 1) - pad_spike_extract_windows(1), data_replace_windows(:, 2) + pad_spike_extract_windows(2)];

            if p.Results.debug_ignore_replace_data
                data_replace = {};
                data_replace_windows = zeros(0, 2);
            end
            rez = reextractSpikesWithFixedTemplates(ks, ...
                    'data_replace', data_replace, ...
                    'data_replace_windows', data_replace_windows, ...
                    'spike_extract_windows', spike_extract_windows);
                
            % now build out the kspr fields
            cluster_offset = -1;
            
            kspr.sort_windows = spike_extract_windows;

            kspr.spike_times = rez.st3(:, 1);
            kspr.spike_templates_preSplit = rez.st3(:, 2);
            kspr.amplitudes = rez.st3(:, 3);
            kspr.spike_templates = rez.st3(:, rez.st3_template_col);
            kspr.spike_clusters = uint32(rez.st3(:, rez.st3_cluster_col) + cluster_offset);
            kspr.template_features = rez.cProj;
            kspr.pc_features = rez.cProjPC;
            
            kspr.cutoff_spike_times = rez.st3_cutoff_invalid(:, 1);
            kspr.cutoff_spike_templates_preSplit = rez.st3_cutoff_invalid(:, 2);
            kspr.cutoff_amplitudes = rez.st3_cutoff_invalid(:, 3);
            kspr.cutoff_spike_templates = rez.st3_cutoff_invalid(:, rez.st3_template_col);
            kspr.cutoff_spike_clusters = uint32(rez.st3_cutoff_invalid(:, rez.st3_cluster_col) + cluster_offset);
            kspr.cutoff_template_features = rez.cProj_cutoff_invalid;
            kspr.cutoff_pc_features = rez.cProjPC_cutoff_invalid;
        end
    end
end