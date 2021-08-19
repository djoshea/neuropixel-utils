classdef KilosortPartialResort < handle & matlab.mixin.Copyable
    
    properties (Transient)
        ks  % npxutils.KilosortDataset
    end
    
    properties
        cluster_ids % stored from ks
        fsAP % stored from ks
        sort_windows (:,2) uint64 % nWindows x 2 list of start, stop samples which were actually resorted
        
        % spike_times.npy - [nSpikes, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        spike_times (:,1) uint64
        
        % amplitudes.npy - [nSpikes, ] double vector with the amplitude scaling factor that was applied to the template when extracting that spike
        amplitudes (:,1) single;
        
        % spike_templates.npy - [nSpikes, ] uint32 vector specifying the identity of the template that was used to extract each spike
        spike_templates (:,1) uint32
        
        % spike_templates.npy - [nSpikes, ] uint32 vector specifying the identity of the template that was originally used to extract each spike, before splitAllClusters
        spike_templates_preSplit (:,1) uint32
        
        % spike_clusters.npy - [nSpikes, ] uint32 vector giving the cluster identity of each spike. This file is optional and
        % if not provided will be automatically created the first time you run the template gui, taking the same values as
        % spike_templates.npy until you do any merging or splitting.
        spike_clusters (:,1) uint32
        
        % [nSpikesCutoff, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        cutoff_spike_times (:,1) uint64
        cutoff_amplitudes (:,1) single
        cutoff_spike_templates (:,1) uint32
        cutoff_spike_templates_preSplit (:,1) uint32
        cutoff_spike_clusters (:,1) uint32
    end
    
    properties (Transient)
        % pc_features.npy - [nSpikes, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike.
        % The channels that those features came from are specified in pc_features_ind.npy. E.g. the value at pc_features[123, 1, 5]
        % is the projection of the 123rd spike onto the 1st PC on the channel given by pc_feature_ind[5].
        pc_features (:,:,:) single
        
        % template_features.npy - [nSpikes, nTemplateRank] single matrix giving the magnitude of the projection of each spike onto nTemplateRank other features.
        % Which other features is specified in template_feature_ind.npy
        template_features (:,:) single
        
        cutoff_template_features (:,:) single % [nSpikesCutoff, nTemplateRank]
        
        % [nSpikesCutoff, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike (from .cProjPC_cutoff_invalid)
        cutoff_pc_features(:,:,:) uint32
    end
    
    properties (Dependent)
        hasFeaturesLoaded
        nSpikes
        nSpikesCutoff
        nSortWindows
        nSortSamples
    end
    
    methods
        function this = KilosortPartialResort()
            
        end
        
        function n = get.nSpikes(this)
            n = numel(this.spike_times);
        end
        
        function n = get.nSpikesCutoff(this)
            n = numel(this.cutoff_spike_times);
        end
        
        function n = get.nSortWindows(this)
            n = size(this.sort_windows, 1);
        end
        
        function n = get.nSortSamples(this)
            durations = this.sort_windows(:, 2) - this.sort_windows(:,1) + uint64(1);
            n = sum(durations);
        end
        
        function tf = get.hasFeaturesLoaded(this)
            tf = ~isempty(this.pc_features);
        end
        
        function splice_into_ks(this, ks)
            arguments
                this
                ks (1,1) npxutils.KilosortDataset
            end
            
            % compute mask of spikes to keep
            mask_keep = true(ks.nSpikes, 1);
            mask_keep_cutoff = true(ks.nSpikesCutoff, 1);
            for iW = 1:this.nSortWindows
                % keep if not in this window
                mask_keep = mask_keep & ~(ks.spike_times >= this.sort_windows(iW, 1) & ks.spike_times <= this.sort_windows(iW, 2));
                mask_keep_cutoff = mask_keep_cutoff & ~(ks.cutoff_spike_times >= this.sort_windows(iW, 1) & ks.cutoff_spike_times <= this.sort_windows(iW, 2));
            end
            
            debug('Removing %d spikes / %d cutoff, inserting %d / %d in %d sort windows\n', nnz(~mask_keep), nnz(~mask_keep_cutoff), this.nSpikes, this.nSpikesCutoff, this.nSortWindows);
            ks.mask_spikes(mask_keep, mask_keep_cutoff);
            ks.append_spikes(this);
            ks.sort_spikes();
        end
        
        function [spike_idx_segmented, cutoff_spike_idx_segmented] ...
                = segment_into_windows_clusters(this, varargin)
            % used mostly for evaluating the response, quickly segments spike times into sort_windows and cluster_ids
            % spike_idx_segmented is a nSortWindows x nClusters
            p = inputParser();
            p.addParameter('cluster_ids', this.cluster_ids, @isvector);
            p.addParameter('elide_padding', [0 0], @isvector); % in samples
            p.parse(varargin{:});
            
            cluster_ids = p.Results.cluster_ids;
            elide_padding = p.Results.elide_padding;
            
            windows = this.sort_windows;
            windows(:,1) = windows(:,1) + elide_padding(1);
            windows(:, 2) = windows(:, 2) - elide_padding(2);
            
            spike_idx_segmented = do_segment(this.spike_times, this.spike_clusters);
            cutoff_spike_idx_segmented = do_segment(this.cutoff_spike_times, this.cutoff_spike_clusters);
            
            function idx_segmented = do_segment(times, clusters)
                nWindows = size(windows, 1);
                nClusters = numel(cluster_ids);
                
                idx = (1:numel(times))';
                [mask_in_cluster, cluster_ind] = ismember(clusters, cluster_ids);
                window_ind = npxutils.internal.discretize_windows(times, windows);
                mask_in_window = ~isnan(window_ind);
                
                mask = mask_in_cluster & mask_in_window;
                subs = [window_ind(mask), cluster_ind(mask)];
                idx_segmented = npxutils.internal.TensorUtils.splitAlongDimensionBySubscripts(idx(mask), 1, [nWindows, nClusters], subs);
            end
        end
        
        function [spike_time_segmented_rel, cutoff_spike_time_segmented_rel] ...
                = segment_align_into_windows_clusters(this, varargin)
            % spike_time_segmented_rel is nSortWindows x nClusters { nTimes } where each time is relative to (start of the sort window plus elide_padding(1))
            % unless convert_to_ms is true, times are in samples, not ms
            p = inputParser();
            p.addParameter('sort_window_mask', true(this.nSortWindows, 1), @isvector);
            p.addParameter('cluster_ids', this.cluster_ids, @isvector);
            p.addParameter('elide_padding', [0 0], @isvector); % in samples
            p.addParameter('convert_to_ms', false, @islogical);
            p.addParameter('align_sample_offset', 0, @isscalar); % for alignment, treat window(:,1) + align_sample_offset as time 0
            p.parse(varargin{:});
            
            % used mostly for evaluating the response, quickly segments spike times into sort_windows and cluster_ids
            % spike_idx_segmented is a nSortWindows x nClusters
            sort_window_mask = npxutils.internal.TensorUtils.vectorIndicesToMask(p.Results.sort_window_mask, this.nSortWindows);
            cluster_ids = p.Results.cluster_ids;
            elide_padding = p.Results.elide_padding;
            convert_to_ms = p.Results.convert_to_ms;
            align_sample_offset = p.Results.align_sample_offset;
            
            windows = this.sort_windows(sort_window_mask, :);
            windows(:,1) = windows(:,1) + elide_padding(1);
            windows(:, 2) = windows(:, 2) - elide_padding(2);
            
            sample0 = int64(windows(:,1) + align_sample_offset);
            
            nWindows = size(windows, 1);
            nClusters = numel(cluster_ids);
            
            if nWindows == 0 || nClusters == 0
                [spike_time_segmented_rel, cutoff_spike_time_segmented_rel] = deal(cell(nWindows, nClusters));
                return;
            end
            
            spike_time_segmented_rel = do_segment(this.spike_times, this.spike_clusters);
            cutoff_spike_time_segmented_rel = do_segment(this.cutoff_spike_times, this.cutoff_spike_clusters);
            
            function time_seg_rel = do_segment(times, clusters)
                if isempty(times)
                    time_seg_rel = cell(nWindows, nClusters);
                    time_seg_rel(:) = {zeros(0, 1, 'like', times)};
                    return;
                end
                
                [mask_in_cluster, cluster_ind] = ismember(clusters, cluster_ids);
                window_ind = npxutils.internal.discretize_windows(times, windows);
                mask_in_window = ~isnan(window_ind);
                
                times_rel = int64(times);
                times_rel(mask_in_window) = times_rel(mask_in_window) - sample0(window_ind(mask_in_window));
                
                if convert_to_ms
                    times_rel = single(times_rel) ./ single(this.fsAP / 1000);
                end
                
                mask = mask_in_cluster & mask_in_window;
                subs = [window_ind(mask), cluster_ind(mask)];
                time_seg_rel = npxutils.internal.TensorUtils.splitAlongDimensionBySubscripts(times_rel(mask), 1, [nWindows, nClusters], subs);
            end
        end
        
        function [spike_counts_segmented, cutoff_spike_counts_segmented] = count_by_window_cluster(this, cluster_ids)
            windows = this.sort_windows;
            nWindows = size(windows, 1);
            nClusters = numel(cluster_ids);
            
            spike_counts_segmented = do_count(this.spike_times, this.spike_clusters);
            cutoff_spike_counts_segmented = do_count(this.cutoff_spike_times, this.cutoff_spike_clusters);
            
            function counts_segmented = do_count(times, clusters)
                [mask_in_cluster, cluster_ind] = ismember(clusters, cluster_ids);
                window_ind = npxutils.internal.discretize_windows(times, windows);
                mask_in_window = ~isnan(window_ind);
                
                mask = mask_in_cluster & mask_in_window;
                subs = [window_ind(mask), cluster_ind(mask)];
                counts_segmented = accumarray(subs, 1, [nWindows, nClusters]);
            end
        end
    end
    
    methods % post-sort modifications matching those in KilosortDataset
        function accept_cutoff_spikes(this, ratings_or_cluster_ids)
            if isempty(this.cutoff_spike_times)
                return;
            end
            
            if isa(ratings_or_cluster_ids, 'npxutils.ClusterRatingInfo')
                cluster_ids = ratings_or_cluster_ids.cluster_ids(ratings_or_cluster_ids.includeCutoffSpikes); %#ok<*PROPLC>
            elseif islogical(ratings_or_cluster_ids)
                assert(numel(ratings_or_cluster_ids) == this.nClusters);
                cluster_ids = this.cluster_ids(ratings_or_cluster_ids);
            else
                cluster_ids = ratings_or_cluster_ids;
            end
            
            accept_cutoff_mask = ismember(this.cutoff_spike_clusters, cluster_ids);
            nCurrent = this.nSpikes;
            nAccepted = nnz(accept_cutoff_mask);
            nTotal = nAccepted + nCurrent;
            [this.spike_times, sortIdx] = sort(cat(1, this.spike_times, this.cutoff_spike_times(accept_cutoff_mask)));
            this.cutoff_spike_times = this.cutoff_spike_times(~accept_cutoff_mask);
            [~, insertOrigAt] = ismember((1:nCurrent)', sortIdx);
            [~, insertCutoffAt] = ismember((nCurrent+1:nTotal)', sortIdx);
            
            function [combined, cutoff] = combineAndSort(orig, cutoff)
                sz = size(orig);
                sz(1) = nTotal;
                combined = zeros(sz, 'like', orig);
                combined(insertOrigAt, :, :, :) = orig;
                combined(insertCutoffAt, :, :, :) = cutoff(accept_cutoff_mask, :, :, :);
                cutoff = cutoff(~accept_cutoff_mask, :, :, :);
            end
            
            [this.spike_templates, this.cutoff_spike_templates] = combineAndSort(this.spike_templates, this.cutoff_spike_templates);
            [this.spike_templates_preSplit, this.cutoff_spike_templates_preSplit] = combineAndSort(this.spike_templates_preSplit, this.cutoff_spike_templates_preSplit);
            [this.amplitudes, this.cutoff_amplitudes] = combineAndSort(this.amplitudes, this.cutoff_amplitudes);
            [this.spike_clusters, this.cutoff_spike_clusters] = combineAndSort(this.spike_clusters, this.cutoff_spike_clusters);
            if this.hasFeaturesLoaded
                [this.pc_features, this.cutoff_pc_features] = combineAndSort(this.pc_features, this.cutoff_pc_features);
                [this.template_features, this.cutoff_template_features] = combineAndSort(this.template_features, this.cutoff_template_features);
            end
        end
        
        function drop_cutoff_spikes(this)
            this.cutoff_spike_times = [];
            this.cutoff_spike_templates = [];
            this.cutoff_spike_templates_preSplit = [];
            this.cutoff_amplitudes = [];
            this.cutoff_spike_clusters = [];
            if this.hasFeaturesLoaded
                this.cutoff_pc_features = [];
                this.cutoff_template_features = [];
            end
        end
        
        function apply_cluster_merge(this, mergeInfo)
            % apply the merges in clusterMergeInfo
            assert(isa(mergeInfo, 'npxutils.ClusterMergeInfo'));
            
            spike_clusters = this.spike_clusters;
            cutoff_spike_clusters = this.cutoff_spike_clusters;
            for iM = 1:mergeInfo.nMerges
                spike_clusters = apply_single_merge(spike_clusters, mergeInfo.new_cluster_ids(iM), mergeInfo.merges{iM});
                cutoff_spike_clusters = apply_single_merge(cutoff_spike_clusters, mergeInfo.new_cluster_ids(iM), mergeInfo.merges{iM});
            end
            this.spike_clusters = spike_clusters;
            this.cutoff_spike_clusters = cutoff_spike_clusters;
            
            function spike_clusters = apply_single_merge(spike_clusters, dst_cluster_id, src_cluster_ids)
                mask_assign_to_dst = ismember(spike_clusters, src_cluster_ids);
                spike_clusters(mask_assign_to_dst) = dst_cluster_id;
            end
        end
        
        function [clusterInds, cluster_ids] = lookup_clusterIds(this, cluster_ids)
            if islogical(cluster_ids)
                cluster_ids = this.cluster_ids(cluster_ids);
            end
            [tf, clusterInds] = ismember(cluster_ids, this.cluster_ids);
            assert(all(tf, 'all'), 'Some cluster ids were not found in kspr.clusterids');
        end
        
        function mask_clusters(this, cluster_ids)
            [~, cluster_ids] = this.lookup_clusterIds(cluster_ids);
            
            mask = ismember(this.spike_clusters, cluster_ids);
            cutoff_mask = ismember(this.cutoff_spike_clusters, cluster_ids);
            this.mask_spikes(mask, cutoff_mask);
        end
        
        function mask_spikes(this, mask, mask_cutoff)
            assert(islogical(mask) && numel(mask) == this.nSpikes);
            assert(islogical(mask_cutoff) && numel(mask_cutoff) == this.nSpikesCutoff);
            
            this.spike_times = this.spike_times(mask);
            this.spike_templates = this.spike_templates(mask);
            if ~isempty(this.spike_templates_preSplit)
                this.spike_templates_preSplit = this.spike_templates_preSplit(mask);
            end
            this.amplitudes = this.amplitudes(mask);
            this.spike_clusters = this.spike_clusters(mask);
            
            this.cutoff_spike_times = this.cutoff_spike_times(mask_cutoff);
            this.cutoff_spike_templates = this.cutoff_spike_templates(mask_cutoff);
            this.cutoff_spike_templates_preSplit = this.cutoff_spike_templates_preSplit(mask_cutoff);
            this.cutoff_amplitudes = this.cutoff_amplitudes(mask_cutoff);
            this.cutoff_spike_clusters = this.cutoff_spike_clusters(mask_cutoff);
            
            if this.hasFeaturesLoaded
                this.pc_features = this.pc_features(mask, :, :);
                this.template_features = this.template_features(mask, :);
                this.cutoff_pc_features = this.cutoff_pc_features(mask_cutoff, :, :);
                this.cutoff_template_features = this.cutoff_template_features(mask_cutoff, :);
            end
        end
        
        function cluster_ids_retained = apply_subgroup_cutoff_merges_selection(this, ...
                subgroup, cluster_rating_info, cluster_merge_info, varargin)
            p = inputParser();
            p.addParameter('ignore_cluster_ratings', false, @islogical); % used primarily for debugging, allowing this to be tested even when no cluster is rated
            p.addParameter('ignore_cluster_merges', false, @islogical); % used primarily for debugging, allowing this to be tested even when no cluster is rated
            p.addParameter('ignore_cutoff_spikes', false, @islogical); % used primarily for waveform extraction, only effective when ignore_cluster_ratings is false
            p.addParameter('remove_duplicate_cluster_ids', [], @(x) isempty(x) || isvector(x));
            p.addParameter('ratings_accept', ["good", "unstable"], @isstring)
            p.addParameter('verbose', true, @islogical);
            p.parse(varargin{:});
            
            ratings_accept = string(p.Results.ratings_accept);
            ignore_cluster_ratings = p.Results.ignore_cluster_ratings;
            ignore_cluster_merges = p.Results.ignore_cluster_merges;
            ignore_cutoff_spikes = p.Results.ignore_cutoff_spikes;
            remove_duplicate_cluster_ids = p.Results.remove_duplicate_cluster_ids;
            verbose = p.Results.verbose;
            assert(ignore_cluster_ratings || isa(cluster_rating_info, 'npxutils.ClusterRatingInfo'));
            assert(ignore_cluster_merges || isa(cluster_merge_info, 'npxutils.ClusterMergeInfo'));
            
            if verbose
                print = @(varargin) debug(varargin{:});
            else
                print = @(varargin) 1;
            end
            
            if ~ignore_cluster_ratings && ~ignore_cutoff_spikes
                % do post-hoc modifications to ks according to user ratings in cluster_rating_info'
                % 1. accept cutoff spikes from selected units
                print('Accepting cutoff spikes from %d / %d clusters\n', nnz(cluster_rating_info.includeCutoffSpikes), cluster_rating_info.nClusters);
                this.accept_cutoff_spikes(cluster_rating_info);
            end
            this.drop_cutoff_spikes();
            
            if ~ignore_cluster_merges
                % 2. apply merges
                print('Performing %d cluster merges\n', cluster_merge_info.nMerges);
                this.apply_cluster_merge(cluster_merge_info);
                
                if ~ignore_cluster_ratings
                    % 3. apply the cluster merges to the ratings list, allowing it to pick the best rating for each cluster
                    cluster_rating_info_merged = copy(cluster_rating_info);
                    cluster_rating_info_merged.apply_cluster_merge(cluster_merge_info);
                end
            elseif ~ignore_cluster_ratings
                cluster_rating_info_merged = cluster_rating_info;
            end
            
            if ~isempty(remove_duplicate_cluster_ids)
                print('Removing %d / %d clusters as detected duplicates\n', numel(remove_duplicate_cluster_ids), numel(this.cluster_ids));
            else
                remove_duplicate_cluster_ids = zeros(0, 1, 'uint32');
            end
            
            if ~ignore_cluster_ratings
                % 4. subselect units acceptable for this subgroup / these subgroups
                if ~ignore_cluster_merges
                    [num_clusters_unrated, num_clusters_postMerge] = cluster_rating_info.computeClusterUnratedCountAfterApplyingMerges(cluster_merge_info);
                else
                    num_clusters_unrated = cluster_rating_info.nClustersUnrated;
                    num_clusters_postMerge = cluster_rating_info.nClusters;
                end
                if num_clusters_unrated > 0
                    warning('%d / %d clusters post-merge are unrated in cluster_rating_info', num_clusters_unrated, num_clusters_postMerge);
                end
                
                cluster_ids_retained = cluster_rating_info_merged.listClusterIdsUsableAcrossSubgroupsWithRating(subgroup, ratings_accept);
                print('Subselecting %d/ %d %s clusters post-merge for subgroup %s\n', numel(cluster_ids_retained), num_clusters_postMerge, strjoin(ratings_accept, "+"), subgroup);
                
                % remove duplicates if requested
                cluster_ids_retained = setdiff(cluster_ids_retained, remove_duplicate_cluster_ids);
                
                this.mask_clusters(cluster_ids_retained);
                
                cluster_rating_info_merged.lookupClusterRatings(cluster_ids_retained);
            else
                cluster_ids_retained = this.cluster_ids;
                
                % remove duplicates if requested
                cluster_ids_retained = setdiff(cluster_ids_retained, remove_duplicate_cluster_ids);
                
                this.mask_clusters(cluster_ids_retained);
            end
        end
    end
    
    methods (Static)
        function kspr = construct_reextractSpikesWithFixedTemplates(ks, varargin)
            p = inputParser();
            % these are used to replace specific time windows in the raw data with
            p.addParameter('data_replace_windows', zeros(0, 2), @(x) ismatrix(x) && size(x, 2) == 2);
            p.addParameter('data_replace', {}, @iscell);
            p.addParameter('debug_ignore_replace_data', false, @islogical);
            p.addParameter('pad_spike_extract_windows', [0 0], @isvector); % pad spike-sorted windows by this amount [pre post] relative to the data_replace_windows
            
            % use this to simply reextract spikes from specific windows, useful when not replacing data in situ
            p.addParameter('spike_extract_windows', zeros(0, 2), @(x) ismatrix(x) && size(x, 2) == 2);
            p.addParameter('suppressEmptySpikeWindowWarning', false, @islogical);
            p.parse(varargin{:});
            
            assert(isa(ks, 'npxutils.KilosortDataset'));
            kspr = npxutils.KilosortPartialResort();
            kspr.ks = ks;
            kspr.fsAP = ks.fsAP;
            kspr.cluster_ids = ks.cluster_ids;
            
            data_replace = p.Results.data_replace;
            data_replace_windows = p.Results.data_replace_windows;
            pad_spike_extract_windows = p.Results.pad_spike_extract_windows;
            spike_extract_windows = p.Results.spike_extract_windows;
            if ismember('spike_extract_windows', p.UsingDefaults)
                if ismember('data_replace_windows', p.UsingDefaults)
                    error('Neither data_replace_windows nor spike_extract_windows specified, not sure what to resort');
                end
                spike_extract_windows = [data_replace_windows(:,1) - pad_spike_extract_windows(1), data_replace_windows(:, 2) + pad_spike_extract_windows(2)];
            else
                spike_extract_windows = [spike_extract_windows(:,1) - pad_spike_extract_windows(1), spike_extract_windows(:, 2) + pad_spike_extract_windows(2)];
            end
            
            if p.Results.debug_ignore_replace_data
                data_replace = {};
                data_replace_windows = zeros(0, 2);
            end
            
            if ~isempty(spike_extract_windows)
                rez = reextractSpikesWithFixedTemplates(ks, ...
                    'data_replace', data_replace, ...
                    'data_replace_windows', data_replace_windows, ...
                    'spike_extract_windows', spike_extract_windows);
                
                % now build out the kspr fields
                cluster_offset = -1;
                
                kspr.sort_windows = spike_extract_windows;
                
                kspr.spike_times = rez.st3(:,1);
                kspr.spike_templates_preSplit = rez.st3(:, 2);
                kspr.amplitudes = rez.st3(:, 3);
                kspr.spike_templates = rez.st3(:, rez.st3_template_col);
                kspr.spike_clusters = uint32(rez.st3(:, rez.st3_cluster_col) + cluster_offset);
                kspr.template_features = rez.cProj;
                kspr.pc_features = rez.cProjPC;
                
                kspr.cutoff_spike_times = rez.st3_cutoff_invalid(:,1);
                kspr.cutoff_spike_templates_preSplit = rez.st3_cutoff_invalid(:, 2);
                kspr.cutoff_amplitudes = rez.st3_cutoff_invalid(:, 3);
                kspr.cutoff_spike_templates = rez.st3_cutoff_invalid(:, rez.st3_template_col);
                kspr.cutoff_spike_clusters = uint32(rez.st3_cutoff_invalid(:, rez.st3_cluster_col) + cluster_offset);
                kspr.cutoff_template_features = rez.cProj_cutoff_invalid;
                kspr.cutoff_pc_features = rez.cProjPC_cutoff_invalid;
            else
                if ~p.Results.suppressEmptySpikeWindowWarning
                    warning('Empty spike_extract_windows provided, no spikes were reextracted')
                end
                kspr.sort_windows = spike_extract_windows;
            end
        end
    end
end
