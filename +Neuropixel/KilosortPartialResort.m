classdef KilosortPartialResort < handle & matlab.mixin.Copyable
    % container for the results of re-running Kilosort2 on a specific sample window with fixed templates 
    % to conduct a re-extraction / sorting of spikes on altered data from specific windows of time
    
    properties(Transient)
        ks  % Neuropixel.KilosortDataset
    end
    
    properties
        cluster_ids % stored from ks
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
        
        is_deduplicated (1, 1) logical = false;
        deduplication_stats struct;
        deduplicate_spikes logical = false;
        deduplicate_cutoff_spikes logical = false;
        deduplicate_within_samples uint64 = 5;
        deduplicate_within_distance single = 50;
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
        hasFeaturesLoaded
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
        
        function tf = get.hasFeaturesLoaded(kspr)
            tf = ~isempty(kspr.pc_features);
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
        
        function [spike_idx_segmented, cutoff_spike_idx_segmented] = segment_into_windows_clusters(kspr, varargin)
            % used mostly for evaluating the response, quickly segments spike times into sort_windows and cluster_ids
            % spike_idx_segmented is a nSortWindows x nClusters 
            p = inputParser();
            p.addParameter('cluster_ids', kspr.cluster_ids, @isvector);
            p.addParameter('elide_padding', [0 0], @isvector); % in samples
            p.parse(varargin{:});
            
            cluster_ids = p.Results.cluster_ids;
            elide_padding = p.Results.elide_padding;
            
            windows = kspr.sort_windows;
            windows(:, 1) = windows(:, 1) + elide_padding(1);
            windows(:, 2) = windows(:, 2) - elide_padding(2);
            
            spike_idx_segmented = do_segment(kspr.spike_times, kspr.spike_clusters);
            cutoff_spike_idx_segmented = do_segment(kspr.cutoff_spike_times, kspr.cutoff_spike_clusters);
            
            function idx_segmented = do_segment(times, clusters)
                nWindows = size(windows, 1);
                nClusters = numel(cluster_ids);
                
                idx = (1:numel(times))';
                [mask_in_cluster, cluster_ind] = ismember(clusters, cluster_ids);
                window_ind = Neuropixel.Utils.discretize_windows(times, windows);
                mask_in_window = ~isnan(window_ind);
                
                mask = mask_in_cluster & mask_in_window;
                subs = [window_ind(mask), cluster_ind(mask)];
                idx_segmented = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(idx(mask), 1, [nWindows, nClusters], subs);
            end
        end
        
        function [spike_time_segmented_rel, cutoff_spike_time_segmented_rel] = segment_align_into_windows_clusters(kspr, varargin)
            % spike_time_segmented_rel is nSortWindows x nClusters { nTimes } where each time is relative to (start of the sort window plus elide_padding(1)) 
            % unless convert_to_ms is true, times are in samples, not ms
            p = inputParser();
            p.addParameter('sort_window_mask', true(kspr.nSortWindows, 1), @isvector);
            p.addParameter('cluster_ids', kspr.cluster_ids, @isvector);
            p.addParameter('elide_padding', [0 0], @isvector); % in samples
            p.addParameter('convert_to_ms', false, @islogical);
            p.addParameter('align_sample_offset', 0, @isscalar); % for alignment, treat window(:, 1) + align_sample_offset as time 0
            p.parse(varargin{:});
            
            % used mostly for evaluating the response, quickly segments spike times into sort_windows and cluster_ids
            % spike_idx_segmented is a nSortWindows x nClusters 
            sort_window_mask = Neuropixel.Utils.TensorUtils.vectorIndicesToMask(p.Results.sort_window_mask, kspr.nSortWindows);
            cluster_ids = p.Results.cluster_ids;
            elide_padding = p.Results.elide_padding;
            convert_to_ms = p.Results.convert_to_ms;
            align_sample_offset = p.Results.align_sample_offset;
            
            windows = kspr.sort_windows(sort_window_mask, :);
            windows(:, 1) = windows(:, 1) + elide_padding(1);
            windows(:, 2) = windows(:, 2) - elide_padding(2);
            
            sample0 = int64(windows(:, 1) + align_sample_offset);
            
            nWindows = size(windows, 1);
            nClusters = numel(cluster_ids);
                
            if nWindows == 0 || nClusters == 0
                [spike_time_segmented_rel, cutoff_spike_time_segmented_rel] = deal(cell(nWindows, nClusters));
                return;
            end
            
            spike_time_segmented_rel = do_segment(kspr.spike_times, kspr.spike_clusters);
            cutoff_spike_time_segmented_rel = do_segment(kspr.cutoff_spike_times, kspr.cutoff_spike_clusters);
            
            function time_seg_rel = do_segment(times, clusters)
                if isempty(times)
                    time_seg_rel = cell(nWindows, nClusters);
                    time_seg_rel(:) = {zeros(0, 1, 'like', times)};
                    return;
                end
            
                [mask_in_cluster, cluster_ind] = ismember(clusters, cluster_ids);
                window_ind = Neuropixel.Utils.discretize_windows(times, windows);
                mask_in_window = ~isnan(window_ind);
                
                times_rel = int64(times);
                times_rel(mask_in_window) = times_rel(mask_in_window) - sample0(window_ind(mask_in_window));
                
                if convert_to_ms
                    times_rel = single(times_rel) ./ single(kspr.fsAP / 1000);
                end
                
                mask = mask_in_cluster & mask_in_window;
                subs = [window_ind(mask), cluster_ind(mask)];
                time_seg_rel = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(times_rel(mask), 1, [nWindows, nClusters], subs);
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
    
    methods % post-sort modifications matching those in KilosortDataset
        function [mask_dup_spikes, mask_dup_spikes_cutoff, stats] = identify_duplicate_spikes(kspr, varargin)
            p = inputParser();
            p.addParameter('include_cutoff_spikes', kspr.deduplicate_cutoff_spikes, @islogical);
            p.addParameter('withinSamples', kspr.deduplicate_within_samples, @isscalar);
            p.addParameter('withinDistance', kspr.deduplicate_within_distance, @isscalar);
            p.addParameter('progress', true, @islogical);
            p.parse(varargin{:});

            assert(~isempty(kspr.ks), 'Requires .ks field set');
            assert(issorted(kspr.spike_times), 'call .sortSpikes first');
            assert(issorted(kspr.cutoff_spike_times), 'call .sortSpikes first');

            withinSamples = p.Results.withinSamples;
            withinDistance = p.Results.withinDistance;

            % lookup best channels by templates
            ks = kspr.ks;
            m = ks.computeMetrics();

            % determine which templates are withinDistance of each other templates based on centroid
            temptempdist = squareform(pdist(m.template_centroid, 'euclidean'));
            temptempprox = temptempdist < withinDistance;

            % begin by combining spikes with spikes cutoff
            if p.Results.include_cutoff_spikes
                spikes = cat(1, kspr.spike_times, kspr.cutoff_spike_times);
                templates = cat(1, kspr.spike_templates, kspr.cutoff_spike_templates);
                clusters = cat(1, kspr.spike_clusters, kspr.cutoff_spike_clusters);
                [spikes, sort_idx] = sort(spikes);
                templates = templates(sort_idx);
                clusters = clusters(sort_idx);

                mask_from_cutoff = true(size(spikes));
                mask_from_cutoff(1:kspr.nSpikes) = false;
                mask_from_cutoff = mask_from_cutoff(sort_idx);
            else
                spikes = kspr.spike_times;
                templates = kspr.spike_templates;
                clusters = kspr.spike_clusters;
                mask_from_cutoff = false(size(spikes));
            end

            cluster_inds = ks.lookup_clusterIds(clusters);

            mask_dup = false(size(spikes));
            dup_from_template = zeros(size(spikes), 'like', templates);
            dup_from_cluster_ind = zeros(size(spikes), 'like', cluster_inds);
            if p.Results.progress
                prog = Neuropixel.Utils.ProgressBar(numel(spikes), 'Checking for duplicate spikes');
            else
                prog.update = @(varargin) true;
                prog.finish = @(varargin) true;
            end
            
            for iS = 1:numel(spikes)
                if mask_dup(iS), continue, end
                temp1 = templates(iS);
                iS2 = iS+1;
                while iS2 < numel(spikes) && spikes(iS2) < spikes(iS) + withinSamples
                    temp2 = templates(iS2);
                    if temptempprox(temp1, temp2)
                        mask_dup(iS2) = true;
                        dup_from_template(iS2) = temp1;
                        dup_from_cluster_ind(iS2) = cluster_inds(iS);
                    end
                    iS2 = iS2 + 1;
                end
                if mod(iS, 1000) == 0
                    prog.update(iS);
                end
            end
            prog.finish();

            if p.Results.include_cutoff_spikes
                mask_dup_spikes = mask_dup(~mask_from_cutoff);
                mask_dup_spikes_cutoff = mask_dup(mask_from_cutoff);
            else
                mask_dup_spikes = mask_dup;
                mask_dup_spikes_cutoff = false(kspr.nSpikesCutoff, 1);
            end

            % compute stats
            mask_dup_same_template = mask_dup & dup_from_template == templates;
            stats.nSameTemplate = nnz(mask_dup_same_template);
            stats.fracRemoved = mean(mask_dup);
            stats.fracRemovedExcludingSameTemplate = mean(mask_dup & ~mask_dup_same_template);

            mask_dup_same_cluster = mask_dup & dup_from_cluster_ind == cluster_inds;
            stats.fracRemovedExcludingSameCluster = mean(mask_dup & ~mask_dup_same_cluster);

            % accumulate (i, j) is the count of spikes with template i that were marked as duplicates of template j
            stats.templateTemplateDuplicateCounts = accumarray([templates(mask_dup), dup_from_template(mask_dup)], 1, [ks.nTemplates, ks.nTemplates], @sum);
            stats.templateDuplicateCounts = sum(stats.templateTemplateDuplicateCounts, 2);
            template_spike_counts =  accumarray(templates, 1, [ks.nTemplates, 1], @sum);
            stats.templateDuplicateFrac = stats.templateDuplicateCounts ./ template_spike_counts;

            % same stats but for clusters
            stats.clusterClusterDuplicateCounts = accumarray([cluster_inds(mask_dup), dup_from_cluster_ind(mask_dup)], 1, [ks.nClusters, ks.nClusters], @sum);
            stats.clusterDuplicateCounts = sum(stats.clusterClusterDuplicateCounts, 2);
            cluster_spike_counts =  accumarray(cluster_inds, 1, [ks.nClusters, 1], @sum);
            stats.clusterDuplicateFrac = stats.clusterDuplicateCounts ./ cluster_spike_counts;
        end
        
        function stats = remove_duplicate_spikes(kspr, varargin)
            if kspr.is_deduplicated
                return;
            end

            debug('Removing duplicate spikes from KSPR dataset\n');
            [mask_dup_spikes, mask_dup_spikes_cutoff, stats] = kspr.identify_duplicate_spikes(varargin{:});
            kspr.mask_spikes(~mask_dup_spikes, ~mask_dup_spikes_cutoff);
            kspr.is_deduplicated = true;
            kspr.deduplication_stats = stats;
        end
        
        function accept_cutoff_spikes(kspr, ratings_or_cluster_ids)
            if isempty(kspr.cutoff_spike_times)
                return;
            end
            
            if isa(ratings_or_cluster_ids, 'Neuropixel.ClusterRatingInfo')
                cluster_ids = ratings_or_cluster_ids.cluster_ids(ratings_or_cluster_ids.includeCutoffSpikes); %#ok<*PROPLC>
            elseif islogical(ratings_or_cluster_ids)
                assert(numel(ratings_or_cluster_ids) == kspr.nClusters);
                cluster_ids = kspr.cluster_ids(ratings_or_cluster_ids);
            else
                cluster_ids = ratings_or_cluster_ids;
            end

            accept_cutoff_mask = ismember(kspr.cutoff_spike_clusters, cluster_ids);
            nCurrent = kspr.nSpikes;
            nAccepted = nnz(accept_cutoff_mask);
            nTotal = nAccepted + nCurrent;
            [kspr.spike_times, sortIdx] = sort(cat(1, kspr.spike_times, kspr.cutoff_spike_times(accept_cutoff_mask)));
            kspr.cutoff_spike_times = kspr.cutoff_spike_times(~accept_cutoff_mask);
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

            [kspr.spike_templates, kspr.cutoff_spike_templates] = combineAndSort(kspr.spike_templates, kspr.cutoff_spike_templates);
            [kspr.spike_templates_preSplit, kspr.cutoff_spike_templates_preSplit] = combineAndSort(kspr.spike_templates_preSplit, kspr.cutoff_spike_templates_preSplit);
            [kspr.amplitudes, kspr.cutoff_amplitudes] = combineAndSort(kspr.amplitudes, kspr.cutoff_amplitudes);
            [kspr.spike_clusters, kspr.cutoff_spike_clusters] = combineAndSort(kspr.spike_clusters, kspr.cutoff_spike_clusters);
            if kspr.hasFeaturesLoaded
                [kspr.pc_features, kspr.cutoff_pc_features] = combineAndSort(kspr.pc_features, kspr.cutoff_pc_features);
                [kspr.template_features, kspr.cutoff_template_features] = combineAndSort(kspr.template_features, kspr.cutoff_template_features);
            end
        end
        
        function drop_cutoff_spikes(kspr)
            kspr.cutoff_spike_times = [];
            kspr.cutoff_spike_templates = [];
            kspr.cutoff_spike_templates_preSplit = [];
            kspr.cutoff_amplitudes = [];
            kspr.cutoff_spike_clusters = [];
            if kspr.hasFeaturesLoaded
                kspr.cutoff_pc_features = [];
                kspr.cutoff_template_features = [];
            end
        end
        
        function apply_cluster_merge(kspr, mergeInfo)
            % apply the merges in clusterMergeInfo
            assert(isa(mergeInfo, 'Neuropixel.ClusterMergeInfo'));

            spike_clusters = kspr.spike_clusters;
            cutoff_spike_clusters = kspr.cutoff_spike_clusters;
            for iM = 1:mergeInfo.nMerges
                spike_clusters = apply_single_merge(spike_clusters, mergeInfo.new_cluster_ids(iM), mergeInfo.merges{iM});
                cutoff_spike_clusters = apply_single_merge(cutoff_spike_clusters, mergeInfo.new_cluster_ids(iM), mergeInfo.merges{iM});
            end
            kspr.spike_clusters = spike_clusters;
            kspr.cutoff_spike_clusters = cutoff_spike_clusters;

            function spike_clusters = apply_single_merge(spike_clusters, dst_cluster_id, src_cluster_ids)
                mask_assign_to_dst = ismember(spike_clusters, src_cluster_ids);
                spike_clusters(mask_assign_to_dst) = dst_cluster_id;
            end
        end
        
        function [clusterInds, cluster_ids] = lookup_clusterIds(kspr, cluster_ids)
            if islogical(cluster_ids)
                cluster_ids = kspr.cluster_ids(cluster_ids);
             end
            [tf, clusterInds] = ismember(cluster_ids, kspr.cluster_ids);
            assert(all(tf, 'all'), 'Some cluster ids were not found in kspr.clusterids');
        end
        
        function mask_clusters(kspr, cluster_ids)
            [~, cluster_ids] = kspr.lookup_clusterIds(cluster_ids);

            mask = ismember(kspr.spike_clusters, cluster_ids);
            cutoff_mask = ismember(kspr.cutoff_spike_clusters, cluster_ids);
            kspr.mask_spikes(mask, cutoff_mask);
        end
        
        function mask_spikes(kspr, mask, mask_cutoff)
            assert(islogical(mask) && numel(mask) == kspr.nSpikes);
            assert(islogical(mask_cutoff) && numel(mask_cutoff) == kspr.nSpikesCutoff);

            kspr.spike_times = kspr.spike_times(mask);
            kspr.spike_templates = kspr.spike_templates(mask);
            if ~isempty(kspr.spike_templates_preSplit)
                kspr.spike_templates_preSplit = kspr.spike_templates_preSplit(mask);
            end
            kspr.amplitudes = kspr.amplitudes(mask);
            kspr.spike_clusters = kspr.spike_clusters(mask);
            
            kspr.cutoff_spike_times = kspr.cutoff_spike_times(mask_cutoff);
            kspr.cutoff_spike_templates = kspr.cutoff_spike_templates(mask_cutoff);
            kspr.cutoff_spike_templates_preSplit = kspr.cutoff_spike_templates_preSplit(mask_cutoff);
            kspr.cutoff_amplitudes = kspr.cutoff_amplitudes(mask_cutoff);
            kspr.cutoff_spike_clusters = kspr.cutoff_spike_clusters(mask_cutoff);
            
            if kspr.hasFeaturesLoaded
                kspr.pc_features = kspr.pc_features(mask, :, :);
                kspr.template_features = kspr.template_features(mask, :);
                kspr.cutoff_pc_features = kspr.cutoff_pc_features(mask_cutoff, :, :);
                kspr.cutoff_template_features = kspr.cutoff_template_features(mask_cutoff, :);
            end
        end
        
        function cluster_ids_retained = apply_subgroup_cutoff_merges_selection(kspr, subgroup, cluster_rating_info, cluster_merge_info, varargin)
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
            assert(ignore_cluster_ratings || isa(cluster_rating_info, 'Neuropixel.ClusterRatingInfo'));
            assert(ignore_cluster_merges || isa(cluster_merge_info, 'Neuropixel.ClusterMergeInfo'));
            
            if verbose
                print = @(varargin) debug(varargin{:});
            else
                print = @(varargin) 1;
            end
            
            if ~ignore_cluster_ratings && ~ignore_cutoff_spikes
                % do post-hoc modifications to ks according to user ratings in cluster_rating_info'
                % 1. accept cutoff spikes from selected units
                print('Accepting cutoff spikes from %d / %d clusters\n', nnz(cluster_rating_info.includeCutoffSpikes), cluster_rating_info.nClusters);
                kspr.accept_cutoff_spikes(cluster_rating_info);
            end
            kspr.drop_cutoff_spikes();
            
            if ~ignore_cluster_merges
                % 2. apply merges
                print('Performing %d cluster merges\n', cluster_merge_info.nMerges);
                kspr.apply_cluster_merge(cluster_merge_info);
                
                if ~ignore_cluster_ratings
                    % 3. apply the cluster merges to the ratings list, allowing it to pick the best rating for each cluster
                    cluster_rating_info_merged = copy(cluster_rating_info);
                    cluster_rating_info_merged.apply_cluster_merge(cluster_merge_info);
                end
            elseif ~ignore_cluster_ratings
                cluster_rating_info_merged = cluster_rating_info;
            end
            
            if ~isempty(remove_duplicate_cluster_ids)
                print('Removing %d / %d clusters as detected duplicates\n', numel(remove_duplicate_cluster_ids), numel(kspr.cluster_ids));
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
                
                kspr.mask_clusters(cluster_ids_retained);
                
                cluster_rating_info_merged.lookupClusterRatings(cluster_ids_retained);
            else
                cluster_ids_retained = kspr.cluster_ids;
                
                % remove duplicates if requested
                cluster_ids_retained = setdiff(cluster_ids_retained, remove_duplicate_cluster_ids);
                
                kspr.mask_clusters(cluster_ids_retained);
            end
        end
    end
    
    methods(Static) % Construction by re-sorting
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
            
            assert(isa(ks, 'Neuropixel.KilosortDataset'));
            kspr = Neuropixel.KilosortPartialResort();
            kspr.ks = ks;
            kspr.fsAP = ks.fsAP;
            kspr.cluster_ids = ks.cluster_ids;
            kspr.deduplicate_spikes = ks.deduplicate_spikes;
            kspr.deduplicate_cutoff_spikes = ks.deduplicate_cutoff_spikes;
            kspr.deduplicate_within_samples = ks.deduplicate_within_samples;
            kspr.deduplicate_within_distance = ks.deduplicate_within_distance;
            
            data_replace = p.Results.data_replace;
            data_replace_windows = p.Results.data_replace_windows;
            pad_spike_extract_windows = p.Results.pad_spike_extract_windows;
            spike_extract_windows = p.Results.spike_extract_windows;
            if ismember('spike_extract_windows', p.UsingDefaults)
                if ismember('data_replace_windows', p.UsingDefaults)
                    error('Neither data_replace_windows nor spike_extract_windows specified, not sure what to resort');
                end
                spike_extract_windows = [data_replace_windows(:, 1) - pad_spike_extract_windows(1), data_replace_windows(:, 2) + pad_spike_extract_windows(2)];
            else
                spike_extract_windows = [spike_extract_windows(:, 1) - pad_spike_extract_windows(1), spike_extract_windows(:, 2) + pad_spike_extract_windows(2)];
            end

            if p.Results.debug_ignore_replace_data
                data_replace = {};
                data_replace_windows = zeros(0, 2);
            end

            if ~isempty(spike_extract_windows)
                % to deduplicate correctly at the edges, we need a little bit of extra padding that we'll strip here
                if kspr.deduplicate_spikes
                    dedup_padding = kspr.deduplicate_within_samples * 2;
                end
                spike_extract_windows_dedup_padded = [spike_extract_windows(:, 1) - dedup_padding, spike_extract_windows(:, 2) + dedup_padding];
            
                rez = Kilosort2.MainLoop.reextractSpikesWithFixedTemplates(ks, ...
                        'data_replace', data_replace, ...
                        'data_replace_windows', data_replace_windows, ...
                        'spike_extract_windows', spike_extract_windows_dedup_padded);
                
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
                
                if kspr.deduplicate_spikes
                    % deduplicate the spikes (which still include spikes in the padding windows
                    kspr.remove_duplicate_spikes();

                    % and now trim the spikes again to remove the dedup_padding
                    mask_within_windows = any(kspr.spike_times >= spike_extract_windows(:, 1)' & kspr.spike_times <= spike_extract_windows(:, 2)', 2); % nSpikes x nWindows --> nSpikes
                    cutoff_mask_within_windows = any(kspr.cutoff_spike_times >= spike_extract_windows(:, 1)' & kspr.cutoff_spike_times <= spike_extract_windows(:, 2)', 2); % nSpikes x nWindows --> nSpikes
                    kspr.mask_spikes(mask_within_windows, cutoff_mask_within_windows);
                end
            else
                if ~p.Results.suppressEmptySpikeWindowWarning
                    warning('Empty spike_extract_windows provided, no spikes were reextracted')
                end
                kspr.sort_windows = spike_extract_windows;
            end
        end
        
        function kspr = construct_resegmentKilosortDatasetUnchanged(ks, varargin)
            p = inputParser();
            % these are used to replace specific time windows in the raw data with 
            p.addParameter('pad_spike_extract_windows', [0 0], @isvector); % pad spike-sorted windows by this amount [pre post] relative to the data_replace_windows
            % use this to simply reextract spikes from specific windows, useful when not replacing data in situ
            p.addParameter('spike_extract_windows', zeros(0, 2), @(x) ismatrix(x) && size(x, 2) == 2);
            p.addParameter('suppressEmptySpikeWindowWarning', false, @islogical);
            p.parse(varargin{:});
            
            assert(isa(ks, 'Neuropixel.KilosortDataset'));
            kspr = Neuropixel.KilosortPartialResort();
            kspr.ks = ks;
            kspr.fsAP = ks.fsAP;
            kspr.cluster_ids = ks.cluster_ids;
            kspr.deduplicate_spikes = ks.deduplicate_spikes;
            kspr.deduplicate_cutoff_spikes = ks.deduplicate_cutoff_spikes;
            kspr.deduplicate_within_samples = ks.deduplicate_within_samples;
            kspr.deduplicate_within_distance = ks.deduplicate_within_distance;
            
            pad_spike_extract_windows = p.Results.pad_spike_extract_windows;
            spike_extract_windows = p.Results.spike_extract_windows;
            if ismember('spike_extract_windows', p.UsingDefaults)
                error('spike_extract_windows must be specified, not sure what to resort');
            end
            spike_extract_windows = [spike_extract_windows(:, 1) - pad_spike_extract_windows(1), spike_extract_windows(:, 2) + pad_spike_extract_windows(2)];
            
            if ~isempty(spike_extract_windows)
                % we use a segmented dataset to do the extraction for us and then copy the fields over
                nWindows = size(spike_extract_windows, 1);
                tsi = Neuropixel.TrialSegmentationInfo(nWindows, ks.fsAP);
                tsi.trialId = 1:nWindows;
                tsi.conditionId = ones(nWindows, 1);
                tsi.idxStart = spike_extract_windows(:, 1);
                tsi.idxStop = spike_extract_windows(:, 2);
                
                seg = Neuropixel.KilosortTrialSegmentedDataset(ks, tsi, tsi.trialId, ...
                    'segment_by_trials', true, 'segment_by_clusters', false, ... % KSPR keeps the spikes from all clusters together
                    'loadFeatures', true, 'loadSync', false, 'cluster_ids', ks.cluster_ids);

                kspr.sort_windows = spike_extract_windows;

                % seg is naturally segmented by trials, but kspr stores all its resorted windows contiguously so that it resembles
                % a KilosortDataset object internally, so here we re concatenate them. 
                kspr.spike_times = cat(1, seg.spike_times{:});
                kspr.spike_templates_preSplit = cat(1, seg.spike_templates_preSplit{:});
                kspr.amplitudes = cat(1, seg.amplitudes{:});
                kspr.spike_templates = cat(1, seg.spike_templates{:});
                kspr.spike_clusters = cat(1, seg.spike_clusters{:});
                kspr.template_features = cat(1, seg.template_features{:});
                kspr.pc_features = cat(1, seg.pc_features{:});
                
                kspr.cutoff_spike_times = cat(1, seg.cutoff_spike_times{:});
                kspr.cutoff_spike_templates_preSplit = cat(1, seg.cutoff_spike_templates_preSplit{:});
                kspr.cutoff_amplitudes = cat(1, seg.cutoff_amplitudes{:});
                kspr.cutoff_spike_templates = cat(1, seg.cutoff_spike_templates{:});
                kspr.cutoff_spike_clusters = cat(1, seg.cutoff_spike_clusters{:});
                kspr.cutoff_template_features = cat(1, seg.cutoff_template_features{:});
                kspr.cutoff_pc_features = cat(1, seg.cutoff_pc_features{:});
            else
                if ~p.Results.suppressEmptySpikeWindowWarning
                    warning('Empty spike_extract_windows provided, no spikes were reextracted')
                end
                kspr.sort_windows = spike_extract_windows;
            end
        end
    end
end
