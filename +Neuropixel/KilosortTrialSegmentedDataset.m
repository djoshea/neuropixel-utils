classdef KilosortTrialSegmentedDataset < handle & matlab.mixin.Copyable

    % Properties that are copied over but not segmented into trials
    properties
        dataset % KilosortDataset
        
        % optional modes indicating which mode of segmentation is used
        is_segmented_by_trials logical = true;
        is_segmented_by_clusters logical = true;

        trial_ids(:, 1) uint32

        % nTrials x 1
        trial_has_data(:, 1) logical

        % nTrials x 1
        trial_start(:, 1) uint64

        % nTrials x 1
        trial_stop(:, 1) uint64

        % cluster ids corresponding to each column of the {nTrials, nClusters} properties
        cluster_ids(:, 1) uint32

        cluster_groups(:, 1) categorical

        % sync channel segmented by trials
        sync(:, 1) cell

        syncBitNames(:, 1) string
    end

    properties(Dependent)
        raw_dataset
        nTrials
        nTrialsHaveData
        nClusters
        nChannelsSorted
        channel_ids_sorted % which channel_ids were used for sorting

        % nTrials x 1
        trial_duration_ms
        
        spike_counts (:, :) % nTrials x nClusters
        trial_has_nonzero_spikes

        fsAP % sampling rate pass thru
        
        spike_clusters (:,:) cell
        cutoff_spike_clusters (:, :) cell
        
        spike_times_rel_start (:, :) cell % spike time of each spike in samples from trial start
    end

    % Properties that are segmented by trial
    % each of these is nTrials x nTemplates cell with the same inner size as in Dataset (as described)
    properties
         % indices into master dataset: trials x clusters
        spike_idx(:, :) cell

        %  [nSpikes, ] double vector with the amplitude scaling factor that was applied to the template when extracting that spike
        amplitudes(:,:) cell

        % [nSpikes, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike.
        % The channels that those features came from are specified in pc_features_ind.npy. E.g. the value at pc_features[123, 1, 5]
        % is the projection of the 123rd spike onto the 1st PC on the channel given by pc_feature_ind[5].
        pc_features(:,:) cell

        % [nSpikes, ] uint64 vector giving the sample index of each spike
        % in the raw data file
        spike_times(:, :) cell

        % [nSpikes, ] double vector giving the spike time of each spike in ms from trial start
        spike_times_ms_rel_start(:, :) cell

        % [nSpikes, nTempFeatures] single matrix giving the magnitude of the projection of each spike onto nTempFeatures other features.
        % Which other features is specified in template_feature_ind.npybrew
        template_features(:, :) cell

        % [nSpikes, ] uint32 vector specifying the identity of the template that was used to extract each spike
        spike_templates(:, :) cell
        spike_templates_preSplit (:, :) cell

        % these are only useful if segment_by_trials or segment_by_clusters are false
        spike_trial_ind(:, :) cell
        spike_cluster_inds(:, :) cell
        
        % same as the above but for spikes dropped by set_cutoff, where nSpikes --> nSpikesCutoff
        cutoff_spike_idx(:, :) cell
        cutoff_amplitudes(:, :) cell
        cutoff_spike_times(:, :)  cell
        cutoff_spike_times_ms_rel_start(:, :) cell
        cutoff_spike_templates(:, :) cell
        cutoff_spike_templates_preSplit (:, :) cell
        cutoff_template_features(:, :) cell
        cutoff_pc_features (:, :) cell
        
        cutoff_spike_trial_inds(:, :) cell
        cutoff_spike_cluster_inds(:, :) cell  
    end

    methods
        function seg = KilosortTrialSegmentedDataset(ks, tsi, trial_ids, varargin)
            p = inputParser();
            p.addParameter('loadCutoff', true, @islogical);
            p.addParameter('loadFeatures', false, @islogical);
            p.addParameter('loadSync', false, @islogical);
            p.addParameter('loadBatchwise', false, @islogical);
            p.addParameter('loadPreSplit', false, @islogical);
            p.addParameter('cluster_ids', [], @(x) isempty(x) || isvector(x)); % if specified, subselect the clusters
            
            % these may be turned off, in which case everythinng will be 1 x nClusters or nTrials x 1 instead of nTrials x nClusters
            p.addParameter('segment_by_trials', true, @islogical);
            p.addParameter('segment_by_clusters', true, @islogical);
            
            p.parse(varargin{:});

            loadCutoff = p.Results.loadCutoff;
            loadFeatures = p.Results.loadFeatures;
            loadSync = p.Results.loadSync;
            segment_by_trials = p.Results.segment_by_trials;
            segment_by_clusters = p.Results.segment_by_clusters;

            % trial_ids specifies the id of each trial that will appear in
            % this segmented dataset. tsi has its own trialId field, and
            % the data will be copied over where these ids match. But the
            % ultimate size will be set according to trial_ids
            ks.load('loadCutoff', loadCutoff, 'loadFeatures', p.Results.loadFeatures, ...
                'loadBatchwise', p.Results.loadBatchwise, 'loadPreSplit', p.Results.loadPreSplit);

            seg.dataset = ks;
            seg.trial_ids = trial_ids;

            % trials here are over the requested trial_ids, which are different
            % from those in trialInfo (which comes from the sync line and is limited
            % to the trials in the neuropixel file)
            nTrials = numel(trial_ids);

            [mask_trial_in_tsi, tsi_ind_each_local_trial] = ismember(trial_ids, tsi.trialId);
            if ~any(mask_trial_in_tsi)
                warning('No trial ids found in TSI were included');
                return;
            end

            seg.trial_has_data = mask_trial_in_tsi;

            seg.trial_start = nan(nTrials, 1);
            seg.trial_start(mask_trial_in_tsi) = tsi.idxStart(tsi_ind_each_local_trial(mask_trial_in_tsi));

            seg.trial_stop = nan(nTrials, 1);
            seg.trial_stop(mask_trial_in_tsi) = tsi.idxStop(tsi_ind_each_local_trial(mask_trial_in_tsi));

            % lookup from columns of all the nTrials x nUnits back into template ids
            seg.cluster_ids = ks.cluster_ids;
            seg.cluster_groups = ks.cluster_groups;
            if isempty(p.Results.cluster_ids)
                seg.cluster_ids = ks.cluster_ids;
            else
                % use the specified cluster_ids list, but check that they all exist in ks.cluster_ids
                seg.cluster_ids = p.Results.cluster_ids;
                found_in_dataset = ismember(seg.cluster_ids, ks.cluster_ids);
                if ~all(found_in_dataset)
                    warning('%d / %d cluster_ids manually specified for KilosortTrialSegmentedDataset are not found in KilosortDataset.cluster_ids', ...
                        nnz(~found_in_dataset), numel(found_in_dataset));
                end
            end
            seg.cluster_groups = seg.dataset.cluster_groups;

            nUnits = numel(seg.cluster_ids);

            % figure out which trial each spike in spike_times belongs
            % indices into TSI's list of trials, now we need to map into seg.trial_ids list
            local_trial_ind_each_spike = get_local_trial_ind_each_time(tsi, trial_ids, seg.dataset.spike_times);

            % which cluster does each spike belong to
            [spike_mask, unit_idx_each_spike] = ismember(seg.dataset.spike_clusters, seg.cluster_ids);

            subs = [local_trial_ind_each_spike, unit_idx_each_spike];
            
            % same for cutoff spikesks.
            if loadCutoff
                cutoff_local_trial_ind_each_spike = get_local_trial_ind_each_time(tsi, trial_ids, seg.dataset.cutoff_spike_times);
                [cutoff_spike_mask, cutoff_unit_idx_each_spike] = ismember(seg.dataset.cutoff_spike_clusters, seg.cluster_ids);
                cutoff_subs = [cutoff_local_trial_ind_each_spike, cutoff_unit_idx_each_spike];
            else
                cutoff_subs = zeros(2, 0);
            end
            
            if segment_by_trials
                nTrialsSegment = nTrials;
                if ~segment_by_clusters
                    subs(:, 2) = 1;
                    cutoff_subs(:, 2) = 1;
                    nUnitsSegment = 1;
                else
                    nUnitsSegment = nUnits;
                end
            else
                nTrialsSegment = 1;
                if segment_by_clusters
                    subs(:, 1) = 1;
                    cutoff_subs(:, 1) = 1;
                    nUnitsSegment = nUnits;
                else
                    error('Either segment_by_trials or segment_by_clusters must be true');
                end
            end
            seg.is_segmented_by_trials = segment_by_trials;
            seg.is_segmented_by_clusters = segment_by_clusters;

            nToLoad = 3;
            if p.Results.loadFeatures
                nToLoad = nToLoad + 3;
            end
            if p.Results.loadSync
                nToLoad = nToLoad + 1;
            end
            prog = Neuropixel.Utils.ProgressBar(nToLoad, 'Segmenting trials: spike_times');
            prog.increment();
            spike_times_grouped = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                ks.spike_times(spike_mask), 1, [nTrialsSegment, nUnitsSegment], subs(spike_mask, :));

            seg.spike_times = spike_times_grouped;

            if loadCutoff
                cutoff_spike_times_grouped = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    ks.cutoff_spike_times(cutoff_spike_mask), 1, [nTrialsSegment, nUnitsSegment], cutoff_subs(cutoff_spike_mask, :));
                seg.cutoff_spike_times = cutoff_spike_times_grouped;
            end

            % convert samples to ms relative to trial_start
            if segment_by_trials
                trial_start = seg.trial_start;
            else
                trial_start = 0;
            end
            for iT = 1:nTrialsSegment
                if segment_by_trials && ~seg.trial_has_data(iT), continue, end
                for iU = 1:nUnitsSegment
                    if isempty(spike_times_grouped{iT, iU})
                        spike_times_grouped{iT, iU} = nan(0, 1, 'single');
                    else
                        spike_times_grouped{iT, iU} = single(spike_times_grouped{iT, iU} - trial_start(iT)) / single(seg.dataset.sample_rate / 1000);
                    end

                    if loadCutoff
                        if isempty(cutoff_spike_times_grouped{iT, iU})
                            cutoff_spike_times_grouped{iT, iU} = nan(0, 1, 'single');
                        else
                            cutoff_spike_times_grouped{iT, iU} = single(cutoff_spike_times_grouped{iT, iU} - trial_start(iT)) / single(seg.dataset.sample_rate / 1000);
                        end
                    end
                end
            end
            seg.spike_times_ms_rel_start = spike_times_grouped;
            if loadCutoff
                seg.cutoff_spike_times_ms_rel_start = cutoff_spike_times_grouped;
            end

            prog.increment('Segmenting trials: spike_idx');
            spike_idx = (1:seg.dataset.nSpikes)';
            seg.spike_idx = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                spike_idx(spike_mask), 1, [nTrialsSegment, nUnitsSegment], subs(spike_mask, :));
            if loadCutoff
                cutoff_spike_idx = (1:seg.dataset.nSpikesCutoff)';
                seg.cutoff_spike_idx = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    cutoff_spike_idx(cutoff_spike_mask), 1, [nTrialsSegment, nUnitsSegment], cutoff_subs(cutoff_spike_mask, :));
            end
            
            if ~segment_by_trials
                prog.increment('Segmenting trials: spike trial inds');
                seg.spike_trial_ind = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    local_trial_ind_each_spike(spike_mask), 1, [nTrialsSegment, nUnitsSegment], subs(spike_mask, :));
                if loadCutoff
                    seg.cutoff_spike_trial_inds = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                        cutoff_local_trial_ind_each_spike(cutoff_spike_mask), 1, [nTrialsSegment, nUnitsSegment], cutoff_subs(cutoff_spike_mask, :));
                end
            end
            
            if ~segment_by_clusters
                prog.increment('Segmenting trials: spike cluster ids');
                seg.spike_cluster_inds = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    unit_idx_each_spike(spike_mask), 1, [nTrialsSegment, nUnitsSegment], subs(spike_mask, :));
                if loadCutoff
                    seg.cutoff_spike_cluster_inds = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                        cutoff_unit_idx_each_spike(cutoff_spike_mask), 1, [nTrialsSegment, nUnitsSegment], cutoff_subs(cutoff_spike_mask, :));
                end
            end
            
            if loadFeatures
                % slice the other fields into trials x unit:
                prog.increment('Segmenting trials: amplitudes');
                seg.amplitudes = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                ks.amplitudes(spike_mask), 1, [nTrialsSegment, nUnitsSegment], subs(spike_mask, :));
                if loadCutoff
                    seg.cutoff_amplitudes = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                        ks.cutoff_amplitudes(cutoff_spike_mask), 1, [nTrialsSegment, nUnitsSegment], cutoff_subs(cutoff_spike_mask, :));
                end

                prog.increment('Segmenting trials:  pc_features');
                seg.pc_features = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    ks.pc_features(spike_mask, :, :), 1, [nTrialsSegment, nUnitsSegment], subs(spike_mask, :));
                if loadCutoff
                    seg.cutoff_pc_features = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                        ks.cutoff_pc_features(cutoff_spike_mask, :, :), 1, [nTrialsSegment, nUnitsSegment], cutoff_subs(cutoff_spike_mask, :));
                end

                prog.increment('Segmenting trials: template_features');
                seg.template_features = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    ks.template_features(spike_mask, :), 1, [nTrialsSegment, nUnitsSegment], subs(spike_mask, :));
                if loadCutoff
                    seg.cutoff_template_features = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                        ks.cutoff_template_features(cutoff_spike_mask, :), 1, [nTrialsSegment, nUnitsSegment], cutoff_subs(cutoff_spike_mask, :));
                end
            end

            prog.increment('Segmenting trials: spike_clusters');
            seg.spike_templates = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                ks.spike_templates(spike_mask), 1, [nTrialsSegment, nUnitsSegment], subs(spike_mask, :));
            if ~isempty(ks.spike_templates_preSplit)
                seg.spike_templates_preSplit = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    ks.spike_templates_preSplit(spike_mask), 1, [nTrialsSegment, nUnitsSegment], subs(spike_mask, :));
            end
            if loadCutoff
                seg.cutoff_spike_templates = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    ks.cutoff_spike_templates(cutoff_spike_mask), 1, [nTrialsSegment, nUnitsSegment], cutoff_subs(cutoff_spike_mask, :));
                seg.cutoff_spike_templates_preSplit = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    ks.cutoff_spike_templates_preSplit(cutoff_spike_mask), 1, [nTrialsSegment, nUnitsSegment], cutoff_subs(cutoff_spike_mask, :));
            end

            seg.syncBitNames = ks.syncBitNames;

            if loadSync && segment_by_trials
                prog.increment('Segmenting trials: sync');
                sync = ks.readSync();
                sample_vec = uint32(1:numel(sync))';
                local_trial_ind_each_sample = get_local_trial_ind_each_time(tsi, trial_ids, sample_vec);

                seg.sync = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    sync, 1, nTrialsSegment, local_trial_ind_each_sample);
            else
                seg.sync = {};
            end

            prog.finish();
            fprintf('\n');

            function local_trial_ind_each_time = get_local_trial_ind_each_time(tsi, trial_ids, times)
                [~, trial_id_each_time, in_trial] = tsi.segmentTimes(times); % indices into TSI's list of trials, now we need to map into seg.trial_ids list

                [mask_in_trial_ids, local_trial_ind_each_time] = ismember(trial_id_each_time, trial_ids);
                local_trial_ind_each_time = single(local_trial_ind_each_time);
                local_trial_ind_each_time(~in_trial) = NaN;
                local_trial_ind_each_time(~mask_in_trial_ids) = NaN;
            end
        end

        function n = get.nTrials(seg)
            n = numel(seg.trial_ids);
        end

        function tf = get.nTrialsHaveData(seg)
            tf = nnz(seg.trial_has_data);
        end

        function n = get.nClusters(seg)
            n = numel(seg.cluster_ids);
        end

        function n = get.nChannelsSorted(seg)
            n = numel(seg.channel_ids_sorted);
        end

        function n = get.channel_ids_sorted(seg)
            n = seg.dataset.channel_ids_sorted;
        end

        function rd = get.raw_dataset(seg)
            rd = seg.dataset.raw_dataset;
        end

        function d = get.trial_duration_ms(seg)
            % trial_start/stop are in samples at fsAP/1000 samples per second
            d = double(seg.trial_stop - seg.trial_start) / double(seg.fsAP) * 1000;
        end

        function fsAP = get.fsAP(ds)
            fsAP = ds.dataset.fsAP;
        end
        
        function spike_clusters = get.spike_clusters(seg)
            spike_clusters = cellfun(@(inds) seg.cluster_ids(inds), seg.spike_cluster_inds, 'UniformOutput', false);
        end
        
        function spike_clusters = get.cutoff_spike_clusters(seg)
            spike_clusters = cellfun(@(inds) seg.cluster_ids(inds), seg.cutoff_spike_cluster_inds, 'UniformOutput', false);
        end

        function spike_times_rel_start = get.spike_times_rel_start(seg)
            [nT, nU] = size(seg.spike_times);
            spike_times_rel_start = cell([nT nU]);
            for iT = 1:nT
                for iU = 1:nU
                    spike_times_rel_start{iT, iU} = seg.spike_times{iT, iU} - seg.trial_start(iT);
                end
            end
        end

        function idx = lookupSyncBitByName(seg, names)
            if ischar(names)
                names = {names};
            end
            assert(iscellstr(names));

            [tf, idx] = ismember(names, seg.syncBitNames);
            idx(~tf) = NaN;
        end

        function setSyncBitNames(ds, idx, names)
            if ~isempty(ds.dataset)
                ds.dataset.setSyncBitNames(idx, names);
            else
                if isscalar(idx) && ischar(names)
                    ds.syncBitNames{idx} = names;
                else
                    assert(iscellstr(names) || isstring(names))
                    ds.syncBitNames(idx) = names;
                end
            end
        end

        function syncBitNames = get.syncBitNames(ds)
            if isempty(ds.dataset)
                if isempty(ds.syncBitNames)
                    syncBitNames = strings(16, 1);
                else
                    syncBitNames = string(ds.syncBitNames);
                end
            else
                syncBitNames = string(ds.dataset.syncBitNames);
            end
        end

        function [clusterInds, cluster_ids] = lookup_clusterIds(ds, cluster_ids)
            if islogical(cluster_ids)
                cluster_ids = ds.cluster_ids(cluster_ids);
             end
            [tf, clusterInds] = ismember(cluster_ids, ds.cluster_ids);
            assert(all(tf), 'Some cluster ids were not found in ds.clusterids');
        end

        function [channelInds, channelIds] = lookup_channelIds(ds, channelIds)
             if islogical(channelIds)
                channelIds = ds.channel_ids(channelIds);
             end
            [tf, channelInds] = ismember(channelIds, ds.channel_ids);
            assert(all(tf), 'Some channel ids not found');
        end
        
        function counts = get.spike_counts(ds) % nTrials x nClusters
            counts = cellfun(@numel, ds.spike_times);
        end
        
        function tf = get.trial_has_nonzero_spikes(ds)
            tf = sum(ds.spike_counts, 2) > 0;
        end
    end

    % pulling things from raw data
    methods
        function snippet_set = getWaveformsFromRawData(seg, varargin)
            % mask_cell is nTrials x nClusters cell of mask or indices over
            % spike times within that bin, if omitted, all spikes will be
            % included

            p = inputParser();
            % from KilosortDataset.readAPSnippets
            p.addOptional('mask_cell', {}, @(x) isempty(x) || iscell(x));
            p.addParameter('cluster_ids', seg.cluster_ids, @isvector); % which clusters each column of mask_cell correspond to
            p.addParameter('channel_ids_by_cluster', [], @(x) isempty(x) || ismatrix(x)); % specify a subset of channels to extract
            p.addParameter('best_n_channels', NaN, @isscalar); % or take the best n channels based on this clusters template when cluster_id is scalar
            p.addParameter('data_distrust_mask', [], @(x) isempty(x) || islogical(x)); % if provided, these samples will be flagged when included in snippets via ss.data_trust_mask

            % other params:
            p.addParameter('window', [-40 41], @isvector); % Number of samples before and after spiketime to include in waveform
            p.addParameter('car', false, @islogical);
            p.addParameter('centerUsingFirstSamples', 20, @(x) isscalar(x) || islogical(x)); % subtract mean of each waveform's first n samples, don't do if false
            p.addParameter('num_waveforms', Inf, @isscalar); % caution: Inf will request ALL waveforms in order (typically useful if spike_times directly specified)

            p.addParameter('primary_template_only', true, @islogical); % only sample waveforms that use the clusters primary template
            
            % time consuming step to remove the contribution of the other clusters to a given snippet. 
            % "auto" uses a mode that subtracts only if it decreases the waveform's variance
            p.addParameter('subtractOtherClusters', false, @(x) islogical(x) || isstringlike(x)); 
            p.addParameter('excludeClusterFromOwnReconstruction', true, @islogical); % avoid reconstructing a cluster with itself (set false if a cell has many spikes in close proximity
             
            p.addParameter('raw_dataset', seg.raw_dataset, @(x) true);
            p.addParameter('fromSourceDatasets', false, @islogical); % go all the way back to the imecDatasets that were concatenated to form ks.raw_dataset
            p.addParameter('from_cutoff_spikes', false, @islogical);
            p.addParameter('random_seed', 'shuffle');

            p.parse(varargin{:});

            % lookup cluster_ids from
            cluster_ids = p.Results.cluster_ids;
            [cluster_ind, cluster_ids] = seg.lookup_clusterIds(cluster_ids);
            nClu = numel(cluster_ids);

            % specify all spikes of clusters if no mask provided
            mask_cell = p.Results.mask_cell;
            if isempty(mask_cell)
                mask_cell = cellfun(@(times) true(size(times)), seg.spike_times(:, cluster_ind), 'UniformOutput', false);
            end

            % check size
            assert(size(mask_cell, 1) == seg.nTrials, 'mask_cell must be nTrials along dim 1');
            assert(size(mask_cell, 2) == nClu, 'mask_cell must be nClusters along dim 2');

            % check that mask sizes match
            idx = seg.spike_idx(:, cluster_ind);
            for iT = 1:seg.nTrials
                for iC = 1:nClu
                    assert(numel(idx{iT, iC}) == numel(mask_cell{iT, iC}), 'Spike mask_cell{%d, %d} mismatched with seg.spike_idx, trial id %d, cluster id %d', ...
                        iT, iC, seg.trial_ids(iT), cluster_ids(iC));
                end
            end
                
            % assemble spike_times we want to collect, all at once
            if ~p.Results.from_cutoff_spikes
                masked_spike_idx = cellfun(@(spike_idx, mask) spike_idx(mask), seg.spike_idx(:, cluster_ind), mask_cell, 'UniformOutput', false);
            else
                masked_spike_idx = cellfun(@(spike_idx, mask) spike_idx(mask), seg.cutoff_spike_idx(:, cluster_ind), mask_cell, 'UniformOutput', false);
            end
            [masked_spike_idx, whichCell] = Neuropixel.Utils.TensorUtils.catWhich(1, masked_spike_idx{:});
            trialIdMat = repmat(seg.trial_ids, 1, numel(cluster_ind));
            trial_idx = trialIdMat(whichCell);

            args = rmfield(p.Results, 'mask_cell');
            snippet_set = seg.dataset.getWaveformsFromRawData('spike_idx', masked_spike_idx, ...
                'cluster_ids', cluster_ids, 'trial_idx', trial_idx, args);
        end

        function snippet_set = getSnippetsFromRawData(seg, rel_start_ms_each_trial, duration_or_window_ms, varargin)
            % rel_start_ms_each_trial is nTrials x 1  in ms relative to
            % the start of each trial, specifying the window to grab on
            % each trial. rel_start_ms_each_trial should have the same length
            % nTrials as 'trialIdx', which defaults to 1:seg.nTrials. If
            % start_ms_each_trial is NaN, the snippet will be set to NaN
            %
            % if duration_or_window is scalar, the window will be [0
            % duration_or_window], else the window will
            % rel_start_ms_each_trial + (window(1):window(2))

            p = inputParser();
            p.addParameter('trial_idx', 1:seg.nTrials, @isvector);
            p.addParameter('channel_ids', seg.channel_ids, @isvector);
            p.addParameter('raw_dataset', seg.raw_dataset, @(x) true);
            p.parse(varargin{:});

            trial_idx = Neuropixel.Utils.TensorUtils.vectorMaskToIndices(p.Results.trial_idx);
            nTrials = numel(trial_idx); %#ok<*PROPLC>
            ms_to_samples =  seg.raw_dataset.fsAP / 1000;

            assert(numel(rel_start_ms_each_trial) == nTrials);

            if isscalar(duration_or_window_ms )
                window_ms = [0 duration_or_window_ms];
            else
                window_ms = duration_or_window_ms;
            end
            window_samples = round(window_ms * ms_to_samples);

            % deal with only non nans with raw dataset
            rel_start_ms_each_trial(~seg.trial_has_data(trial_idx)) = NaN;
            mask_non_nan = ~isnan(rel_start_ms_each_trial);
            rel_start_ms_each_trial = rel_start_ms_each_trial(mask_non_nan);

            % zero time of each snippet in samples
            trial_starts = seg.trial_start(trial_idx);
            trial_starts = trial_starts(mask_non_nan);
            req_zero = uint64(round(rel_start_ms_each_trial * ms_to_samples)) + uint64(trial_starts);

            channel_ids = p.Results.channel_ids;

            % inflate back to full trials
            snippet_set = p.Results.raw_dataset.readAPSnippetSet(req_zero, window_samples, channel_ids);
            snippet_set.data = Neuropixel.Utils.TensorUtils.inflateMaskedTensor(...
                snippet_set.data, 3, mask_non_nan, 0);

            snippet_set.valid = mask_non_nan;
            snippet_set.trial_idx = trial_idx;
        end
    end

    methods
        function tRelStart = convertIdxEachTrialToMsRelStart(seg, idxEachTrial)
            tRelStart = double(idxEachTrial - seg.trial_start) / double(seg.fsAP) * 1000;
        end

        function idx = convertMsRelStartToIdx(seg, tRelStart)
            idx = uint64(round(tRelStart / 1000 * seg.fsAP)) + seg.trial_start;
        end

        function [rates, duration, tStart, tStop] = computeFiringRateEachTrial(seg, varargin)
            p = inputParser();
            p.addParameter('tStart', zeros(seg.nTrials, 1), @isvector);
            p.addParameter('tStop', inf(seg.nTrials, 1), @isvector);
            p.parse(varargin{:});

            % keep times within trial boundaries
            tStart = p.Results.tStart;
            tStop = p.Results.tStop;
            mask = ~isnan(tStart) & ~isnan(tStop);

            tStart = max(0, tStart);
            tStop = min(seg.trial_duration_ms, p.Results.tStop);
            duration = tStop - tStart;

            rates = nan(seg.nTrials, seg.nClusters);
            for iT = 1:seg.nTrials
                if ~mask(iT), continue, end
                for iC = 1:seg.nClusters
                    spikes = seg.spike_times_ms_rel_start{iT, iC};
                    count = nnz(spikes >= tStart(iT) & spikes <= tStop(iT));
                    rates(iT, iC) = count / duration(iT) * 1000;
                end
            end
        end

        function css = buildClusterStabilitySummary(seg, varargin)
            p = inputParser();
            p.addParameter('tStart', zeros(seg.nTrials, 1), @isvector);
            p.addParameter('tBinWidth', 500, @isscalar);
            p.addParameter('condition_ids', zeros(seg.nTrials, 1), @isvector);
            p.parse(varargin{:});

            tBinWidth = p.Results.tBinWidth;
            tStart = p.Results.tStart;
            tStop = tStart + tBinWidth;
            [rates, ~, tStart] = seg.computeFiringRateEachTrial('tStart', tStart, 'tStop', tStop);

            css = Neuropixel.ClusterStabilitySummary();
            css.trial_ids = seg.trial_ids;
            css.trial_has_data = seg.trial_has_data & ~isnan(tStart) & ~isnan(tStop);
            css.condition_ids = p.Results.condition_ids;
            css.tBinWidth = tBinWidth;
            css.cluster_ids = seg.cluster_ids;
            css.rates = rates;
            css.idxStart = seg.convertMsRelStartToIdx(tStart);
            css.fs = seg.fsAP;
        end
        
        function [nSpikes, nCutoffSpikes] = computeTotalSpikeCounts(seg)
            nSpikes = sum(cellfun(@numel, seg.spike_idx), 'all');
            nCutoffSpikes = sum(cellfun(@numel, seg.cutoff_spike_idx), 'all');
        end
    end

%     methods(Static)
%         function seg = emptyWithSize(nTrials, nUnits)
%             seg = Kilosort.TrialSegmentedDataset;
%
%             seg.trial_ids = nan(nTrials, 1);
%             seg.trial_has_data = false(nTrials, 1);
%             seg.trial_start = zeros(nTrials, 1, 'uint64');
%             seg.cluster_ids = zeros(nTrials, 1, 'uint32');
%
%             c = cell(nTrials, nUnits);
%             seg.spike_idx = c;
%             seg.amplitudes = c;
%             seg.pc_features = c;
%             seg.spike_times_ms_rel_start = c;
%             seg.template_features = c;
%             seg.spike_cluster_inds = c;
%         end
%
%         % merge over trials
%         function dsmerge = mergeDatasets(varargin)
%             if numel(varargin) == 1
%                 dsmerge = varargin{1};
%                 return;
%             end
%
%             nTrials = cellfun(@(seg) seg.nTrials, varargin);
%             assert(all(nTrials == nTrials(1)), 'nTrials do not match');
%             nTrials = nTrials(1);
%
%             nUnits = cellfun(@(seg) seg.nUnits, varargin);
%             assert(all(nUnits == nUnits(1)), 'nUnits do not match');
%             nUnits = nUnits(1);
%
%             dsmerge = Kilosort.TrialSegmentedDataset.emptyWithSize(nTrials, nUnits);
%             dsmerge.dataset = varargin{1}.dataset;
%             nMerge = numel(varargin);
%
%             occmat = nan(nTrials, nMerge);
%             for iM = 1:nMerge
%                 occmat(:, iM) = varargin{iM}.trial_has_data;
%             end
%             nOcc = sum(occmat, 2);
%             assert(max(nOcc) < 2, 'At least one trial has data in multiple segmented datasets');
%
%             for iM = 1:nMerge
%                 seg = varargin{iM};
%                 idx = occmat(:, iM);
%                 dsmerge.trial_ids(idx) = seg.trial_ids(idx);
%                 dsmerge.trial_has_data = seg.trial_has_data(idx);
%                 dsmerge.trial_start(idx) = seg.trial_start(idx);
%                 dsmerge.cluster_ids(idx) = seg.cluster_ids(idx);
%
%                 dsmerge.spike_idx(idx, :) = seg.spike_idx(idx, :);
%                 dsmerge.amplitudes(idx, :) = seg.amplitudes(idx, :);
%                 dsmerge.pc_features(idx, :) = pc_features(idx, :);
%                 dsmerge.spike_times_ms_rel_start(idx, :) = spike_times_rel_start(idx, :);
%                 dsmerge.template_features(idx, :) = template_features(idx, :);
%                 dsmerge.spike_clusters(idx, :) = spike_clusters(idx, :);
%             end
%         end
%     end

end
