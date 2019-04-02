classdef KiloSortTrialSegmentedDataset < handle & matlab.mixin.Copyable
    
    % Properties that are copied over but not segmented into trials
    properties
        dataset

        trial_ids(:, 1) uint32

        % nTrials x 1
        trial_has_data(:, 1) logical

        % nTrials x 1
        trial_start(:, 1) uint64

        % nTrials x 1
        trial_stop(:, 1) uint64

        % indices into master dataset: trials x clusters
        spike_idx(:, :) cell

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
        channel_ids % which channel_ids were used for sorting
        
        % nTrials x 1
        trial_duration_ms
        
        fsAP % sampling rate pass thru
    end

    % Properties that are segmented by trial
    % each of these is nTrials x nTemplates cell with the same inner size as in Dataset (as described)
    properties
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
    end

    methods
        function seg = KiloSortTrialSegmentedDataset(dataset, tsi, trial_ids, varargin)
            p = inputParser();
            p.addParameter('loadFeatures', true, @islogical);
            p.addParameter('loadSync', true, @islogical);
            p.parse(varargin{:});
            
            % trial_ids specifies the id of each trial that will appear in
            % this segmented dataset. tsi has its own trialId field, and
            % the data will be copied over where these ids match. But the
            % ultimate size will be set according to trial_ids
            dataset.load();
            
            seg.dataset = dataset;
            seg.trial_ids = trial_ids;

            % trials here are over the requested trial_ids, which are different
            % from those in trialInfo (which comes from the sync line and is limited
            % to the trials in the neuropixel file)
            nTrials = numel(trial_ids);
            
            % filter trialInfo included in trial_ids
            tsi_trial_ids = tsi.trialId;
            tsi_start_idx = tsi.idxStart;
            tsi_stop_idx = tsi.idxStop;

            [trial_info_included, trial_idx_each_trial_info] = ismember(tsi_trial_ids, trial_ids);
            if ~any(trial_info_included)
                warning('No trial ids found in trialInfo were included');
                return;
            end
            tsi_trial_ids = tsi_trial_ids(trial_info_included);
            tsi_start_idx = tsi_start_idx(trial_info_included);
            tsi_stop_idx = tsi_stop_idx(trial_info_included);
            trial_idx_each_trial_info = trial_idx_each_trial_info(trial_info_included);

            seg.trial_start = nan(nTrials, 1);
            seg.trial_start(trial_idx_each_trial_info) = tsi_start_idx;

            seg.trial_stop = nan(nTrials, 1);
            seg.trial_stop(trial_idx_each_trial_info) = tsi_stop_idx;

            seg.trial_has_data = false(nTrials, 1);
            seg.trial_has_data(trial_idx_each_trial_info) = true;

            % lookup from columns of all the nTrials x nUnits back into template ids
            seg.cluster_ids = dataset.cluster_ids;
            seg.cluster_groups = dataset.cluster_groups;
            
            nUnits = numel(seg.cluster_ids);

            % figure out which trial each spike in spike_times belongs
            % to do discard data that occurs after trial stop, here we assume that each
            % trial ends at the next's start and the last trial ends at EOF
            edges = [tsi_start_idx; tsi_stop_idx(end)];
            trial_info_idx_each_spike = discretize(seg.dataset.spike_times, edges);
            trial_info_trial_id_each_spike = Neuropixel.Utils.TensorUtils.selectAlongDimensionWithNaNs(tsi_trial_ids, 1, trial_info_idx_each_spike);

            % convert trial info idx into trial_ids idx
            [~, trial_idx_each_spike] = ismember(trial_info_trial_id_each_spike, trial_ids);

            % which cluster does each spike belong to
            [~, unit_idx_each_spike] = ismember(seg.dataset.spike_clusters, seg.cluster_ids);

            subs = [trial_idx_each_spike, unit_idx_each_spike];

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
                dataset.spike_times, 1, [nTrials, nUnits], subs);

            seg.spike_times = spike_times_grouped;
            
            % convert samples to ms relative to trial_start
            for iT = 1:nTrials
                if ~seg.trial_has_data(iT), continue, end
                for iU = 1:nUnits
                    if isempty(spike_times_grouped{iT, iU})
                        spike_times_grouped{iT, iU} = nan(0, 1, 'single');
                    else
                        spike_times_grouped{iT, iU} = single(spike_times_grouped{iT, iU} - seg.trial_start(iT)) / single(seg.dataset.sample_rate / 1000);
                    end
                end
            end
            seg.spike_times_ms_rel_start = spike_times_grouped;

            prog.increment('Segmenting trials: spike_idx');
            seg.spike_idx = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                (1:seg.dataset.nSpikes)', 1, [nTrials, nUnits], subs);

            if p.Results.loadFeatures
                % slice the other fields into trials x unit:
                prog.increment('Segmenting trials: amplitudes');
                seg.amplitudes = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    dataset.amplitudes, 1, [nTrials, nUnits], subs);

                prog.increment('Segmenting trials: pc_features');
                seg.pc_features = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    dataset.amplitudes, 1, [nTrials, nUnits], subs);

                prog.increment('Segmenting trials: template_features');
                seg.template_features = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    dataset.template_features, 1, [nTrials, nUnits], subs);
            end
            
            prog.increment('Segmenting trials: spike_clusters');
            seg.spike_templates = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                dataset.spike_templates, 1, [nTrials, nUnits], subs);

            seg.syncBitNames = dataset.syncBitNames;
            
            
            if p.Results.loadSync
                prog.increment('Segmenting trials: sync');
                sync = dataset.readSync();
                trial_info_idx_each_sample = discretize(1:numel(sync), edges);
                trial_info_trial_id_each_sample = Neuropixel.Utils.TensorUtils.selectAlongDimensionWithNaNs(tsi_trial_ids, 1, trial_info_idx_each_sample);
                % convert trial info idx into trial_ids idx
                [~, trial_idx_each_sample] = ismember(trial_info_trial_id_each_sample, trial_ids);
                seg.sync = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(...
                    sync, 1, nTrials, trial_idx_each_sample);
            else
                seg.sync = {};
            end
            
            prog.finish();
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
            n = numel(seg.channel_ids);
        end
        
        function n = get.channel_ids(seg)
            n = seg.dataset.channel_ids;
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
    end

    % pulling things from raw data
    methods
        function snippet_set = getWaveformsFromRawData(seg, varargin)
            % mask_cell is nTrials x nClusters cell of mask or indices over
            % spike times within that bin, if omitted, all spikes will be
            % included
            
            p = inputParser();
            % from KiloSortDataset.readAPSnippets
            p.addOptional('mask_cell', {}, @(x) isempty(x) || iscell(x));
            p.addParameter('cluster_ids', seg.cluster_ids, @isvector); % which clusters each column of mask_cell correspond to
            p.addParameter('channel_ids_by_cluster', [], @(x) isempty(x) || ismatrix(x)); % specify a subset of channels to extract
            p.addParameter('best_n_channels', NaN, @isscalar); % or take the best n channels based on this clusters template when cluster_id is scalar 

            % other params:
            p.addParameter('window', [-40 41], @isvector); % Number of samples before and after spiketime to include in waveform
            p.addParameter('car', false, @islogical);
            p.addParameter('centerUsingFirstSamples', 20, @(x) isscalar(x) || islogical(x)); % subtract mean of each waveform's first n samples, don't do if false
            p.addParameter('num_waveforms', Inf, @isscalar); % caution: Inf will request ALL waveforms in order (typically useful if spike_times directly specified)
             
            p.addParameter('subtractOtherClusters', false, @islogical); % time consuming step to remove the contribution of the other clusters to a given snippet
            
            p.addParameter('raw_dataset', seg.raw_dataset, @(x) true); 
            
            p.parse(varargin{:});
            
            % lookup cluster_ids from
            cluster_ids = p.Results.cluster_ids;
            [tf, cluster_ind] = ismember(cluster_ids, seg.cluster_ids);
            if any(~tf)
                error('Some cluster_ids were not found in dataset');
            end

            nClu = numel(cluster_ids);
            
            % specify all spikes of clusters if no mask provided
            mask_cell = p.Results.mask_cell;
            if isempty(mask_cell)
                mask_cell = cellfun(@(times) true(size(times)), seg.spike_times(:, cluster_ind), 'UniformOutput', false);
            end

            % check size
            assert(size(mask_cell, 1) == seg.nTrials, 'mask_cell must be nTrials along dim 1');
            assert(size(mask_cell, 2) == nClu, 'mask_cell must be nClusters along dim 2');

%             % assemble spike_times we want to collect, all at once
            masked_spike_idx = cellfun(@(spike_idx, mask) spike_idx(mask), seg.spike_idx(:, cluster_ind), mask_cell, 'UniformOutput', false);
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
            
    end

%     methods(Static)
%         function seg = emptyWithSize(nTrials, nUnits)
%             seg = KiloSort.TrialSegmentedDataset;
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
%             seg.spike_clusters = c;
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
%             dsmerge = KiloSort.TrialSegmentedDataset.emptyWithSize(nTrials, nUnits);
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
