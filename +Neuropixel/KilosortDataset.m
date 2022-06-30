classdef KilosortDataset < handle & matlab.mixin.Copyable
    % wrapper around a Kilosort dataset
    % todo - load cluster ratings from cluster_groups.tsv
    % Note 1: in the context of this file, time refers to samples, 1-indexed
    % Note 2: this file will differ from raw_dataset in nChannelsSorted. Here, nChannelsSorted means the number of channels
    %   in .channel_ids_sorted (which will match the other properties)

    properties
        path(1, :) char

        dataset_variant char;

        raw_dataset % Neuropixel.ImecDataset instance

        channelMap % Neuropixel.ChannelMap

        fsAP % sampling rate pass thru to raw_dataset or specified during construction

        apGain
        apScaleToUv % multiply raw samples by this to get uV

        meta % ap metadata loaded

        concatenationInfo

        isLoaded logical = false;
        isLoadedBatchwise logical = false;
        isLoadedFeatures logical = false;
        isLoadedCutoff logical = false;
        isLoadedPreSplit logical = false;

        trimToNumSamples logical = true;
        modifiedInMemory logical = false;
    end

    % Computed properties
    properties(Hidden)
        metrics
    end

    properties(Dependent)
        pathLeaf
        hasRawDataset
        hasFeaturesLoaded
        hasCutoffLoaded
        hasBatchwiseLoaded
        hasPreSplitLoaded

        nChannelsSorted
        nSpikes
        nClusters
        nTemplates
        nTemplatesPreSplit
        nTemplateRank
        nPCFeatures
        nFeaturesPerChannel
        nTemplateTimepoints
        nTemplateChannels
        templateTimeRelative

        nBatches
        nSpikesCutoff
        
        isLoadedAll
    end

    properties(Constant)
        CUTOFF_SPIKES_EXPORT_OFFSET = 10000;
    end

    properties(SetAccess=protected)
        syncBitNames(:, 1) string
    end

    properties % Primary data properties
        % dat_path - location of raw data file as specified in params.py
        dat_path(1, :) char

        % n_channels_dat - total number of rows in the data file (not just those that have your neural data on them. This is for loading the file)
        n_channels_dat(1, 1) uint64

        % dtype - data type to read, e.g. 'int16'
        dtype(1, :) char

        % offset - number of bytes at the beginning of the file to skip
        offset uint64

        % sample_rate - in Hz
        sample_rate double

        % hp_filtered - True/False, whether the data have already been filtered
        hp_filtered(1,1) logical

        % amplitudes.npy - [nSpikes, ] double vector with the amplitude scaling factor that was applied to the template when extracting that spike
        amplitudes(:,1) double;

        % channel_map.npy - [nChannelsSorted, ] uint32 vector with the channel map, i.e. which row of the data file to look in for the channel in question
        channel_ids_sorted(:, 1) uint32

        % channel_positions.npy - [nChannelsSorted, 2] double matrix with each row giving the x and y coordinates of that channel. Together with the channel map, this determines how waveforms will be plotted in WaveformView (see below).
        channel_positions_sorted(:, :) double

        % pc_features.npy - [nSpikes, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike.
        % The channels that those features came from are specified in pc_features_ind.npy. E.g. the value at pc_features[123, 1, 5]
        % is the projection of the 123rd spike onto the 1st PC on the channel given by pc_feature_ind[5].
        pc_features(:, :, :) single

        % pc_feature_ind.npy - [nTemplates, nPCFeatures] uint32 matrix specifying which channels contribute to each entry in dim 3 of the pc_features matrix
        % e.g. if pc_feature(iT
        pc_feature_ind(:, :) uint32

        % similar_templates.npy - [nTemplates, nTemplates] single matrix giving the similarity score (larger is more similar) between each pair of templates
        similar_templates(:, :) single

        % spike_templates.npy - [nSpikes, ] uint32 vector specifying the identity of the template that was used to extract each spike
        spike_templates(:, 1) uint32

        % spike_templates.npy - [nSpikes, ] uint32 vector specifying the identity of the template that was originally used to extract each spike, before splitAllClusters
        spike_templates_preSplit(:, 1) uint32

        % spike_times.npy - [nSpikes, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        spike_times(:, 1) uint64

        % template_features.npy - [nSpikes, nTemplateRank] single matrix giving the magnitude of the projection of each spike onto nTemplateRank other features.
        % Which other features is specified in template_feature_ind.npy
        template_features(:, :) single

        % template_feature_ind.npy - [nTemplates, nTemplateRank] uint32 matrix specifying which templateFeatures are included in the template_features matrix.
        template_feature_ind(:, :) uint32

        template_sample_offset (1,1) uint64 % scalar value indicating what index into templates (time dim 2) the spike time occurs at

        % templates.npy - [nTemplates, nTemplateTimepoints, nTemplateChannels] single matrix giving the template shapes on the channels given in templates_ind.npy
        templates(:, :, :) single

        % templates_ind.npy - [nTemplates, nTempChannels] double matrix specifying the channels on which each template is defined.
        % In the case of Kilosort templates_ind is just the integers from 0 to nChannelsSorted-1, since templates are defined on all channels.
        templates_ind(:, :) double

        % whitening_mat.npy - [nChannelsSorted, nChannelsSorted] double whitening matrix applied to the data during automatic spike sorting
        whitening_mat(:, :) double

        % whitening_mat_inv.npy - [nChannelsSorted, nChannelsSorted] double, the inverse of the whitening matrix.
        whitening_mat_inv(:, :) double

        % spike_clusters.npy - [nSpikes, ] uint32 vector giving the cluster identity of each spike. This file is optional and
        % if not provided will be automatically created the first time you run the template gui, taking the same values as
        % spike_templates.npy until you do any merging or splitting.
        spike_clusters(:, 1) uint32

        % original of spike_clusters before being adjusted by Phy or
        spike_clusters_ks2orig(:, 1) uint32

        % cluster_groups - "cluster group" of each cluster as set in Phy
        cluster_groups(:, 1) categorical

        % unique clusters in spike_clusters [nClusters]
        cluster_ids (:, 1) uint32

        % has the list of spikes already been deduplicated
        is_deduplicated (1, 1) logical = false;
        deduplication_stats struct;

        deduplicate_spikes logical = false;
        deduplicate_cutoff_spikes logical = false;
        deduplicate_within_samples uint64 = 5;
        deduplicate_within_distance single = 50;
    end

    properties % Secondary data, loaded from rez.mat for Kilosort2 only
        ops struct % options used by Kilosort

        kilosort_version (1, :) uint32 % scalar or vector (leading digit is 1 or 2 for Kilosort1 or Kilosort2)

        % low-rank decomposition of templates into spatial and temporal components
        W (:, :, :) single % [nTemplateTimepoints, nTemplates, nTemplateRank] - temporal decomposition of templates
        U (:, :, :) single % [nChannelsSorted, nTemplates, nTemplateRank] - spatial decomposition of templates
        mu (:, 1) single % [nTemplates, 1] - original template scale factor determined in learnAndSolve8b

        W_preSplit (:, :, :) double % [nTemplateTimepoints, nTemplatesPreSplit, nTemplateRank] - temporal decomposition of templates
        U_preSplit (:, :, :) double % [nChannelsSorted, nTemplatesPreSplit, nTemplateRank] - spatial decomposition of templates
        mu_preSplit (:, 1) double % [nTemplatesPreSplit, 1] - original template scale factor determined in learnAndSolve8b
        iW_preSplit(:, 1) int32 % nTemplatesPreSplit - channel indices where the biggest deflection is (at time nt0min), used mainly for spike reextraction with the same settings

        % batch information used in drift correction
        batch_starts (:, 1) uint64 % [nBatches] timesteps where each batch begins (unsorted)

        batchwise_cc (:, :) single % (rez.ccb) [nBatches, nBatches] batch-wise cross-correlation in original data ordering

        batch_sort_order (:, 1) uint32 % (rez.iorig) [nBatches] sorted batch indices used in KS2

        % batchwise template data
        W_batch (:, :, :, :) single % [nTemplateTimepoints, nTemplates, nTemplateRank, nBatches] - full spatial templates by batch

        W_batch_preSplit (:, :, :, :) single % [nTemplateTimepoints, nTemplatesPreSplit, nTemplateRank, nBatches] - full spatial templates by batch, for original pre-split templates

        % per-template svd of W_batch U*S --> W_batch_a, V --> W_batch_b
        W_batch_US (:, :, :, :) single % [nTemplateTimepoints, nTemplateRank, nTemplatePCs, nTemplates] - reshaped from rez.W_a
        W_batch_V (:, :, :) single % [nBatches, nTemplatePCs, nTemplates] % rez.W_b

        U_batch (:, :, :, :) single  % [nChannelsSorted, nTemplates, nTemplateRank, nBatches] - full temporal templates by batch

        U_batch_preSplit (:, :, :, :) single % [nTemplateTimepoints, nTemplatesPreSplit, nTemplateRank, nBatches] - full temporal templates by batch, for original pre-split templates

        % per-template svd of U_batch U*S --> U_batch_a, V --> U_batch_b
        U_batch_US (:, :, :, :) single % [nChannelsSorted, nTemplateRank, nTemplatePCs, nTemplates]
        U_batch_V (:, :, :) single % [nBatches, nTemplatePCs, nTemplates]

        mu_batch (:, :) single % [nTemplates, nBatches] - original template scale factors determined in learnAndSolve8b

        mu_batch_preSplit (:, :) single % [nTemplatesPreSplit, nBatches] - original template scale factors determined in learnAndSolve8b

        % good vs. mua label as assigned by Kilosort2 in cluster_KSLabel.tsv (good, mua)
        cluster_ks_label(:, 1) categorical % [nClusters]

        cluster_est_contam_rate (:, 1) double % [nClusters]

        % optional split and merge meta information (only on djoshea branch of Kilosort2)
        % each is nTemplates x 1 corresponding and points to cluster_ids (meaning are 0 indexed)
        cluster_merge_count (:, 1) single
        cluster_merge_dst(:, 1) single % cluster_merge_dst(i) == j implies that spikes from template i, cluster i-1 were merged into cluster j-1
        cluster_split_src (:, 1) single
        cluster_split_dst (:, 1) cell
        cluster_split_auc (:, 1) cell
        cluster_split_candidate (:, 1) logical
        cluster_orig_template (:, 1) uint32
        cluster_split_projections(:, 1) cell

        % [nSpikesInvalid, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        cutoff_spike_times(:, 1) uint64

        cutoff_thresholds(:, 1) double % cutoffs for spike_amplitude decided in vexp, these are in RMS units

        cutoff_amplitudes(:, 1) double

        % [nSpikesCutoff, ] uint32 vector specifying the identity of the template that was used to extract each spike
        cutoff_spike_templates(:, 1) uint32
        cutoff_spike_templates_preSplit(:, 1) uint32
        cutoff_spike_clusters(:, 1) uint32
        cutoff_spike_clusters_ks2orig(:, 1) uint32

         % [nSpikesCutoff, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike (from .cProjPC_cutoff_invalid)
        cutoff_pc_features(:, :, :) uint32

        % [nTemplates, nTemplateTimepoints, nTemplateChannels] single matrix giving the template shapes on the channels given in templates_ind (from .cProj_cutoff_invalid)
        cutoff_template_features(:, :) single
    end

    properties(Dependent)
        clusters_good
        clusters_mua
        clusters_noise
        clusters_unsorted

        cluster_spike_counts % (:, 1) uint32 % nClusters x 1 number of spikes assigned to each cluster
        cutoff_cluster_spike_counts % (:, 1) uint32 % nClusters x 1 number of spikes assigned to each cluster

        cluster_ids_have_spikes % list of cluster_ids with spikes (typically useful after spikes have been masked)
    end

    methods % Dependent properties
        function leaf = get.pathLeaf(ks)
            [~, leaf] = fileparts(ks.path);
        end

        function tf = get.hasRawDataset(ks)
            tf = ~isempty(ks.raw_dataset);
        end

        function tf = get.hasFeaturesLoaded(ks)
            tf = ~isempty(ks.pc_features);
        end

        function tf = get.hasCutoffLoaded(ks)
            tf = ~isempty(ks.cutoff_spike_clusters);
        end

        function tf = get.hasBatchwiseLoaded(ks)
            tf = ~isempty(ks.W_batch);
        end

        function tf = get.hasPreSplitLoaded(ks)
            tf = ~isempty(ks.W_preSplit);
        end

        function n = get.nSpikes(ks)
            if isempty(ks.spike_times)
                n = NaN;
            else
                n = size(ks.spike_times, 1);
            end
        end

        function n = get.nSpikesCutoff(ks)
            if isempty(ks.cutoff_spike_times)
                n = 0;
            else
                n = size(ks.cutoff_spike_times, 1);
            end
        end

        function n = get.nClusters(ks)
            if isempty(ks.cluster_ids)
                n = NaN;
            else
                n = numel(ks.cluster_ids);
            end
        end

        function n = get.nTemplates(ks)
            if isempty(ks.templates)
                n = NaN;
            else
                n = size(ks.templates, 1);
            end
        end

        function n = get.nTemplatesPreSplit(ks)
            if isempty(ks.W_preSplit)
                if isempty(ks.spike_templates_preSplit)
                    n = NaN;
                else
                    n = max(ks.spike_templates_preSplit);
                end
            else
                n = size(ks.W_preSplit, 2);
            end
        end

        function n = get.nTemplateRank(ks)
            if isempty(ks.template_features)
                n = NaN;
            else
                n = size(ks.template_features, 2);
            end
        end

        function n = get.nPCFeatures(ks)
            if isempty(ks.pc_features)
                n = NaN;
            else
                n = size(ks.pc_features, 3);
            end
        end

        function n = get.nFeaturesPerChannel(ks)
            if isempty(ks.pc_features)
                n = NaN;
            else
                n = size(ks.pc_features, 2);
            end
        end

        function n = get.nTemplateTimepoints(ks)
            if isempty(ks.templates)
                n = NaN;
            else
                n = size(ks.templates, 2);
            end
        end

        function n = get.nTemplateChannels(ks)
            if isempty(ks.templates)
                n = NaN;
            else
                n = size(ks.templates, 3);
            end
        end

        function tvec = get.templateTimeRelative(ks)
            T = int64(ks.nTemplateTimepoints);
            off = int64(ks.template_sample_offset);
            start = -off + int64(1);
            stop = start + T - int64(1);
            tvec = (start:stop)';
        end

        function n = get.nChannelsSorted(ks)
            if isempty(ks.channel_ids_sorted)
                n = NaN;
            else
                n = numel(ks.channel_ids_sorted);
            end
        end

        function map = get.channelMap(ks)
            if ~isempty(ks.raw_dataset) && ~isempty(ks.raw_dataset.channelMap)
                % pass thru to raw dataset
                map = ks.raw_dataset.channelMap;
            else
                map = ks.channelMap;
            end
        end

        function nBatches = get.nBatches(ks)
            nBatches = numel(ks.batch_sort_order);
        end

        function fsAP = get.fsAP(ks)
            if isempty(ks.fsAP)
                if isempty(ks.raw_dataset)
                    fsAP = NaN;
                else
                    fsAP = ks.raw_dataset.fsAP;
                end
            else
                fsAP = ks.fsAP;
            end
        end

        function apGain = get.apGain(ks)
            if isempty(ks.apGain)
                if isempty(ks.raw_dataset)
                    apGain = NaN;
                else
                    apGain = ks.raw_dataset.apGain;
                end
            else
                apGain = ks.apGain;
            end
        end

        function apScaleToUv = get.apScaleToUv(ks)
            if isempty(ks.apScaleToUv)
                if isempty(ks.raw_dataset)
                    apScaleToUv = NaN;
                else
                    apScaleToUv = ks.raw_dataset.apScaleToUv;
                end
            else
                apScaleToUv = ks.apScaleToUv;
            end
        end

        function syncBitNames = get.syncBitNames(ks)
            if isempty(ks.raw_dataset)
                if isempty(ks.syncBitNames)
                    syncBitNames = strings(16, 1);
                else
                    syncBitNames = string(ks.syncBitNames);
                end
            else
                syncBitNames = string(ks.raw_dataset.syncBitNames);
            end
        end

        function c = get.clusters_good(ks)
            c = ks.cluster_ids(ks.cluster_groups == "good");
        end

        function c = get.clusters_mua(ks)
            c = ks.cluster_ids(ks.cluster_groups == "mua");
        end

        function c = get.clusters_noise(ks)
            c = ks.cluster_ids(ks.cluster_groups == "noise");
        end

        function c = get.clusters_unsorted(ks)
            c = ks.cluster_ids(ks.cluster_groups == "unsorted");
        end

        function n = get.cluster_spike_counts(ks)
            if isempty(ks.cluster_ids)
                n = zeros(0, 1, 'uint32');
            else
                n = uint32(Neuropixel.Utils.discrete_histcounts(ks.spike_clusters, ks.cluster_ids));
            end
        end

        function n = get.cutoff_cluster_spike_counts(ks)
            if isempty(ks.cluster_ids)
                n = zeros(0, 1, 'uint32');
            else
                n = uint32(Neuropixel.Utils.discrete_histcounts(ks.cutoff_spike_clusters, ks.cluster_ids));
            end
        end

        function ids = get.cluster_ids_have_spikes(ks)
            mask = ks.cluster_spike_counts > 0 | ks.cutoff_cluster_spike_counts > 0;
            ids = ks.cluster_ids(mask);
        end
    end

    methods
        function ks = KilosortDataset(path, varargin)
            if nargin == 0
                return;
            end

            p = inputParser();
            p.addParameter('channelMap', [], @(x) true);
            p.addParameter('imecDataset', [], @(x) true);
            p.addParameter('fsAP', [], @(x) isempty(x) || isscalar(x));
            p.addParameter('apScaleToUv', [], @(x) isempty(x) || isscalar(x));

            p.addParameter('trimToNumSamples', true, @islogical); % a bug in my modified Kilosort code

            % these will affect how load() will be performed later on, but if no saved deduplicated spikes are found we'll do it now
            p.addParameter('deduplicate_spikes', true, @islogical);
            p.addParameter('deduplicate_cutoff_spikes', true, @islogical);
            p.addParameter('deduplicate_within_samples', 5, @isscalar);
            p.addParameter('deduplicate_within_distance', 50, @isscalar);

            p.parse(varargin{:});

            if isa(path, 'Neuropixel.ImecDataset')
                ks.raw_dataset = path;
                path = ks.raw_dataset.pathRoot;
            end

            assert(exist(path, 'dir') == 7, 'Path %s not found', path);
            ks.path = path;

            if isempty(ks.raw_dataset)
                if isa(p.Results.imecDataset, 'Neuropixel.ImecDataset')
                    ks.raw_dataset = p.Results.imecDataset;
                elseif ischar(p.Results.imecDataset) || isempty(p.Results.imecDataset)
                    raw_path = p.Results.imecDataset;
                    if isempty(raw_path), raw_path = path; end
                    if Neuropixel.ImecDataset.folderContainsDataset(raw_path)
                        try
                            ks.raw_dataset = Neuropixel.ImecDataset(raw_path, 'channelMap', p.Results.channelMap);
                        catch exc
                            warning('KilosortDataset could not be loaded with imec_cleaned, exception was: %s', exc.message());
                        end
                    end
                end
            end

            if isempty(ks.raw_dataset)
                warning('No ImecDataset found in Kilosort path, specify imecDataset parameter directly');
            end

            % these will pass thru to raw_dataset if provided
            ks.fsAP = p.Results.fsAP;
            ks.apScaleToUv = p.Results.apScaleToUv;
            if isempty(ks.fsAP) || isnan(ks.fsAP) || isempty(ks.apScaleToUv) || isnan(ks.apScaleToUv)
                % will pass thru to raw_dataset if found
                % try reading sample_rate from params.py
                ks.readParamsPy();
                ks.fsAP = ks.sample_rate;
            end

            % manually specify some additional props
            channelMap = p.Results.channelMap;
            if isempty(channelMap)
                if ~isempty(ks.raw_dataset)
                    channelMap = ks.raw_dataset.channelMap;
                end
                if isempty(channelMap)
                    % try default paths with ks.path
                    channelMap = Neuropixel.Utils.searchForChannelMapInDirectory(ks.path);
                    if channelMap == ""
                        channelMap = Neuropixel.Utils.getDefaultChannelMapFile(true);
                    end
                end
            end

            if ischar(channelMap) || isstring(channelMap)
                channelMap = Neuropixel.ChannelMap(channelMap);
            end
            ks.channelMap = channelMap;

            if ~isempty(ks.raw_dataset)
                ks.meta = ks.raw_dataset.readAPMeta();
                ks.concatenationInfo = ks.raw_dataset.concatenationInfoAP;
            end

            ks.trimToNumSamples = p.Results.trimToNumSamples;
            ks.deduplicate_spikes = p.Results.deduplicate_spikes;
            ks.deduplicate_cutoff_spikes = p.Results.deduplicate_cutoff_spikes;
            ks.deduplicate_within_samples = p.Results.deduplicate_within_samples;
            ks.deduplicate_within_distance = p.Results.deduplicate_within_distance;
        end

        function s = computeBasicStats(ks, varargin)
            % a way of computing basic statistics without fully loading
            % from disk. see printBasicStats

            p = inputParser();
            p.addParameter('frThresholdHz', 3, @isscalar);
            p.parse(varargin{:});

            % can work without loading raw data, much faster
            if ks.isLoaded
                s.spike_clusters = ks.spike_clusters; %#ok<*PROP>
                s.cluster_ids = ks.cluster_ids;
                s.offset = ks.offset;
                s.sample_rate = ks.sample_rate;
                s.spike_times = ks.spike_times;

                if ks.hasCutoffLoaded
                    s.cutoff_spike_times = ks.cutoff_spike_times;
                    s.cutoff_spike_clusters = ks.cutoff_spike_clusters;
                else
                    s.cutoff_spike_times = [];
                    s.cutoff_spike_clusters = [];
                end
            else
                % partial load
                s.spike_clusters = read('spike_clusters');
                mask = s.spike_clusters < ks.CUTOFF_SPIKES_EXPORT_OFFSET; % eliminate cutoff spikes if included
                s.spike_clusters = s.spike_clusters(mask);
                s.cluster_ids = unique(s.spike_clusters);
                params = Neuropixel.readINI(fullfile(ks.path, 'params.py'));
                s.sample_rate = params.sample_rate;
                s.offset = params.offset;
                s.spike_times = read('spike_times');
                s.spike_times = s.spike_times(mask);

                if existp('cutoff_spike_times')
                    s.cutoff_spike_times = read('cutoff_spike_times');
                    s.cutoff_spike_clusters = read('cutoff_spike_clusters');
                else
                    s.cutoff_spike_times = [];
                    s.cutoff_spike_clusters = [];
                end
            end

            s.nSpikes = numel(s.spike_clusters);
            s.nClusters = numel(s.cluster_ids);
            s.nSec = double(max(s.spike_times) - s.offset) / double(s.sample_rate);

            counts = histcounts(s.spike_clusters, sort(s.cluster_ids));
            if ~isempty(s.cutoff_spike_clusters)
                cutoff_counts = histcounts(s.cutoff_spike_clusters, sort(s.cluster_ids));
            else
                cutoff_counts = zeros(size(counts));
            end
            s.fr = counts / double(s.nSec);
            s.fr_cutoff_only = cutoff_counts / double(s.nSec);
            s.fr_plus_cutoff = s.fr + s.fr_cutoff_only;
            s.thresh = p.Results.frThresholdHz;

            s.clusterMask = s.fr > s.thresh;
            s.nClustersAboveThresh = nnz(s.clusterMask);
            s.nSpikesAboveThresh = nnz(ismember(s.spike_clusters, s.cluster_ids(s.clusterMask)));

            function tf = existp(file)
                tf = exist(fullfile(ks.path, [file '.npy']), 'file') == 2;
            end

            function out = read(file)
                out = Neuropixel.readNPY(fullfile(ks.path, [file '.npy']));
            end
        end

        function s = printBasicStats(ks, varargin)
            p = inputParser();
            p.addParameter('color', 'k', @(x) true);
            p.KeepUnmatched = true;
            p.parse(varargin{:})
            s = ks.computeBasicStats(p.Unmatched);

            fprintf('%s: %.1f sec, %d (%d) spikes, %d (%d) clusters (with fr > %g Hz)\n', ks.pathLeaf, s.nSec, ...
                s.nSpikes, s.nSpikesAboveThresh, ...
                s.nClusters, s.nClustersAboveThresh, s.thresh);

            h = plot(sort(s.fr, 'descend'), '-', 'Color', p.Results.color);
            Neuropixel.Utils.showFirstInLegend(h, ks.pathLeaf);
            if any(s.fr_cutoff_only)
                hold on
                h = plot(sort(s.fr_plus_cutoff, 'descend'), '--', 'Color', p.Results.color);
                Neuropixel.Utils.hideInLegend(h);
            end
            xlabel('Cluster');
            ylabel('# spikes');
            hold on
            h = yline(s.thresh, 'Color', 'r');
            Neuropixel.Utils.hideInLegend(h);
            hold off;
            box off;
            grid on;
        end

        function readParamsPy(ks)
            params = Neuropixel.readINI(fullfile(ks.path, 'params.py'));
            ks.dat_path = params.dat_path;
            ks.n_channels_dat = params.n_channels_dat;
            ks.dtype = params.dtype;
            ks.offset = params.offset;
            ks.sample_rate = params.sample_rate;
            ks.hp_filtered = params.hp_filtered;
            if isfield(params, 'scale_to_uv')
                ks.apScaleToUv = params.scale_to_uv;
            end
        end
        
        function tf = get.isLoadedAll(ks)
            tf = ks.isLoadedBatchwise && ks.isLoadedFeatures && ks.isLoadedCutoff && ks.isLoadedPreSplit;
        end

        function load(ks, varargin)
            p = inputParser();
            p.addParameter('reload', false, @islogical);
            p.addParameter('loadBatchwise', true, @islogical);
            p.addParameter('loadFeatures', true, @islogical);
            p.addParameter('loadCutoff', true, @islogical);
            p.addParameter('loadPreSplit', true, @islogical);
            p.addParameter('progressInitializeFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(nUpdates) to print update
            p.addParameter('progressIncrementFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(updateString) to print update
            p.parse(varargin{:});

            reload = p.Results.reload;

            loadFeatures = p.Results.loadFeatures;
            loadBatchwise = p.Results.loadBatchwise;
            loadCutoff = p.Results.loadCutoff;
            loadPreSplit = p.Results.loadPreSplit;

            if ks.isLoaded && ~reload
                % check whether we're missing any requested things to load
                if (~p.Results.loadBatchwise || ks.isLoadedBatchwise) && ...
                   (~p.Results.loadFeatures || ks.isLoadedFeatures) && ...
                   (~p.Results.loadCutoff || ks.isLoadedCutoff) && ...
                   (~p.Results.loadPreSplit || ks.isLoadedPreSplit)
                    return;
                end
            end

            if ks.modifiedInMemory
                error('Cannot reload aspects of KilosortDataset after in-memory modifications have been made');
            end

            ks.isLoaded = false;

            ks.readParamsPy();

            path = ks.path;
            existp = @(file) exist(fullfile(path, file), 'file') > 0;

            nProg = 15;
            has_tsv = existp('cluster_KSLabel.tsv');
            if has_tsv
                ks.kilosort_version = 2;
                nProg = nProg + 24;
            else
                ks.kilosort_version = 1;
            end

            initStr = sprintf('Loading Kilosort %d dataset', ks.kilosort_version);
            if isempty(p.Results.progressInitializeFn) && isempty(p.Results.progressIncrementFn)
                prog = Neuropixel.Utils.ProgressBar(nProg, initStr);
                progIncrFn = @(text) prog.increment(text);
            else
                if ~isempty(p.Results.progressInitializeFn)
                    p.Results.progressInitializeFn(nProg, initStr);
                end
                progIncrFn = p.Results.progressIncrementFn;
            end

            if existp('ops.mat')
                progIncrFn('Loading ops.mat');
                ld = load(fullfile(path, 'ops.mat'), 'ops');
                ks.ops = ld.ops;

                % replace absolute paths in case ks.path has changed
                ks.ops.root = ks.path;
                ks.ops.saveDir = ks.path;
                ks.ops.fbinary = fullfile(ks.path, ks.dat_path);
            end

            ks.amplitudes = read('amplitudes');
            ks.channel_ids_sorted = read('channel_map');
            ks.channel_ids_sorted = ks.channel_ids_sorted + ones(1, 'like', ks.channel_ids_sorted); % make channel map 1 indexed
            ks.channel_positions_sorted = read('channel_positions');
            if loadFeatures
                ks.pc_features = read('pc_features');
                ks.pc_feature_ind = read('pc_feature_ind');
                ks.pc_feature_ind = ks.pc_feature_ind + ones(1, 'like', ks.pc_feature_ind); % convert plus zero indexed to 1 indexed channels
            end
            ks.similar_templates = read('similar_templates');
            ks.spike_templates = read('spike_templates');
            ks.spike_templates = ks.spike_templates + ones(1, 'like', ks.spike_templates);
            if existp('spike_templates_preSplit.npy')
                ks.spike_templates_preSplit = read('spike_templates_preSplit');
                ks.spike_templates_preSplit = ks.spike_templates_preSplit + ones(1, 'like', ks.spike_templates_preSplit);
            end
            ks.spike_times = read('spike_times');
            if loadFeatures
                ks.template_features = read('template_features');
                ks.template_feature_ind = read('template_feature_ind');
                ks.template_feature_ind = ks.template_feature_ind + ones(1, 'like', ks.template_feature_ind); % 0 indexed templates to 1 indexed templates
            end

            ks.templates = read('templates');
            T = size(ks.templates, 2);
            ks.template_sample_offset = uint64(floor(T/2)); % phy templates are centered on the spike time, with one extra sample post (e.g. -41 : 41)

            ks.templates_ind = read('templates_ind');
            ks.templates_ind = ks.templates_ind + ones(1, 'like', ks.templates_ind); % convert plus zero indexed to 1 indexed channels

            ks.whitening_mat = read('whitening_mat');
            ks.whitening_mat_inv = read('whitening_mat_inv');

            % ensure that spike_clusters wont change if load_cutoff is false
            ks.spike_clusters = read('spike_clusters');
            cutoff_spike_clusters = readOr('cutoff_spike_clusters');

            unique_cluster_ids = unique(cat(1, ks.spike_clusters, cutoff_spike_clusters));
            if existp('unique_cluster_ids.npy')
                % if unique_cluster_ids is hardcoded
                ks.cluster_ids = read('unique_cluster_ids');
                % check that the loaded cluster_ids is a superset of the cluster_ids seen in ks.spike_clusters
                assert(all(ismember(unique_cluster_ids, ks.cluster_ids)), ...
                    'cluster_ids list loaded from unique_cluster_ids.npy file must be a superset of all cluster ids encountered in spike_clusters.npy');
            else
                ks.cluster_ids = unique_cluster_ids;
            end

            if loadCutoff
                ks.cutoff_spike_clusters = cutoff_spike_clusters;
            end

            assert(~isempty(ks.cluster_ids));

            if existp('spike_clusters_ks2orig.npy')
                ks.spike_clusters_ks2orig = read('spike_clusters_ks2orig');
                if numel(ks.spike_clusters_ks2orig) ~= numel(ks.spike_times)
                    warning('Saved spike_clusters_ks2orig size does not match other fields, ignoring');
                    ks.spike_clusters_ks2orig = ks.spike_clusters;
                end
            else
                ks.spike_clusters_ks2orig = ks.spike_clusters;
            end

            % filter out clusters > 10000, typically used for cutoff clusters during export so that Phy can see them
            mask = ks.spike_clusters < ks.CUTOFF_SPIKES_EXPORT_OFFSET;
            ks.amplitudes = ks.amplitudes(mask);
            ks.spike_times = ks.spike_times(mask);
            ks.spike_templates = ks.spike_templates(mask);
            ks.spike_clusters = ks.spike_clusters(mask);
            ks.spike_clusters_ks2orig = ks.spike_clusters_ks2orig(mask);
            if loadFeatures
                ks.pc_features = ks.pc_features(mask, :, :);
                ks.template_features = ks.template_features(mask, :);
            end

            % load KS cluster labels
            ks.cluster_ks_label = repmat(categorical("unsorted"), numel(ks.cluster_ids), 1);
            if existp('cluster_KSLabel.tsv')
                tbl = readClusterMetaTSV('cluster_KSLabel.tsv', 'KSLabel', 'categorical');
                [tf, ind] = ismember(tbl.cluster_id, ks.cluster_ids);
                ks.cluster_ks_label(ind(tf)) = tbl{tf, 2};
            end

            % load cluster groups (mapping from cluster_ids to {good, mua, noise, unsorted})
            % important to lookup the cluster inds since we've uniquified spike_clusters
            fnameSearch = {'cluster_groups.csv', 'cluster_group.tsv'};
            ks.cluster_groups = repmat(categorical("unsorted"), numel(ks.cluster_ids), 1);
            for iF = 1:numel(fnameSearch)
                file = fnameSearch{iF};
                if existp(file)
                    tbl = readClusterMetaTSV(file, 'group', 'categorical');

                    [tf, ind] = ismember(tbl.cluster_id, ks.cluster_ids);
                    ks.cluster_groups(ind(tf)) = tbl{tf, 2};
                end
            end

%             % load cluster ratings from disk
%             ks.cluster_ratings = repmat(categorical("unrated"), numel(ks.cluster_ids), 1);
%             if existp('cluster_Rating.tsv')
%                 tbl = readClusterMetaTSV('cluster_Rating.tsv', 'Rating', 'categorical');
%                 [tf, ind] = ismember(tbl.cluster_id, ks.cluster_ids);
%                 ks.cluster_rating(ind(tf)) = tbl{tf, 2};
%             end

            if ks.kilosort_version == 2
                if loadFeatures
                    ks.W = readOr('template_W');
                    ks.U = readOr('template_U');
                    ks.mu = readOr('template_mu');

                    if loadPreSplit
                        ks.W_preSplit = readOr('template_W_presplit');
                        ks.U_preSplit = readOr('template_U_presplit');
                        ks.mu_preSplit = readOr('template_mu_presplit');
                        ks.iW_preSplit = readOr('template_iW_presplit');
                    end
                end

                % strip leading zeros off of ks.templates based on size of W
                if ~isempty(ks.W)
                    nTimeTemp = size(ks.W, 1);
                    assert(nTimeTemp ~= 0);
                    nStrip = size(ks.templates, 2) - nTimeTemp;
                elseif isfield(ks.ops, 'nt0')
                    nTimeTemp = ks.ops.nt0;
                    assert(nTimeTemp ~= 0);
                    nStrip = size(ks.templates, 2) - nTimeTemp;
                else
                    W = readOr('template_W');
                    if isempty(W)
                        % first 21 tend to be zero in kilosort2
                        firstNonZero = find(any(ks.templates ~= 0, [1 3]), 1, 'first');
                        nStrip = firstNonZero - 1;
                    else
                        nTimeTemp = size(W, 1);
                        assert(nTimeTemp ~= 0);
                        nStrip = size(ks.templates, 2) - nTimeTemp;
                    end
                end

                if nStrip > 0
                    ks.templates = ks.templates(:, nStrip+1:end, :);
                    ks.template_sample_offset = ks.template_sample_offset - uint64(nStrip);
                end

                if loadBatchwise
                    ks.batchwise_cc = readOr('batchwise_ccb');
                    ks.batch_sort_order = readOr('batch_sort_order');
                    ks.batch_starts = readOr('batch_starts');

                    ks.W_batch = readOr('template_W_batch');

                    ks.W_batch_US = readOr('template_W_batch_US');
                    ks.W_batch_V = readOr('template_W_batch_V');
                    ks.U_batch = readOr('template_U_batch');
                    ks.U_batch_US = readOr('template_U_batch_US');
                    ks.U_batch_V = readOr('template_U_batch_V');
                    ks.mu_batch = readOr('template_mu_batch');

                    if loadPreSplit % mostly used for reextracting spikes
                        ks.W_batch_preSplit = readOr('template_W_batch_presplit');
                        ks.U_batch_preSplit = readOr('template_U_batch_presplit');
                        ks.mu_batch_preSplit = readOr('template_mu_batch_presplit');
                    end
                end

                ks.cluster_est_contam_rate = readOr('cluster_est_contam_rate');

                if existp('splitMergeInfo.mat')
                    progIncrFn('Loading splitMergeInfo.mat');
                    ld = load(fullfile(path, 'splitMergeInfo.mat'), 'splitMergeInfo');
                    splitMergeInfo = ld.splitMergeInfo;
                    ks.cluster_merge_count = splitMergeInfo.mergecount;
                    ks.cluster_merge_dst = splitMergeInfo.mergedst;
                    ks.cluster_split_src = splitMergeInfo.splitsrc;
                    ks.cluster_split_dst = splitMergeInfo.splitdst;
                    ks.cluster_split_auc = splitMergeInfo.splitauc;
                    ks.cluster_split_candidate = splitMergeInfo.split_candidate;
                    ks.cluster_orig_template = splitMergeInfo.split_orig_template;
                    ks.cluster_split_projections = splitMergeInfo.split_projections;
                end

                if loadCutoff
                    ks.cutoff_thresholds = readOr('cutoff_thresholds');

                    ks.cutoff_spike_times = readOr('cutoff_spike_times');
                    cutoff_temps = readOr('cutoff_spike_templates');
                    if ~isempty(cutoff_temps)
                        ks.cutoff_spike_templates = cast(cutoff_temps, 'like', ks.spike_templates) + ones(1, 'like', ks.spike_templates); % 0 indexed templates to 1 indexed templates
                        cutoff_temps_preSplit = readOr('cutoff_spike_templates_preSplit');
                        ks.cutoff_spike_templates_preSplit = cast(cutoff_temps_preSplit, 'like', ks.spike_templates) + ones(1, 'like', ks.spike_templates); % 0 indexed templates to 1 indexed templates
                    end

                    ks.cutoff_amplitudes = readOr('cutoff_amplitudes');
                    % now loaded above
%                     ks.cutoff_spike_clusters = readOr('cutoff_spike_clusters'); % rez.st3_cutoff is 1 indexed, cluster ids are 0 indexed
                    if existp('cutoff_spike_clusters_ks2orig.npy')
                        ks.cutoff_spike_clusters_ks2orig = read('cutoff_spike_clusters_ks2orig');
                        if numel(ks.cutoff_spike_clusters_ks2orig) ~= numel(ks.cutoff_spike_times)
                            warning('Saved cutoff_spike_clusters_ks2orig size does not match other fields, ignoring');
                            ks.cutoff_spike_clusters_ks2orig = ks.cutoff_spike_clusters;
                        end
                    else
                        ks.cutoff_spike_clusters_ks2orig = ks.cutoff_spike_clusters;
                    end

                    if loadFeatures
                        ks.cutoff_pc_features = readOr('cutoff_pc_features');
                        ks.cutoff_template_features = readOr('cutoff_template_features');
                    end
                end
            end

            ks.metrics = []; % old metrics no longer valid
            ks.isLoaded = true;
            if loadBatchwise
                ks.isLoadedBatchwise = true;
            end
            if loadFeatures
                ks.isLoadedFeatures = true;
            end
            if loadCutoff
                ks.isLoadedCutoff = true;
            end
            if loadPreSplit
                ks.isLoadedPreSplit = true;
            end

            if ks.deduplicate_spikes || ks.deduplicate_cutoff_spikes
                need_deduplication = true;
                if existp('spike_deduplication_mask.mat')
                    progIncrFn('Applying saved spike_deduplication_mask.mat');
                    temp = load(fullfile(path, 'spike_deduplication_mask.mat'), 'spike_deduplication_mask', 'cutoff_spike_deduplication_mask', 'deduplication_stats');
                    if isfield(temp, 'spike_deduplication_mask') && isfield(temp, 'deduplication_stats')
                        mask = temp.spike_deduplication_mask;
                        if isfield(temp, 'cutoff_spike_deduplication_mask')
                            cutoff_mask = temp.cutoff_spike_deduplication_mask;
                        else
                            cutoff_mask = false(0, 1);
                        end
                        ks.mask_spikes(mask, cutoff_mask);
                        ks.is_deduplicated = true;
                        ks.deduplication_stats = temp.deduplication_stats;

                        need_deduplication = false;
                    else
                        warning('Spike dedpulication file spike_deduplication_mask.mat missing required fields, redoing duplication')
                    end
                end

                if need_deduplication
                    progIncrFn('Deduplicating spikes');
                    ks.remove_duplicate_spikes();
                end
            else
                ks.is_deduplicated = false;
                ks.deduplication_stats = struct();
            end

            if ks.trimToNumSamples
                if ~isempty(ks.raw_dataset)
                    nSamples = ks.raw_dataset.nSamplesAP;
                elseif ~isempty(ks.ops) && isfield(ks.ops, 'tend')
                    nSamples = ks.ops.tend;
                else
                    nSamples = NaN;
                end

                if isnan(nSamples)
                    warning('Cannot trim spikes to nSamples, nSamples not found in raw_dataset or in ops');
                else
                    mask = ks.spike_times <= nSamples;
                    mask_cutoff = ks.cutoff_spike_times <= nSamples;
                    ks.mask_spikes(mask, mask_cutoff, 'setModifiedInMemory', false);
                end
            end

            if exist('prog', 'var')
                prog.finish();
            end

            function out = readOr(file, default)
                ffile = fullfile(path, [file '.npy']);
                if exist(ffile, 'file') > 0
                    out = read(file);
                elseif nargin > 2
                    out = default;
                else
                    out = [];
                end
            end

            function out = read(file)
                progIncrFn(sprintf('Loading %s', file));
                out = Neuropixel.readNPY(fullfile(path, [file '.npy']));
            end

            function tbl = readClusterMetaTSV(file, field, type)
                opts = delimitedTextImportOptions("NumVariables", 2);
                    % Specify range and delimiter
                    opts.DataLines = [2, Inf];
                    opts.Delimiter = "\t";

                    % Specify column names and types
                    opts.VariableNames = ["cluster_id", string(field)];
                    opts.VariableTypes = ["uint32", string(type)];
                    opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");
                    opts.ExtraColumnsRule = "ignore";
                    opts.EmptyLineRule = "read";

                    % Import the data
                    tbl = readtable(fullfile(path,file), opts);
            end

        end

        function create_spike_clusters_ks2orig_if_missing(ks)
            % copies the file spike_clusters.npy to spike_clusters_ks2orig.npy if the latter is nonexistent
            dest_path = fullfile(ks.path, 'spike_clusters_ks2orig.npy');
            if ~exist(dest_path, 'file')
                src_path = fullfile(ks.path, 'spike_clusters.npy');
                copyfile(src_path, dest_path);
            end

            dest_path = fullfile(ks.path, 'cutoff_spike_clusters_ks2orig.npy');
            if ~exist(dest_path, 'file')
                src_path = fullfile(ks.path, 'cutoff_spike_clusters.npy');
                if exist(src_path, 'file')
                    copyfile(src_path, dest_path);
                end
            end
        end

        function save_spike_clusters_to_disk(ks)
            writeNPY_local = @(v, fname) writeNPY(v, fullfile(ks.path, fname));
            writeNPY_local(ks.spike_clusters, 'spike_clusters.npy');
            writeNPY_local(ks.cutoff_spike_clusters, 'cutoff_spike_clusters.npy');
        end
        
        function markModifiedInMemory(ks)
            ks.modifiedInMemory = true;
            ks.metrics = [];
        end

        function apply_cluster_merge(ks, mergeInfo)
            if ~ks.isLoaded
                ks.load('loadFeatures', false, 'loadBatchwise', false);
            end
            ks.create_spike_clusters_ks2orig_if_missing();

            % apply the merges in clusterMergeInfo
            assert(isa(mergeInfo, 'Neuropixel.ClusterMergeInfo'));

            spike_clusters = ks.spike_clusters;
            cutoff_spike_clusters = ks.cutoff_spike_clusters;
            for iM = 1:mergeInfo.nMerges
                spike_clusters = apply_single_merge(spike_clusters, mergeInfo.new_cluster_ids(iM), mergeInfo.merges{iM});
                cutoff_spike_clusters = apply_single_merge(cutoff_spike_clusters, mergeInfo.new_cluster_ids(iM), mergeInfo.merges{iM});
            end
            ks.spike_clusters = spike_clusters;
            ks.cutoff_spike_clusters = cutoff_spike_clusters;

            ks.markModifiedInMemory();

            function spike_clusters = apply_single_merge(spike_clusters, dst_cluster_id, src_cluster_ids)
                mask_assign_to_dst = ismember(spike_clusters, src_cluster_ids);
                spike_clusters(mask_assign_to_dst) = dst_cluster_id;
            end
        end

        function accept_cutoff_spikes(ks, ratings_or_cluster_ids)
            % assumes that cutoff spikes are already loaded
            assert(ks.isLoadedCutoff, 'Must have loaded cutoff spikes');

            % THIS WILL NEED TO BE UDPATED IF ADDITIONAL PROPS ARE ADDED
            if nargin < 2
                cluster_ids = ks.cluster_ids;
            elseif isa(ratings_or_cluster_ids, 'Neuropixel.ClusterRatingInfo')
                cluster_ids = ratings_or_cluster_ids.cluster_ids(ratings_or_cluster_ids.includeCutoffSpikes);
            elseif islogical(ratings_or_cluster_ids)
                assert(numel(ratings_or_cluster_ids) == ks.nClusters);
                cluster_ids = ks.cluster_ids(ratings_or_cluster_ids);
            else
                cluster_ids = ratings_or_cluster_ids;
            end

            if isempty(ks.cutoff_spike_times)
                return;
            end

            accept_cutoff_mask = ismember(ks.cutoff_spike_clusters, cluster_ids);
            nCurrent = ks.nSpikes;
            nAccepted = nnz(accept_cutoff_mask);
            nTotal = nAccepted + nCurrent;
            [ks.spike_times, sortIdx] = sort(cat(1, ks.spike_times, ks.cutoff_spike_times(accept_cutoff_mask)));
            ks.cutoff_spike_times = ks.cutoff_spike_times(~accept_cutoff_mask);
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

            [ks.spike_templates, ks.cutoff_spike_templates] = combineAndSort(ks.spike_templates, ks.cutoff_spike_templates);
            [ks.spike_templates_preSplit, ks.cutoff_spike_templates_preSplit] = combineAndSort(ks.spike_templates_preSplit, ks.cutoff_spike_templates_preSplit);
            [ks.amplitudes, ks.cutoff_amplitudes] = combineAndSort(ks.amplitudes, ks.cutoff_amplitudes);
            [ks.spike_clusters, ks.cutoff_spike_clusters] = combineAndSort(ks.spike_clusters, ks.cutoff_spike_clusters);
            [ks.spike_clusters_ks2orig, ks.cutoff_spike_clusters_ks2orig] = combineAndSort(ks.spike_clusters_ks2orig, ks.cutoff_spike_clusters_ks2orig);

            if ks.hasFeaturesLoaded
                [ks.pc_features, ks.cutoff_pc_features] = combineAndSort(ks.pc_features, ks.cutoff_pc_features);
                [ks.template_features, ks.cutoff_template_features] = combineAndSort(ks.template_features, ks.cutoff_template_features);
            end

            ks.markModifiedInMemory();
        end

%         function remove_zero_spike_clusters(ks)
%             ks.cluster_ids = ks.cluster_ids(ks.cluster_spike_counts > 0 | ks.cutoff_cluster_spike_counts > 0);
%         end

        function sort_spikes(ks)
            ks.metrics = [];
            [ks.spike_times, sortIdx] = sort(ks.spike_times);
            ks.spike_templates = ks.spike_templates(sortIdx);
            ks.spike_templates_preSplit = ks.spike_templates_preSplit(sortIdx);
            ks.amplitudes = ks.amplitudes(sortIdx);
            ks.spike_clusters = ks.spike_clusters(sortIdx);
            ks.spike_clusters_ks2orig = ks.spike_clusters_ks2orig(sortIdx);
            if ks.hasFeaturesLoaded
                ks.pc_features = ks.pc_features(sortIdx, :, :);
                ks.template_features = ks.template_features(sortIdx, :);
            end

            [ks.cutoff_spike_times, sortIdx] = sort(ks.cutoff_spike_times);
            ks.cutoff_spike_templates = ks.cutoff_spike_templates(sortIdx);
            ks.cutoff_spike_templates_preSplit = ks.cutoff_spike_templates_preSplit(sortIdx);
            ks.cutoff_amplitudes = ks.cutoff_amplitudes(sortIdx);
            ks.cutoff_spike_clusters = ks.cutoff_spike_clusters(sortIdx);
            ks.cutoff_spike_clusters_ks2orig = ks.cutoff_spike_clusters_ks2orig(sortIdx);
            if ks.hasFeaturesLoaded
                ks.cutoff_pc_features = ks.cutoff_pc_features(sortIdx, :, :);
                ks.cutoff_template_features = ks.cutoff_template_features(sortIdx, :);
            end

            ks.markModifiedInMemory();
        end

        function mask_spikes(ks, mask, mask_cutoff, varargin)
            p = inputParser();
            p.addParameter('setModifiedInMemory', true, @islogical);
            p.parse(varargin{:});

            assert(islogical(mask) && numel(mask) == ks.nSpikes);
            assert(islogical(mask_cutoff) && numel(mask_cutoff) == ks.nSpikesCutoff);

            ks.metrics = [];

            ks.spike_times = ks.spike_times(mask);
            ks.spike_templates = ks.spike_templates(mask);
            if ~isempty(ks.spike_templates_preSplit)
                ks.spike_templates_preSplit = ks.spike_templates_preSplit(mask);
            end
            ks.amplitudes = ks.amplitudes(mask);
            ks.spike_clusters = ks.spike_clusters(mask);
            ks.spike_clusters_ks2orig = ks.spike_clusters_ks2orig(mask);
            if ks.hasFeaturesLoaded
                ks.pc_features = ks.pc_features(mask, :, :);
                ks.template_features = ks.template_features(mask, :);
            end

            ks.cutoff_spike_times = ks.cutoff_spike_times(mask_cutoff);
            ks.cutoff_spike_templates = ks.cutoff_spike_templates(mask_cutoff);
            ks.cutoff_spike_templates_preSplit = ks.cutoff_spike_templates_preSplit(mask_cutoff);
            ks.cutoff_amplitudes = ks.cutoff_amplitudes(mask_cutoff);
            ks.cutoff_spike_clusters = ks.cutoff_spike_clusters(mask_cutoff);
            ks.cutoff_spike_clusters_ks2orig = ks.cutoff_spike_clusters_ks2orig(mask_cutoff);
            if ks.hasFeaturesLoaded
                ks.cutoff_pc_features = ks.cutoff_pc_features(mask_cutoff, :, :);
                ks.cutoff_template_features = ks.cutoff_template_features(mask_cutoff, :);
            end

            if p.Results.setModifiedInMemory
                ks.markModifiedInMemory();
            end
        end

        function mask_clusters(ks, cluster_ids)
            [~, cluster_ids] = ks.lookup_clusterIds(cluster_ids);

            mask = ismember(ks.spike_clusters, cluster_ids);
            cutoff_mask = ismember(ks.cutoff_spike_clusters, cluster_ids);
            ks.mask_spikes(mask, cutoff_mask);
        end

        function drop_cutoff_spikes(ks)
            ks.metrics = [];
            ks.cutoff_spike_times = [];
            ks.cutoff_spike_templates = [];
            ks.cutoff_spike_templates_preSplit = [];
            ks.cutoff_amplitudes = [];
            ks.cutoff_spike_clusters = [];
            if ks.hasFeaturesLoaded
                ks.cutoff_pc_features = [];
                ks.cutoff_template_features = [];
            end

            ks.markModifiedInMemory();
        end

        function append_spikes(ks, append)
            ks.metrics = [];

            nNew = numel(append.spike_times);
            ks.spike_times = cat(1, ks.spike_times, append.spike_times);
            ks.spike_templates = cat(1, ks.spike_templates, append.spike_templates);
            ks.spike_templates_preSplit = cat(1, ks.spike_templates_preSplit, append.spike_templates_preSplit);
            ks.amplitudes = cat(1, ks.amplitudes, append.amplitudes);
            ks.spike_clusters = cat(1, ks.spike_clusters, append.spike_clusters);
            ks.spike_clusters_ks2orig = cat(1, ks.spike_clusters_ks2orig, append.spike_clusters);
            if ks.hasFeaturesLoaded
                if isempty(append.pc_features)
                    f = zeros([nNew, size(ks.pc_features, [2 3])], 'like', ks.pc_features);
                else
                    f = append.pc_features;
                end
                ks.pc_features = cat(1, ks.pc_features, f);

                if isempty(append.template_features)
                    f = zeros([nNew, size(ks.template_features, 2)], 'like', ks.template_features);
                else
                    f = append.template_features;
                end
                ks.template_features = cat(1, ks.template_features, f);
            end

            nNewCutoff = numel(append.cutoff_spike_times);
            ks.cutoff_spike_times = cat(1, ks.cutoff_spike_times, append.cutoff_spike_times);
            ks.cutoff_spike_templates = cat(1, ks.cutoff_spike_templates, append.cutoff_spike_templates);
            ks.cutoff_spike_templates_preSplit = cat(1, ks.cutoff_spike_templates_preSplit, append.cutoff_spike_templates_preSplit);
            ks.cutoff_amplitudes = cat(1, ks.cutoff_amplitudes, append.cutoff_amplitudes);
            ks.cutoff_spike_clusters = cat(1, ks.cutoff_spike_clusters, append.cutoff_spike_clusters);
            ks.cutoff_spike_clusters_ks2orig = cat(1, ks.cutoff_spike_clusters_ks2orig, append.cutoff_spike_clusters);
            if ks.hasFeaturesLoaded
                if isempty(append.pc_features)
                    f = zeros([nNewCutoff, size(ks.cutoff_pc_features, [2 3])], 'like', ks.cutoff_pc_features);
                else
                    f = append.cutoff_pc_features;
                end
                ks.cutoff_pc_features = cat(1, ks.cutoff_pc_features, f);

                if isempty(append.template_features)
                    f = zeros([nNewCutoff, size(ks.cutoff_template_features, 2)], 'like', ks.cutoff_template_features);
                else
                    f = append.cutoff_template_features;
                end
                ks.cutoff_template_features = cat(1, ks.cutoff_template_features, f);
            end

            ks.markModifiedInMemory();
        end

%         function loadFromRez(ks, rez)
%             % NOTE: not completed, currently just holding useful code
%             error('not yet functional');
%              % load rez.mat
%             if hasRez
%                 prog.increment('Loading rez.mat');
%                 d = load(fullfile(p, 'rez.mat'));
%                 rez = d.rez;
%
%                 if isfield(rez, 'est_contam_rate')
%                     ks.kilosort_version = 2;
%                 else
%                     ks.kilosort_version = 1;
%                 end
%
%                 if ks.kilosort_version == 2
%                     ops = rez.ops;
%                     ks.W = rez.W;
%                     ks.U = rez.U;
%
%                     % strip leading zeros off of ks.templates based on size of W
%                     nTimeTempW = size(ks.W, 1);
%                     nStrip = size(ks.templates, 2) - nTimeTempW;
%                     if nStrip > 0
%                         ks.templates = ks.templates(:, nStrip+1:end, :);
%                         ks.template_sample_offset = ks.template_sample_offset - uint64(nStrip);
%                     end
%
%                     nBatches  = ops.Nbatch;
%                     NT  	= uint64(ops.NT);
%                     ks.batch_starts = (uint64(1) : NT : (NT*uint64(nBatches) + uint64(1)))';
%                     ks.batchwise_cc = rez.ccb;
%                     ks.batch_sort_order = rez.iorig;
%
%                     % avoiding referencing the dependent properties in the constructor
%                     % so we don't need to worry about what properties they reference
%                     % that might not have been set yet
%                     nTemplateRank = size(rez.WA, 3);
%                     nTemplatePCs = size(rez.W_a, 2);
%                     nTemplates = size(ks.templates, 1);
%                     nTemplateTimepoints = size(ks.templates, 2);
%                     nChannelsSorted = numel(ks.channel_ids_sorted);
%
%                     ks.W_batch = single(rez.WA); % [nTemplateTimepoints, nTemplates, nTemplateRank, nBatches]
%                     ks.W_batch_US = reshape(single(rez.W_a), [nTemplateTimepoints, nTemplateRank, nTemplatePCs, nTemplates]);
%                     ks.W_batch_V = single(rez.W_b);
%
%                     ks.U_batch = single(rez.UA);
%                     ks.U_batch_US = reshape(single(rez.U_a), [nChannelsSorted, nTemplateRank, nTemplatePCs, nTemplates]);
%                     ks.U_batch_V = single(rez.U_b);
%
%                     ks.cluster_est_contam_rate = getOr(rez, 'est_contam_rate');
%                     ks.cluster_merge_count = getOr(rez, 'mergecount');
%                     ks.cluster_split_src = getOr(rez, 'splitsrc');
%                     ks.cluster_split_dst = getOr(rez, 'splitdst');
%                     ks.cluster_split_auc = getOr(rez, 'splitauc');
%                     ks.cluster_split_candidate = getOr(rez, 'split_candidate');
%                     ks.cluster_orig_template = getOr(rez, 'cluster_split_orig_template');
%
%                     if isfield(rez, 'st3_cutoff_invalid')
%                         % include invalid spikeas
%                         cutoff_spike_times = rez.st3_cutoff_invalid(:, 1);
%                         [~, isort] = sort(cutoff_spike_times);
%                         cols = size(rez.st3_cutoff_invalid, 2);
%                         if cols > 5
%                             col = cols;
%                         else
%                             col = 2;
%                         end
%
%                         ks.cutoff_spike_times = cutoff_spike_times(isort);
%                         ks.cutoff_spike_templates = rez.st3_cutoff_invalid(isort, 2);
%                         ks.cutoff_spike_templates = ks.cutoff_spike_templates + ones(1, 'like', ks.spike_templates); % 0 indexed templates to 1 indexed templates
%                         ks.cutoff_amplitudes = rez.st3_cutoff_invalid(isort, 3);
%                         ks.cutoff_spike_clusters = uint32(rez.st3_cutoff_invalid(isort, col) - 1); % rez.st3_cutoff is 1 indexed, cluster ids are 0 indexed
%
%                         ks.cutoff_pc_features = read('cutoff_pc_features');
%                         ks.cutoff_template_features = read('cutoff_template_features');
%                     end
%                 end
%             end
%
%                 function v = getOr(s, fld, def)
%                     if isfield(s, fld)
%                         v = s.(fld);
%                     elseif nargin > 2
%                         v = def;
%                     else
%                         v = [];
%                     end
%                 end
%             end
%         end

        function checkLoaded(ks)
            if isempty(ks.amplitudes)
                ks.load();
            end
        end

        function sync = readSync(ks)
            if isempty(ks.raw_dataset)
                sync = zeros(0, 1, 'uint16');
            else
                sync = ks.raw_dataset.readSync();
            end
        end

        function sync = loadSync(ks, varargin)
            if isempty(ks.raw_dataset)
                sync = zeros(0, 1, 'uint16');
            else
                sync = ks.raw_dataset.readSync(varargin{:});
            end
        end

        function setSyncBitNames(ks, idx, names)
            if ~isempty(ks.raw_dataset)
                ks.raw_dataset.setSyncBitNames(idx, names);
            else
                if isscalar(idx) && ischar(names)
                    ks.syncBitNames{idx} = names;
                else
                    assert(iscellstr(names)) %#ok<ISCLSTR>
                    ks.syncBitNames(idx) = names;
                end
            end
        end

        function idx = lookupSyncBitByName(ks, names)
            if ischar(names)
                names = {names};
            end
            assert(iscellstr(names));

            [tf, idx] = ismember(names, ks.syncBitNames);
            idx(~tf) = NaN;
        end

        function seg = segmentIntoClusters(ks)
            % generates a KilosortTrialSegmentedDataset with only 1 trial.
            % Segments as though there is one all-inclusive trials, so that
            % clusters are split but not trials
            tsi = Neuropixel.TrialSegmentationInfo(1);

            tsi.trialId = 0;
            tsi.conditionId = 0;
            tsi.idxStart = 1;

            if ks.hasRawDataset
                tsi.idxStop = ks.raw_dataset.nSamplesAP;
            else
                tsi.idxStop = max(ks.spike_times) + 1;
            end

            seg = ks.segmentIntoTrials(tsi, 0);
        end

        function seg = segmentIntoTrials(ks, tsi, trial_ids)

            % trialInfo hs fields:
            %   trialId
            %   conditionId
            %   idxStart
            %   idxStop

            seg = Neuropixel.KilosortTrialSegmentedDataset(ks, tsi, trial_ids);
        end

        function m = computeMetrics(ks, recompute, varargin)
            p = inputParser();
            p.addParameter('loadBatchwise', false, @islogical);
            p.addParameter('loadFeatures', false, @islogical);
            p.addParameter('loadCutoff', true, @islogical);
            p.addParameter('loadPreSplit', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            if isempty(ks.metrics) || ~isvalid(ks.metrics) || (nargin >= 2 && recompute)
                ks.load(p.Results);
                ks.metrics = Neuropixel.KilosortMetrics(ks, p.Unmatched);
            end
            m = ks.metrics;
        end

        function [clusterInds, cluster_ids] = lookup_clusterIds(ks, cluster_ids)
            if islogical(cluster_ids)
                cluster_ids = ks.cluster_ids(cluster_ids);
             end
            [tf, clusterInds] = ismember(cluster_ids, ks.cluster_ids);
            assert(all(tf, 'all'), 'Some cluster ids were not found in ks.clusterids');
        end

        function  [sortedChannelInds, channelIds] = lookup_sortedChannelIds(ks, channelIds)
             if islogical(channelIds)
                channelIds = ks.channel_ids_sorted(channelIds);
             end
            [tf, sortedChannelInds] = ismember(channelIds, ks.channel_ids_sorted);
            assert(all(tf, 'all'), 'Some channel ids not found');
        end

        function [channelInds, channelIds] = lookup_channelIds(ks, channelIds)
             [channelInds, channelIds] = ks.channelMap.lookup_channelIds(channelIds);
        end

        function [fileInds, origSampleInds] = lookup_sampleIndexInConcatenatedFile(ks, inds)
           [fileInds, origSampleInds] = ks.concatenationInfo.lookup_sampleIndexInSourceFiles(inds);
        end
    end

    methods % Snippet set construction
        function [snippetSet, reconstructionFromOtherClusters] = getWaveformsFromRawData(ks, varargin)
             % Extracts individual spike waveforms from the raw datafile, for multiple
             % clusters. Returns the waveforms and their means within clusters.
             %
             % Based on original getWaveForms function by C. Schoonover and A. Fink
             %
             % OUTPUT
             % wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
             % wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
             % wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
             % wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
             %                                          % nClu: number of different clusters in .spikeClusters
             %                                          % nSWf: number of samples per waveform
             %
             % % USAGE
             % wf = getWaveForms(gwfparams);

             % Load .dat and Kilosort/Phy output

             p = inputParser();
             % specify (spike_idx or spike_times) and/or cluster_ids
             % if both are specified, spikes in the list of times with belonging to those cluster_ids will be kept
             p.addParameter('spike_idx', [], @(x) isempty(x) || isvector(x)); % manually specify which idx into spike_times
             p.addParameter('spike_times', [], @(x) isempty(x) || isvector(x)); % manually specify which times directly to extract
             p.addParameter('cluster_ids', [], @(x) isempty(x) || isvector(x)); % manually specify all spikes from specific cluster_ids
             
             % for averaging by cluster
             p.addParameter('average_by_cluster_id', false, @islogical); % returns a snippet set with one snippet (mean) for each cluster_ids
             p.addParameter('average_by_cluster_id_batchwise', false, @islogical); % returns a snippet set with one snippet (mean) for each cluster_id for each batch (stored in ss.group
             p.addParameter('average_by_group_id', false, @islogical);
             p.addParameter('group_ids', [], @(x) isempty(x) || isvector(x)); % for manually indicating groups to average, must match numel(spike_idx) or numel(spike_times)
             
             % and ONE OR NONE of these to pick channels (or channels for each cluster)
             p.addParameter('channel_ids', [], @(x) isempty(x) || isvector(x));
             p.addParameter('channel_ids_by_cluster', [], @(x) isempty(x) || ismatrix(x));
             p.addParameter('best_n_channels', NaN, @isscalar); % or take the best n channels based on this clusters template when cluster_id is scalar
    
             % if specified, spikes times will be filtered within this window
             p.addParameter('filter_window', [], @(x) isempty(x) || isvector(x));

             % other params:
             p.addParameter('num_waveforms', Inf, @isscalar); % caution: Inf will request ALL waveforms in order (typically useful if spike_times directly specified)
             p.addParameter('subselect_waveforms_mode', "random", @isstringlike); % modes include "random", "first". When num_waveforms is infinite, how to pick the ones we keep
             p.addParameter('random_seed', 'shuffle'); 
             p.addParameter('primary_template_only', true, @islogical); % only sample waveforms that use each cluster's primary template
            
             p.addParameter('window', [-40 41], @isvector); % Number of samples before and after spiketime to include in waveform
             p.addParameter('car', false, @islogical);
             p.addParameter('centerUsingFirstSamples', 20, @(x) isscalar(x) || islogical(x)); % subtract mean of each waveform's first n samples, don't do if false
             p.addParameter('subtractOtherClusters', false, @(x) islogical(x) || isstringlike(x)); % time consuming step to remove the contribution of the other clusters to a given snippet. Auto is a mode that subtracts only if it decreases the wavefroms variance
             p.addParameter('excludeClusterFromOwnReconstruction', true, @islogical); % avoid reconstructing a cluster with itself (set false if a cell has many spikes in close proximity
             p.addParameter('applyScaling', false, @islogical); % scale to uV
             
             p.addParameter('band', 'ap', @ischar); % 'ap' or 'lf'
             p.addParameter('raw_dataset', ks.raw_dataset, @(x) true);
             p.addParameter('fromSourceDatasets', false, @islogical); % go all the way back to the imecDatasets that were concatenated to form ks.raw_dataset

             % other metadata set in snippetSet
             p.addParameter('trial_idx', [], @isvector);

             p.addParameter('data_distrust_mask', [], @(x) isempty(x) || islogical(x)); % if provided, these samples will be flagged when included in snippets via ss.data_trust_mask

             p.addParameter('from_cutoff_spikes', false, @islogical);

             p.parse(varargin{:});

             raw_dataset = p.Results.raw_dataset;
             if isempty(raw_dataset)
                 error('KilosortDataset has no raw_dataset and raw_dataset not provided');
             end
             
             band = p.Results.band;
             fromSource = p.Results.fromSourceDatasets;
             switch band
                 case 'ap'
                     if fromSource
                         assert(raw_dataset.hasSourceAP);
                     else
                         assert(raw_dataset.hasAP);
                     end
                 case 'lf'
                     assert(~p.Results.subtractOtherClusters, ...
                         'Subtract other clusters not supported for lf band, since templates are in ap band');
                     if fromSource
                         assert(raw_dataset.hasSourceLF);
                     else
                         assert(raw_dataset.hasLF);
                     end
                 otherwise
                     error('Unknown band %s', band);
             end

             if p.Results.average_by_cluster_id && p.Results.subtractOtherClusters
                 error('Cannot average and subtractOtherClusters');
             end

             ks.checkLoaded();

             from_cutoff_spikes = p.Results.from_cutoff_spikes;

             if ~isempty(p.Results.spike_times)
                 spike_times = p.Results.spike_times; %#ok<*PROPLC>
                 if ~from_cutoff_spikes
                    [tf, spike_idx] = ismember(spike_times, ks.spike_times);
                 else
                     [tf, spike_idx] = ismember(spike_times, ks.cutoff_spike_times);
                 end
                 if any(~tf)
                     error('Not all spike times were found in KilosortDataset');
                 end


             elseif ~isempty(p.Results.spike_idx)
                 spike_idx = p.Results.spike_idx;
                 if ~from_cutoff_spikes
                    spike_times = ks.spike_times(spike_idx);
                 else
                     spike_times = ks.cutoff_spike_times(spike_idx);
                 end

             elseif ~isempty(p.Results.cluster_ids)
                 clu = p.Results.cluster_ids;

                 if ~from_cutoff_spikes
                     if isscalar(clu)
                         spike_idx = find(ks.spike_clusters == clu);
                         spike_times = ks.spike_times(spike_idx);
                     else
                         mask = ismember(ks.spike_clusters, clu);
                         spike_times = ks.spike_times(mask);
                         spike_idx = find(mask);
                     end
                 else
                     if isscalar(clu)
                         spike_idx = find(ks.cutoff_spike_clusters == clu);
                         spike_times = ks.cutoff_spike_times(spike_idx);
                     else
                         mask = ismember(ks.cutoff_spike_clusters, clu);
                         spike_times = ks.cutoff_spike_times(mask);
                         spike_idx = find(mask);
                     end
                 end
             else
                 error('Must specify one of spike_times, spike_idx, or cluster_id');
             end

             mask = true(numel(spike_idx), 1);

             % filter by cluster_ids if specified
             unique_cluster_ids = p.Results.cluster_ids; % allow caller to specify directly

             if ~isempty(unique_cluster_ids)
                 if ~from_cutoff_spikes
                     mask = mask & ismember(ks.spike_clusters(spike_idx), unique_cluster_ids);
                 else
                     mask = mask & ismember(ks.cutoff_spike_clusters(spike_idx), unique_cluster_ids);
                 end
             end

             % filter by sample_window range
             filter_window = p.Results.filter_window;
             if ~isempty(filter_window)
                 mask = mask & spike_times >= uint64(filter_window(1)) & spike_times <= uint64(filter_window(2));
             end

             if ~isempty(p.Results.primary_template_only)
                 % filter only spikes that use its cluster's primary template
                 which_cluster_id = ks.spike_clusters(spike_idx);
                 which_template = ks.spike_templates(spike_idx);

                 metrics = ks.computeMetrics();
                 which_cluster_ind = metrics.lookup_clusterIds(which_cluster_id);
                 which_template_primary = metrics.cluster_template_mostUsed(which_cluster_ind);
                 
                 mask = mask & which_template == which_template_primary;
             end

             % apply mask
             spike_idx = spike_idx(mask);
             spike_times = spike_times(mask);

             % fetch other info
             if ~from_cutoff_spikes
                cluster_ids = ks.spike_clusters(spike_idx);
                spike_amp = ks.amplitudes(spike_idx);
             else
                cluster_ids = ks.cutoff_spike_clusters(spike_idx);
                spike_amp = ks.cutoff_amplitudes(spike_idx);
             end
             if isempty(unique_cluster_ids)
                 unique_cluster_ids = unique(cluster_ids);
             end
             
             % handle clusterwise / batchwise / groupwise averaging
             average_by_cluster_id = p.Results.average_by_cluster_id;
             average_by_group_id = p.Results.average_by_group_id;
             group_ids = p.Results.group_ids;
             if ~isempty(group_ids)
                % manual grouping
                assert(numel(group_ids) == numel(mask), 'Manually specified group_ids must match number of spikes');
                group_ids = group_ids(mask);
                
             elseif p.Results.average_by_cluster_id_batchwise
                % grouping by ks batch
                group_ids = ks.compute_which_batch(spike_times);
                assert(~average_by_group_id, 'average_by_group_id ignored when average_by_cluster_id_batchwise is set');
                average_by_group_id = true;
                average_by_cluster_id = true;
                
             else
                 assert(~average_by_group_id, 'group_ids must be specified for average_by_group_id == true');
             end
             
             trial_idx = p.Results.trial_idx;
             if ~isempty(trial_idx)
                 trial_idx = trial_idx(mask);
                 assert(numel(trial_idx) == numel(spike_idx));
             end

             % take max of num_waveforms from each cluster
             subselect_waveforms_mode = string(p.Results.subselect_waveforms_mode);
             rs = RandStream('mt19937ar','Seed', p.Results.random_seed);
             if isfinite(p.Results.num_waveforms)
                 mask = false(numel(spike_idx), 1);
                 nSample = p.Results.num_waveforms;
                 [~, uclust_ind] = ismember(cluster_ids, unique_cluster_ids);
                 nClu = numel(unique_cluster_ids);
                 for iC = 1:nClu
                     thisC = find(uclust_ind == iC);
                     if numel(thisC) <= nSample
                         mask(thisC) = true;
                     else
                         % TODO: might want to add more intelligent subselection modes here, e.g. for sampling evenly with respect to batches
                         switch subselect_waveforms_mode
                             case "random"
                                 mask(thisC(randsample(rs, numel(thisC), nSample, false))) = true;
                             case "first"
                                 mask(thisC(1:nSample)) = true;
                             case "largest"
                                 % sort by amplitude
                                 [~, inds] = maxk(spike_amp(thisC), nSample);
                                 mask(thisC(inds)) = true;
                             otherwise
                                 error('Unknown subselect_waveforms_mode %s', subselect_waveforms_mode);
                         end
                     end
                 end

                 spike_idx = spike_idx(mask); %#ok<NASGU>
                 spike_times = spike_times(mask);
                 cluster_ids = cluster_ids(mask);
                 if ~isempty(trial_idx)
                     trial_idx = trial_idx(mask);
                 end
             end

             if nnz(mask) == 0
                 error('No spikes collected');
             end

             % figure out channels requested, one set for each cluster
             if ~isnan(p.Results.best_n_channels)
                 metrics = ks.computeMetrics();
                 cluster_best_template_channels = metrics.cluster_best_channels;

                 % okay to have multiple clusters
                 [~, cluster_ind] = ismember(unique_cluster_ids, metrics.cluster_ids);
                 if any(cluster_ind == 0)
                     error('Some cluster idx not found in cluster_ids');
                 end
                 channel_ids_by_cluster = cluster_best_template_channels(cluster_ind, 1:p.Results.best_n_channels)';

                 channel_id_args = {'channel_ids_by_cluster', channel_ids_by_cluster};
             elseif ~isempty(p.Results.channel_ids_by_cluster)
                 channel_ids_by_cluster = p.Results.channel_ids_by_cluster;
                 channel_id_args = {'channel_ids_by_cluster', channel_ids_by_cluster};
             elseif ~isempty(p.Results.channel_ids)
                 % same channel ids for each cluster
                 channel_id_args = {'channel_ids', p.Results.channel_ids};
             else
                 channel_ids = ks.channel_ids_sorted;
                 channel_id_args = {'channel_ids', channel_ids};
             end
             
             switch band
                 case 'ap'
                     spike_times_for_band = spike_times;
                 case 'lf'
                     spike_times_for_band = raw_dataset.closestSampleLFForAP(spike_times);
                 otherwise
                     error('Unknown band %s', band);
             end

             % channel_ids is provided since raw data often has additional channels that we're not interested in
             window = p.Results.window;
             snippetSet = raw_dataset.readSnippetSet(band, spike_times_for_band, ...
                 window, channel_id_args{:}, ...
                 'unique_cluster_ids', unique_cluster_ids, 'cluster_ids_by_snippet', cluster_ids, ...
                 'car', p.Results.car, 'fromSourceDatasets', p.Results.fromSourceDatasets, ...
                 'average_by_cluster_id', average_by_cluster_id, ...
                 'average_by_group_id', average_by_group_id, 'group_ids', group_ids, ...
                 'applyScaling', p.Results.applyScaling, ...
                 'data_distrust_mask', p.Results.data_distrust_mask);
             snippetSet.trial_idx = trial_idx;
             snippetSet.ks = ks;

             subtractOtherClusters = p.Results.subtractOtherClusters;
             subtractOtherClustersAuto = false;
             if ~islogical(subtractOtherClusters)
                 if strcmpi(subtractOtherClusters, "auto")
                     subtractOtherClusters = true;
                     subtractOtherClustersAuto = true;
                 else
                     error('Unknown value for subtractOtherClusters. Supported values are true, false and "auto"')
                 end
             end
             if subtractOtherClusters
                 reconstructionFromOtherClusters = ks.reconstructSnippetSetFromTemplates(snippetSet, ...
                     'excludeClusterFromOwnReconstruction', p.Results.excludeClusterFromOwnReconstruction);
                 % nChannels x time x nSnippets
                 candidate = snippetSet.data - reconstructionFromOtherClusters.data;
                 if subtractOtherClustersAuto
                     % apply only if the subtraction reduces a given snippet's variance
                     mask_apply = var(single(snippetSet.data), 0, [1 2]) >= var(single(candidate), 0, [1 2]);
                     snippetSet.data(:, :, mask_apply) = candidate(:, :, mask_apply);
                 else
                    snippetSet.data = candidate;
                 end
             end

             if p.Results.centerUsingFirstSamples
                 snippetSet.data = snippetSet.data - mean(snippetSet.data(:, 1:p.Results.centerUsingFirstSamples, :), 2, 'native');
             end
        end

        function snippetSet = readAPSnippetSet(ks, times, window, varargin)
            snippetSet = ks.raw_dataset.readAPSnippetSet(times, window, ...
                'channel_ids', ks.channel_ids_sorted, varargin{:});
            snippetSet.ks = ks;
        end

        function [spike_inds, template_start_ind, scaled_templates] = findSpikesOverlappingWithWindow(ks, sample_window, varargin)
            % spike_inds is nSpikes x 1 index into ks.spike_times indicating which spikes overlap with this window
            % template_start_ind indicates where each spike would begin wihtin sample_window
            % templates is nCh x nTemplateTimepoints x nSpikes
            p = inputParser();
            p.addParameter('cluster_ids', ks.cluster_ids, @(x) isempty(x) || isvector(x));
            p.addParameter('channel_ids', ks.channel_ids_sorted, @(x) isempty(x) || isvector(x)); % which channel_ids to construct templates into
            p.parse(varargin{:});

             % find spikes that would lie within this window (with padding),
            relTvec_template = ks.templateTimeRelative;
            minT = sample_window(1) - int64(relTvec_template(end));
            maxT = sample_window(2) - int64(relTvec_template(1));
            spike_inds = find(ks.spike_times >= minT & ks.spike_times <= maxT);

            % only those from clusters we wish to include
            cluster_ids = p.Results.cluster_ids;
            spike_inds(~ismember(ks.spike_clusters(spike_inds), cluster_ids)) = [];

            template_offset = find(relTvec_template == 0);
            % spike time - template_offset + int64(1) is sample index where start of template belongs
            % this above - sample_window(1) + 1 gives the index where the template should start in the sample_window
            template_start_ind = int64(ks.spike_times(spike_inds)) - int64(template_offset) + int64(1) - int64(sample_window(1)) + int64(1);

            channel_ids = p.Results.channel_ids;
            nChannels = numel(channel_ids);
            sorted_channel_inds = ks.lookup_sortedChannelIds(channel_ids);
            scaled_templates = nan(nChannels, ks.nTemplateTimepoints, numel(spike_inds));

            metrics = ks.computeMetrics();
            templates_unw =  metrics.template_unw; % unwhitened templates, but not scaled and still in quantized units (not uV)

            for iS = 1:numel(spike_inds)
                ind = spike_inds(iS);
                amp = ks.amplitudes(ind) * ks.apScaleToUv;

                % figure out time overlap and add to reconstruction
                scaled_templates(:, :, iS)  = amp .* permute(templates_unw(ks.spike_templates(ind), :, sorted_channel_inds), [3 2 1]); % 1 x T x C --> C x T x 1
            end
        end

%         function reconstruction = reconstructRawSnippetsFromTemplates_slow(ks, times, window, varargin)
%             % generate the best reconstruction of each snippet using the amplitude-scaled templates
%             % this is the internal workhorse for reconstructSnippetSetFromTemplates
%
%             p = inputParser();
%             % these define the lookup table of channels for each cluster
%             p.addParameter('channel_ids_by_snippet', [], @(x) isempty(x) || ismatrix(x));
%             p.addParameter('unique_cluster_ids', 1, @isvector);
%             % and this defines the cluster corresponding to each spike
%             p.addParameter('cluster_ids', ones(numel(times), 1), @isvector);
%
%             p.addParameter('exclude_cluster_ids_each_snippet', [], @(x) isempty(x) || isvector(x) || iscell(x));
%             p.addParameter('exclude_cluster_ids_all_snippets', [], @(x) isempty(x) || isvector(x));
%             p.addParameter('showPlots', false, @islogical);
%             p.addParameter('rawData', [], @isnumeric); % for plotting only
%
%             p.addParameter('use_batchwise_templates', ~isempty(ks.W_batch), @islogical);
%             p.parse(varargin{:});
%             showPlots = p.Results.showPlots;
%
%             % check sizes of everything
%             nTimes = numel(times);
%             channel_ids_by_snippet = p.Results.channel_ids_by_snippet;
%             assert(~isempty(channel_ids_by_snippet));
%             if size(channel_ids_by_snippet, 2) == 1
%                 channel_ids_by_snippet = repmat(channel_ids_by_snippet, 1, numel(times));
%             end
%             assert(numel(times) <= size(channel_ids_by_snippet, 2));
%
%             nChannelsSorted = size(channel_ids_by_snippet, 1);
%             unique_cluster_ids = p.Results.unique_cluster_ids;
%
%             cluster_ids = p.Results.cluster_ids;
%             assert(numel(cluster_ids) >= nTimes);
%
%             exclude_cluster_ids_each_snippet = p.Results.exclude_cluster_ids_each_snippet;
%             exclude_cluster_ids_all_snippets = p.Results.exclude_cluster_ids_all_snippets;
%
%             % templates post-whitening is nTemplates x nTimepoints x nChannelsFull
%             use_batchwise_templates = p.Results.use_batchwise_templates;
%             if use_batchwise_templates && (isempty(ks.W_batch) || isempty(ks.U_batch))
%                 warning('Cannot use batchwise tempaltes as W_batch or U_batch was not loaded successfully from rez.mat');
%                 use_batchwise_templates = false;
%             end
%
%             if ~use_batchwise_templates
%                 metrics = ks.computeMetrics();
%                 templates =  metrics.template_unw; % unwhitened templates, but not scaled and still in quantized units (not uV)
%             end
%
%             relTvec_template = ks.templateTimeRelative;
%             relTvec_snippet = int64(window(1):window(2));
%             reconstruction = zeros(nChannelsSorted, numel(relTvec_snippet), nTimes, 'int16');
%
%             if exist('ProgressBar', 'class') == 8
%                 prog = ProgressBar(nTimes, 'Reconstructing templates around snippet times');
%             else
%                 prog = Neuropixel.Utils.ProgressBar(nTimes, 'Reconstructing templates around snippet times');
%             end
%
%             for iT = 1:nTimes
%                 prog.update(iT);
%
%                 reconstruction_this = zeros(nChannelsSorted, numel(relTvec_snippet), 'single');
%
%                 % find spikes that would lie within this window (with padding),
%                 % excluding those from clusters we wish to exclude
%                 t = int64(times(iT));
%                 minT = t + window(1) - int64(relTvec_template(end));
%                 maxT = t + window(2) - int64(relTvec_template(1));
%                 if iscell(exclude_cluster_ids_each_snippet)
%                     exclude_this = exclude_cluster_ids_each_snippet{iT};
%                 elseif ~isempty(exclude_cluster_ids_each_snippet)
%                     exclude_this = exclude_cluster_ids_each_snippet(iT);
%                 else
%                     exclude_this = [];
%                 end
%                 exclude_this = union(exclude_this, exclude_cluster_ids_all_snippets);
%
%                 nearby_spike_inds = find(ks.spike_times >= minT & ks.spike_times <= maxT);
%                 %nearby_spike_inds = find(ks.spike_times == t);
%
%                 nearby_spike_inds(ismember(ks.spike_clusters(nearby_spike_inds), exclude_this)) = [];
%
%                 if use_batchwise_templates
%                     spike_batches = ks.compute_which_batch(ks.spike_times(nearby_spike_inds));
%                 end
%
%                 % figure out what channels we need to reconstruct onto
%                 cluster_ids_this = cluster_ids(iT);
%                 [~, cluster_ind_this] = ismember(cluster_ids_this, unique_cluster_ids);
%                 assert(cluster_ind_this > 0, 'Cluster for times(%d) not found in unique_cluster_ids', iT);
%                 channel_ids_this = channel_ids_by_snippet(:, iT);
%                 sorted_channel_inds_this = ks.lookup_sortedChannelIds(channel_ids_this);
%
%                 if showPlots
%                     clf;
%                     plot(relTvec_snippet, p.Results.rawData(1, :, iT), 'k-', 'LineWidth', 2);
%                     hold on;
%                 end
%
%                 % loop over the enarby sp
%                 for iS = 1:numel(nearby_spike_inds)
%                     ind = nearby_spike_inds(iS);
%                     amp = ks.amplitudes(ind);
%
%                     % figure out time overlap and add to reconstruction
%                     tprime = int64(ks.spike_times(ind));
%                     indFromTemplate = relTvec_template + tprime >= t + relTvec_snippet(1) & relTvec_template + tprime <= t + relTvec_snippet(end);
%                     indInsert = relTvec_snippet + t >= relTvec_template(1) + tprime & relTvec_snippet + t <= relTvec_template(end) + tprime;
%
%                     if use_batchwise_templates
%                         template_this = ks.construct_batchwise_templates(ks.spike_templates(ind), 'batches', spike_batches(iS));
%                         template_this = template_this(:, indFromTemplate, sorted_channel_inds_this);
%                     else
%                         template_this = templates(ks.spike_templates(ind), indFromTemplate, sorted_channel_inds_this);
%                     end
%                     insert = amp .* permute(template_this, [3 2 1]);
%                     reconstruction_this(:, indInsert) = reconstruction_this(:, indInsert) + insert;
%
%                     if showPlots
%                         if tprime == t
%                             args = {'LineWidth', 1, 'Color', 'g'};
%                         else
%                             args = {};
%                         end
%                         plot(relTvec_template(indInsert), amp .* template_this(:, :, 1), 'Color', [0.7 0.7 0.7], args{:});
%                     end
%                 end
%
%                 if showPlots
%                     plot(relTvec_snippet, reconstruction(1, :, iT), 'r--', 'LineWidth', 2);
%                     plot(relTvec_snippet, p.Results.rawData(1, :, iT) - int16(reconstruction(1, :, iT)), 'b--');
%                     pause;
%                 end
%
%                 reconstruction(:, :, iT) = int16(reconstruction_this);
%             end
%             prog.finish();
%         end

        function reconstruction = reconstructRawSnippetsFromTemplates(ks, times, window, varargin)
            % generate the best reconstruction of each snippet using the amplitude-scaled templates
            % this is the internal workhorse for reconstructSnippetSetFromTemplates

            p = inputParser();
            % these define the lookup table of channels for each cluster
            p.addParameter('channel_ids_by_snippet', [], @(x) isempty(x) || ismatrix(x));
            p.addParameter('unique_cluster_ids', 1, @isvector);
            % and this defines the cluster corresponding to each spike
            p.addParameter('cluster_ids', ones(numel(times), 1), @isvector);

            p.addParameter('exclude_cluster_ids_each_snippet', [], @(x) isempty(x) || isvector(x) || iscell(x));
            p.addParameter('exclude_cluster_ids_all_snippets', [], @(x) isempty(x) || isvector(x));
            p.addParameter('rawData', [], @isnumeric); % for plotting only
            p.addParameter('showPlots', false, @islogical); % not supported
            p.addParameter('use_batchwise_templates', ~isempty(ks.W_batch), @islogical);
            p.addParameter('parallel', false, @islogical);
            p.parse(varargin{:});

            if p.Results.showPlots
                error('No longer supported');
            end

            % check sizes of everything
            nTimes = numel(times);
            channel_ids_by_snippet = p.Results.channel_ids_by_snippet;
            assert(~isempty(channel_ids_by_snippet));
            if size(channel_ids_by_snippet, 2) == 1
                channel_ids_by_snippet = repmat(channel_ids_by_snippet, 1, numel(times));
            end
            assert(numel(times) <= size(channel_ids_by_snippet, 2)); % TODO CHANGE TO ==

            debug('Preparing to reconstruct templates around snippet times\n');

            nChannelsSorted = size(channel_ids_by_snippet, 1);
            unique_cluster_ids = p.Results.unique_cluster_ids;
            sorted_channel_inds_by_snippet = ks.lookup_sortedChannelIds(channel_ids_by_snippet);

            % lookup cluster_ids
            cluster_ids = p.Results.cluster_ids;
            assert(numel(cluster_ids) >= nTimes); % TODO CHANGE TO ==
            [~, cluster_inds_in_unique] = ismember(cluster_ids, unique_cluster_ids);
            assert(all(cluster_inds_in_unique > 0), 'Some cluster_ids not found in unique_cluster_ids');

            exclude_cluster_ids_each_snippet = p.Results.exclude_cluster_ids_each_snippet;
            if isempty(exclude_cluster_ids_each_snippet)
                exclude_cluster_ids_each_snippet = cell(nTimes, 1);
            elseif isnumeric(exclude_cluster_ids_each_snippet)
                exclude_cluster_ids_each_snippet = num2cell(exclude_cluster_ids_each_snippet);
            end

            exclude_cluster_ids_all_snippets = p.Results.exclude_cluster_ids_all_snippets;

            % templates post-whitening is nTemplates x nTimepoints x nChannelsFull
            use_batchwise_templates = p.Results.use_batchwise_templates;
            if use_batchwise_templates && (isempty(ks.W_batch) || isempty(ks.U_batch))
                warning('Cannot use batchwise tempaltes as W_batch or U_batch was not loaded successfully from rez.mat');
                use_batchwise_templates = false;
            end

            % decide on templates we're using
            if ~use_batchwise_templates
                metrics = ks.computeMetrics();
                templates =  metrics.template_unw; % unwhitened templates, but not scaled and still in quantized units (not uV)
            else
                if isempty(ks.W_batch) || isempty(ks.U_batch)
                    error('No batchwise information present in KilosortDataset');
                end
                templates = [];
            end

            relTvec_template = ks.templateTimeRelative;
            relTvec_snippet = int64(window(1):window(2));
            reconstruction = zeros(nChannelsSorted, numel(relTvec_snippet), nTimes, 'int16');

            % pre-extract everything in case we switch to parfor
            spike_times = ks.spike_times;
            spike_templates = ks.spike_templates;
            spike_clusters = ks.spike_clusters;
            amplitudes = ks.amplitudes;
            
            if use_batchwise_templates
                % use batchwise templates for reconstructing each spike more accurately
                W_batch = ks.W_batch;
                U_batch = ks.U_batch;
                spike_batches = ks.compute_which_batch(ks.spike_times);
                whitening_mat_inv = ks.whitening_mat_inv;
                assert(size(W_batch, 3) == 3);
            end
    
            % find spikes that would lie within this window (with padding),
            % excluding those from clusters we wish to exclude
            minT = int64(window(1)) - int64(relTvec_template(end));
            maxT = int64(window(2)) - int64(relTvec_template(1));
%             t_center_offset = (minT + maxT) / 2;
%             search_half_width = ceil(maxT - t_center_offset);
%             c_nearby_spike_inds = rangesearch(double(ks.spike_times), double(times) + double(t_center_offset), search_half_width, 'SortIndices', false);
%
            c_nearby_spike_inds = Neuropixel.Utils.simple_rangesearch(ks.spike_times, times, [minT maxT]);

            prog = Neuropixel.Utils.ProgressBar(nTimes, 'Reconstructing templates around snippet times');
            for iT = 1:nTimes
                t = int64(times(iT));
                reconstruction_this = zeros(nChannelsSorted, numel(relTvec_snippet), 'single');

                % exclude specified clusters from nearby spikes
                nearby_spike_inds = c_nearby_spike_inds{iT};
                nearby_spike_clusters = spike_clusters(nearby_spike_inds);
                mask_nearby_okay = ~ismember(nearby_spike_clusters, exclude_cluster_ids_each_snippet{iT}) & ...
                                    ~ismember(nearby_spike_clusters, exclude_cluster_ids_all_snippets);

                % figure out what channels we need to reconstruct onto
                sorted_channel_inds_this = sorted_channel_inds_by_snippet(:, iT);

                % loop over the enarby sp
                for iS = 1:numel(nearby_spike_inds)
                    if ~mask_nearby_okay(iS)
                        continue;
                    end
                    ind = nearby_spike_inds(iS);
                    amp = amplitudes(ind);

                    % figure out time overlap and add to reconstruction
                    tprime = int64(spike_times(ind));
                    maskFromTemplate = relTvec_template + tprime >= t + relTvec_snippet(1) & relTvec_template + tprime <= t + relTvec_snippet(end);
                    maskInsert = relTvec_snippet + t >= relTvec_template(1) + tprime & relTvec_snippet + t <= relTvec_template(end) + tprime;

                    if use_batchwise_templates
                        template_ind = spike_templates(ind);
                        batch_ind = spike_batches(ind);
                        wmi_this = whitening_mat_inv(:, sorted_channel_inds_this);

                        W = W_batch(maskFromTemplate, template_ind, :, batch_ind); % nTT x 1 x 3 x 1
                        U = U_batch(:, template_ind, :, batch_ind); % nCh x 1 x 3 x 1

                        % manual loop unrolling to make this faster
                        template_this = shiftdim(W(:,:,1,:)*(U(:,:,1,:)'*wmi_this) + W(:,:,2,:)*(U(:,:,2,:)'*wmi_this) + W(:,:,3,:)*(U(:,:,3,:)'*wmi_this), -1);
                    else
                        template_this = templates(spike_templates(ind), maskFromTemplate, sorted_channel_inds_this);
                    end
                    insert = amp .* permute(template_this, [3 2 1]);
                    reconstruction_this(:, maskInsert) = reconstruction_this(:, maskInsert) + insert;
                end

                reconstruction(:, :, iT) = int16(reconstruction_this);
                prog.update(iT);
            end
            prog.finish();
        end

        function ssReconstruct = reconstructSnippetSetFromTemplates(ks, ss, varargin)
            % used to build a snippet set consisting of all the templates inserted at all of the spike times
            % that overlap with the time windows present in ss (the source snippet set)
            % this is typically used to clean the snippet set by subtracting the inferred spike templates for a subset
            % of clusters (those not of interest)
            p = inputParser();
            p.addParameter('excludeClusterFromOwnReconstruction', false, @islogical);
            p.addParameter('exclude_cluster_ids', [], @(x) isempty(x) || isvector(x));
            p.addParameter('showPlots', false, @islogical);
            p.parse(varargin{:});

            if p.Results.excludeClusterFromOwnReconstruction
                exclude_cluster_ids_each_snippet = ss.cluster_ids;
            else
                exclude_cluster_ids_each_snippet = [];
            end

            reconstruction = ks.reconstructRawSnippetsFromTemplates(ss.sample_idx, ss.window, ...
                'channel_ids_by_snippet', ss.channel_ids_by_snippet, 'unique_cluster_ids', ss.unique_cluster_ids, ...
                'cluster_ids', ss.cluster_ids, 'exclude_cluster_ids_each_snippet', exclude_cluster_ids_each_snippet, ...
                'exclude_cluster_ids_all_snippets', p.Results.exclude_cluster_ids, ...
                'showPlots', p.Results.showPlots, 'rawData', ss.data);

            ssReconstruct = Neuropixel.SnippetSet(ks);
            ssReconstruct.data = int16(reconstruction);
            ssReconstruct.sample_idx = ss.sample_idx;
            ssReconstruct.channel_ids_by_snippet = ss.channel_ids_by_snippet;
            ssReconstruct.cluster_ids = ss.cluster_ids;
            ssReconstruct.window = ss.window;
        end

        function snippetSet = readAPSnippsetSet_clusterIdSubset(ks, times, window, cluster_ids, varargin)
            % a utility used by several methods in KilosortMetrics. Extracts snippets at times in times, but sets the snippet set to
            % overlay waveforms for specific cluster_ids of interest, and also selects only the channels for those clusters
             p = inputParser();

             % and ONE OR NONE of these to pick channels (or channels for each cluster) that will be highlighted
             p.addParameter('channel_ids', [], @(x) isempty(x) || ismatrix(x));
             p.addParameter('best_n_channels', [], @(x) isempty(x) || isscalar(x)); % if specified, the best n channels for each cluster will be chosen

             p.addParameter('car', false, @islogical);
             p.addParameter('centerUsingFirstSamples', 20, @(x) isscalar(x) || islogical(x)); % subtract mean of each waveform's first n samples, don't do if false

             p.addParameter('subtractOtherClusters', false, @islogical); % time consuming step to remove the contribution of the other clusters (those not in cluster_ids) to a given snippet
             p.addParameter('raw_dataset', ks.raw_dataset, @(x) true);

             p.parse(varargin{:});

             raw_dataset = p.Results.raw_dataset;
             channel_ids = p.Results.channel_ids;


             % automatically pick clusters
             if isempty(p.Results.channel_ids)
                 if isempty(p.Results.best_n_channels)
                    channel_ids = raw_dataset.goodChannels;
                 else
                     m = ks.computeMetrics();
                     channel_ids = m.gather_best_channels_multiple_clusters(cluster_ids, p.Results.best_n_channels);
                 end
             end

             snippetSet = raw_dataset.readAPSnippetSet(times, window, 'channel_ids', channel_ids, 'car', p.Results.car);
             snippetSet.ks = ks;

             if p.Results.subtractOtherClusters
                 reconstructionFromOtherClusters = ks.reconstructSnippetSetFromTemplates(snippetSet, ...
                     'exclude_cluster_ids', cluster_ids, ...
                     'excludeClusterFromOwnReconstruction', false);
                 snippetSet.data = snippetSet.data - reconstructionFromOtherClusters.data;
             end

             if p.Results.centerUsingFirstSamples
                 snippetSet.data = snippetSet.data - mean(snippetSet.data(:, 1:p.Results.centerUsingFirstSamples, :), 2, 'native');
             end

             % specify these cluster_ids for overlaying waveforms
             snippetSet.overlay_cluster_ids = cluster_ids;

        end
    end

    methods
        function [batch_ind, sorted_batch_ind] = compute_which_batch(ks, sample_idx)
            assert(~isempty(ks.batch_starts), 'batch_starts is empty. Load KS dataset with loadBatchwise=true');
            batch_ind = discretize(sample_idx, [ks.batch_starts; Inf]);
            [~, batch_sorted_location] = ismember((1:ks.nBatches)', ks.batch_sort_order);
            sorted_batch_ind = batch_sorted_location(batch_ind);
        end

        function data = construct_batchwise_templates(ks, template_inds, varargin)
            p = inputParser();
            p.addParameter('batches', [], @isvector);
            p.addParameter('every_n_batches', 1, @isscalar);
            p.addParameter('average_skipped_batches', false, @islogical);
            p.addParameter('whitened', true, @islogical);
            % not sure this is correct, plus we haven't saved muA to disk anyway
%             p.addParameter('scaleByBatchAmplitude', false, @islogical); % multiply by muA to scale according to batch, dont use if you plan to scale by amplitudes
            p.parse(varargin{:});
            batch_inds = p.Results.batches;
            average_skipped_batches = p.Results.average_skipped_batches;

            if isempty(batch_inds)
                batch_inds = 1:p.Results.every_n_batches:ks.nBatches;
            end
%             scaleByBatchAmplitude = p.Results.scaleByBatchAmplitude;

            % returns batch-wise templates for a set of templates
            % nTemplates x nTimePoints x nTemplateChannels x nBatches
            nTemplates = numel(template_inds);
            nBatches = numel(batch_inds);

            if isempty(ks.W_batch) || isempty(ks.U_batch)
                error('No batchwise information present in KilosortDataset');
            end

            nTT = ks.nTemplateTimepoints;
            nCh = ks.nChannelsSorted;
            data = zeros(nTemplates, nTT, nCh, nBatches);
            if p.Results.whitened
                wmi = ks.whitening_mat_inv;
            else
                wmi = eye(size(ks.whitening_mat_inv));
            end

            assert(size(ks.W_batch, 3) == 3);

            for iT = 1:nTemplates
                for iB = 1:nBatches
                    if average_skipped_batches
                        if iB == nBatches
                            batch_inds_this = batch_inds(iB):ks.nBatches;
                        else
                            batch_inds_this = batch_inds(iB):batch_inds(iB+1);
                        end
                    else
                        batch_inds_this = batch_inds(iB);
                    end

                    for iiB = 1:numel(batch_inds_this)
                        W = ks.W_batch(:, template_inds(iT), :, batch_inds_this(iiB)); % nTT x 1 x 3 x 1
                        if all(W==0, 'all')
                            continue;
                        end
                        U = ks.U_batch(:, template_inds(iT), :, batch_inds_this(iiB)); % nCh x 1 x 3 x 1

                        % manual loop unrolling to make this faster

                        data(iT, :, :, iB) = data(iT, :, :, iB) + shiftdim(W(:,:,1,:)*(U(:,:,1,:)'*wmi) + W(:,:,2,:)*(U(:,:,2,:)'*wmi) + W(:,:,3,:)*(U(:,:,3,:)'*wmi), -1);

%                         for iP = 1:size(U, 3)
%                             this_pc = W(:, :, iP, :) * (U(:, :, iP, :)' * wmi); % nTT x nCh
%                             data(iT, :, :, iB) = data(iT, :, :, iB) + reshape(this_pc, [1, nTT, nCh]);
%                         end
                    end
                    data(iT, :, :, iB) = data(iT, :, :, iB) ./ numel(batch_inds_this);
                end
            end
        end

        function [mask_dup_spikes, mask_dup_spikes_cutoff, stats] = identify_duplicate_spikes(ks, varargin)
            p = inputParser();
            p.addParameter('include_cutoff_spikes', ks.deduplicate_cutoff_spikes, @islogical);
            p.addParameter('withinSamples', ks.deduplicate_within_samples, @isscalar);
            p.addParameter('withinDistance', ks.deduplicate_within_distance, @isscalar);
            p.parse(varargin{:});

            assert(issorted(ks.spike_times), 'call .sortSpikes first');
            assert(issorted(ks.cutoff_spike_times), 'call .sortSpikes first');

            withinSamples = p.Results.withinSamples;
            withinDistance = p.Results.withinDistance;

            % lookup best channels by templates
            m = ks.computeMetrics();

            % determine which templates are withinDistance of each other templates based on centroid
            temptempdist = squareform(pdist(m.template_centroid, 'euclidean'));
            temptempprox = temptempdist < withinDistance;

            % mechanism originally described in Siegle et al. (2019) using template best channel instead of centroid
%             template_best_sorted_channel_ind = m.lookup_sortedChannelIds(m.template_best_channels(:, 1));
%             chchprox = squareform(pdist(ks.channel_positions_sorted, 'euclidean')) <= withinDistance;
%             [R, C] = ndgrid(template_best_sorted_channel_ind, template_best_sorted_channel_ind);
%             idx = sub2ind(size(chchprox), R, C);
%             temptempprox = chchprox(idx);

            % begin by combining spikes with spikes cutoff
            if p.Results.include_cutoff_spikes
                spikes = cat(1, ks.spike_times, ks.cutoff_spike_times);
                templates = cat(1, ks.spike_templates, ks.cutoff_spike_templates);
                clusters = cat(1, ks.spike_clusters, ks.cutoff_spike_clusters);
                [spikes, sort_idx] = sort(spikes);
                templates = templates(sort_idx);
                clusters = clusters(sort_idx);

                mask_from_cutoff = true(size(spikes));
                mask_from_cutoff(1:ks.nSpikes) = false;
                mask_from_cutoff = mask_from_cutoff(sort_idx);
            else
                spikes = ks.spike_times;
                templates = ks.spike_templates;
                clusters = ks.spike_clusters;
                mask_from_cutoff = false(size(spikes));
            end

            cluster_inds = ks.lookup_clusterIds(clusters);

            mask_dup = false(size(spikes));
            dup_from_template = zeros(size(spikes), 'like', templates);
            dup_from_cluster_ind = zeros(size(spikes), 'like', cluster_inds);
            prog = Neuropixel.Utils.ProgressBar(numel(spikes), 'Checking for duplicate spikes');
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
                mask_dup_spikes_cutoff = false(ks.nSpikesCutoff, 1);
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

        function stats = remove_duplicate_spikes(ks, varargin)
            p = inputParser();
            p.addParameter('saveToDisk', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            saveToDisk = p.Results.saveToDisk;

            if ks.modifiedInMemory
                error('Calling ks.remove_duplicate_spikes after modified in memory');
            end

            if ks.is_deduplicated
                if saveToDisk
                    error('Calling ks.remove_duplicate_spikes with saveToDisk true when ks.is_deduplicated is true');
                else
                    warning('Calling ks.remove_duplicate_spikes when ks.is_deduplicated is true');
                end
            end

            fprintf('Removing duplicate spikes from KS dataset\n');
            [mask_dup_spikes, mask_dup_spikes_cutoff, stats] = ks.identify_duplicate_spikes(p.Unmatched);
            ks.mask_spikes(~mask_dup_spikes, ~mask_dup_spikes_cutoff);
            ks.is_deduplicated = true;
            ks.deduplication_stats = stats;
            ks.modifiedInMemory = true;

            if saveToDisk
                toSave.spike_deduplication_mask = ~mask_dup_spikes;
                toSave.cutoff_spike_deduplication_mask = ~mask_dup_spikes_cutoff;
                toSave.deduplicate_cutoff_spikes = ks.deduplicate_cutoff_spikes;
                toSave.deduplicate_within_samples = ks.deduplicate_within_samples;
                toSave.deduplicate_within_distance = ks.deduplicate_within_distance;
                toSave.deduplication_stats = stats;

                save(fullfile(ks.path, 'spike_deduplication_mask.mat'), '-v7.3', '-struct', 'toSave');
            end
        end
    end

    methods
        function writeChannelMap(ks, outfile)
            map = ks.channelMap;
            chanMap = map.channelIdsMapped;
            xcoords = map.xcoords;
            ycoords = map.ycoords;
            % this is a mask over mapped channels that is true if channel is good
            connected = ismember(1:map.nChannelsMapped, ks.channel_ids_sorted);
            kcoords = map.shankInd;
            fs = ks.fsAP;

            if nargin < 2
                outfile = fullfile(ks.path, 'chanMap.mat');
            end
            save(outfile, 'chanMap', 'xcoords', 'ycoords', 'connected', 'kcoords', 'fs');
        end

        function writeToDisk(ks, outpath, varargin)
            p = inputParser();
            p.addParameter('progressInitializeFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(nUpdates) to print update
            p.addParameter('progressIncrementFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(updateString) to print update
            p.addParameter('allowOverwriteSelf', false, @islogical);
            p.addParameter('updateParamsDotPyDatPath', "", @isstringlike);
            p.parse(varargin{:});

            overwriteSelf = strcmp(outpath, ks.path);
            assert(p.Results.allowOverwriteSelf || ~overwriteSelf, 'Refusing to overwrite on disk, pass allowOverwriteSelf true to allow this');

            if exist(outpath, 'dir') ~= 7
                mkdirRecursive(outpath);
            end

            oldp = @(file) fullfile(ks.path, file);
            newp = @(file) fullfile(outpath, file);
            existnew = @(file) exist(newp(file), 'file') > 0;

            if existnew('spike_deduplication_mask.mat')
                % this will cause problems if it doesn't match the spikes being written down
                delete(newp('spike_deduplication_mask.mat'));
            end

            nProg = 17;
            if ks.kilosort_version == 2
                nProg = nProg + 25;
            end

            initStr = sprintf('Writing Kilosort %d dataset to %s', ks.kilosort_version, outpath);
            if isempty(p.Results.progressInitializeFn) && isempty(p.Results.progressIncrementFn)
                prog = Neuropixel.Utils.ProgressBar(nProg, initStr);
                progIncrFn = @(text) prog.increment(text);
            else
                if ~isempty(p.Results.progressInitializeFn)
                    p.Results.progressInitializeFn(nProg, initStr);
                end
                progIncrFn = p.Results.progressIncrementFn;
            end

            if ~overwriteSelf
                progIncrFn('Copying params.py');
                new_dat_path = string(p.Results.updateParamsDotPyDatPath);
                if new_dat_path == ""
                    % leave dat_path alone in params.py
                    copyfile(oldp('params.py'), newp('params.py'));
                else
                    lines = readlines(oldp('params.py'));
                    dat_path_line = find(startsWith(lines, "dat_path"), 1, 'last');
                    if isempty(dat_path_line)
                        error('Could not find line beginning with dat_path in params.py');
                    end
                    lines(dat_path_line) = sprintf("dat_path = '%s'", new_dat_path);
                    writelines(lines, newp('params.py'));
                end
            end

            if ~isempty(ks.ops)
                progIncrFn('Copying ops.mat');
                ops = ks.ops;
                save(newp('ops.mat'), 'ops');
            end

            % write channelMap
            progIncrFn('Writing chanMap.mat');
            chanMapFile = fullfile(outpath, 'chanMap.mat');
            ks.writeChannelMap(chanMapFile);

            write(ks.amplitudes, 'amplitudes');
            write(ks.channel_ids_sorted - ones(1, 'like', ks.channel_ids_sorted), 'channel_map');
            write(ks.channel_positions_sorted, 'channel_positions');
            if ks.hasFeaturesLoaded
                write(ks.pc_features, 'pc_features');
                write(ks.pc_feature_ind - ones(1, 'like', ks.pc_feature_ind), 'pc_feature_ind'); % convert 1 indexed to 0 indexed
            end
            write(ks.similar_templates, 'similar_templates');

            write(ks.spike_templates - ones(1, 'like', ks.spike_templates), 'spike_templates');
            write(ks.spike_templates_preSplit - ones(1, 'like', ks.spike_templates), 'spike_templates_preSplit');

            write(ks.spike_times, 'spike_times');
            if ks.hasFeaturesLoaded
                write(ks.template_features, 'template_features');
                write(ks.template_feature_ind - ones(1, 'like', ks.template_feature_ind), 'template_feature_ind'); % 1 indexed to 0 indeded
            end

            % add back leading zeros stripped off of ks.templates based on size of W
            prepad_templates = ks.ops.nt0 - 2*ks.ops.nt0min - 1;
            szPad = [size(ks.templates, 1), prepad_templates, size(ks.templates, 3)];
            templates_padded = cat(2, zeros(szPad, 'like', ks.templates), ks.templates);
            write(templates_padded, 'templates');
            write(ks.templates_ind - ones(1, 'like', ks.templates_ind), 'templates_ind');
            write(ks.whitening_mat, 'whitening_mat');
            write(ks.whitening_mat_inv, 'whitening_mat_inv');
            write(ks.spike_clusters, 'spike_clusters');
            write(ks.spike_clusters_ks2orig, 'spike_clusters_ks2orig');

            % write the unique cluster ids as well
            write(ks.cluster_ids, 'unique_cluster_ids');

            if ks.kilosort_version == 2
                if ks.hasFeaturesLoaded
                    write(ks.W, 'template_W');
                    write(ks.U, 'template_U');
                    write(ks.mu, 'template_mu');

                    if ks.hasPreSplitLoaded
                        write(ks.W_preSplit, 'template_W_presplit');
                        write(ks.U_preSplit, 'template_U_presplit');
                        write(ks.mu_preSplit, 'template_mu_presplit');
                        write(ks.iW_preSplit, 'template_iW_presplit');
                    end
                end

                if ks.hasBatchwiseLoaded
                    write(ks.batchwise_cc , 'batchwise_ccb');
                    write(ks.batch_sort_order, 'batch_sort_order');
                    write(ks.batch_starts, 'batch_starts');

                    write(ks.W_batch, 'template_W_batch');

                    write(ks.W_batch_US, 'template_W_batch_US');
                    write(ks.W_batch_V, 'template_W_batch_V');
                    write(ks.U_batch, 'template_U_batch');
                    write(ks.U_batch_US, 'template_U_batch_US');
                    write(ks.U_batch_V, 'template_U_batch_V');
                    write(ks.mu_batch, 'template_mu_batch');

                    if ks.hasPreSplitLoaded % mostly used for reextracting spikes
                        write(ks.W_batch_preSplit, 'template_W_batch_presplit');
                        write(ks.U_batch_preSplit, 'template_U_batch_presplit');
                        write(ks.mu_batch_preSplit, 'template_mu_batch_presplit');
                    end
                end

                write(ks.cluster_est_contam_rate, 'cluster_est_contam_rate');

                if ~isempty(ks.cluster_merge_count) || ~isempty(ks.cluster_split_src)
                    progIncrFn('Writing splitMergeInfo.mat');
                    splitMergeInfo.mergecount = ks.cluster_merge_count;
                    splitMergeInfo.mergedst = ks.cluster_merge_dst;
                    splitMergeInfo.splitsrc = ks.cluster_split_src;
                    splitMergeInfo.splitdst = ks.cluster_split_dst;
                    splitMergeInfo.splitauc = ks.cluster_split_auc;
                    splitMergeInfo.split_candidate = ks.cluster_split_candidate;
                    splitMergeInfo.split_orig_template = ks.cluster_orig_template;
                    splitMergeInfo.split_projections = ks.cluster_split_projections;
                    save(newp('splitMergeInfo.mat'), 'splitMergeInfo');
                end

                if ks.hasCutoffLoaded
                    write(ks.cutoff_thresholds, 'cutoff_thresholds');
                    write(ks.cutoff_spike_times, 'cutoff_spike_times');
                    write(ks.cutoff_spike_templates - ones(1, 'like', ks.cutoff_spike_templates), 'cutoff_spike_templates');
                    write(ks.cutoff_spike_templates_preSplit - ones(1, 'like', ks.cutoff_spike_templates_preSplit), 'cutoff_spike_templates_preSplit');
                    write(ks.cutoff_amplitudes, 'cutoff_amplitudes');

                    write(ks.cutoff_spike_clusters, 'cutoff_spike_clusters');
                    write(ks.cutoff_spike_clusters_ks2orig, 'cutoff_spike_clusters_ks2orig');
                    if ks.hasFeaturesLoaded
                        write(ks.cutoff_pc_features, 'cutoff_pc_features');
                        write(ks.cutoff_template_features, 'cutoff_template_features');
                    end
                end

                % write cluster_KSLabel.tsv file - important so it is recognized as KS2
                writeClusterMetaTSV('cluster_KSLabel', 'KSLabel', ks.cluster_ks_label);
            end

            if exist('prog', 'var')
                prog.finish();
            end

            function write(data, file)
                if isempty(data)
                    return
                end
                progIncrFn(sprintf('Writing %s', file));
                ffile = fullfile(outpath, [file '.npy']);
                Neuropixel.writeNPY(data, ffile);
            end

            function writeClusterMetaTSV(file_noext, field, values)
                file = fullfile(outpath, [file_noext, '.tsv']);
                progIncrFn(sprintf('Writing %s', file));
                fid = fopen(file, 'w');
                assert(fid > 0, 'Error opening %s for writing', file);

                fprintf(fid, 'cluster_id\t%s\n', field);

                cluster_ids = ks.cluster_ids;
                values = string(values);
                for iC = 1:numel(cluster_ids)
                    fprintf(fid, '%d\t%s\n', cluster_ids(iC), values(iC));
                end

                fclose(fid);
            end
        end

        function [ks, rez] = copyReextractSpikesWithFixedTemplates(ks, varargin)
            p = inputParser();
            p.addParameter('rez_reextract', [], @(x) isstruct(x) || isempty(x));
            p.KeepUnmatched = true;
            p.parse(varargin{:});

            ks = copy(ks);

            % rerun Kilosort using ks's templates as a given
            if ~isempty(p.Results.rez_reextract)
                rez = p.Results.rez_reextract;
            else
                rez = Kilosort2.MainLoop.reextractSpikesWithFixedTemplates(ks, p.Unmatched);
            end

            % copy the updated fields from rez back into ks
            cluster_col = rez.st3_cluster_col;
            template_col = rez.st3_template_col;

            ks.spike_times = uint64(rez.st3(:, 1));
            ks.spike_templates_preSplit = uint32(rez.st3(:, 2));
            ks.amplitudes = rez.st3(:, 3);
            ks.spike_templates = uint32(rez.st3(:, template_col));

            % ensure that we subtract one to replicate what rez2phy does
            cluster_offset = -1;
            ks.spike_clusters = uint32(rez.st3(:, cluster_col) + cluster_offset);
            ks.spike_clusters_ks2orig = ks.spike_clusters;

            ks.template_features = rez.cProj;
            ks.pc_features = rez.cProjPC;

            ks.cutoff_spike_times = uint64(rez.st3_cutoff_invalid(:, 1));
            ks.cutoff_spike_templates_preSplit = uint32(rez.st3_cutoff_invalid(:, 2));
            ks.cutoff_amplitudes = rez.st3_cutoff_invalid(:, 3);
            ks.cutoff_spike_templates = uint32(rez.st3_cutoff_invalid(:, template_col));

            ks.cutoff_spike_clusters = uint32(rez.st3_cutoff_invalid(:, cluster_col) + cluster_offset);
            ks.cutoff_spike_clusters_ks2orig = ks.cutoff_spike_clusters;

            ks.cutoff_template_features = rez.cProj_cutoff_invalid;
            ks.cutoff_pc_features = rez.cProjPC_cutoff_invalid;
        end
    end
end
