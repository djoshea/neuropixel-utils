classdef KilosortDataset < handle
    % wrapper around a Kilosort dataset
    % todo - load cluster ratings from cluster_groups.tsv
    % Note 1: in the context of this file, time refers to samples, 1-indexed
    % Note 2: this file will differ from raw_dataset in nChannelsSorted. Here, nChannelsSorted means the number of channels
    %   in .channel_ids_sorted (which will match the other properties)

    properties
        path(1, :) char

        raw_dataset % Neuropixel.ImecDataset instance

        channelMap % Neuropixel.ChannelMap

        fsAP % sampling rate pass thru to raw_dataset or specified during construction

        apGain
        apScaleToUv % multiply raw samples by this to get uV

        meta % ap metadata loaded

        concatenationInfo

        isLoaded logical = false;
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

        nChannelsSorted
        nSpikes
        nClusters
        nTemplates
        nTemplateRank
        nPCFeatures
        nFeaturesPerChannel
        nTemplateTimepoints
        templateTimeRelative

        nBatches
        nSpikesCutoff
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
        channel_positions_sorted(:, 2) double

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

        cluster_spike_counts(:, 1) uint32 % nClusters x 1 number of spikes assigned to each cluster

        % cluster_groups - "cluster group" of each cluster as set in Phy
        cluster_groups(:, 1) categorical

        % custom cluster_ratings completely divorced from Phy
        cluster_ratings

        % unique clusters in spike_clusters [nClusters]
        cluster_ids (:, 1) uint32
    end

    properties % Secondary data, loaded from rez.mat for Kilosort2 only
        ops struct % options used by Kilosort

        kilosort_version (1, :) uint32 % scalar or vector (leading digit is 1 or 2 for Kilosort1 or Kilosort2)

        % low-rank decomposition of templates into spatial and temporal components
        W (:, :, :) double % [nTemplateTimepoints, nTemplates, nTemplateRank] - temporal decomposition of templates

        U (:, :, :) double % [nChannelsSorted, nTemplates, nTemplateRank] - spatial decomposition of templates

        % batch information used in drift correction
        batch_starts (:, 1) uint64 % [nBatches] timesteps where each batch begins (unsorted)

        batchwise_cc (:, :) single % (rez.ccb) [nBatches, nBatches] batch-wise cross-correlation in original data ordering

        batch_sort_order (:, 1) uint32 % (rez.iorig) [nBatches] sorted batch indices used in KS2

        % batchwise template data
        W_batch (:, :, :, :) single % [nTemplateTimepoints, nTemplates, nTemplateRank, nBatches] - full spatial templates by batch

        % per-template svd of W_batch U*S --> W_batch_a, V --> W_batch_b
        W_batch_US (:, :, :, :) single % [nTemplateTimepoints, nTemplateRank, nTemplatePCs, nTemplates] - reshaped from rez.W_a
        W_batch_V (:, :, :) single % [nBatches, nTemplatePCs, nTemplates] % rez.W_b

        U_batch (:, :, :, :) single  % [nChannelsSorted, nTemplates, nTemplateRank, nBatches] - full temporal templates

        % per-template svd of U_batch U*S --> U_batch_a, V --> U_batch_b
        U_batch_US (:, :, :, :) single % [nChannelsSorted, nTemplateRank, nTemplatePCs, nTemplates]
        U_batch_V (:, :, :) single % [nBatches, nTemplatePCs, nTemplates]

        % good vs. mua label as assigned by Kilosort2 in cluster_KSLabel.tsv (good, mua)
        cluster_ks_label(:, 1) categorical % [nClusters]

        cluster_est_contam_rate (:, 1) double % [nClusters]

        % optional split and merge meta information (only on djoshea branch of Kilosort2)
        cluster_merge_count (:, 1) uint32
        cluster_split_src (:, 1) uint32
        cluster_split_dst (:, 1) uint32
        cluster_split_auc (:, 1) double
        cluster_split_candidate (:, 1) logical
        cluster_orig_template (:, 1) uint32

        % [nSpikesInvalid, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        cutoff_spike_times(:, 1) uint64

        cutoff_amplitudes(:, 1) double

        % [nSpikesCutoff, ] uint32 vector specifying the identity of the template that was used to extract each spike
        cutoff_spike_templates(:, 1) uint32
        cutoff_spike_clusters(:, 1) uint32
        cutoff_spike_clusters_ks2orig(:, 1) uint32

        cutoff_cluster_spike_counts(:, 1) uint32 % nClusters x 1 number of spikes assigned to each cluster

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
            nBatches = numel(ks.batch_starts);
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
                apScaleToUv = ks.apScaleToUv;k
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
                        ks.raw_dataset = Neuropixel.ImecDataset(raw_path, 'channelMap', p.Results.channelMap);
                    end
                end
            end

            if isempty(ks.raw_dataset)
                warning('No ImecDataset found in Kilosort path, specify imecDataset parameter directly');
            end

            % these will pass thru to raw_dataset if provided
            ks.fsAP = p.Results.fsAP;
            if isempty(ks.fsAP) || isnan(ks.fsAP) % will pass thru to raw_dataset if found
                % try reading sample_rate from params.py
                ks.readParamsPy();
                ks.fsAP = ks.sample_rate;
            end

            ks.apScaleToUv = p.Results.apScaleToUv;

            % manually specify some additional props
            channelMap = p.Results.channelMap;
            if isempty(channelMap)
                if ~isempty(ks.raw_dataset)
                    channelMap = ks.raw_dataset.channelMap;
                end
                if isempty(channelMap)
                    channelMap = Neuropixel.Utils.getDefaultChannelMapFile(true);
                end
            end

            if ischar(channelMap)
                channelMap = Neuropixel.ChannelMap(channelMap);
            end
            ks.channelMap = channelMap;

            if ~isempty(ks.raw_dataset)
                ks.meta = ks.raw_dataset.readAPMeta();
                ks.concatenationInfo = ks.raw_dataset.concatenationInfoAP;
            end
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
            if any(s.fr_cutoff_only) > 0
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
        end

        function load(ks, varargin)
            p = inputParser();
            p.addOptional('reload', false, @islogical);
            p.addParameter('loadBatchwise', true, @islogical);
            p.addParameter('loadFeatures', true, @islogical);
            p.addParameter('loadCutoff', true, @islogical);
            p.addParameter('progressInitializeFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(nUpdates) to print update
            p.addParameter('progressIncrementFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(updateString) to print update
            p.parse(varargin{:});

            reload = p.Results.reload;
            if ks.isLoaded && ~reload
                return;
            end
            ks.isLoaded = false;

            ks.readParamsPy();

            loadFeatures = p.Results.loadFeatures;
            loadBatchwise = p.Results.loadBatchwise;
            loadCutoff = p.Results.loadCutoff;

            path = ks.path;
            existp = @(file) exist(fullfile(path, file), 'file') > 0;

            nProg = 15;
            has_ccb = existp('batchwise_ccb.npy');
            if has_ccb
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
            ks.spike_clusters = read('spike_clusters');

            if existp('spike_clusters_ks2orig.npy')
                ks.spike_clusters_ks2orig = read('spike_clusters_ks2orig');
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
            ks.cluster_ids = unique(ks.spike_clusters);

            progIncrFn('Computing cluster spike counts');
            assert(~isempty(ks.cluster_ids));
            ks.cluster_spike_counts = uint32(Neuropixel.Utils.discrete_histcounts(ks.spike_clusters, ks.cluster_ids));

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

            % load cluster ratings from disk
            ks.cluster_ratings = repmat(categorical("unrated"), numel(ks.cluster_ids), 1);
            if existp('cluster_Rating.tsv')
                tbl = readClusterMetaTSV('cluster_Rating.tsv', 'Rating', 'categorical');
                [tf, ind] = ismember(tbl.cluster_id, ks.cluster_ids);
                ks.cluster_rating(ind(tf)) = tbl{tf, 2};
            end

            if existp('ops.mat')
                temp = load(fullfile(p, 'ops.mat'), 'ops');
                ks.ops = temp.ops;
            end

            if ks.kilosort_version == 2
                ks.W = readOr('template_W');
                ks.U = readOr('template_U');

                % strip leading zeros off of ks.templates based on size of W
                nTimeTempW = size(ks.W, 1);
                nStrip = size(ks.templates, 2) - nTimeTempW;
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
                end

                ks.cluster_est_contam_rate = readOr('cluster_est_contam_rate');
                ks.cluster_merge_count = readOr('cluster_mergecount');
                ks.cluster_split_src = readOr('cluster_splitsrc');
                ks.cluster_split_dst = readOr('cluster_splitdst');
                ks.cluster_split_auc = read('cluster_splitauc');
                ks.cluster_split_candidate = read('cluster_split_candidate');
                ks.cluster_orig_template = read('cluster_split_orig_template');

                if loadCutoff
                    ks.cutoff_spike_times = readOr('cutoff_spike_times');
                    ks.cutoff_spike_templates = readOr('cutoff_spike_templates') + ones(1, 'like', ks.spike_templates); % 0 indexed templates to 1 indexed templates
                    ks.cutoff_amplitudes = readOr('cutoff_amplitudes');
                    ks.cutoff_spike_clusters = readOr('cutoff_spike_clusters'); % rez.st3_cutoff is 1 indexed, cluster ids are 0 indexed
                    if existp('cutoff_spike_clusters_ks2orig.npy')
                        ks.cutoff_spike_clusters_ks2orig = read('cutoff_spike_clusters_ks2orig');
                    else
                        ks.cutoff_spike_clusters_ks2orig = ks.cutoff_spike_clusters;
                    end

                    progIncrFn('Computing cluster cutoff spike counts');
                    ks.cutoff_cluster_spike_counts = uint32(Neuropixel.Utils.discrete_histcounts(ks.cutoff_spike_clusters, ks.cluster_ids));

                    if loadFeatures
                        ks.cutoff_pc_features = readOr('cutoff_pc_features');
                        ks.cutoff_template_features = readOr('cutoff_template_features');
                    end
                end
            end

            if exist('prog', 'var')
                prog.finish();
            end
            ks.metrics = []; % old metrics no longer valid
            ks.isLoaded = true;

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
                copyfile(src_path, dest_path);
            end
        end

        function save_spike_clusters_to_disk(ks)
            writeNPY_local = @(v, fname) writeNPY(v, fullfile(ks.path, fname));
            writeNPY_local(ks.spike_clusters, 'spike_clusters.npy');
            writeNPY_local(ks.cutoff_spike_clusters, 'cutoff_spike_clusters.npy');
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
            if isempty(ks.metrics) || ~isvalid(ks.metrics) || (nargin >= 2 && recompute)
                ks.load();
                ks.metrics = Neuropixel.KilosortMetrics(ks, varargin{:});
            end
            m = ks.metrics;
        end

        function [clusterInds, cluster_ids] = lookup_clusterIds(ks, cluster_ids)
            if islogical(cluster_ids)
                cluster_ids = ks.cluster_ids(cluster_ids);
             end
            [tf, clusterInds] = ismember(cluster_ids, ks.cluster_ids);
            assert(all(tf), 'Some cluster ids were not found in ks.clusterids');
        end

        function  [sortedChannelInds, channelIds] = lookup_sortedChannelIds(ks, channelIds)
             if islogical(channelIds)
                channelIds = ks.channel_ids_sorted(channelIds);
             end
            [tf, sortedChannelInds] = ismember(channelIds, ks.channel_ids_sorted);
            assert(all(tf), 'Some channel ids not found');
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
             p.addParameter('spike_idx', [], @isvector); % manually specify which idx into spike_times
             p.addParameter('spike_times', [], @isvector); % manually specify which times directly to extract
             p.addParameter('cluster_ids', [], @isvector); % manually specify all spikes from specific cluster_ids

             % and ONE OR NONE of these to pick channels (or channels for each cluster)
             p.addParameter('channel_ids_by_cluster', [], @(x) isempty(x) || ismatrix(x));
             p.addParameter('best_n_channels', NaN, @isscalar); % or take the best n channels based on this clusters template when cluster_id is scalar

             % if specified, spikes times will be filtered within this window
             p.addParameter('filter_window', [], @(x) isempty(x) || isvector(x));

             % other params:
             p.addParameter('num_waveforms', Inf, @isscalar); % caution: Inf will request ALL waveforms in order (typically useful if spike_times directly specified)
             p.addParameter('window', [-40 41], @isvector); % Number of samples before and after spiketime to include in waveform
             p.addParameter('car', false, @islogical);
             p.addParameter('centerUsingFirstSamples', 20, @(x) isscalar(x) || islogical(x)); % subtract mean of each waveform's first n samples, don't do if false

             p.addParameter('subtractOtherClusters', false, @islogical); % time consuming step to remove the contribution of the other clusters to a given snippet

             p.addParameter('raw_dataset', ks.raw_dataset, @(x) true);
             p.addParameter('fromSourceDatasets', false, @islogical); % go all the way back to the imecDatasets that were concatenated to form ks.raw_dataset

             % other metadata set in snippetSet
             p.addParameter('trial_idx', [], @isvector);

             p.addParameter('from_cutoff_spikes', false, @islogical);

             p.parse(varargin{:});

             assert(ks.hasRawDataset, 'KilosortDataset has no raw ImecDataset');

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
                     end
                 else
                     if isscalar(clu)
                         spike_idx = find(ks.cutoff_spike_clusters == clu);
                         spike_times = ks.cutoff_spike_times(spike_idx);
                     else
                         mask = ismember(ks.cutoff_spike_clusters, clu);
                         spike_times = ks.cutoff_spike_times(mask);
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

             % apply mask
             spike_idx = spike_idx(mask);
             spike_times = spike_times(mask);

             % fetch other info
             if ~from_cutoff_spikes
                cluster_ids = ks.spike_clusters(spike_idx);
             else
                cluster_ids = ks.cutoff_spike_clusters(spike_idx);
             end
             if isempty(unique_cluster_ids)
                 unique_cluster_ids = unique(cluster_ids);
             end

             trial_idx = p.Results.trial_idx;
             if ~isempty(trial_idx)
                 trial_idx = trial_idx(mask);
                 assert(numel(trial_idx) == numel(spike_idx));
             end

             % take max of num_waveforms from each cluster
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
                         mask(thisC(randsample(numel(thisC), nSample, false))) = true;
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
             elseif ~isempty(p.Results.channel_ids_by_cluster)
                 channel_ids_by_cluster = p.Results.channel_ids_by_cluster;
             else
                 channel_ids_by_cluster = ks.channel_ids_sorted;
             end

             % channel_ids is provided since raw data often has additional channels that we're not interested in
             window = p.Results.window;
             snippetSet = p.Results.raw_dataset.readAPSnippetSet(spike_times, ...
                 window, 'channel_ids_by_cluster', channel_ids_by_cluster, ...
                 'unique_cluster_ids', unique_cluster_ids, 'cluster_ids_by_snippet', cluster_ids, ...
                 'car', p.Results.car, 'fromSourceDatasets', p.Results.fromSourceDatasets);
             snippetSet.trial_idx = trial_idx;
             snippetSet.ks = ks;

             if p.Results.subtractOtherClusters
                 reconstructionFromOtherClusters = ks.reconstructSnippetSetFromTemplates(snippetSet, ...
                     'excludeClusterFromOwnReconstruction', true);
                 snippetSet.data = snippetSet.data - reconstructionFromOtherClusters.data;
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
            p.addParameter('showPlots', false, @islogical);
            p.addParameter('rawData', [], @isnumeric); % for plotting only

            p.addParameter('use_batchwise_templates', ~isempty(ks.W_batch), @islogical);
            p.parse(varargin{:});
            showPlots = p.Results.showPlots;

            % check sizes of everything
            nTimes = numel(times);
            channel_ids_by_snippet = p.Results.channel_ids_by_snippet;
            assert(~isempty(channel_ids_by_snippet));
            if size(channel_ids_by_snippet, 2) == 1
                channel_ids_by_snippet = repmat(channel_ids_by_snippet, 1, numel(times));
            end
            assert(numel(times) == size(channel_ids_by_snippet, 2));

            nChannelsSorted = size(channel_ids_by_snippet, 1);
            unique_cluster_ids = p.Results.unique_cluster_ids;

            cluster_ids = p.Results.cluster_ids;
            assert(numel(cluster_ids) == nTimes);

            exclude_cluster_ids_each_snippet = p.Results.exclude_cluster_ids_each_snippet;
            exclude_cluster_ids_all_snippets = p.Results.exclude_cluster_ids_all_snippets;

            % templates post-whitening is nTemplates x nTimepoints x nChannelsFull
            use_batchwise_templates = p.Results.use_batchwise_templates;
            if use_batchwise_templates && (isempty(ks.W_batch) || isempty(ks.U_batch))
                warning('Cannot use batchwise tempaltes as W_batch or U_batch was not loaded successfully from rez.mat');
                use_batchwise_templates = false;
            end

            if ~use_batchwise_templates
                metrics = ks.computeMetrics();
                templates =  metrics.template_unw; % unwhitened templates, but not scaled and still in quantized units (not uV)
            end

            relTvec_template = ks.templateTimeRelative;
            relTvec_snippet = int64(window(1):window(2));
            reconstruction = zeros(nChannelsSorted, numel(relTvec_snippet), nTimes, 'int16');

            if exist('ProgressBar', 'class') == 8
                prog = ProgressBar(nTimes, 'Reconstructing templates around snippet times');
            else
                prog = Neuropixel.Utils.ProgressBar(nTimes, 'Reconstructing templates around snippet times');
            end

            for iT = 1:nTimes
                prog.update(iT);

                reconstruction_this = zeros(nChannelsSorted, numel(relTvec_snippet), 'single');

                % find spikes that would lie within this window (with padding),
                % excluding those from clusters we wish to exclude
                t = int64(times(iT));
                minT = t + window(1) - int64(relTvec_template(end));
                maxT = t + window(2) - int64(relTvec_template(1));
                if iscell(exclude_cluster_ids_each_snippet)
                    exclude_this = exclude_cluster_ids_each_snippet{iT};
                elseif ~isempty(exclude_cluster_ids_each_snippet)
                    exclude_this = exclude_cluster_ids_each_snippet(iT);
                else
                    exclude_this = [];
                end
                exclude_this = union(exclude_this, exclude_cluster_ids_all_snippets);

                nearby_spike_inds = find(ks.spike_times >= minT & ks.spike_times <= maxT);
                %nearby_spike_inds = find(ks.spike_times == t);

                nearby_spike_inds(ismember(ks.spike_clusters(nearby_spike_inds), exclude_this)) = [];

                if use_batchwise_templates
                    spike_batches = ks.compute_which_batch(ks.spike_times(nearby_spike_inds));
                end

                % figure out what channels we need to reconstruct onto
                cluster_ids_this = cluster_ids(iT);
                [~, cluster_ind_this] = ismember(cluster_ids_this, unique_cluster_ids);
                assert(cluster_ind_this > 0, 'Cluster for times(%d) not found in unique_cluster_ids', iT);
                channel_ids_this = channel_ids_by_snippet(:, iT);
                sorted_channel_inds_this = ks.lookup_sortedChannelIds(channel_ids_this);

                if showPlots
                    clf;
                    plot(relTvec_snippet, p.Results.rawData(1, :, iT), 'k-', 'LineWidth', 2);
                    hold on;
                end

                % loop over the enarby sp
                for iS = 1:numel(nearby_spike_inds)
                    ind = nearby_spike_inds(iS);
                    amp = ks.amplitudes(ind);

                    % figure out time overlap and add to reconstruction
                    tprime = int64(ks.spike_times(ind));
                    indFromTemplate = relTvec_template + tprime >= t + relTvec_snippet(1) & relTvec_template + tprime <= t + relTvec_snippet(end);
                    indInsert = relTvec_snippet + t >= relTvec_template(1) + tprime & relTvec_snippet + t <= relTvec_template(end) + tprime;

                    if use_batchwise_templates
                        template_this = ks.construct_batchwise_templates(ks.spike_templates(ind), 'batches', spike_batches(iS));
                        template_this = template_this(:, indFromTemplate, sorted_channel_inds_this);
                    else
                        template_this = templates(ks.spike_templates(ind), indFromTemplate, sorted_channel_inds_this);
                    end
                    insert = amp .* permute(template_this, [3 2 1]);
                    reconstruction_this(:, indInsert) = reconstruction_this(:, indInsert) + insert;

                    if showPlots
                        if tprime == t
                            args = {'LineWidth', 1, 'Color', 'g'};
                        else
                            args = {};
                        end
                        plot(relTvec_template(indInsert), amp .* template_this(:, :, 1), 'Color', [0.7 0.7 0.7], args{:});
                    end
                end

                if showPlots
                    plot(relTvec_snippet, reconstruction(1, :, iT), 'r--', 'LineWidth', 2);
                    plot(relTvec_snippet, p.Results.rawData(1, :, iT) - int16(reconstruction(1, :, iT)), 'b--');
                    pause;
                end

                reconstruction(:, :, iT) = int16(reconstruction_this);
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
            batch_ind = discretize(sample_idx, [ks.batch_starts; Inf]);
            [~, batch_sorted_location] = ismember((1:ks.nBatches)', ks.batch_sort_order);
            sorted_batch_ind = batch_sorted_location(batch_ind);
        end

        function data = construct_batchwise_templates(ks, template_inds, varargin)
            p = inputParser();
            p.addParameter('batches', 1:ks.nBatches, @isvector);
            p.addParameter('whitened', true, @islogical);
            p.addParameter('scaleByBatchAmplitude', false, @islogical); % multiply by muA to scale according to batch, dont use if you plan to scale by amplitudes
            p.parse(varargin{:});
            batch_inds = p.Results.batches;
            scaleByBatchAmplitude = p.Results.scaleByBatchAmplitude;

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

            for iT = 1:nTemplates
                for iB = 1:nBatches
                    if scaleByBatchAmplitude
                        amp = ks.muA(template_inds(iT), batch_inds(iB));
                    else
                        amp = 1;
                    end

                	W = ks.W_batch(:, template_inds(iT), :, batch_inds(iB)); % nTT x 1 x 3 x 1
                    U = ks.U_batch(:, template_inds(iT), :, batch_inds(iB)); % nCh x 1 x 3 x 1

                    for iP = 1:size(U, 3)
                        this_pc = amp * W(:, :, iP, :) * (U(:, :, iP, :)' * wmi); % nTT x nCh
                        data(iT, :, :, iB) = data(iT, :, :, iB) + reshape(this_pc, [1, nTT, nCh]);
                    end
                end
            end
        end

    end
end
