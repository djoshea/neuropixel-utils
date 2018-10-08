classdef KiloSortDataset < handle
    % wrapper around a KiloSort dataset
    % todo - load cluster ratings from cluster_groups.tsv
    % Note 1: in the context of this file, time refers to samples, 1-indexed
    % Note 2: this file will differ from raw_dataset in nChannels. Here, nChannels means the number of channels
    %   in .channel_map (which will match the other properties)

    properties
        path(1, :) char

        raw_dataset % Neuropixel.ImecDataset instance

        cluster_best_template_channels % nClusters x nClosest, indexed as in data file
        
        loadSync logical = false; % if false, won't load sync channel on load()
    end

    properties(Dependent)
        pathLeaf
        isLoaded
        hasRawDataset
        nSpikes
        nChannelsSorted % number of channels in channel_map
        nClusters
        nTemplates
        nPCFeatures
        nFeaturesPerChannel
        nChannels
        channelMap
        
        % pass through to sync information in raw_dataset if present
        sync(:, 1) uint16 % nSamples sync channel contents
    end
    
    properties(Dependent, SetAccess=protected)
        syncBitNames(:, 1) cell
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

        % channel_map.npy - [nChannels, ] uint32 vector with the channel map, i.e. which row of the data file to look in for the channel in question
        channel_map(:, 1) uint32

        % channel_positions.npy - [nChannels, 2] double matrix with each row giving the x and y coordinates of that channel. Together with the channel map, this determines how waveforms will be plotted in WaveformView (see below).
        channel_positions(:, 2) double

        % pc_features.npy - [nSpikes, nFeaturesPerChannel, nPCFeatures] single matrix giving the PC values for each spike.
        % The channels that those features came from are specified in pc_features_ind.npy. E.g. the value at pc_features[123, 1, 5]
        % is the projection of the 123rd spike onto the 1st PC on the channel given by pc_feature_ind[5].
        pc_features(:, :, :) single

        % pc_feature_ind.npy - [nTemplates, nPCFeatures] uint32 matrix specifying which pcFeatures are included in the pc_features matrix.
        pc_feature_ind(:, :) uint32

        % similar_templates.npy - [nTemplates, nTemplates] single matrix giving the similarity score (larger is more similar) between each pair of templates
        similar_templates(:, :) single

        % spike_templates.npy - [nSpikes, ] uint32 vector specifying the identity of the template that was used to extract each spike
        spike_templates(:, 1) uint32

        % spike_times.npy - [nSpikes, ] uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.
        spike_times(:, 1) uint64

        % template_features.npy - [nSpikes, nTempFeatures] single matrix giving the magnitude of the projection of each spike onto nTempFeatures other features.
        % Which other features is specified in template_feature_ind.npy
        template_features(:, :) single

        % template_feature_ind.npy - [nTemplates, nTempFeatures] uint32 matrix specifying which templateFeatures are included in the template_features matrix.
        template_feature_ind(:, :) uint32

        % templates.npy - [nTemplates, nTimePoints, nTempChannels] single matrix giving the template shapes on the channels given in templates_ind.npy
        templates(:, :, :) single

        % templates_ind.npy - [nTemplates, nTempChannels] double matrix specifying the channels on which each template is defined.
        % In the case of Kilosort templates_ind is just the integers from 0 to nChannels-1, since templates are defined on all channels.
        templates_ind(:, :) double

        % whitening_mat.npy - [nChannels, nChannels] double whitening matrix applied to the data during automatic spike sorting
        whitening_mat(:, :) double

        % whitening_mat_inv.npy - [nChannels, nChannels] double, the inverse of the whitening matrix.
        whitening_mat_inv(:, :) double

        % spike_clusters.npy - [nSpikes, ] int32 vector giving the cluster identity of each spike. This file is optional and
        % if not provided will be automatically created the first time you run the template gui, taking the same values as
        % spike_templates.npy until you do any merging or splitting.
        spike_clusters(:, 1) int32

        % cluster_groups - comma-separated value text file giving the "cluster group" of each cluster (0=noise, 1=MUA, 2=Good, 3=unsorted)
        cluster_groups
        
        % unique clusters in spike_clusters [nClusters
        cluster_ids (:, 1) int32
    end

    methods % Dependent properties
        function leaf = get.pathLeaf(ds)
            [~, leaf] = fileparts(ds.path);
        end
        
        function tf = get.isLoaded(ds)
            tf = ~isempty(ds.spike_times);
        end
        
        function tf = get.hasRawDataset(ds)
            tf = ~isempty(ds.raw_dataset);
        end
        
        function n = get.nSpikes(ds)
            if isempty(ds.spike_times)
                n = NaN;
            else
                n = size(ds.spike_times, 1);
            end
        end

        function n = get.nChannelsSorted(ds)
            if isempty(ds.channel_map)
                n = NaN;
            else
                n = size(ds.channel_map, 1);
            end
        end

        function n = get.nClusters(ds)
            if isempty(ds.cluster_ids)
                n = NaN;xi
            else
                n = numel(ds.cluster_ids);
            end
        end

        function n = get.nTemplates(ds)
            if isempty(ds.pc_feature_ind)
                n = NaN;
            else
                n = size(ds.pc_feature_ind, 1);
            end
        end

        function n = get.nPCFeatures(ds)
            if isempty(ds.pc_features)
                n = NaN;
            else
                n = size(ds.pc_features, 3);
            end
        end

        function n = get.nFeaturesPerChannel(ds)
            if isempty(ds.pc_features)
                n = NaN;
            else
                n = size(ds.pc_features, 2);
            end
        end
        
        function n = get.nChannels(ds)
            if isempty(ds.channel_map)
                n = NaN;
            else
                n = numel(ds.channel_map);
            end
        end
        
        function map = get.channelMap(ds)
            map = ds.raw_dataset.channelMap;
        end
        
        function sync = get.sync(ds)
            if isempty(ds.raw_dataset) || ~ds.isLoaded || ~ds.loadSync
                % don't autoload on print out make sure this only gets
                % called by load
                sync = zeros(0, 1, 'uint16');
            else
                sync = ds.raw_dataset.readSyncChannel();
            end
        end
        
        function syncBitNames = get.syncBitNames(ds)
            if isempty(ds.raw_dataset)
                if isempty(ds.syncBitNames)
                    syncBitNames = repmat({''}, 16, 1);
                else
                    syncBitNames = ds.syncBitNames;
                end
            else
                syncBitNames = ds.raw_dataset.syncBitNames;
            end
        end
    end

    methods
        function ds = KiloSortDataset(path, varargin)
            if nargin == 0
                return;
            end
            
            p = inputParser();
            p.addParameter('channelMap', Neuropixel.Utils.getDefaultChannelMapFile(), @(x) true);
            p.addParameter('imecDataset', [], @(x) true);
            p.addParameter('loadSync', false, @islogical);
            p.parse(varargin{:});
            
            if isa(path, 'Neuropixel.ImecDataset')
                ds.raw_dataset = path;
                path = ds.raw_dataset.pathRoot;
            end
            
            assert(exist(path, 'dir') == 7, 'Path %s not found', path);
            ds.path = path;
            
            ds.loadSync = p.Results.loadSync; % if false, won't load sync channel on load()
            
            if isempty(ds.raw_dataset)
                if isa(p.Results.imecDataset, 'Neuropixel.ImecDataset')
                    ds.raw_dataset = p.Results.imecDataset;
                elseif ischar(p.Results.imecDataset) || isempty(p.Results.imecDataset)
                    raw_path = p.Results.imecDataset;
                    if isempty(raw_path), raw_path = path; end
                    if Neuropixel.ImecDataset.folderContainsDataset(raw_path)
                        ds.raw_dataset = Neuropixel.ImecDataset(raw_path, 'channelMap', p.Results.channelMap);
                    end
                end
            end
        end
        
        function s = computeBasicStats(ds, varargin)
            % a way of computing basic statistics without fully loading
            % from disk. see printBasicStats
            
            p = inputParser();
            p.addParameter('frThresholdHz', 3, @isscalar);
            p.parse(varargin{:});

            % can work without loading raw data, much faster
            if ds.isLoaded
                s.spike_clusters = ds.spike_clusters; %#ok<*PROP>
                s.cluster_ids = ds.cluster_ids; 
                s.offset = ds.offset;
                s.sample_rate = ds.sample_rate;
                s.spike_times = ds.spike_times;
            else
                % partial load
                s.spike_clusters = read('spike_clusters');
                s.cluster_ids = unique(s.spike_clusters);
                params = Neuropixel.readINI(fullfile(ds.path, 'params.py'));
                s.sample_rate = params.sample_rate;
                s.offset = params.offset;
                s.spike_times = read('spike_times');
            end
            
            s.nSpikes = numel(s.spike_clusters);
            s.nClusters = numel(s.cluster_ids);
            s.nSec = double(max(s.spike_times) - s.offset) / double(s.sample_rate);
            
            counts = histcounts(s.spike_clusters, sort(s.cluster_ids));
            s.fr = counts / double(s.nSec);
            s.thresh = p.Results.frThresholdHz;
            
            s.clusterMask = s.fr > s.thresh;
            s.nClustersAboveThresh = nnz(s.clusterMask);
            s.nSpikesAboveThresh = nnz(ismember(s.spike_clusters, s.cluster_ids(s.clusterMask)));
            
            function out = read(file)
                out = Neuropixel.readNPY(fullfile(ds.path, [file '.npy']));
            end
        end
        
        function printBasicStats(ds, varargin)
            s = ds.computeBasicStats(varargin{:});
            
            fprintf('%s: %.1f sec, %d (%d) spikes, %d (%d) clusters (with fr > %g Hz)\n', ds.pathLeaf, s.nSec, ...
                s.nSpikes, s.nSpikesAboveThresh, ...
                s.nClusters, s.nClustersAboveThresh, s.thresh);
            
            plot(sort(s.fr, 'descend'), 'k-')
            xlabel('Cluster');
            ylabel('# spikes');
            hold on
            horzLine(s.thresh, 'Color', 'r');
            hold off;
            box off;
            grid on;
        end

        function load(ds, reload)
            if nargin < 2
                reload = false;
            end
            if ds.isLoaded && ~reload
                return;
            end
                
            params = Neuropixel.readINI(fullfile(ds.path, 'params.py'));
            ds.dat_path = params.dat_path;
            ds.n_channels_dat = params.n_channels_dat;
            ds.dtype = params.dtype;
            ds.offset = params.offset;
            ds.sample_rate = params.sample_rate;
            ds.hp_filtered = params.hp_filtered;

            p = ds.path;
            prog = ProgressBar(15, 'Loading KiloSort dataset: ');

            ds.amplitudes = read('amplitudes');
            ds.channel_map = read('channel_map');
            ds.channel_map = ds.channel_map + ones(1, 'like', ds.channel_map); % make channel map 1 indexed
            ds.channel_positions = read('channel_positions');
            ds.pc_features = read('pc_features');
            ds.pc_feature_ind = read('pc_feature_ind');
            ds.similar_templates = read('similar_templates');
            ds.spike_templates = read('spike_templates');
            ds.spike_templates = ds.spike_templates + ones(1, 'like', ds.spike_templates);
            ds.spike_times = read('spike_times');
            ds.template_features = read('template_features');
            ds.template_feature_ind = read('template_feature_ind');
            ds.templates = read('templates');

            ds.templates_ind = read('templates_ind');
            ds.templates_ind = ds.templates_ind + ones(1, 'like', ds.templates_ind); % convert plus zero indexed to 1 indexed

            ds.whitening_mat = read('whitening_mat');
            ds.whitening_mat_inv = read('whitening_mat_inv');
            ds.spike_clusters = read('spike_clusters');
            ds.cluster_ids = unique(ds.spike_clusters);
            
            % do auto loading now
            ds.sync;

            prog.finish()

            function out = read(file)
                prog.increment('Loading KiloSort dataset: %s', file);
                out = Neuropixel.readNPY(fullfile(p, [file '.npy']));
            end
        end

        function checkLoaded(ds)
            if isempty(ds.amplitudes)
                ds.load();
            end
        end
        
        function setSyncBitNames(ds, idx, names)
            if ~isempty(ds.raw_dataset)
                ds.raw_dataset.setSyncBitNames(idx, names);
            else
                if isscalar(idx) && ischar(names)
                    ds.syncBitNames{idx} = names;
                else
                    assert(iscellstr(names)) %#ok<ISCLSTR>
                    ds.syncBitNames(idx) = names;
                end
            end
        end
        
        function idx = lookupSyncBitByName(ds, names)
            if ischar(names)
                names = {names};
            end
            assert(iscellstr(names));

            [tf, idx] = ismember(names, ds.syncBitNames);
            idx(~tf) = NaN;
        end

        function seg = segmentIntoClusters(ds)
            % generates a KiloSortTrialSegmentedDataset with only 1 trial.
            % Segments as though there is one all-inclusive trials, so that
            % clusters are split but not trials
            tsi = Neuropixel.TrialSegmentationInfo(1);

            tsi.trialId = 0;
            tsi.conditionId = 0;
            tsi.idxStart = 1;
            
            if ds.hasRawDataset
                tsi.idxStop = ds.raw_dataset.nSamplesAP;
            else
                tsi.idxStop = max(ds.spike_times) + 1;
            end

            seg = ds.segmentIntoTrials(tsi, 0);
        end

        function seg = segmentIntoTrials(ds, tsi, trial_ids)

            % trialInfo hs fields:
            %   trialId
            %   conditionId
            %   idxStart
            %   idxStop

            seg = Neuropixel.KiloSortTrialSegmentedDataset(ds, tsi, trial_ids);
        end

        function computeBestChannelsForClusters(ds)
            % set nClosest to nChannelsSorted so that we have all the
            % distances sorted ahead of time

            % split spike_templates by cluster idspike_templates_by_unit
            ds.checkLoaded();
            [~, unit_idx_each_spike] = ismember(ds.spike_clusters, ds.cluster_ids);
            spike_templates_by_unit = TensorUtils.splitAlongDimensionBySubscripts(...
                ds.spike_templates, 1, ds.nClusters, unit_idx_each_spike);

            % nChannelsSorted x nChannelsSorted, include each channel in its own
            % closest list
            closest_lookup = [ds.channel_map, ds.channelMap.getClosestConnectedChannels(ds.nChannelsSorted-1, ds.channel_map)];

            best_template_channels = nan(ds.nClusters, ds.nChannelsSorted);
            for iU = 1:ds.nClusters
                template_id = mode(spike_templates_by_unit{iU});

                % numel_template_ind
                template = ds.templates(template_id, :, :);
                template = ds.unwhiten_templates(template);
                [~, bestChannelSortedInd] = max(range(template, 2));

                best_template_channels(iU, :) = closest_lookup(bestChannelSortedInd, :);
            end
            ds.cluster_best_template_channels = best_template_channels;
        end

        function templates = unwhiten_templates(ds, templates)
            % templates is [nTemplates, nTimePoints, nChannelsSorted]
            % TODO if templates not on all channels as they are with kilosort, this function must be revised

            wmi = ds.whitening_mat_inv;
            assert(size(wmi, 1) == size(templates, 3), 'dim 3 of templates must match whitening matrix inverse');
            templates = TensorUtils.linearCombinationAlongDimension(templates, 3, wmi');
        end

        function snippetSet = getWaveformsFromRawData(ds, varargin)
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

             % Load .dat and KiloSort/Phy output

             p = inputParser();
             % specify EXACTLY one of these:
             p.addParameter('cluster_id', [], @isscalar);
             p.addParameter('spike_idx', [], @isvector); % manually specify which idx into spike_times
             p.addParameter('spike_times', [], @isvector); % manually specify which times directly to extract

             % and ONE OR NONE of these
             p.addParameter('channel_idx', ds.channel_map, @isvector); % specify a subset of channels to extract
             p.addParameter('best_n_channels', NaN, @isscalar); % or take the best n channels based on this clusters template when cluster_id is scalar

             % other params:
             p.addParameter('num_waveforms', 2000, @isscalar); % caution: Inf will request ALL waveforms in order (typically useful if spike_times directly specified)
             p.addParameter('window', [-40 41], @isvector); % Number of samples before and after spiketime to include in waveform
             p.addParameter('car', false, @islogical);
             p.addParameter('centerUsingFirstSamples', 20, @(x) isscalar(x) || islogical(x)); % subtract mean of each waveform's first n samples, don't do if false

             p.addParameter('raw_dataset', ds.raw_dataset, @(x) true);

             p.parse(varargin{:});

             assert(ds.hasRawDataset, 'KiloSortDataset has no raw ImecDataset');
             
             ds.checkLoaded();

             cluster_idx = [];
             if ~isempty(p.Results.spike_times)
                 spike_times = p.Results.spike_times; %#ok<*PROPLC>
             elseif ~isempty(p.Results.spike_idx)
                 spike_times = ds.spike_times(p.Results.spike_idx);
             elseif ~isempty(p.Results.cluster_id)
                 clu = p.Results.cluster_id;
                 if isscalar(clu)
                     spike_times = ds.spike_times(ds.spike_clusters == clu);
                     cluster_idx = repmat(clu, size(spike_times));
                 else
                     [mask, which] = ismember(ds.spike_cluster, clu);
                     spike_times = ds.spike_times(mask);
                     cluster_idx = clu(which(mask));
                 end
             end

             % figure out actual times requested
             channel_idx = p.Results.channel_idx;
             if ~isnan(p.Results.best_n_channels)
                 % okay to have multiple clusters, just use first cluster
                 % to pick channels
%                  assert(isscalar(p.Results.cluster_id), 'best_n_channels_for_cluster only valid for scalar cluster_id');
                 cluster_ind = find(ds.cluster_ids == p.Results.cluster_id(1), 1);
                 if isempty(cluster_ind)
                     error('Cluster %d not found in cluster_ids', p.Results.cluster_id);
                 end
                 channel_idx = ds.cluster_best_template_channels(cluster_ind, 1:p.Results.best_n_channels);
             end

             if isfinite(p.Results.num_waveforms)
                 spike_times = randsample(spike_times, p.Results.num_waveforms);
             end

             % channel_map is provided since raw data often has additional channels that we're not interested in
             window = p.Results.window;
             snippetSet = p.Results.raw_dataset.readAPSnippetSet(spike_times, ...
                 window, channel_idx, 'car', p.Results.car);

             if p.Results.centerUsingFirstSamples
                 snippetSet.data = snippetSet.data - mean(snippetSet.data(:, 1:p.Results.centerUsingFirstSamples, :), 2, 'native');
             end
             snippetSet.cluster_idx = cluster_idx;
        end
         
        function snippetSet = readAPSnippetSet(ds, times, window, channel_idx, varargin)
            if nargin < 4
                channel_idx = 1:ds.nChannels;
            end
            snippetSet = ds.raw_dataset.readAPSnippetSet(times, window, channel_idx, varargin{:});
        end
            
    end
end
