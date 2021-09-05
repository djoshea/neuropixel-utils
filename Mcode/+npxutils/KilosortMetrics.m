classdef KilosortMetrics < handle
    % this class stores a set of computed information about a KilosortDataset,
    % particularly regarding template shapes, depths, and spiking drift
    % Note that this code borrows heavily from the math in github.com/cortex-lab/spikes
    
    properties (Transient) % not saved with metrics
        ks % npxutils.KilosortDataset
    end
    
    properties (Dependent)
        nSpikes
        nChannelsSorted
        nTemplates
        nTemplateTimepoints
        nTemplateChannels
        nClusters
        nConcatenatedFiles
        maxTemplatesPerCluster
        
        templateTimeRelative
        templateTimeRelativeMs
        
        nBatches
        nSpatialDims
    end
    
    properties
        % applied during construction
        ampRelativeThreshCentroid (1,1) single
        centroidMethod char
        extThreshLocalizedTemplate (1,1) single
        distThreshLocalizedTemplate (1,1) single
        
        % copied over from ks
        fs
        channelMap
        channel_ids_sorted
        concatenationInfo
        
        % nChannelsSorted x nChannelsSorted
        whitening_mat_inv (:,:) single
        
        % per template properties
        
        % nTemplates x nTimePoints x nTemplateChannels
        template_unw single % unwhitened, unscaled templates (i.e. in an arbitrary scale determined by raw * whitening_mat_inv)
        template_scale_factor single % nTemplates x 1 (mult factor from template_unw into template_scaled (uV)
        template_scaled single % unwhitened, scaled templates
        
        % nTemplates x nSpatialCoord
        template_centroid % single x,y,z, coords for each template's center of mass
        
        % nTemplates x 1
        template_is_localized (:,1) logical
        
        % nTemplates x nTimepoints
        template_waveform single % single, maximum amplitude waveform
        
        % nTemplates x 1
        template_waveform_ch (:,1) uint32 % channel where the maximum amplitude part of the template lives
        template_amplitude (:, 1) single % in uV
        template_ttp (:,1) single % trough to peak
        
        % nTemplates x nChannelsSorted
        template_best_channels uint32 % nTemplates x nChannelsSorted matrix indicating the closest channels to the max (typically take first 20 cols)
        
        % [nTemplates, nTemplates] = same as  similar_templates but updated using the final post-split, post-merge templates
        similar_templates_recomputed (:,:) single
        similar_templates_best_lag (:,:) int16
        
        % per spike properties
        
        % nSpikes x 1
        spike_times (:,1) uint64
        
        % nSpikes x 1
        spike_amplitude (:,1) single % in Uv,  product of template's largest amplitude over channels and scaling of each template (by ks.amplitudes)
        
        % nSpikes x nSpatialCoord
        spike_centroid single % single x,y,z coords for each spike's center of mass based on private PCs
        
        % nSpikes x 1
        spike_templates (:,1) uint32
        spike_clusters (:,1) uint32
        
        % per template properties
        
        % unique cluster ids corresponding to each slot in the other cluster properties
        cluster_ids (:,1) uint32
        
        % nClusters x 1
        cluster_template_mostUsed (:,1) uint32 % nClusters x 1 which template is most used by each cluster
        
        % nClusters x maxTemplatesPerCluster
        cluster_template_list cell
        
        % nClusters x nTemplates
        cluster_template_useCount uint64 % nClusters x nTemplates number of spikes in cluster i using template j
        cluster_template_weight single % nClusters x nTemplates fraction of spikes using template j of all cluster i's spikes
        cluster_num_templates uint32 % nClusters x 1 number of templates
        
        % nClusters x nChannelsSorted channel ids
        cluster_best_channels uint32
        
        % nClusters x nSpatialCoord
        cluster_centroid single
        
        % nClusters x 1 logical
        cluster_is_localized (:,1) logical
        
        % nClusters x nTimePoint x maxTemplatesPerCluster
        cluster_waveform single
        
        % [nClusters, nClusters] single matrix giving the similarity score (larger is more similar) between each pair of clusters, based on similar_templates
        similar_clusters (:,:) single
        
        % this version is based on similar_templates_recomputed and factors in adjusted
        % templates after splits are made
        similar_clusters_recomputed (:,:) single
        similar_clusters_best_lag (:,:) int16
        
        % nClusters x maxTemplatesPerCluster
        cluster_waveform_ch uint32
        cluster_amplitude single
        cluster_ttp single
    end
    
    properties % batch-wise versions
        % per template properties
        batchesComputed (:,1) uint32
        
        % nTemplates x nBatchesComputed
        template_useCount_batchwise (:,:) single
        
        % nTemplates x nSpatialCoord x nBatchesComputed
        template_centroid_batchwise (:,:,:) single % x,y,z, coords for each template's center of mass
        
        % nTemplates x nTemplateTime x nBatchesComputed
        template_waveform_batchwise (:,:,:) single % template waveform on channel template_waveform_ch
        
        % nTemplates x nBatchesComputed
        template_amplitude_batchwise (:,:) single % in uV, per-batch template amplitude (max - min) on ALL channels, not just template_waveform_ch
        
        % per cluster properties
        
        % nClusters x nTemplates x nBatchesComputed
        cluster_template_useCount_batchwise (:,:,:) uint32 % nClusters x nTemplates x nBatches number of spikes in cluster i using template j in batch k
        
        % nClusters x nSpatialCoord x nBatchesComputed
        cluster_centroid_batchwise (:,:,:) single
        
        % nClusters x nTimePoint x maxTemplatesPerCluster x nBatchesComputed
        cluster_waveform_batchwise (:,:,:) single
        
        % nClusters x maxTemplatesPerCluster x nBatchesComputed
        cluster_amplitude_batchwise (:,:) single
    end
    
    properties(Dependent)
        has_computed_batchwise
        nBatchesComputed
        
        template_depth % column 2 (y) of center of mass
        spike_depth % column 2 (y) of center of mass
        spike_is_localized % logical
        cluster_depth
        
        cluster_signed_peak % nClusters x 1
        
        cluster_signed_extr_ratio % ratio of trough to peak ratio, sign matches which one is bigger
    end
    
    methods % Constructor / metrics and batchwise metrics computation
        function this = KilosortMetrics(ks, varargin)
            p = inputParser();
            p.addParameter('ampRelativeThreshCentroid', 0.3, @isscalar);
            p.addParameter('centroidMethod', 'com_best_timepoint', @ischar);
            p.addParameter('extThreshLocalizedTemplate', 0.5, @isscalar);
            p.addParameter('distThreshLocalizedTemplate', 200, @isscalar);
            
            p.addParameter('progressInitializeFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(nUpdates) to print update
            p.addParameter('progressIncrementFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(updateString) to print update
            p.parse(varargin{:});
            
            assert(isa(ks, 'npxutils.KilosortDataset'));
            assert(ks.isLoaded, 'KilosortDataset is not loaded');
            
            nSteps = 8;
            initStr = 'Computing KilosortMetrics for dataset';
            if isempty(p.Results.progressInitializeFn) && isempty(p.Results.progressIncrementFn)
                prog = npxutils.internal.ProgressBar(nSteps, initStr);
                progIncrFn = @(text) prog.increment(text);
            else
                if ~isempty(p.Results.progressInitializeFn)
                    p.Results.progressInitializeFn(nSteps, initStr);
                end
                progIncrFn = p.Results.progressIncrementFn;
            end
            
            % store configurable settings
            this.centroidMethod = p.Results.centroidMethod;
            this.ampRelativeThreshCentroid = p.Results.ampRelativeThreshCentroid;
            this.extThreshLocalizedTemplate = p.Results.extThreshLocalizedTemplate;
            this.distThreshLocalizedTemplate = p.Results.distThreshLocalizedTemplate;
            
            this.ks = ks;
            this.fs = ks.fsAP;
            this.spike_templates = ks.spike_templates;
            this.spike_clusters = ks.spike_clusters;
            this.cluster_ids = ks.cluster_ids;
            this.channelMap = ks.channelMap;
            this.channel_ids_sorted = ks.channel_ids_sorted;
            
            this.concatenationInfo = ks.concatenationInfo;
            
            this.spike_times = ks.spike_times;
            
            % A. compute unwhitened templates
            % ks.templates is potentially only specified on a subset of channels which may differ across templates.
            % although for Kilosort this is every channel
            ks = this.ks; ks.checkLoaded(); %#ok<*PROP>
            templates = ks.templates;
            template_unw = zeros(size(templates, 1), size(templates, 2), ks.nChannelsSorted, 'like', templates);
            wmi = single(ks.whitening_mat_inv);
            this.whitening_mat_inv = wmi;
            assert(size(wmi, 1) == size(template_unw, 3), 'dim 3 of templates must match whitening matrix inverse');
            for iT = 1:size(templates, 1)
                whichChannels = ks.templates_ind(iT, :);
                template_unw(iT, :, whichChannels) = templates(iT, :, :);
            end
            sz = size(template_unw);
            this.template_unw = reshape(reshape(template_unw, sz(1)*sz(2), sz(3)) * wmi, sz);
            
            % for each channel i, list other channels in order of spatial proximity
            % nChannelsSorted x nChannelsSorted, include each channel in its own closest list
            closest_lookup = [ks.channel_ids_sorted, ks.channelMap.getClosestChannels(ks.nChannelsSorted-1, ks.channel_ids_sorted, ks.channel_ids_sorted)];
            this.template_best_channels = nan(ks.nClusters, ks.nChannelsSorted);
            for iT = 1:ks.nTemplates
                [~, bestChannelInd] = max(range(this.template_unw(iT, :, :), 2));
                this.template_best_channels(iT, :) = closest_lookup(bestChannelInd, :)';
            end
            
            progIncrFn('Template center of mass');
            % B. determine template center of mass
            %   1. compute template amp on each channel, zeroing out small (< 0.3 max) faraway channels
            templateUnscaledAmps = squeeze(max(this.template_unw,[],2) - min(this.template_unw,[],2)); % nTemplates x nTemplateChannels
            [templateUnscaledMaxAmp, templateMaxChInd] = max(templateUnscaledAmps, [], 2);
            %             threshAmp = templateUnscaledMaxAmp * ampRelativeThreshCentroid;
            %             templateUnscaledAmps(templateUnscaledAmps < threshAmp) = 0;
            
            %   2. compute template channel Centroids (nTemplates x nCh x nSpatialDim)
            this.template_centroid = npxutils.internal.computeWaveformImageCentroid(this.template_unw, ...
                1:ks.nChannelsSorted, ks.channel_positions_sorted, ...
                'method', this.centroidMethod, 'relativeThreshold', this.ampRelativeThreshCentroid);
            
            %             templateChannelPos = reshape(ks.channel_positions_sorted(ks.templates_ind(:), :), [size(ks.templates_ind), size(ks.channel_positions_sorted, 2)]);
            
            %   3. compute template center of mass
            %             m.template_centroid = npxutils.internal.TensorUtils.squeezeDims(sum(templateUnscaledAmps .* templateChannelPos, 2) ./ ...
            %                 sum(templateUnscaledAmps, 2), 2);
            
            % C. determine which templates are localized
            this.template_is_localized = false(ks.nTemplates, 1);
            for iT = 1:ks.nTemplates
                % find channels where the template has weight at least 0.5 of max
                threshExt = max(abs(this.template_unw(iT, :))) * this.extThreshLocalizedTemplate;
                maskCh = squeeze(max(abs(this.template_unw(iT, :, :)), [], 2)) >= threshExt; % nChannelsSorted x 1
                
                if nnz(maskCh) < 2
                    this.template_is_localized(iT) = true;
                else
                    distMat = pdist(ks.channel_positions_sorted(maskCh, :));
                    this.template_is_localized(iT) = max(distMat(:)) < this.distThreshLocalizedTemplate;
                end
            end
            
            % D. scale the spikes by both templates's largest intrinsic amplitude and then ks.amplitudes
            progIncrFn('Spike / template amplitude scaling');
            this.spike_amplitude = ks.amplitudes .* templateUnscaledMaxAmp(ks.spike_templates) * ks.apScaleToUv;
            
            % E. determine template scale by averaging ks.amplitudes that use that template
            assert(~isnan(ks.apScaleToUv));
            mean_template_amplitude = accumarray(ks.spike_templates, ks.amplitudes, [ks.nTemplates, 1], @mean);
            this.template_scale_factor = mean_template_amplitude * ks.apScaleToUv;
            this.template_scaled = this.template_unw .* this.template_scale_factor;
            
            % F. determine template largest scaled waveform
            progIncrFn('Template waveforms');
            template_waveform_ch = nan(ks.nTemplates, 1);
            template_waveform = nan(ks.nTemplates, size(ks.templates, 2));
            for iT = 1:ks.nTemplates
                template_waveform_ch(iT) = ks.templates_ind(iT, templateMaxChInd(iT));
                template_waveform(iT, :) = this.template_scaled(iT, :, template_waveform_ch(iT));
            end
            this.template_waveform_ch = template_waveform_ch;
            this.template_waveform = template_waveform;
            this.template_amplitude = templateUnscaledMaxAmp .* this.template_scale_factor;
            
            % G. determine trough to peak time of each waveform, this should be refined
            [~, time_template_trough] = min(this.template_waveform, [], 2);
            [~, time_template_peak] = max(this.template_waveform, [], 2);
            this.template_ttp = (time_template_peak - time_template_trough) / ks.fsAP * 1000;
            this.template_ttp(this.template_ttp <= 0) = NaN;
            
            if ks.hasFeaturesLoaded
                progIncrFn('Individual spike center of mass');
                % H. use private pcs to determine spike center of mass a la driftmap
                %   1. weight channels based on squared projection onto 1st pc
                pc1proj = squeeze(ks.pc_features(:, 1, :)); % nSpikes x nPCFeatures
                pc1proj(pc1proj < 0) = 0; % only positive entries contribute
                pc1weight = pc1proj.^2;
                
                %   2. compute center of mass. (spikeFeatureChannelPos is nSpikes x nCh x nSpatialDim)
                spike_pcfeat_chind = ks.pc_feature_ind(ks.spike_templates, :);
                this.spike_centroid = npxutils.internal.computeWaveformImageCentroid(reshape(pc1weight, [size(pc1weight, 1), 1, size(pc1proj, 2)]), ...
                    spike_pcfeat_chind, ks.channel_positions_sorted, ...
                    'method', 'com_via_provided_amplitudes', 'relativeThreshold', this.ampRelativeThreshCentroid);
                %                 spikeFeatureChannelPos = reshape(ks.channel_positions_sorted(spike_pcfeat_chind(:), :), [size(spike_pcfeat_chind), size(ks.channel_positions_sorted, 2)]);
                %                 m.spike_centroid = npxutils.internal.TensorUtils.squeezeDims(sum(pc1weight .* spikeFeatureChannelPos, 2) ./ ...
                %                     sum(pc1weight, 2), 2);
            end
            
            % recompute template similarity from W and U in order to get best lags
            if ks.hasFeaturesLoaded && ~isempty(ks.W)
                progIncrFn('Recomputing template similarity scores');
                % recompute cluster similarity from W and U
                [this.similar_templates_recomputed, this.similar_templates_best_lag] = ...
                    npxutils.internal.computeTemplateSimilarity(ks.W, ks.U);
            else
                this.similar_templates_recomputed = [];
                this.similar_templates_best_lag = zeros(size(ks.similar_templates), 'int16');
            end
            % I. compute cluster weighting over templates and list of templates used by each cluster, sorted by number of uses
            progIncrFn('Cluster weighting over templates');
            [mask_spike_in_cluster, cluster_ind_each_spike] = ismember(ks.spike_clusters, ks.cluster_ids);
            this.cluster_template_useCount = accumarray([cluster_ind_each_spike(mask_spike_in_cluster), ks.spike_templates(mask_spike_in_cluster)], ...
                ones(ks.nSpikes, 1), [ks.nClusters, ks.nTemplates]);
            this.cluster_template_weight = single(this.cluster_template_useCount) ./ single(ks.cluster_spike_counts);
            this.cluster_template_weight(isnan(this.cluster_template_weight)) = 0;
            
            [~, this.cluster_template_mostUsed] = max(this.cluster_template_useCount, [], 2);
            this.cluster_best_channels = this.template_best_channels(this.cluster_template_mostUsed, :);
            
            this.cluster_num_templates = sum(this.cluster_template_useCount > 0, 2);
            maxTemplatesPerCluster = max(this.cluster_num_templates);
            this.cluster_template_list = cell(ks.nClusters, 1);
            for iC = 1:ks.nClusters
                inds = find(this.cluster_template_useCount(iC, :));
                [~, sortIdx] = sort(this.cluster_template_useCount(iC, inds), 'descend');
                this.cluster_template_list{iC} = inds(sortIdx)';
            end
            
            % J. cluster is localized if all templates used are localized
            this.cluster_is_localized = arrayfun(@(iC) all(this.template_is_localized(this.cluster_template_list{iC})), 1:ks.nClusters);
            
            progIncrFn('Cluster center of mass');
            % K. find cluster center of mass
            this.cluster_centroid = npxutils.internal.TensorUtils.linearCombinationAlongDimension(this.template_centroid, 1, ...
                this.cluster_template_weight, 'replaceNaNWithZero', true); % needed in case any templates have nan centroid
            
            % L. cluster waveforms
            progIncrFn('Cluster waveform');
            cluster_waveform = nan(ks.nClusters, size(this.template_waveform, 2), maxTemplatesPerCluster);
            cluster_waveform_ch = nan(ks.nClusters);
            for iC = 1:ks.nClusters
                cluster_waveform(iC, :, 1:this.cluster_num_templates(iC)) = permute(this.template_waveform(this.cluster_template_list{iC}, :), [3 2 1]);
                cluster_waveform_ch(iC, 1:this.cluster_num_templates(iC)) = this.template_waveform_ch(this.cluster_template_list{iC})';
            end
            this.cluster_waveform = cluster_waveform;
            this.cluster_waveform_ch = cluster_waveform_ch;
            % J. cluster_amplitudes and ttp - weighted mean of template amplitude / ttp
            progIncrFn('Cluster amplitudes');
            this.cluster_amplitude = accumarray(cluster_ind_each_spike(mask_spike_in_cluster), this.spike_amplitude(mask_spike_in_cluster), [ks.nClusters, 1], @mean);
            template_mask = ~isnan(this.template_ttp);
            this.cluster_ttp = (double(this.cluster_template_useCount(:, template_mask)) * this.template_ttp(template_mask)) ./ sum(this.cluster_template_useCount(:, template_mask), 2);
            
            progIncrFn('Cluster similarity');
            % computed template_use_count weighted similarity to other templates
            % similar_clusters_i,j = sum_t1 sum_t2 [ template_weight(i, t1) * template_weight(j, t2) * similar_templates(t1, t2) ]
            this.similar_clusters = this.cluster_template_weight * (this.cluster_template_weight * ks.similar_templates)';
            if ~isempty(this.similar_templates_recomputed)
                tw = this.cluster_template_weight;
                tsim = this.similar_templates_recomputed;
                tlag = this.similar_templates_best_lag;
                
                this.similar_clusters_recomputed = tw * tsim * tw';
                this.similar_clusters_best_lag = round(tw * single(tlag) * tw') .* (this.similar_clusters_recomputed >= 0.1);
                
            else
                this.similar_clusters_recomputed = [];
                this.similar_clusters_best_lag = zeros(size(this.similar_clusters), 'int16');
            end
            if exist('prog', 'var')
                prog.finish();
            end
        end
        
        function computeBatchwiseMetrics(this, varargin)
            p = inputParser();
            p.addParameter('recompute', false, @islogical);
            p.addParameter('every_n_batches', 10, @islogical);
            p.addParameter('average_skipped_batches', false, @islogical);
            p.addParameter('progressInitializeFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(nUpdates) to print update
            p.addParameter('progressIncrementFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % f(updateString) to print update
            p.parse(varargin{:});
            
            % decide if we need to recompute anything
            ks = this.ks;
            recompute = p.Results.recompute;
            batchesComputed = uint32(1:p.Results.every_n_batches:ks.nBatches)';
            if ~isequal(batchesComputed, this.batchesComputed)
                recompute = true; % ensure we recompute everything if batch indices have changed
            end
            if this.has_computed_batchwise && ~p.Results.recompute
                return;
            end
            this.batchesComputed = batchesComputed;
            nBatchesComputed = numel(this.batchesComputed);
            
            average_skipped_batches = p.Results.average_skipped_batches;
            
            nSteps = 1 + this.nTemplates + this.nClusters;
            initStr = 'Computing batchwise KilosortMetrics for dataset';
            if isempty(p.Results.progressInitializeFn) && isempty(p.Results.progressIncrementFn)
                prog = npxutils.internal.ProgressBar(nSteps, initStr);
                progIncrFn = @(varargin) prog.increment(varargin{:});
            else
                if ~isempty(p.Results.progressInitializeFn)
                    p.Results.progressInitializeFn(nSteps, initStr);
                end
                progIncrFn = p.Results.progressIncrementFn;
            end
            
            if recompute || isempty(this.cluster_template_useCount_batchwise)
                progIncrFn('Computing batchwise cluster weighting over templates');
                batch_each_spike = ks.compute_which_batch(ks.spike_times);
                % compute batch membership to include all spikes until the next counted batch
                batch_ind_each_spike = discretize(batch_each_spike, [batchesComputed; ks.nBatches]);
                mask_spike_in_batch = ~isnan(batch_ind_each_spike);
                %                 [mask_spike_in_batch, batch_ind_each_spike] = ismember(batch_each_spike, m.batchesComputed);
                
                mask_spike = mask_spike_in_batch;
                this.template_useCount_batchwise = accumarray(...
                    [ks.spike_templates(mask_spike), batch_ind_each_spike(mask_spike)], ...
                    ones(nnz(mask_spike), 1, 'uint64'), [this.nTemplates, nBatchesComputed]);
                
                [mask_spike_in_cluster, cluster_ind_each_spike] = ismember(ks.spike_clusters, ks.cluster_ids);
                mask_spike = mask_spike_in_batch & mask_spike_in_cluster;
                this.cluster_template_useCount_batchwise = accumarray(...
                    [cluster_ind_each_spike(mask_spike), ks.spike_templates(mask_spike), batch_ind_each_spike(mask_spike)], ...
                    ones(nnz(mask_spike), 1, 'uint64'), [this.nClusters, this.nTemplates, nBatchesComputed]);
            end
            
            if recompute || isempty(this.template_centroid_batchwise)
                template_centroid_batchwise = nan(this.nTemplates, nBatchesComputed, this.nSpatialDims, 'single');
                template_waveform_batchwise = nan(this.nTemplates, this.nTemplateTimepoints, nBatchesComputed, 'single');
                template_amplitude_batchwise =  nan(this.nTemplates, nBatchesComputed, 'single');
                
                for iT = 1:this.nTemplates
                    progIncrFn();
                    
                    % 1 x time x channels x batch --> batch x time x channels
                    % these templates are all on all channel_ids_sorted
                    template_batch_time_channel = permute(this.construct_scaled_template_batchwise(iT, 'batches', this.batchesComputed, ...
                        'average_skipped_batches', average_skipped_batches), [4 2 3 1]);
                    
                    % some of these can be nan if the template is entirely zero
                    template_centroid_batchwise(iT, :, :) = npxutils.internal.computeWaveformImageCentroid(...
                        template_batch_time_channel, 1:this.nChannelsSorted, this.ks.channel_positions_sorted);
                    
                    template_waveform_batchwise(iT, :, :) = template_batch_time_channel(:, :, this.template_waveform_ch(iT))';
                    
                    % for amplitude, we take the global max range over all channels, not just the best tempalte_waveform_ch
                    template_amplitude_batchwise(iT, :) = max(max(template_batch_time_channel, [], 2) - min(template_batch_time_channel, [], 2), [], 3);
                end
                
                this.template_centroid_batchwise = template_centroid_batchwise;
                this.template_waveform_batchwise = template_waveform_batchwise;
                this.template_amplitude_batchwise = template_amplitude_batchwise;
            end
            
            if recompute || isempty(this.cluster_centroid_batchwise)
                cluster_centroid_batchwise  = nan(this.nClusters, this.nBatchesComputed, this.nSpatialDims, 'single');
                cluster_waveform_batchwise = nan(this.nClusters, this.nTemplateTimepoints, this.nBatches, 'single');
                cluster_amplitude_batchwise =  nan(this.nClusters, this.nBatchesComputed, 'single');
                
                % nClusters x nSpatialCoord x nBatches
                for iC = 1:this.nClusters
                    progIncrFn();
                    for iB = 1:this.nBatchesComputed
                        if any(this.cluster_template_useCount_batchwise(iC, :, iB) > 0)
                            % some weights will be nan where no spikes found for cluster, that's okay as we don't know where the cluster is then anyway
                            % this is a row vector of weights over templates
                            weights_temp_batch = single(this.cluster_template_useCount_batchwise(iC, :, iB)) ./ sum(this.cluster_template_useCount_batchwise(iC, :, iB), 2);
                            
                            % we have to mask over templates because unusued templates will have NaN centroids and the product will be NaN even where weight is 0
                            mask_template = weights_temp_batch ~= 0;
                            cluster_centroid_batchwise(iC, iB, :) =  weights_temp_batch(mask_template) * permute(this.template_centroid_batchwise(mask_template, iB, :), [1 3 2]);
                            
                            cluster_waveform_batchwise(iC, :, iB) = weights_temp_batch(mask_template) * this.template_waveform_batchwise(mask_template, :, iB);
                            cluster_amplitude_batchwise(iC, iB) = weights_temp_batch(mask_template) * this.template_amplitude_batchwise(mask_template, iB);
                        end
                    end
                end
                
                this.cluster_centroid_batchwise = cluster_centroid_batchwise;
                this.cluster_waveform_batchwise = cluster_waveform_batchwise;
                this.cluster_amplitude_batchwise = cluster_amplitude_batchwise;
            end
            
            if exist('prog', 'var')
                prog.finish();
            end
        end
        
        function templates = construct_scaled_template_batchwise(this, template_inds, varargin)
            templates = this.ks.construct_batchwise_templates(template_inds, varargin{:});
            templates = templates .* this.template_scale_factor(template_inds);
        end
    end
    
    methods % Dependent prop get impl
        function n = get.nSpikes(this)
            n = size(this.spike_times, 1);
        end
        
        function n = get.nChannelsSorted(this)
            n = size(this.template_unw, 3);
        end
        
        function n = get.nTemplates(this)
            n = size(this.template_unw, 1);
        end
        
        function n = get.nTemplateTimepoints(this)
            n = size(this.template_unw, 2);
        end
        
        function n = get.nTemplateChannels(this)
            n = size(this.template_unw, 3);
        end
        
        function n = get.nBatches(this)
            n = this.ks.nBatches;
        end
        
        function n = get.nBatchesComputed(this)
            n = numel(this.batchesComputed);
        end
        
        function n = get.nSpatialDims(this)
            n = this.channelMap.nSpatialDims;
            if size(this.ks.channel_positions_sorted, 2) < 2
                n = size(this.ks.channel_positions_sorted, 2);
            end
        end
        
        function n = get.nClusters(this)
            n = size(this.cluster_ids, 1);
        end
        
        function n = get.maxTemplatesPerCluster(this)
            n = size(this.cluster_waveform, 3);
        end
        
        function d = get.spike_depth(this)
            d = this.spike_centroid(:, 2);
        end
        
        function d = get.cluster_depth(this)
            d = this.cluster_centroid(:, 2);
        end
        
        function amp = get.cluster_signed_peak(this)
            % take a min and max over time dimension of the most used cluster waveform
            % and return either the min or the max, depending on which absolute value is larger
            hi = max(this.cluster_waveform(:,:,1), [], 2, 'omitnan');
            lo = min(this.cluster_waveform(:,:,1), [], 2, 'omitnan');
            
            amp = hi;
            amp(-lo > hi) = lo(-lo > hi);
        end
        
        function ratio = get.cluster_signed_extr_ratio(this)
            % take a min and max over time dimension of the most used cluster waveform
            % and return either the min or the max, depending on which absolute value is larger
            hi = max(this.cluster_waveform(:,:,1), [], 2, 'omitnan');
            lo = min(this.cluster_waveform(:,:,1), [], 2, 'omitnan');
            
            ratio = hi ./ -lo; % will be positive (typically axonal spikes)
            this = -lo > hi;
            ratio(this) = lo(this) ./ hi(this); % will be negative
        end
        
        function tvec = get.templateTimeRelative(this)
            T = int64(this.nTemplateTimepoints);
            off = int64(this.ks.template_sample_offset);
            start = -off + int64(1);
            stop = start + T - int64(1);
            tvec = (start:stop)';
        end
        
        function tvec = get.templateTimeRelativeMs(this)
            tvec = single(this.templateTimeRelative) / this.fs * 1000;
        end
        
        function tf = get.spike_is_localized(this)
            tf = this.template_is_localized(this.spike_templates);
        end
        
        function tf = get.has_computed_batchwise(this)
            tf = ~isempty(this.template_centroid_batchwise) && ~isempty(this.cluster_centroid_batchwise);
        end
        
        function assertHasKs(this)
            assert(~isempty(this.ks), 'Must set .ks to KilosortDataset');
        end
    end
    
    methods % Convenience methods
        function [channel_ids_unique, channel_ids_by_cluster] = ...
                gather_best_channels_multiple_clusters(this, cluster_ids, n_best_each)
            if nargin < 3
                n_best_each = 24;
            end
            
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            
            channel_ids_by_cluster = this.cluster_best_channels(cluster_inds, 1:n_best_each);
            channel_ids_unique = unique(channel_ids_by_cluster(:));
        end
    end
    
    methods % Cluster distance stats
        function cc_dist = computeClusterDistanceMatrix(this, varargin)
            p = inputParser();
            p.addParameter('cluster_ids', this.cluster_ids, @isvector);
            p.parse(varargin{:});
            
            cluster_inds = this.lookup_clusterIds(p.Results.cluster_ids);
            
            pos = this.cluster_centroid(cluster_inds, :);
            cc_dist = squareform(pdist(pos));
        end
    end
    
    methods % Cluster similarity statistics and CCGs
        function [similar_cluster_ids, similarity, best_lag] = computeSimilarClusters(this, cluster_ids, varargin)
            p = inputParser();
            p.addParameter('max_similar', this.nClusters-1, @isscalar);
            p.addParameter('min_similarity', 0, @isscalar);
            p.parse(varargin{:});
            
            cluster_ids = npxutils.internal.makecol(cluster_ids);
            [cluster_inds, ~] = this.lookup_clusterIds(cluster_ids);
            if isempty(this.similar_clusters_recomputed)
                sim = this.similar_clusters;
            else
                sim = this.similar_clusters_recomputed;
            end
            lags = this.similar_clusters_best_lag;
            
            [similar_clusters, which_input_clusters] = max(sim(:, cluster_inds), [], 2);
            similar_clusters(cluster_inds) = -1; % clear the self-similarity terms within group so they aren't considered
            
            max_similar = p.Results.max_similar;
            min_similarity = max(0, p.Results.min_similarity);
            
            [similarity,  similar_cluster_inds] = maxk(similar_clusters, max_similar);
            mask = similarity >= min_similarity;
            similar_cluster_inds = similar_cluster_inds(mask);
            similar_cluster_ids = this.cluster_ids(similar_cluster_inds);
            similarity = similarity(mask);
            
            % figure out best lag using the
            if isempty(lags)
                best_lag = zeros(size(similarity), 'int16');
            else
                which_input_clusters = which_input_clusters(similar_cluster_inds);
                lag_idx = sub2ind(size(lags), similar_cluster_inds, cluster_inds(which_input_clusters));
                best_lag = lags(lag_idx);
            end
        end
        
        function best_lag = computeBestLagForSimilarClusters(this, cluster_ids, similar_cluster_ids)
            if isempty(this.similar_clusters_best_lag) || isempty(cluster_ids) || isempty(similar_cluster_ids)
                best_lag = zeros(size(similar_cluster_ids), 'int16');
                return;
            end
            
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            sim_cluster_inds = this.lookup_clusterIds(similar_cluster_ids);
            
            if isempty(this.similar_clusters_recomputed)
                sim = this.similar_clusters(sim_cluster_inds, cluster_inds);
            else
                sim = this.similar_clusters_recomputed(sim_cluster_inds, cluster_inds);
            end
            lags = this.similar_clusters_best_lag(sim_cluster_inds, cluster_inds);
            
            [~, which_primary_cluster] = max(sim, [], 2);
            lag_idx = sub2ind(size(lags), (1:numel(sim_cluster_inds))', which_primary_cluster);
            best_lag = lags(lag_idx);
        end
        
        function [K, bins] = computeCCG(this, cluster_ids1, cluster_ids2, varargin)
            p = inputParser();
            p.addParameter('windowMs', 25, @isscalar);
            p.addParameter('binMs', 1, @isscalar);
            p.addParameter('include_cutoff_spikes', false, @islogical);
            p.parse(varargin{:});
            
            windowMs = p.Results.windowMs;
            binMs = p.Results.binMs;
            include_cutoff_spikes = p.Results.include_cutoff_spikes;
            
            st1 = this.internal_getSpikes(cluster_ids1, include_cutoff_spikes);
            st2 = this.internal_getSpikes(cluster_ids2, include_cutoff_spikes);
            
            [K, bins] = this.internal_computeCCG(st1, st2, windowMs, binMs);
        end
        
        function [K, bins] = computeACG(this, cluster_ids, varargin)
            p = inputParser();
            p.addParameter('windowMs', 25, @isscalar);
            p.addParameter('binMs', 1, @isscalar);
            p.addParameter('include_cutoff_spikes', false, @islogical);
            p.parse(varargin{:});
            
            windowMs = p.Results.windowMs;
            binMs = p.Results.binMs;
            include_cutoff_spikes = p.Results.include_cutoff_spikes;
            
            st = this.internal_getSpikes(cluster_ids, include_cutoff_spikes);
            [K, bins] = this.internal_computeCCG(st, [], windowMs, binMs);
        end
        
        function st = internal_getSpikes(this, cluster_ids, include_cutoff_spikes)
            if isscalar(cluster_ids)
                mask = this.ks.spike_clusters == cluster_ids;
            else
                mask = ismember(this.ks.spike_clusters, cluster_ids);
            end
            st = this.ks.spike_times(mask);
            
            if include_cutoff_spikes
                if isscalar(cluster_ids)
                    mask = this.ks.cutoff_spike_clusters == cluster_ids;
                else
                    mask = ismember(this.ks.cutoff_spike_clusters, cluster_ids);
                end
                st = sort(cat(1, st, this.ks.cutoff_spike_times(mask)));
            end
        end
        
        function [K, bins] = internal_computeCCG(this, st1, st2, windowMs, binWidthMs)
            % this is based on Marius's ccg() method
            
            if isempty(st2)
                st2 = st1;
                isACG = true;
            else
                isACG = false;
            end
            
            nbins = ceil(windowMs / binWidthMs);
            % st1 and st2 are both in ks samples defined by ks.fsAP
            
            dt_ms = nbins*binWidthMs;
            bins = (-dt_ms : binWidthMs : dt_ms)';
            
            % in samples
            dt = dt_ms * this.fs / 1000;
            tbin = binWidthMs * this.fs / 1000;
            
            ilow = 1;
            ihigh = 1;
            j = 1;
            
            K = zeros(2*nbins+1, 1);
            
            while j<=numel(st2)
                while (ihigh<=numel(st1)) && (double(st1(ihigh)) < double(st2(j))+dt)
                    ihigh = ihigh + 1;
                end
                while (ilow<=numel(st1)) && double(st1(ilow)) <= double(st2(j))-dt
                    ilow = ilow + 1;
                end
                if ilow>numel(st1)
                    break;
                end
                if double(st1(ilow)) > double(st2(j))+dt
                    j = j+1;
                    continue;
                end
                for k = ilow:(ihigh-1)
                    if isACG && k == j, continue, end
                    ibin = round((double(st2(j)) - double(st1(k)))/tbin);
                    K(ibin+ nbins+1) = K(ibin + nbins+1) + 1;
                end
                j = j+1;
            end
        end
    end
    
    methods % Computing / plotting stability over time
        function driftDistanceByCluster = computeClusterDriftDistance(this)
            driftDistanceByCluster = nan(this.nClusters, 1);
            for iC = 1:this.nClusters
                centroids = permute(this.cluster_centroid_batchwise(iC, :, :), [2 3 1]); % nBatch x spatial dim
                driftDistanceByCluster(iC) = max(pdist2(centroids, centroids, 'euclidean', 'Largest', 1));
            end
        end
        
        function [fracBinsValid, fracBinsValidExcludingEdges] = computeClusterValidTimeSpan(this, varargin)
            % fracBinsValid is the total fraction of bins that had at least threshSpikes
            % Then we divide up the day into any zeros at the beginning and end of the file, and the middle region between these bands of zeros
            
            p = inputParser();
            p.addParameter('binsClose', 2, @isscalar);
            p.addParameter('threshSpikes', 5, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            [counts, ~] = this.computeClusterBinnedCounts(p.Unmatched);
            valid = imclose(counts >= p.Results.threshSpikes, ones(1, p.Results.binsClose));
            
            fracBinsValid = mean(valid, 2);
            
            fracBinsValidExcludingEdges = zeros(size(fracBinsValid));
            for iC = 1:size(counts, 1)
                v = valid(iC, :);
                if ~any(v), continue, end
                idxStart = find(v, 1, 'first');
                idxStop = find(v, 1, 'last');
                
                fracBinsValidExcludingEdges(iC) = nnz(v) / (idxStop - idxStart + 1);
            end
        end
        
        function [counts, tvec] = computeClusterBinnedCounts(this, varargin)
            p = inputParser();
            p.addParameter('cluster_ids', [], @isvector);
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'npxutils.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('maskRegionsOutsideTrials', true, @islogical);
            p.addParameter('binWidth', 10, @isscalar);
            p.parse(varargin{:});
            
            tsi = p.Results.tsi;
            if ~isempty(tsi) && p.Results.maskRegionsOutsideTrials
                mask = this.computeSpikeMaskWithinTrials(tsi);
            else
                mask = true(numel(this.spike_times), 1);
            end
            if ~isempty(p.Results.cluster_ids)
                cluster_ids = p.Results.cluster_ids;
                mask = mask & ismember(this.spike_clusters, p.Results.cluster_ids);
            else
                cluster_ids = this.cluster_ids;
            end
            
            spikeTimes = double(this.spike_times(mask)) / this.ks.fsAP; % convert to seconds
            spikeClu = this.spike_clusters(mask); %#ok<*PROPLC>
            
            nClu = numel(cluster_ids);
            
            edges = (0:p.Results.binWidth:max(spikeTimes))';
            tvec = edges(1:end-1);
            counts = zeros(nClu, numel(tvec));
            for iC = 1:nClu
                thisClu = spikeClu == cluster_ids(iC);
                
                times = spikeTimes(thisClu);
                
                nSp = histcounts(times, edges);
                nSp = nSp(1:end-1);
                
                counts(iC, :) = nSp;
            end
            
            % keep only bins within trial regions
            if ~isempty(tsi)
                mask = false(numel(tvec), 1);
                tvec_samples = tvec * this.fs;
                bin_width_samples = p.Results.binWidth * this.fs;
                [idxStart, idxStop, ~] = tsi.computeActiveRegions();
                for iR = 1:numel(idxStart)
                    mask(tvec_samples >= idxStart(iR) & tvec_samples+bin_width_samples <= idxStop(iR)) = true;
                end
                counts = counts(:, mask);
                tvec = tvec(mask);
            end
        end
        
        function meanFR = computeClusterMeanFR(this, varargin)
            p = inputParser();
            p.addParameter('binWidth', 10, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            counts = this.computeClusterBinnedCounts('binWidth', p.Results.binWidth, varargin{:});
            meanFR = mean(counts, 2) / p.Results.binWidth;
        end
        
        
        function mask = computeSpikeMaskWithinTrials(this, tsi)
            mask = false(numel(this.spike_times), 1);
            [idxStart, idxStop, ~] = tsi.computeActiveRegions();
            
            for iR = 1:numel(idxStart)
                mask(this.spike_times >= idxStart(iR) & this.spike_times <= idxStop(iR)) = true;
            end
        end
    end
    
    methods % Plotting drift maps
        function [probe_offsets, time_samples] = computeProbeDrift(this, varargin)
            % compute a single estimate of probe dist along a specific dimension (typically 2)
            p = inputParser();
            p.addParameter('spatial_dimension', 2, @isscalar);
            p.addParameter('spikeAmpQuantile', 0.5, @isscalar); % consider only spikes larger than quantile fo amplitude
            p.addParameter('depth_bin', 2, @isscalar); % um for building histogram depth
            p.addParameter('time_bin_sec', 20, @isscalar); % seconds for producing depth estimates
            p.addParameter('weight_by_amp', true, @islogical); % instead of counting spikes, use amplitude-weighted counts for alignment
            p.addParameter('max_drift', 200, @isscalar) % max allowable drift in um
            
            % for plotting
            p.addParameter('show_plot', false, @islogical);
            p.addParameter('markerSize', 4, @isscalar);
            p.addParameter('timeInSeconds', true, @islogical);
            p.addParameter('show_all_spikes', true, @islogical);
            p.parse(varargin{:});
            
            spikeAmpQuantile = p.Results.spikeAmpQuantile;
            spikeAmpStdThresh = erfinv(2*(spikeAmpQuantile-0.5));
            drift_dim = p.Results.spatial_dimension;
            D = p.Results.depth_bin; % um
            T = p.Results.time_bin_sec * this.ks.fsAP; % samples
            weight_by_amp = p.Results.weight_by_amp;
            max_drift = p.Results.max_drift;
            
            time_samples = uint64(0:T:this.spike_times(end)-T)';
            n_time_bins = numel(time_samples);
            
            spikeTimes = this.spike_times;
            spikeAmps = this.spike_amplitude;
            spikePos = this.spike_centroid(:, drift_dim);
            
            % mask out spikes whose amplitude is >= mean + spikeAmpStdThresh * std
            amp_thresh = mean(spikeAmps) + spikeAmpStdThresh*std(spikeAmps);
            mask_spikes = this.spike_is_localized & spikeAmps > amp_thresh; % large spikes in current segment
            spikeAmps = spikeAmps(mask_spikes);
            spikePos = spikePos(mask_spikes);
            spikeTimes = spikeTimes(mask_spikes);
            
            bins = min(spikePos)-D:D:max(spikePos)+D;
            bin_pos = bins(1:end - 1) + D/2;
            n_depth_bins = numel(bin_pos);
            
            histmat = zeros(n_depth_bins, n_time_bins);
            for iT = 1:n_time_bins
                if iT == n_time_bins
                    mask_time = spikeTimes > time_samples(iT);
                else
                    mask_time = spikeTimes > time_samples(iT) & spikeTimes <= time_samples(iT+1);
                end
                
                if weight_by_amp
                    [~, ~, which_bin] = histcounts(spikePos(mask_time), bins);
                    h = accumarray(which_bin, spikeAmps(mask_time), [n_depth_bins, 1]);
                else
                    h = histcounts(spikePos(mask_time), bins);
                end
                histmat(:, iT) = h;
            end
            
            ref_col = histmat(:, floor(n_time_bins/2));
            delays = finddelay(ref_col, histmat, floor(max_drift / D));
            
            probe_offsets = single((delays - delays(1)) * D)';
            
            if p.Results.show_plot
                cla
                if p.Results.show_all_spikes
                    spike_mask_args = {};
                else
                    spike_mask_args = {'spike_mask', mask_spikes};
                end
                this.plotSpikesByAmplitude('spatial_dimension', drift_dim, 'timeInSeconds', p.Results.timeInSeconds, ...
                    spike_mask_args{:}, 'localizedOnly',true, 'markerSize', p.Results.markerSize);
                
                if p.Results.timeInSeconds
                    time_plot = [time_samples / this.ks.fsAP; (time_samples(end) + T) / this.ks.fsAP];
                else
                    time_plot = [time_samples; time_samples(end)+T];
                end
                probe_offsets_plot = [probe_offsets; probe_offsets(end)] + bin_pos(floor(n_depth_bins/2));
                hold on;
                stairs(time_plot, probe_offsets_plot, 'Color', 'r', 'LineWidth', 2);
                hold off;
                ylabel('position');
                xlabel('time')
            end
            
        end
        
        function plotClusterDriftmap(this, varargin)
            p = inputParser();
            p.addParameter('cluster_ids', [], @isvector);
            p.addParameter('spike_mask', [], @(x) isempty(x) || isvector(x));
            p.addParameter('localizedOnly', true, @islogical);
            p.addParameter('minAmpQuantile', 0, @isscalar);
            p.addParameter('showIndividual', false, @islogical);
            p.addParameter('highlightGaps', true, @islogical);
            p.addParameter('onlyGaps', false, @islogical),
            p.addParameter('minGapSeconds', 300, @isscalar);
            
            p.addParameter('showSmooth', true, @islogical); % smooth amplitudes over spikes by this width
            p.addParameter('smoothWidthSeconds', 50, @isscalar); % smoothing kernel width in seconds
            p.addParameter('scaleWithAmp', false, @islogical);
            p.addParameter('colorByAmp', false, @islogical);
            p.addParameter('alpha', 1, @isscalar);
            p.addParameter('zShuffleClusters', true, @islogical);
            p.addParameter('maxClustersPlot', Inf, @isscalar);
            p.addParameter('markerSize', 4, @isscalar);
            
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'npxutils.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('maskRegionsOutsideTrials', true, @islogical);
            p.addParameter('exciseRegionsOutsideTrials', false, @islogical);
            p.parse(varargin{:});
            
            mask = p.Results.spike_mask;
            if isempty(mask), mask = true(this.nSpikes, 1); end
            tsi = p.Results.tsi;
            if ~isempty(tsi) && (p.Results.maskRegionsOutsideTrials || p.Results.exciseRegionsOutsideTrials)
                mask = mask & this.computeSpikeMaskWithinTrials(tsi);
            end
            if p.Results.localizedOnly
                mask = mask & this.spike_is_localized;
            end
            if ~isempty(p.Results.cluster_ids)
                mask = mask & ismember(this.spike_clusters, p.Results.cluster_ids);
            end
            if p.Results.minAmpQuantile > 0
                thresh = quantile(this.spike_amplitude, p.Results.minAmpQuantile);
                mask = mask & this.spike_amplitude >= thresh;
            end
            
            showSmooth = p.Results.showSmooth;
            showIndividual = p.Results.showIndividual;
            highlightGaps = p.Results.highlightGaps;
            onlyGaps = p.Results.onlyGaps;
            alpha = p.Results.alpha;
            if alpha == 1 && showSmooth
                alpha = 0.3;
            end
            colorByAmplitude = p.Results.colorByAmp;
            
            spikeTimes = this.spike_times(mask);
            if p.Results.exciseRegionsOutsideTrials
                timeShifts = tsi.computeShiftsExciseRegionsOutsideTrials();
                spikeTimes = timeShifts.shiftTimes(spikeTimes);
            end
            spikeTimes = double(spikeTimes) / this.ks.fsAP; % convert to seconds
            spikeAmps = this.spike_amplitude(mask);
            spikeYpos = this.spike_depth(mask);
            spikeClu = this.spike_clusters(mask); %#ok<*PROPLC>
            ampMax = quantile(spikeAmps, 0.95);
            
            uClu = unique(spikeClu);
            clusterInds = this.lookup_clusterIds(uClu);
            nClu = numel(uClu);
            clusterAmps = this.cluster_amplitude(clusterInds);
            [~, clusterAmpsSortOrder] = sort(clusterAmps, 'ascend');
            
            bgcolor = [0.92 0.92 0.95];
            
            if colorByAmplitude
                cluAmpMax = quantile(clusterAmps, 0.95);
                cluAmpNormalized = min(1, clusterAmps / cluAmpMax);
                %
                %                 cmap = amplitudeCmap(clusterAmpsSortOrder / numel(clusterAmps));
                [cmap, cmap_base] = amplitudeCmap(cluAmpNormalized);
                cmap_lims = [0 cluAmpMax];
            else
                cmap = npxutils.internal.graphics.distinguishable_colors(nClu, [1 1 1; 0 1 0]);
                cmap_base = cmap;
                if nClu == 1
                    cmap_lims = [1 2];
                else
                    cmap_lims = [1 nClu];
                end
            end
            
            rng(1);
            
            if ~isinf(p.Results.maxClustersPlot) && p.Results.maxClustersPlot < nClu
                uClu = uClu(ismember((1:nClu)', randperm(nClu, p.Results.maxClustersPlot)));
                nClu = numel(uClu);
            end
            if p.Results.zShuffleClusters
                clusterAmpsSortOrder = clusterAmpsSortOrder(randperm(nClu));
            end
            
            cla;
            smoothWidthSeconds = p.Results.smoothWidthSeconds;
            smoothBinWidth = smoothWidthSeconds / 5;
            minGapBins = ceil(p.Results.minGapSeconds / smoothBinWidth);
            spikeTimeBins = (0:smoothBinWidth:max(spikeTimes))';
            for iC_sorted = 1:nClu
                iC = clusterAmpsSortOrder(iC_sorted);
                thisClu = spikeClu == uClu(iC);
                
                x = spikeTimes(thisClu);
                y = spikeYpos(thisClu);
                
                if p.Results.scaleWithAmp
                    sz = spikeAmps(thisClu) / ampMax * 5*p.Results.markerSize;
                else
                    sz = p.Results.markerSize;
                end
                
                ud = struct('cluster_id', uClu(iC), 'cluster_amplitude', sprintf('%.1f uV', clusterAmps(iC)), ...
                    'cluster_is_localized', this.cluster_is_localized(clusterInds(iC)), ...
                    'xname', 'Time', 'yname', 'Depth', 'xunits', 'sec', 'yunits', 'um');
                
                if showSmooth || onlyGaps
                    [xsm, ysm, gaps] = binsmooth(x, y, spikeTimeBins, 5, smoothBinWidth, minGapBins); % min 10 spikes in 10 second increments
                    hasGaps = ~isempty(gaps);
                end
                
                if showIndividual && (~onlyGaps || hasGaps)
                    h = plot(x, y, '.', 'Color', cmap(iC,:), 'MarkerSize', sz, 'UserData', ud);
                    TrialDataUtilities.Plotting.setMarkerOpacity(h, alpha);
                    
                    hold on;
                end
                if showSmooth && (~onlyGaps || hasGaps)
                    %                     if showIndividual, width = 2; else, width = 1; end
                    width = 2;
                    plot(xsm, ysm, '-', 'Color', cmap(iC,:), 'LineWidth', width, 'UserData', ud);
                    hold on;
                    
                    if (highlightGaps || onlyGaps) && hasGaps
                        plot(xsm(gaps), ysm(gaps), 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', cmap(iC,:), 'UserData', ud);
                    end
                end
            end
            
            xlabel('time (sec)')
            ylabel('y position (um)')
            hold off;
            box off;
            axh = gca;
            axh.Color = bgcolor;
            axh.GridColor = [1 1 1];
            axh.GridAlpha = 1;
            axh.MinorGridColor = [0.96 0.96 0.96];
            axh.MinorGridAlpha = 1;
            axh.MinorGridLineStyle = '-';
            axh.TickDir = 'out';
            axh.ColorSpace.Colormap = cmap_base;
            axh.CLim = cmap_lims;
            set(gcf, 'InvertHardcopy', 'off');
            axis tight;
            
            if colorByAmplitude
                hbar = colorbar(axh);
                ylabel(hbar, 'Amplitude (uV)');
            end
            
            npxutils.internal.configureDataTipsFromUserData(gcf);
            
            function [xb, yb, gaps] = binsmooth(x, y, edges, smoothBy, minSpikesPerBin, minGap)
                [nSp, ~, whichBin] = histcounts(x, edges);
                xb = 0.5 * (edges(1:end-1) + edges(2:end));
                nSp = nSp(1:end-1);
                
                mask = whichBin > 0;
                yb = accumarray(whichBin(mask), y(mask), [numel(xb), 1], @mean, single(NaN));
                
                % don't trust estimates with fewer than nSpikes
                mask_invalidate = nSp < minSpikesPerBin;
                yb(mask_invalidate) = NaN;
                
                % find large gaps
                mask_gap_eligible = nSp < minSpikesPerBin;
                mask_gap = imopen(mask_gap_eligible, ones(1, minGap));
                first_gap = [false, diff(mask_gap) == 1];
                
                % drop invalid samples (so we plot over them) unless inside large gap (so that large gaps still appear)
                mask_keep = ~mask_invalidate | first_gap;
                xb = xb(mask_keep);
                yb = yb(mask_keep);
                
                gaps = find(first_gap(mask_keep));
                gaps = gaps(gaps > 1) - 1; % ensure the index points to before the gap starts
                
                if smoothBy > 0
                    wasnan = isnan(yb);
                    yb = smooth(yb, smoothBy, 'lowess', 2);
                    yb(wasnan) = NaN;
                end
            end
            
            function [cmap, cmap_base] = amplitudeCmap(amp)
                cmap_base = npxutils.internal.graphics.cmocean('thermal');
                N = size(cmap_base, 1);
                lerp = @(x, a, b) x*(b-a) + a; % map x from [0 1] to [a b]
                cmap_ind = round(lerp(amp, 1, N));
                cmap = cmap_base(cmap_ind, :);
                %
                %                 cmap_hsl = cat(2, rand(N, 1), repmat(0.7, N, 1), repmat(0.5, N, 1));
                %
                %                 lum = lerp(amp, 0, 0.9);
                % %                 sat = lerp(amp, 0, 0.9);
                % %                 cmap_hsl(:, 2) = sat;
                %                 cmap_hsl(:, 3) = lum;
                %                 cmap = npxutils.internal.graphics.hsl2rgb(cmap_hsl);
            end
        end
        
        % Inputs: spikeTimes, spikeAmps, spikeYpos - names self explanatory
        %         opt - optional, empty by default; 'mark' - will mark detected drifts, 'show' - will generate a different plot,
        %               where only large spikes are used, and the detection of drift locations is demonstrated
        function info = plotDriftmap(this, varargin)
            % Note: this is pretty much copied wholesale from cortex-lab/spikes, with some nice annotations added on
            
            p = inputParser();
            p.addParameter('mode', 'mark', @ischar); % or 'show' or 'mark'
            p.addParameter('spikeAmpQuantile', 0.966, @isscalar); % consider only spikes larger than quantile fo amplitude
            p.addParameter('segmentDepth', 800, @isscalar); % um to segment probe into
            p.addParameter('nAmpBins', 20, @isscalar);
            
            p.addParameter('spike_mask', [], @(x) isempty(x) || isvector(x));
            p.addParameter('cluster_ids', [], @isvector);
            p.addParameter('localizedOnly', true, @islogical);
            
            p.addParameter('driftThreshold', 6, @isscalar); % um
            p.addParameter('driftTimeWindow', 10, @isscalar); % in seconds
            p.addParameter('minSpikesDrift', 50, @isscalar);
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'npxutils.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('trialTickColorIndex', [], @(x) isempty(x) || isvector(x));
            p.addParameter('maskRegionsOutsideTrials', false, @islogical);
            p.addParameter('markBoundaries', true, @islogical);
            p.addParameter('exciseRegionsOutsideTrials', false, @islogical);
            
            p.addParameter('markSpecificTrialIds', [], @(x) isempty(x) || isvector(x));
            p.addParameter('markSpecificTrialLabels', [], @(x) isempty(x) || isvector(x));
            
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('ampRange', [], @(x) isvector(x) || isempty(x));
            p.addParameter('markerSize', 4, @isscalar);
            
            p.addParameter('timeInSeconds', true, @islogical);
            p.addParameter('sample_window', [], @(x) isempty(x) || isvector(x)); % subselect a window to plot
            
            p.addParameter('plotDriftEventsInline', false, @islogical);
            p.addParameter('plotDriftEventTicks', true, @islogical);
            p.addParameter('configureDataTips', true, @islogical);
            
            p.parse(varargin{:});
            opt = p.Results.mode;
            segDepth = p.Results.segmentDepth;
            tWindow = p.Results.driftTimeWindow;
            driftThreshold = p.Results.driftThreshold;
            minSpikesDrift = p.Results.minSpikesDrift;
            tsi = p.Results.tsi;
            spikeAmpQuantile = p.Results.spikeAmpQuantile;
            spikeAmpStdThresh = erfinv(spikeAmpQuantile);
            xOffset = double(p.Results.xOffset);
            sample_window = p.Results.sample_window;
            timeInSeconds = p.Results.timeInSeconds;
            
            % mask
            mask = p.Results.spike_mask;
            if isempty(mask), mask = true(this.nSpikes, 1); end
            if ~isempty(tsi) && (p.Results.maskRegionsOutsideTrials || p.Results.exciseRegionsOutsideTrials)
                mask = mask & this.computeSpikeMaskWithinTrials(tsi);
            end
            if p.Results.localizedOnly
                mask = mask & this.spike_is_localized;
            end
            if ~isempty(p.Results.cluster_ids)
                mask = mask & ismember(this.spike_clusters, p.Results.cluster_ids);
            end
            
            % this would be present if we're looking at the cleaned file where times have already
            % been shifted around to excise gaps
            timeShiftsFromRaw = [];
            if ~isempty(this.concatenationInfo)
                timeShiftsAlreadyApplied = this.concatenationInfo.timeShifts;
            else
                timeShiftsAlreadyApplied = [];
            end
            
            spikeTimes = this.spike_times(mask);
            
            % this allows us to apply new time shifts to the data to see what the effect would be
            if p.Results.exciseRegionsOutsideTrials
                assert(isempty(timeShiftsAlreadyApplied), 'Cannot do exciseRegionsOutsideTrials when time shifts already applied given in concatenationInfo');
                timeShiftsAppliedHere = tsi.computeShiftsExciseRegionsOutsideTrials();
                spikeTimes = timeShiftsAppliedHere.shiftTimes(spikeTimes);
                timeShiftsFromRaw = timeShiftsAppliedHere;
            else
                timeShiftsAppliedHere = [];
            end
            
            if ~isempty(sample_window)
                mask_subselect = spikeTimes >= sample_window(1) & spikeTimes <= sample_window(2);
                spikeTimes = spikeTimes(mask_subselect);
                mask(mask) = mask_subselect;
            end
            
            spikeTimesSeconds = double(spikeTimes) / this.fs;
            if timeInSeconds
                spikeTimes = spikeTimesSeconds; % convert to seconds
            end
            spikeAmps = this.spike_amplitude(mask);
            spikeYpos = this.spike_depth(mask);
            
            this.plotSpikesByAmplitude('spike_mask', mask, 'time_shifts', timeShiftsAppliedHere, 'nAmpBins', p.Results.nAmpBins, ...
                'ampRange', p.Results.ampRange, 'timeInSeconds', timeInSeconds, ...
                'localizedOnly', p.Results.localizedOnly, 'tsi', p.Results.tsi, 'xOffset', xOffset, ...
                'markerSize', p.Results.markerSize, 'configureDataTips', p.Results.configureDataTips);
            
            ylim(this.channelMap.ylim);
            set(gca, 'XLimSpec', 'tight');
            box off;
            
            if p.Results.plotDriftEventsInline || p.Results.plotDriftEventTicks
                nD = floor(max(spikeYpos) / segDepth);
                driftEventsAll = cell(nD, 1);
                for iD = 1:nD % break the recording into 800 um segments
                    d = segDepth*(iD-1);
                    tmp = spikeAmps(spikeYpos >= d & spikeYpos < d+segDepth);
                    I = spikeAmps > mean(tmp) + spikeAmpStdThresh*std(tmp) & spikeYpos >= d & spikeYpos < d+segDepth; % large spikes in current segment
                    driftEvents = detectDriftEvents(spikeTimesSeconds(I), spikeYpos(I), strcmp(opt, 'show'));
                    if ~timeInSeconds && ~isempty(driftEvents)
                        driftEvents(:, 1) = driftEvents(:,1) * this.fs;
                    end
                    driftEventsAll{iD} = driftEvents;
                    if p.Results.plotDriftEventsInline && ~isempty(driftEvents)
                        plot(driftEvents(:,1) + xOffset, driftEvents(:,2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r')
                        % text(driftEvents(:,1)+1, driftEvents(:,2), num2str(round(driftEvents(:,3))), 'Color', 'r') % the magnitude of the drift
                    end
                end
                
                driftEvents = cat(1, driftEventsAll{:});
                if ~isempty(driftEvents)
                    driftTimes = driftEvents(:, 1);
                    if p.Results.plotDriftEventTicks
                        hold on;
                        npxutils.internal.graphics.rugplot(driftTimes + xOffset, 'side', 'top', 'Color', [1 0.2 0.2]);
                    end
                else
                    driftTimes = [];
                end
            end
            
            if ~isempty(tsi)
                % make ticks at the bottom of the plot
                hold on;
                if isempty(p.Results.trialTickColorIndex)
                    % single color for all trials
                    tsi.markTrialTicks('sample_window', sample_window, ...
                        'timeInSeconds', timeInSeconds, 'trialTickColorIndex', p.Results.trialTickColorIndex, ...
                        'time_shifts', timeShiftsAppliedHere, 'xOffset', xOffset, 'side', 'bottom');
                else
                    cIndex = p.Results.trialTickColorIndex;
                    uinds = unique(cIndex);
                    uinds = setdiff(uinds, [NaN 0]);
                    N = numel(uinds);
                    cmap = npxutils.internal.graphics.colorcet('C2', 'N', N);
                    for iC = 1:N
                        mask_this = cIndex == uinds(iC);
                        if ~any(mask_this), continue; end
                        tsi.markTrialTicks('sample_window', sample_window, ...
                            'timeInSeconds', timeInSeconds, 'maskTrials', mask_this, ...
                            'time_shifts', timeShiftsAppliedHere, 'xOffset', xOffset, 'side', 'bottom', 'Color', cmap(iC, :));
                    end
                end
                
                % optionally mark specific trial starts
                if ~isempty(p.Results.markSpecificTrialIds)
                    markTrialIds = npxutils.internal.makecol(p.Results.markSpecificTrialIds);
                    markTrialLabels = p.Results.markSpecificTrialLabels;
                    tsi.markSpecificTrials(markTrialIds, 'labels', markTrialLabels, 'Color', cmap, 'sample_window', sample_window, ...
                        'timeInSeconds', timeInSeconds, 'time_shifts', timeShiftsAppliedHere, 'xOffset', xOffset);
                end
            end
            
            if p.Results.markBoundaries
                if isempty(this.concatenationInfo)
                    warning('KilosortDataset had empty concatenationInfo, cannot mark excision boundaries');
                else
                    % the concatenation info struct might contain info on the applied time shifts
                    % but we want them to be shown in the already time shifted time axis
                    this.concatenationInfo.markExcisionBoundaries('sample_window', sample_window, ...
                        'timeInSeconds', timeInSeconds, ...
                        'time_shifts', timeShiftsFromRaw, 'xOffset', xOffset);
                end
            end
            
            if p.Results.exciseRegionsOutsideTrials && p.Results.markBoundaries
                this.markExcisionBoundaries(timeShiftsAppliedHere, 'sample_window', sample_window, ...
                    'timeInSeconds', timeInSeconds, 'time_shifts', timeShiftsAppliedHere, ...
                    'xOffset', xOffset);
            end
            if isempty(this.concatenationInfo)
                warning('KilosortDataset had empty concatenationInfo, cannot mark concatenation boundaries');
            else
                this.concatenationInfo.markConcatenatedFileBoundaries('sample_window', sample_window, ...
                    'timeInSeconds', timeInSeconds, ...
                    'time_shifts', timeShiftsAppliedHere, 'xOffset', xOffset);
            end
            
            set(gca, 'LooseInset', [0 0 0 0])
            hold off;
            
            info.driftTimes = driftTimes;
            info.xMax = max(spikeTimes) + xOffset;
            
            % driftEvents will contain a column of times, a column of depths, and a column of drift magnitudes
            function driftEvents = detectDriftEvents(spikeTimes, spikeDepths, doPlot)
                % spikeTimes must be in seconds
                if nargin < 3
                    doPlot = false;
                end
                driftEvents = [];
                if isempty(spikeTimes)
                    return;
                end
                D = 2; % um
                bins = min(spikeDepths)-D:D:max(spikeDepths)+D;
                h = histcounts(spikeDepths, bins);
                
                h = h(1:end-1); % last bin represents the scalar value bins(end), not an interval
                bins = bins(1:end - 1) + D/2; % now it's the centre of each interval
                
                [~, locs] = findpeaks(h);
                
                if doPlot
                    ax(1) = subplot(1, 5, 1); hold on;
                    plot(h, bins, 'k')
                    box off
                    ylabel('y position')
                    ax(2) = subplot(1, 5, 2:5); hold on;
                    plot(spikeTimes, spikeDepths, '.', 'Color', 0.5*[1 1 1])
                    linkaxes(ax, 'y')
                    xlim([0 spikeTimes(end)+1])
                    xlabel('time')
                end
                
                for iP = 1:numel(locs)
                    if h(locs(iP)) < 0.3*spikeTimes(end)
                        continue
                        % we want the peaks to correspond to some minimal firing rate (otherwise peaks by very few spikes will be considered as well...)
                    end
                    if doPlot
                        subplot(1, 5, 1); hold on;
                    end
                    
                    posBegin = find(h(1:locs(iP)) < 0.05*h(locs(iP)), 1, 'last');
                    if isempty(posBegin)
                        posBegin = 1;
                    end
                    posEnd   = find(h(locs(iP):end) < 0.05*h(locs(iP)), 1, 'first') + locs(iP) - 1;
                    if isempty(posEnd)
                        posEnd = numel(bins);
                    end
                    if (iP > 1 && posBegin < locs(iP-1)) || (iP < numel(locs) && posEnd > locs(iP+1))
                        if doPlot
                            plot(h(locs(iP)), bins(locs(iP)), 'bo')
                        end
                        continue % no clean enough separation from neighbour peak(s
                    elseif doPlot
                        plot(h(locs(iP)), bins(locs(iP)), 'ro')
                        plot(xlim, bins(posBegin)*[1 1], '--', 'Color', 0.5*[1 1 1])
                        plot(xlim, bins(posEnd)*[1 1], '--', 'Color', 0.5*[1 1 1])
                    end
                    
                    I = spikeDepths > bins(posBegin) & spikeDepths < bins(posEnd);
                    
                    currentspikeDepths = spikeDepths(I);
                    currentspikeTimes  = spikeTimes(I);
                    for t = 0:tWindow:spikeTimes(end)-tWindow
                        I = currentspikeTimes >= t & currentspikeTimes <= t+tWindow;
                        driftSize = bins(locs(iP)) - median(currentspikeDepths(I));
                        if abs(driftSize) > driftThreshold && sum(I) > minSpikesDrift % 6 um is the hardcoded threshold for drift, and we want at least 10 spikes for the median calculation
                            driftEvents(end+1,:) = [t+5, bins(locs(iP)), driftSize]; %#ok<AGROW>
                        end
                    end
                    if doPlot && ~isempty(driftEvents)
                        subplot(1, 5, 2:5); hold on;
                        plot(driftEvents(:,1), driftEvents(:,2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r')
                        text(driftEvents(:,1)+1, driftEvents(:,2), num2str(round(driftEvents(:,3))), 'Color', 'r') % the magnitude of the drift
                    end
                end % loop on peak locations
            end
        end
        
        function plotSpikesByAmplitude(this, varargin)
            % this plots each clusters' spikes in their own color, either as markers or as a line with smoothing
            p = inputParser();
            p.addParameter('spatial_dimension', 2, @isscalar);
            p.addParameter('spike_mask', [], @(x) isempty(x) || isvector(x));
            p.addParameter('localizedOnly', true, @islogical);
            p.addParameter('nAmpBins', 20, @isscalar);
            p.addParameter('ampRange', [], @(x) isvector(x) || isempty(x));
            p.addParameter('markerSize', 2, @isscalar);
            p.addParameter('timeInSeconds', true, @islogical);
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'npxutils.TimeShiftSpec'));
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'npxutils.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('configureDataTips', true, @islogical);
            p.parse(varargin{:});
            
            mask = p.Results.spike_mask;
            if isempty(mask)
                mask = true(numel(this.spike_times), 1);
            end
            if p.Results.localizedOnly
                mask = mask & this.spike_is_localized;
            end
            timeShifts =p.Results.time_shifts;
            spikeTimes = this.spike_times(mask);
            spikeTimesOrig = spikeTimes;
            timeInSeconds = p.Results.timeInSeconds;
            spatial_dim = p.Results.spatial_dimension;
            
            if ~isempty(timeShifts)
                spikeTimes = timeShifts.shiftTimes(spikeTimes);
            end
            
            if timeInSeconds
                spikeTimes = double(spikeTimes) / this.ks.fsAP; % convert to seconds
            end
            spikeClusters = this.spike_clusters(mask);
            spikeAmps = this.spike_amplitude(mask);
            spikeYpos = this.spike_centroid(mask, spatial_dim);
            
            nAmpBins = p.Results.nAmpBins;
            if ~isempty(p.Results.ampRange)
                ampRange = p.Results.ampRange;
            else
                ampRange = quantile(spikeAmps, [0.1 0.9]);
            end
            colorBins = linspace(ampRange(1), ampRange(2), nAmpBins);
            xOffset = p.Results.xOffset;
            
            colors = gray(nAmpBins+1);
            colors = colors(end-1:-1:1, :); % first bin is smalles spikes, starts white
            for b = 1:nAmpBins-1
                theseSpikes = spikeAmps>=colorBins(b) & spikeAmps<=colorBins(b+1);
                
                h = plot(spikeTimes(theseSpikes) + xOffset, spikeYpos(theseSpikes), '.', 'Color', colors(b,:), 'MarkerSize', p.Results.markerSize);
                
                if ~verLessThan('matlab', '9.6.0') && p.Results.configureDataTips % R2019a
                    h.DataTipTemplate.DataTipRows(1).Label = 'Shifted Time';
                    h.DataTipTemplate.DataTipRows(2).Format = '%g sec';
                    h.DataTipTemplate.DataTipRows(2).Label = 'Probe Y';
                    h.DataTipTemplate.DataTipRows(2).Format = '%g um';
                    row = dataTipTextRow('Raw sample', double(spikeTimesOrig(theseSpikes)'), '%d');
                    h.DataTipTemplate.DataTipRows(end+1) = row;
                    row = dataTipTextRow('Cluster', double(spikeClusters(theseSpikes)'), '%d');
                    h.DataTipTemplate.DataTipRows(end+1) = row;
                    row = dataTipTextRow('Amplitude', double(spikeAmps(theseSpikes)'), '%.2f mV');
                    h.DataTipTemplate.DataTipRows(end+1) = row;
                    
                    if ~isempty(p.Results.tsi)
                        [~, spikeTrialIds] = p.Results.tsi.segmentTimes(spikeTimesOrig(theseSpikes));
                        row = dataTipTextRow('Trial Id', double(spikeTrialIds), '%d');
                        h.DataTipTemplate.DataTipRows(end+1) = row;
                    end
                end
                
                if ~isempty(this.concatenationInfo) && this.concatenationInfo.nDatasets > 1
                    % show original file name and ind
                    [fileInd, origInd] = this.lookup_sampleIndexInConcatenatedFile(spikeTimesOrig(theseSpikes));
                    
                    % grab associated file name (too slow)
                    %                     fileNames = strings(numel(spikeTimesOrig), 1);
                    %                     mask = ~isnan(fileInd);
                    %                     fileNames(mask) = m.concatenatedNames(fileInd(mask));
                    
                    if ~verLessThan('matlab', '9.6.0')  && p.Results.configureDataTips % R2019a
                        row = dataTipTextRow('Orig File Ind', double(fileInd), '%d');
                        h.DataTipTemplate.DataTipRows(end+1) = row;
                        row = dataTipTextRow('Orig Sample Ind', double(origInd), '%d');
                        h.DataTipTemplate.DataTipRows(end+1) = row;
                    end
                end
                
                hold on;
            end
            xlabel('time (sec)')
            ylabel('y position (um)')
            
            h = zoom;
            h.Motion = 'horizontal';
        end
        
        function [fileInds, origSampleInds] = lookup_sampleIndexInConcatenatedFile(this, inds)
            if ~isempty(this.concatenationInfo)
                [fileInds, origSampleInds] = this.concatenationInfo.lookup_sampleIndexInSourceFiles(inds);
            else
                fileInds = ones(size(inds));
                origSampleInds = inds;
            end
        end
    end
    
    methods % Cluster amplitudes vs. time
        function h = plotClusterAmplitudeVsTime(this, cluster_ids, varargin)
            [tf, which_cluster] = ismember(this.spike_clusters, cluster_ids);
            spike_inds = find(tf);
            which_cluster = which_cluster(tf);
            value = this.spike_amplitude(tf);
            
            h = this.internal_plotSpikeVsTime(spike_inds, value, 'color', which_cluster, 'valueLabel', 'amplitude', ...
                'valueFormat', '%.01f', 'valueUnits', 'uV', varargin{:}); %#ok<FNDSB>
        end
        
        function h = internal_plotSpikeVsTime(this, spike_inds, spike_value, varargin)
            p = inputParser();
            p.addParameter('color', [], @(x) true); % can be N x 1 or single color / value to map into linear colormap
            p.addParameter('size', 4^2, @isvector); % can be N x 1 or single size
            p.addParameter('valueLabel', 'value', @ischar);
            p.addParameter('valueFormat', '', @ischar);
            p.addParameter('valueUnits', '', @ischar);
            p.addParameter('timeInSeconds', false, @islogical);
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'npxutils.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('maskRegionsOutsideTrials', true, @islogical);
            p.addParameter('exciseRegionsOutsideTrials', false, @islogical);
            p.addParameter('showDataTips', true, @islogical);
            p.parse(varargin{:});
            
            tsi = p.Results.tsi;
            if ~isempty(tsi) && (p.Results.maskRegionsOutsideTrials || p.Results.exciseRegionsOutsideTrials)
                mask = this.computeSpikeMaskWithinTrials(tsi);
                spike_inds = spike_inds(mask);
                spike_value = spike_value(mask);
            end
            spike_clusters =this.spike_clusters(spike_inds);
            spike_times = this.spike_times(spike_inds);
            
            if p.Results.exciseRegionsOutsideTrials
                timeShifts = tsi.computeShiftsExciseRegionsOutsideTrials();
                spike_times = timeShifts.shiftTimes(spike_times);
            end
            if p.Results.timeInSeconds
                spike_times = double(spike_times) / this.ks.fsAP; % convert to seconds
            end
            
            color = p.Results.color;
            size = p.Results.size;
            
            if isempty(color)
                colorArg  = {};
            else
                colorArg = {color};
            end
            h = scatter(spike_times, spike_value, size, colorArg{:}, 'filled');
            h.MarkerFaceAlpha = 0.4;
            h.MarkerEdgeAlpha = 0.7;
            if ~isempty(color)
                colormap(npxutils.internal.graphics.phy_cluster_colors());
            end
            
            valueLabel = p.Results.valueLabel;
            valueFormat = p.Results.valueFormat;
            valueUnits = p.Results.valueUnits;
            
            useDataTips = p.Results.showDataTips && ~verLessThan('matlab', '9.6.0');
            if useDataTips
                h.DataTipTemplate.DataTipRows(2).Label = valueLabel;
                if ~isempty(valueFormat) || ~isempty(valueUnits)
                    if isempty(valueUnits)
                        fmat = valueFormat;
                    else
                        fmat = [valueFormat ' ' valueUnits];
                    end
                    h.DataTipTemplate.DataTipRows(2).Format = fmat;
                end
                
                if p.Results.timeInSeconds
                    h.DataTipTemplate.DataTipRows(1).Label = 'Time';
                    h.DataTipTemplate.DataTipRows(1).Format = '%.03f sec';
                else
                    h.DataTipTemplate.DataTipRows(1).Label = 'Sample';
                    h.DataTipTemplate.DataTipRows(1).Format = '%d';
                end
                
                h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Cluster', double(spike_clusters), '%d');
                h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Spike Ind', double(spike_inds), '%d');
            end
            
            if p.Results.timeInSeconds
                xlabel('time (sec)')
            else
                xlabel('time (samples)')
            end
            if isempty(valueUnits)
                ylabel(valueLabel);
            else
                ylabel(sprintf('%s (%s)', valueLabel, valueUnits));
            end
            hold off;
            box off;
            axh = gca;
            axh.Color = [0.92 0.92 0.95];
            axh.GridColor = [1 1 1];
            axh.GridAlpha = 1;
            axh.MinorGridColor = [0.96 0.96 0.96];
            axh.MinorGridAlpha = 1;
            axh.MinorGridLineStyle = '-';
            axh.XGrid = 'on';
            axh.YGrid = 'on';
            axh.TickDir = 'out';
            set(gcf, 'InvertHardcopy', 'off');
            axis tight;
        end
    end
    
    methods % Plotting cluster waveforms
        
        function [clusterInds, cluster_ids] = lookup_clusterIds(this, cluster_ids)
            if islogical(cluster_ids)
                cluster_ids = this.cluster_ids(cluster_ids);
            end
            [tf, clusterInds] = ismember(cluster_ids, this.cluster_ids);
            assert(all(tf), 'Some cluster ids were not found in ds.clusterids');
        end
        
        function sortedChannelInds = lookup_sortedChannelIds(this, channel_ids)
            % lookup channel inds into the subset of sorted channels passed to kilosort
            [tf, sortedChannelInds] = ismember(channel_ids, this.channel_ids_sorted);
            assert(all(tf(:)), 'Some cluster ids were not found in m.cluster_ids');
        end
        
        function [channelInds, channel_ids] = lookup_channelIds(this, channel_ids)
            % lookup channel Inds in the full channelMap, not just sorted channels
            [channelInds, channel_ids] = this.channelMap.lookup_channelIds(channel_ids);
        end
        
        function cluster_templates = buildClusterTemplateScaled(this, varargin)
            p = inputParser();
            p.addParameter('cluster_ids', this.cluster_ids, @isvector);
            p.parse(varargin{:});
            
            clusterInds = this.lookup_clusterIds(p.Results.cluster_ids);
            templateInds = this.cluster_template_mostUsed(clusterInds);
            cluster_templates = this.template_scaled(templateInds, :, :);
        end
        
        function plotClusterImage(this, cluster_ids, varargin)
            clusterInds = this.lookup_clusterIds(cluster_ids);
            templateLists = this.cluster_template_list(clusterInds);
            templateInds = cat(1, templateLists{:});
            this.plotTemplateImage(templateInds,  varargin{:});
        end
        
        function channel_ids_by_template = plotTemplateImage(this, templateInds, varargin)
            channel_ids_by_template = this.plotTemplateImageInternal(templateInds, varargin{:});
        end
        
        function channel_ids_by_template = plotTemplateImageBatchwise(this, templateInds, varargin)
            channel_ids_by_template = this.plotTemplateImageInternal(templateInds, 'batchwise', true, varargin{:});
        end
        
        function cmap = getDefaultLinearColormap(this, N) %#ok<INUSL>
            %             cmap = npxutils.internal.graphics.colorcet('CBL2', 'N', N);
            cmap = npxutils.internal.graphics.cmocean('thermal', N);
        end
        
        function cmap = getDefaultBatchColormap(this, nBatch)
            if nargin < 2
                nBatch = this.nBatchesComputed;
            end
            cmap = npxutils.internal.graphics.colorcet('d6', 'N', nBatch);
            %             cmap = npxutils.internal.graphics.cmocean('haline', nBatch);
        end
        
        function cmap = getDefaultCategoricalColormap(this, nItems) %#ok<INUSL>
            if nargin < 2
                nItems = 10;
            end
            if nItems <= 10
                cmap = npxutils.internal.graphics.seaborn_color_palette('colorblind');
                cmap = cmap(1:nItems, :);
            else
                cmap = npxutils.internal.graphics.distinguishable_colors(nItems);
            end
        end
        
        function channel_ids_by_template = plotTemplateImageInternal(this, templateInds, varargin)
            p = inputParser();
            p.addParameter('xmag', 1.5, @isscalar);
            p.addParameter('ymag', 1.5, @isscalar);
            
            p.addParameter('cluster_ids', this.cluster_ids, @isvector);
            p.addParameter('template_colormap', [], @(x) true);
            
            p.addParameter('batchwise', false, @islogical);
            p.addParameter('batch_colormap', [], @(x) true);
            p.addParameter('colorbar', false, @islogical);
            
            p.addParameter('templateLabels', {}, @iscell);
            
            p.addParameter('plotCentroids', true, @islogical);
            p.addParameter('centroidSize', 8, @isscalar);
            
            % and ONE OR NONE of these to pick channels (or channels for each cluster)
            p.addParameter('channel_ids_by_template', [], @(x) isempty(x) || ismatrix(x));
            p.addParameter('best_n_channels', NaN, @isscalar); % or take the best n channels based on this clusters template when cluster_id is scalar
            p.addParameter('axh', [], @(x) true);
            p.parse(varargin{:});
            
            batchwise = p.Results.batchwise;
            plotCentroids = p.Results.plotCentroids;
            
            axh = p.Results.axh;
            if isempty(axh)
                axh = gca;
            end
            isholding = ishold(axh);
            
            % figure out actual channels requested
            if isfinite(p.Results.best_n_channels)
                channel_ids_by_template = this.template_best_channels(templateInds, 1:p.Results.best_n_channels);
            elseif ~isempty(p.Results.channel_ids_by_template)
                channel_ids_by_template = p.Results.channel_ids_by_template;
            else
                channel_ids_by_template = this.template_best_channels(templateInds, :);
            end
            channel_ind_by_template = this.lookup_sortedChannelIds(channel_ids_by_template);
            
            yspacing = this.channelMap.yspacing;
            xspacing = this.channelMap.xspacing;
            xmag = p.Results.xmag;
            ymag = p.Results.ymag;
            
            % plot relative time vector
            tvec = linspace(0, xspacing * xmag, size(this.template_scaled, 2));
            tvec_shift = tvec - mean(tvec);
            
            nTemp = numel(templateInds);
            nChannelsSorted = size(channel_ind_by_template, 2);
            
            % gather data
            if batchwise
                assert(this.has_computed_batchwise, 'Call m.computeBatchwiseMetrics() first');
                batches = this.batchesComputed;
                nBatchesPlotted = numel(batches);
                data = nan(nTemp, size(this.template_scaled, 2), nChannelsSorted, nBatchesPlotted, 'single');
                for iT = 1:nTemp
                    full_temp = this.construct_scaled_template_batchwise(templateInds(iT), 'batches', batches);
                    data(iT, :, :, :) = full_temp(1, :, channel_ind_by_template(iT, :), :);
                end
            else
                data = nan(nTemp, size(this.template_scaled, 2), nChannelsSorted, 'single');
                for iT = 1:nTemp
                    data(iT, :, :) = this.template_scaled(templateInds(iT), :, channel_ind_by_template(iT, :));
                end
            end
            data = data - mean(data(:, :), 2); % center over time
            data = data ./ (max(data(:)) - min(data(:))) * yspacing * ymag; % normalize amplitudes
            
            if p.Results.batchwise
                batch_cmap = p.Results.batch_colormap;
                if isempty(batch_cmap)
                    batch_cmap = this.getDefaultBatchColormap();
                elseif isa(batch_cmap, 'function_handle')
                    batch_cmap = batch_cmap(nBatchesPlotted);
                end
            else
                template_cmap = p.Results.template_colormap;
                if isempty(template_cmap)
                    template_cmap = this.getDefaultCategoricalColormap(nTemp);
                elseif isa(template_cmap, 'function_handle')
                    template_cmap = template_cmap(nTemp);
                end
            end
            
            templateLabels = p.Results.templateLabels;
            if isempty(templateLabels)
                templateLabels = arrayfun(@(ind) sprintf("template %d", ind), templateInds);
            end
            
            %             if batchwise
            %                 if ~isholding
            %                     cla(axh);
            %                 end
            %                 axh.ColorOrder = batch_cmap;
            %                 axh.ColorOrderIndex = 1;
            %             end
            %
            for iT = 1:nTemp
                this_channel_ind = channel_ind_by_template(iT, :);
                %                 this_channel_ids = channel_ids_by_template(iT, :);
                xc = this.ks.channel_positions_sorted(this_channel_ind, 1);
                yc = this.ks.channel_positions_sorted(this_channel_ind, 2);
                
                for iC = 1:nChannelsSorted
                    if batchwise
                        waves = npxutils.internal.TensorUtils.squeezeDims(data(iT, :, iC, :), [1 3]) + yc(iC);
                        h = plot(axh, tvec_shift + xc(iC), waves, 'LineWidth', 0.5);
                        for iH = 1:numel(h)
                            h(iH).Color = batch_cmap(iH, :);
                        end
                    else
                        wave = npxutils.internal.TensorUtils.squeezeDims(data(iT, :, iC), 1) + yc(iC);
                        
                        %                     ud = struct('template_ind', template_inds(iT), 'template_amplitude', sprintf('%.1f uV', m.template_amplitude(template_inds(iT)),
                        %                     'channel_id', this_channel_ids(iC), 'template_is_localized', m.template_is_localized(template_inds(iT)), ...
                        %                     'xname', 'Time', 'yname', 'Voltage', 'xoffset', xc(iC) - mean(tvec), 'yoffset', yc(iC) + dataOffset(iT, 1, iC), 'xscale', 1, 'yscale', waveScalingFactor_umtouV, 'xunits', 'ms', 'yunits', 'uV');
                        %
                        h = plot(axh, tvec_shift + xc(iC), wave, 'Color', template_cmap(iT, :), 'LineWidth', 0.5);
                        
                    end
                    if iC == 1
                        npxutils.internal.graphics.showFirstInLegend(h, templateLabels{iT});
                    else
                        npxutils.internal.graphics.hideInLegend(h);
                    end
                    hold(axh, 'on');
                end
                
                if plotCentroids
                    if batchwise
                        centroid = npxutils.internal.TensorUtils.squeezeDims(this.template_centroid_batchwise(templateInds(iT), :, [1 2]), 1);
                        plot(axh, centroid(:, 1), centroid(:, 2), '-', 'Color', [0.8 0.8 0.8]);
                        scatter(axh, centroid(:, 1), centroid(:, 2), p.Results.centroidSize^2, batch_cmap, '+');
                        
                    else
                        centroid = this.template_centroid(templateInds(iT), [1 2]);
                        plot(axh, centroid(1), centroid(2), '+', 'MarkerSize', p.Results.centroidSize, 'Color', template_cmap(iT, :));
                    end
                end
                
            end
            
            axis(axh, 'off', 'tight');
            if batchwise
                colormap(axh, batch_cmap);
                caxis(axh, [1 this.nBatches]);
                if p.Results.colorbar
                    h = colorbar(axh);
                    h.YDir = 'reverse';
                    h.Ticks = [1 this.nBatches];
                    h.TickLabels = ["Early", "Late"];
                    ylabel(h, 'Batches');
                end
            end
            %             axh.Color = [0.92 0.92 0.95];
            box(axh, 'off');
            if ~isholding, hold(axh, 'off'), end
        end
        
        function plotTemplateCentroids(this, template_inds, varargin)
            p = inputParser();
            p.addParameter('batchwise', false, @islogical);
            p.addParameter('batch_colormap', [], @(x) true);
            p.addParameter('template_colormap', [], @(x) true);
            p.addParameter('colorbar', false, @islogical);
            p.addParameter('centroidSize', 8, @isscalar);
            p.addParameter('batchMinSpikes', 1, @isscalar)
            p.addParameter('axh', [], @(x) true);
            p.addParameter('ignoreDriftX', false, @islogical);
            p.parse(varargin{:})
            
            axh = p.Results.axh;
            if isempty(axh), axh = gca; end
            isholding = ishold(axh);
            
            if nargin < 2 || isempty(template_inds)
                template_inds = find(this.template_is_localized);
            end
            
            if p.Results.batchwise
                batch_cmap = p.Results.batch_colormap;
                if isempty(batch_cmap)
                    batch_cmap = this.getDefaultBatchColormap();
                elseif isa(batch_cmap, 'function_handle')
                    batch_cmap = batch_cmap(nBatchesPlotted);
                end
            else
                template_cmap = p.Results.template_colormap;
                if isempty(template_cmap)
                    template_cmap = this.getDefaultCategoricalColormap(numel(template_inds));
                elseif isa(template_cmap, 'function_handle')
                    template_cmap = template_cmap(nTemp);
                end
            end
            
            if p.Results.batchwise
                for iT = 1:numel(template_inds)
                    mask_valid = this.template_useCount_batchwise(template_inds(iT), :) >= p.Results.batchMinSpikes;
                    centroid = npxutils.internal.TensorUtils.squeezeDims(this.template_centroid_batchwise(template_inds(iT), :, [1 2]), 1);
                    
                    xv = centroid(mask_valid, 1);
                    if p.Results.ignoreDriftX && ~isempty(xv)
                        xv(2:end) = xv(1);
                    end
                    plot(axh, xv, centroid(mask_valid, 2), '-', 'Color', [0.8 0.8 0.8]);
                    scatter(axh, xv, centroid(mask_valid, 2), p.Results.centroidSize^2, batch_cmap(mask_valid, :), 'filled');
                    hold(axh, 'on');
                end
            else
                for iT = 1:numel(template_inds)
                    centroid = this.template_centroid(template_inds(iT), [1 2]);
                    plot(axh, centroid(1), centroid(2), '+', 'MarkerSize', p.Results.centroidSize, 'Color', template_cmap(iT, :));
                    hold(axh, 'on');
                end
            end
            
            if ~isholding, hold(axh, 'off'); end
            
        end
        
        function plotClusterCentroids(this, cluster_ids, varargin)
            p = inputParser();
            p.addParameter('batchwise', false, @islogical);
            p.addParameter('batch_colormap', [], @(x) true);
            p.addParameter('cluster_colormap', [], @(x) true);
            p.addParameter('colorbar', false, @islogical);
            p.addParameter('centroidSize', 8, @isscalar);
            p.addParameter('batchMinSpikes', 1, @isscalar)
            p.addParameter('axh', [], @(x) true);
            p.addParameter('ignoreDriftX', false, @islogical);
            p.parse(varargin{:})
            
            axh = p.Results.axh;
            if isempty(axh), axh = gca; end
            isholding = ishold(axh);
            
            cluster_inds = this.lookup_clusterIds(cluster_ids);
            
            if p.Results.batchwise
                batch_cmap = p.Results.batch_colormap;
                if isempty(batch_cmap)
                    batch_cmap = this.getDefaultBatchColormap();
                elseif isa(batch_cmap, 'function_handle')
                    batch_cmap = batch_cmap(nBatchesPlotted);
                end
            else
                cluster_cmap = p.Results.cluster_colormap;
                if isempty(cluster_cmap)
                    cluster_cmap = this.getDefaultCategoricalColormap(numel(cluster_inds));
                elseif isa(cluster_cmap, 'function_handle')
                    cluster_cmap = cluster_cmap(nTemp);
                end
            end
            
            cluster_count_batchwise = sum(this.cluster_template_useCount_batchwise, 2);
            if p.Results.batchwise
                for iT = 1:numel(cluster_inds)
                    mask_valid = cluster_count_batchwise(cluster_inds(iT), 1, :) >= p.Results.batchMinSpikes;
                    centroid = npxutils.internal.TensorUtils.squeezeDims(this.cluster_centroid_batchwise(cluster_inds(iT), :, [1 2]), 1);
                    
                    xv = centroid(mask_valid, 1);
                    if p.Results.ignoreDriftX && ~isempty(xv)
                        xv(2:end) = xv(1);
                    end
                    plot(axh, xv, centroid(mask_valid, 2), '-', 'Color', [0.8 0.8 0.8]);
                    scatter(axh, xv, centroid(mask_valid, 2), p.Results.centroidSize^2, batch_cmap(mask_valid, :), 'filled');
                    hold(axh, 'on');
                end
            else
                for iT = 1:numel(cluster_inds)
                    centroid = this.cluster_centroid(cluster_inds(iT), [1 2]);
                    plot(axh, centroid(1), centroid(2), '+', 'MarkerSize', p.Results.centroidSize, 'Color', cluster_cmap(iT, :));
                    hold(axh, 'on');
                end
            end
            
            if ~isholding, hold(axh, 'off'); end
        end
        
        function plotClusterCentroidsBatchwiseReasonablyLocalized(this, varargin)
            p = inputParser();
            p.addParameter('maxDrift', 80, @isscalar);
            p.parse(varargin{:});
            
            cluster_mask = this.cluster_is_localized;
            
            % clusters x batches x spatial dims
            maxDriftByCluster = this.computeClusterDriftDistance();
            cluster_mask = cluster_mask & maxDriftByCluster <= p.Results.maxDrift;
            
            this.plotClusterCentroids(this.cluster_ids(cluster_mask), 'batchwise', true, varargin{:});
        end
        
        function plotClusterDriftSummary(this, varargin)
            p = inputParser();
            p.addParameter('spatialDim', 2, @isscalar); % 1 is x, 2 is y, 3 is z
            p.addParameter('referenceBatchInd', 1, @isscalar);
            p.addParameter('maxDrift', 80, @isscalar);
            p.addParameter('binWidth', 200, @isscalar);
            p.addParameter('colormap', [], @isscalar);
            p.addParameter('smoothBy', 10, @iscalar);
            p.addParameter('compressInitialOffsetsBy', 10, @isscalar);
            p.parse(varargin{:});
            
            this.computeBatchwiseMetrics(); % ensure batchwise metrics are computed
            
            sdim = p.Results.spatialDim;
            referenceBatchInd = p.Results.referenceBatchInd;
            
            cluster_mask = this.cluster_is_localized;
            
            % clusters x batches x spatial dims
            maxDriftByCluster = this.computeClusterDriftDistance();
            cluster_mask = cluster_mask & maxDriftByCluster <= p.Results.maxDrift;
            
            % fill in missing batches with NaNs but don't extrapolate edges
            centroids = fillmissing(this.cluster_centroid_batchwise, 'previous', 2, 'EndValues', 'none');
            
            % grab initial position at reference batch to define the binning
            cluster_mask = cluster_mask & ~isnan(centroids(:, referenceBatchInd, sdim));
            ref_pos = centroids(cluster_mask, referenceBatchInd, sdim);
            
            [~, ~, bins] = histcounts(ref_pos, 'BinWidth', p.Results.binWidth);
            nBins = max(bins);
            binmedian_ref_pos = accumarray(bins, ref_pos, [nBins, 1], @median);
            
            compressed_offsets = (binmedian_ref_pos - mean(binmedian_ref_pos)) / p.Results.compressInitialOffsetsBy;
            compressed_offsets = compressed_offsets(bins);
            
            pos_shift = centroids(cluster_mask, :, sdim) - ref_pos + compressed_offsets;
            
            % now histogram each time slice by bin
            
            
            binnedQuantiles = nan(nBins, size(centroids, 2), 3);
            for iB = 1:nBins
                binnedQuantiles(iB, :, :) = quantile(pos_shift(bins == iB, :), [0.25 0.5 0.75], 1)';
            end
            
            cmap = p.Results.colormap;
            if isempty(cmap)
                cmap = this.getDefaultLinearColormap(nBins);
            end
            
            tvals = this.ks.batch_starts(this.batchesComputed) / this.fs;
            for iB = 1:nBins
                %                  TrialDataUtilities.Plotting.errorshadeInterval(m.batchesComputed, binnedQuantiles(iB, :, 1), binnedQuantiles(iB, :, 3), cmap(iB, :));
                %                hold on
                yvals = binnedQuantiles(iB, :, 2);
                if p.Results.smoothBy > 1
                    yvals = smooth(yvals, p.Results.smoothBy);
                end
                plot(tvals, yvals, '-', 'Color', cmap(iB, :), 'LineWidth', 2);
                hold on
            end
            
            axis tight
            aa = AutoAxis.replaceScaleBars('xUnits', 'sec', 'yUnits', 'um');
            aa.axisMarginLeft = aa.axisMarginRight;
            grid on
            aa.update();
        end
        
        function plotRecordingSites(this, varargin)
            p = inputParser();
            p.addParameter('channel_ids', this.channel_ids_sorted, @isvector)
            p.addParameter('showChannelLabels', false, @islogical);
            p.addParameter('labelArgs', {}, @iscell);
            p.addParameter('color', [0.8 0.8 0.8], @(x) true);
            p.addParameter('markerSize', 25, @isscalar);
            p.parse(varargin{:});
            
            [channelInds, channel_ids] = this.lookup_channelIds(p.Results.channel_ids);
            
            xc = this.channelMap.xcoords(channelInds);
            yc = this.channelMap.ycoords(channelInds);
            plot(xc, yc, '.', 'Marker', 's', 'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', p.Results.color, 'MarkerSize', p.Results.markerSize);
            if p.Results.showChannelLabels
                for iC = 1:numel(channelInds)
                    text(xc(iC), yc(iC), sprintf('ch %d', channel_ids(iC)), ...
                        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                        'Background', 'none', ...
                        p.Results.labelArgs{:});
                end
            end
        end
        
        function plotClusterWaveformAtCentroid(this, varargin)
            p = inputParser();
            p.addParameter('waveformScale', 10, @isscalar);
            p.addParameter('waveformWidth', this.channelMap.xspacing/20, @isscalar);
            p.addParameter('waveformHeight', this.channelMap.yspacing*3, @isscalar);
            p.addParameter('cluster_ids', this.cluster_ids, @isvector);
            p.addParameter('colormap', [], @(x) isempty(x) || ismatrix(x));
            p.addParameter('plotRecordingSites', false, @islogical);
            p.addParameter('useAutoAxis', false, @islogical);
            p.addParameter('LineWidth', 1, @isscalar);
            p.addParameter('recordingSitesMarkerSize', 5, @isscalar);
            p.addParameter('recordingSitesColor', [0.5 0.5 0.5], @(x) true);
            p.parse(varargin{:});
            
            cluster_ids = p.Results.cluster_ids;
            clusterInds = this.lookup_clusterIds(cluster_ids);
            waves = this.cluster_waveform(clusterInds, :, 1); % nClusters x nTimepoints
            
            xvec = linspace(-p.Results.waveformWidth/2, p.Results.waveformWidth/2, size(waves, 2)) * p.Results.waveformScale;
            timeScaleFactor_umtoms = (size(waves, 2) / this.fs * 1000) / range(xvec);
            
            waveScalingFactor_umtouV = double((max(waves(:)) - min(waves(:))) / p.Results.waveformHeight / p.Results.waveformScale);
            waves = waves ./ waveScalingFactor_umtouV;
            
            colormap = p.Results.colormap;
            if isempty(colormap)
                %                 if numel(cluster_ids) < 20
                if exist('turbo', 'file')
                    colormap = turbo(numel(clusterInds));
                else
                    colormap = npxutils.internal.graphics.turbomap(numel(clusterInds));
                end
                %                     colormap = cbrewer('div', 'Spectral', numel(clusterInds));
                %                 % color by cluster amplitude
                %                 colormap = npxutils.internal.graphics.cmocean('haline', numel(clusterInds));
                % %                 colormap = npxutils.internal.evalColorMapAt(colormap, linspace(0, 1, numel(clusterInds)));
                %                 [~, sortIdx] = sort(m.cluster_amplitude(clusterInds), 'ascend');
                %                 [~, cmapSort] = ismember(1:numel(clusterInds), sortIdx);
                %                 colormap = colormap(cmapSort, :);
            end
            
            % nClusters x 2 or 3
            com = this.cluster_centroid(clusterInds, :);
            
            newplot
            holding = ishold;
            
            % plot recording sites
            if p.Results.plotRecordingSites
                this.plotRecordingSites('color', p.Results.recordingSitesColor, 'markerSize', p.Results.recordingSitesMarkerSize);
            end
            ax = gca;
            %ax.Color = [0.92 0.92 0.95];
            ax.TickDir = 'out';
            ax.YDir = 'normal';
            
            hold on
            
            xlim(double([min(com(:,1)) - this.channelMap.xspacing/2, max(com(:,1)) + this.channelMap.xspacing/2]));
            ylim(double([min(com(:,2)) - this.channelMap.yspacing*5, max(com(:,2)) + this.channelMap.yspacing*5]));
            
            for iC = 1:size(waves, 1)
                color = colormap(iC, :);
                
                ud = struct('cluster_id', cluster_ids(iC), 'cluster_amplitude', sprintf('%.1f uV', this.cluster_amplitude(clusterInds(iC))), ...
                    'cluster_is_localized', this.cluster_is_localized(clusterInds(iC)), ...
                    'xname', 'Time', 'yname', 'Voltage', 'xoffset', xvec(1) + com(iC, 1), 'yoffset', com(iC, 2), 'xscale', timeScaleFactor_umtoms, 'yscale', waveScalingFactor_umtouV, 'xunits', 'ms', 'yunits', 'uV');
                plot(xvec + com(iC, 1), waves(iC, :) + com(iC, 2), '-', 'Color', color, 'UserData', ud, 'LineWidth', p.Results.LineWidth);
                hold on;
            end
            
            if ~holding, hold off; end
            
            if p.Results.useAutoAxis
                aa = AutoAxis.replaceScaleBars('xUnits', 'ms', 'xLength', 2, 'xScaleFactor', timeScaleFactor_umtoms, ...
                    'yUnits', 'uV', 'yLength', 200, 'yScaleFactor', waveScalingFactor_umtouV);
                xlabel('x (um)', 'BackgroundColor', 'none');
                ylabel('y (um)', 'BackgroundColor', 'none');
                aa.addAutoAxisX();
                aa.addAutoAxisY();
                aa.axisPaddingRight = 0.4;
                aa.axisMarginLeft = 1;
                aa.axisMarginRight = 1;
            else
                xlabel('x (μm)', 'BackgroundColor', 'none');
                ylabel('y (μm)', 'BackgroundColor', 'none');
            end
            %
            %             bgcolor = [0.9 0.9 0.9];
            %             axh = gca;
            %             axh.Color = bgcolor;
            box off;
            
            npxutils.internal.configureDataTipsFromUserData(gcf);
            if p.Results.useAutoAxis
                aa.update();
            end
            hold off;
        end
        
        function ch_ids_sorted = plotClusterHeatmap(this, cluster_id, varargin)
            cluster_ind = this.lookup_clusterIds(cluster_id);
            
            template_inds = cat(1, this.cluster_template_list{cluster_ind});
            ch_ids_sorted = this.plotTemplateHeatmap(template_inds, varargin{:});
        end
        
        function ch_ids_sorted = plotTemplateHeatmap(this, template_inds, varargin)
            % plots one snippet as a heatmap with each of the clusters that occur during that period
            
            p = inputParser();
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.addParameter('reconstruct_best_channels', 24, @isscalar);
            p.addParameter('successive_residuals', false, @islogical);
            p.addParameter('label_templates', [], @(x) isempty(x) || islogical(x));
            p.parse(varargin{:});
            
            if p.Results.timeInMilliseconds
                time = this.templateTimeRelativeMs;
            else
                time = this.templateTimeRelative;
            end
            
            nChannels = p.Results.reconstruct_best_channels;
            ch_ids_by_proximity = this.template_best_channels(template_inds(1), 1:nChannels);
            ch_ids_sorted = this.channelMap.sortChannelsVertically(ch_ids_by_proximity);
            ch_inds_sorted = this.lookup_sortedChannelIds(ch_ids_sorted);
            data = this.template_scaled(template_inds, :, ch_inds_sorted);
            
            % nTemplates x nTemplates x nBestChannels --> (nBestChannels*nTemplates) x nTime
            data_stacked = npxutils.internal.TensorUtils.reshapeByConcatenatingDims(data, {[1 3], 2});
            
            npxutils.internal.pmatbal(data_stacked, 'x', time);
            
            label_templates = p.Results.label_templates;
            if isempty(label_templates)
                label_templates = numel(template_inds) > 1;
            end
            if label_templates
                for iC = 1:size(data, 1)
                    yline(nChannels * (iC-1) + 0.5, 'k-', sprintf('template %u', template_inds(iC)), 'LabelVerticalAlignment', 'bottom');
                end
            end
            
        end
    end
    
    methods % Pairwise cluster comparison
        function ss = plot_spike_heatmap(this, spike_ind, varargin)
            ss = this.ks.getWaveformsFromRawData('spike_idx', spike_ind, 'best_n_channels', 24, varargin{:});
            ss.plotHeatmapWithTemplates(1);
        end
        
        function ss = plot_spike_overlay(this, spike_ind, varargin)
            p = inputParser();
            p.addParameter('overlay_all_clusters', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:})
            ss = this.ks.getWaveformsFromRawData('spike_idx', spike_ind, 'best_n_channels', 24, p.Unmatched);
            
            if p.Results.overlay_all_clusters
                ss.overlay_cluster_ids = this.ks.cluster_ids;
            end
            
            ss.plotStackedTracesWithOverlays('maskSnippet', 1, 'overlay_templates', true);
        end
        
        function ss = extractSnippets_clusterInWindow(this, cluster_id, filter_window, varargin)
            ss = this.ks.getWaveformsFromRawData('cluster_ids', cluster_id, 'filter_window', filter_window, 'best_n_channels', 24, varargin{:});
        end
        
        function [times1, times2, lags] = pairwiseClusterFindSpikesWithLag(this, cluster_ids, varargin)
            % find spikes with a
            p = inputParser();
            p.addParameter('lagWindowMs', [-5 5], @isvector); % in ms
            p.addParameter('N', Inf, @isscalar);
            p.addParameter('sortClosestToCentralLag', false, @isscalar);
            p.parse(varargin{:});
            lagWindow = p.Results.lagWindowMs / 1000 * this.fs; % convert to samples
            
            assert(numel(cluster_ids) == 2);
            [~, which_cluster] = ismember(this.spike_clusters, cluster_ids);
            
            % now we want to find the set of pairs of spikes where spikes
            times1 = double(this.spike_times(which_cluster == 1));
            times2 = double(this.spike_times(which_cluster == 2));
            
            lagCentral = (lagWindow(2) + lagWindow(1)) / 2;
            lagRadius = (lagWindow(2) - lagWindow(1)) / 2;
            inds2 = rangesearch(times2 - lagCentral, times1, lagRadius);
            
            [inds2, inds1] = npxutils.internal.TensorUtils.catWhich(2, inds2{:});
            inds1 = inds1';
            inds2 = inds2';
            
            times1 = times1(inds1);
            times2 = times2(inds2);
            lags = (times2 - times1) / this.fs * 1000;
            
            if p.Results.sortClosestToCentralLag
                [~, sort_idx] = sort(abs(lags - lagCentral), 'ascend');
                times1 = times1(sort_idx);
                times2 = times2(sort_idx);
                lags = lags(sort_idx);
            end
            
            if ~isinf(p.Results.N)
                select = 1:p.Results.N;
                times1 = times1(select);
                times2 = times2(select);
                lags = lags(select);
            end
        end
        
        function [times1, times2, lags] = findSpikePairsWithAutoLag(this, cluster_id, varargin)
            % find spike pairs that live in a certain bin on an autocorrelogram
            p = inputParser();
            p.addParameter('lagWindowMs', [0 5], @isvector); % in ms
            p.addParameter('N', Inf, @isscalar);
            p.addParameter('sortSmallestLag', false, @isscalar);
            p.parse(varargin{:});
            lagWindow = p.Results.lagWindowMs / 1000 * this.fs; % convert to samples
            if isscalar(lagWindow)
                lagWindow = [0 lagWindow];
            end
            
            assert(numel(cluster_id) == 1);
            [~, which_cluster] = ismember(this.spike_clusters, cluster_id);
            
            % now we want to find the set of pairs of spikes where spikes
            times = double(this.spike_times(which_cluster == 1));
            dt = diff(times);
            inds1 = find(dt >= lagWindow(1) & dt <= lagWindow(2));
            inds2 = inds1 + 1;
            
            times1 = times(inds1);
            times2 = times(inds2);
            lags = (times2 - times1) / this.fs * 1000;
            
            if p.Results.sortSmallestLag
                [~, sort_idx] = sort(abs(lags), 'ascend');
                times1 = times1(sort_idx);
                times2 = times2(sort_idx);
                lags = lags(sort_idx);
            end
            
            if ~isinf(p.Results.N)
                select = 1:p.Results.N;
                times1 = times1(select);
                times2 = times2(select);
                lags = lags(select);
            end
        end
        
        function [snippetSet, lags] = extractSnippets_pairwiseClusterFindSpikesWithLag(this, cluster_ids, varargin)
            p = inputParser();
            p.addParameter('lagWindowMs', [-5 5], @isvector);
            p.addParameter('N', Inf, @isscalar);
            p.addParameter('sortClosestToCentralLag', false, @isscalar);
            p.addParameter('best_n_channels', 24, @iscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            args = rmfield(p.Results, 'best_n_channels');
            [times1, times2, lags] = this.pairwiseClusterFindSpikesWithLag(cluster_ids, args);
            
            window_width = max(int64(times2) - int64(times1)) + int64(this.nTemplateTimepoints);
            window = [-window_width / 2, window_width / 2];
            
            times = (int64(times1) + int64(times2)) / int64(2);
            
            snippetSet = this.ks.readAPSnippsetSet_clusterIdSubset(times, window, cluster_ids, ...
                'best_n_channels', p.Results.best_n_channels, p.Unmatched);
        end
        
        function [snippetSet, lags] = extractSnippets_clusterFindSpikesWithAutoLag(this, cluster_id, varargin)
            p = inputParser();
            p.addParameter('lagWindowMs', [0 5], @isvector);
            p.addParameter('N', Inf, @isscalar);
            p.addParameter('sortSmallestLag', false, @isscalar);
            p.addParameter('best_n_channels', 24, @iscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            args = rmfield(p.Results, 'best_n_channels');
            [times1, times2, lags] = this.findSpikePairsWithAutoLag(cluster_id, args);
            
            window_width = max(int64(times2) - int64(times1)) + int64(this.nTemplateTimepoints);
            window = [-window_width / 2, window_width / 2];
            
            times = (int64(times1) + int64(times2)) / int64(2);
            
            snippetSet = this.ks.readAPSnippsetSet_clusterIdSubset(times, window, cluster_id, ...
                'best_n_channels', p.Results.best_n_channels, p.Unmatched);
        end
        
        function [timesByCluster, indsByCluster, windowWidth] = ...
                multiClusterFindSpikesWithinWindow(this, cluster_ids, varargin)
            p = inputParser();
            p.addParameter('window', 20, @isvector);
            p.addParameter('N', Inf, @isscalar);
            p.addParameter('sortSmallestWindow', true, @isscalar);
            p.parse(varargin{:});
            window = p.Results.window * this.fs; % convert to samples
            nClusters = numel(cluster_ids);
            
            % take times for the first cluster and find times for all other clusters within +/- window
            % which we will narrow down later
            [~, which_cluster] = ismember(this.spike_clusters, cluster_ids);
            times1 = double(this.spike_times(which_cluster == 1));
            mask_other = which_cluster > 1;
            indsOther = find(mask_other);
            timesOther = double(this.spike_times(mask_other));
            indsIntoTimesOtherNearest = rangesearch(timesOther, times1, window, 'SortIndices', false);
            
            valid = false(numel(times1), 1);
            windowWidth = nan(numel(times1), 1);
            indsByCluster = nan(numel(times1), nClusters);
            indsByCluster(:, 1) = find(which_cluster == 1);
            
            for i = 1:numel(times1)
                % find the smallest interval that includes times1(i) and one spike from every other cluster
                % this is the set cover problem on the integers
                
                indsOtherThis = indsOther(indsIntoTimesOtherNearest{i});
                whichClusterThis = which_cluster(indsOtherThis);
                relTimesOtherThis = timesOther(indsIntoTimesOtherNearest{i}) - times1(i);
                
                [valid(i), tempInds, windowWidth(i)] = bestWindowCover(relTimesOtherThis, whichClusterThis);
                indsByCluster(i, 2:end) = indsOtherThis(tempInds);
            end
            
            indsByCluster = indsByCluster(valid, :);
            windowWidth = windowWidth(valid);
            
            if p.Results.sortSmallestWindow
                [windowWidth, sortIdx] = sort(windowWidth, 'ascend');
                indsByCluster = indsByCluster(sortIdx, :);
            end
            timesByCluster = this.spike_times(indsByCluster);
            windowWidth = windowWidth / this.fs;
            
            function [valid, indsByCluster, windowWidth] = bestWindowCover(relativeTimes, whichCluster)
                % given a set of nSpikes relative times between +/- window and cluster identities whichCluster
                % return the ind from each cluster that fall within the smallest max - min window
                
                
                % offsets defines time windows running from offsets(i):offsets(i)+window-1
                % indMat temporary assignment of this spike for a given cluster and a given window start
                % deltaMat is the minimum window required to include that spike (absolute value)
                offsets = [relativeTimes(relativeTimes < 0)', 0];
                [indMat, deltaMat] = deal(nan(nClusters-1, numel(offsets)));
                
                for iS = 1:numel(relativeTimes)
                    % update the indByStartByCluster using this spike if its better than what's there before
                    t = relativeTimes(iS);
                    c = whichCluster(iS) - 1;
                    
                    % this spike lies in the window if t >= offsets (the window start) and t <= offsets + window - 1 (window end)
                    % and the minimum window width to contain it is abs(t)
                    mask_update = t >= offsets & t < offsets + window & ~(deltaMat(c, :) < abs(t));
                    deltaMat(c, mask_update) = abs(t);
                    indMat(c, mask_update) = iS;
                end
                
                % now pick the best window
                [windowWidth, whichOffset] = min(max(deltaMat, [], 1, 'includenan'));
                if isnan(windowWidth)
                    valid = false;
                    indsByCluster = nan(nClusters, 1);
                else
                    valid = true;
                    indsByCluster = indMat(:, whichOffset);
                end
            end
        end
        
        function out = plotSampleIdx_highlightClusterWaveforms(this, sampleIdx, cluster_ids, varargin)
            p = inputParser();
            p.addParameter('channel_ids', this.channelMap.channelIdsMapped, @isvector);
            p.addParameter('channel_ids_by_cluster', [], @ismatrix); % specify this manually
            p.addParameter('best_n_channels', 24, @isscalar); % or take the best n channels based on this clusters template when cluster_id is scalar
            
            p.addParameter('invertChannels', true, @islogical);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 0.95, @isscalar);
            p.addParameter('car', false, @islogical);
            p.addParameter('downsample',1, @isscalar);
            p.addParameter('timeInSeconds', true, @islogical);
            
            p.addParameter('waveform_window', [-40 41], @isvector); % Number of samples before and after spiketime to include in waveform
            p.addParameter('cluster_colormap', npxutils.internal.graphics.distinguishable_colors(numel(cluster_ids)), @(x) ismatrix(x));
            p.parse(varargin{:});
            waveform_window = p.Results.waveform_window;
            
            df = this.ks.raw_dataset;
            sampleIdx = npxutils.internal.makecol(sampleIdx);
            mat = df.readAP_idx(sampleIdx); % C x T
            
            [channelInds, channelIds] = df.lookup_channelIds(p.Results.channel_ids);
            channelsHighlightByCluster = p.Results.channel_ids_by_cluster;
            clusterInds = this.lookup_clusterIds(cluster_ids);
            if isempty(channelsHighlightByCluster)
                n_best = min(p.Results.best_n_channels, size(this.cluster_best_channels, 2));
                channelsHighlightByCluster = this.cluster_best_channels(clusterInds, 1:n_best);
            end
            
            mat = mat(channelInds, :);
            labels = df.channelNamesPadded(channelInds);
            
            if p.Results.downsample > 1
                mat = mat(:, 1:p.Results.downsample:end);
                sampleIdx = sampleIdx(1:p.Results.downsample:end);
            end
            mat = double(mat);
            if p.Results.car
                mat = mat - median(mat, 1);
            end
            
            if ~p.Results.showLabels
                labels = [];
            end
            
            if p.Results.timeInSeconds
                time = double(sampleIdx) / df.fsAP;
            else
                time = sampleIdx;
            end
            
            % build overlay maps for each cluster
            waveformOverlayLabels = zeros(size(mat));
            for iC = 1:numel(cluster_ids)
                channelsIdsHighlight = channelsHighlightByCluster(iC, :);
                [tf, channelIndsHighlight] = ismember(channelsIdsHighlight, channelIds);
                cinds = channelIndsHighlight(tf);
                
                % nTimes x 1
                times = single(this.spike_times(this.spike_clusters == cluster_ids(iC)));
                
                % nSamples x nTimes --> nSamples
                within_window = any(sampleIdx >= times' + waveform_window(1) & sampleIdx <= times' + waveform_window(2), 2);
                
                waveformOverlayLabels(cinds, within_window) = iC;
            end
            
            out = npxutils.internal.plotStackedTraces(time, mat', 'labels', labels, ...
                'gain', p.Results.gain, 'invertChannels', p.Results.invertChannels, 'normalizeEach', false, ...
                'colorOverlayLabels', waveformOverlayLabels', 'colorOverlayMap', p.Results.cluster_colormap);
        end
        
        
    end
    
    methods(Static) % multi plotting to simulate concatenated files
        function multiple_plotDriftmap(mSet, varargin)
            p = inputParser();
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'npxutils.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('xGap', 0, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            xGap = p.Results.xGap;
            
            holding = ishold;
            
            N = numel(mSet);
            tsiSet = p.Results.tsi;
            assert(numel(tsiSet) == N);
            
            % compute global amplitude range
            ampRange = quantile(cat(1, mSet.spike_amplitude), [0.1 0.9]);
            
            xOffset = 0;
            for i = 1:N
                info = mSet(i).plotDriftmap('tsi', tsiSet(i), 'xOffset', xOffset, 'ampRange', ampRange, p.Unmatched);
                hold on;
                xOffset = info.xMax + xGap;
            end
            if ~holding
                hold off
            end
        end
    end
    
end