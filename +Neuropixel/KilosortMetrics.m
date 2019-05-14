classdef KilosortMetrics < handle
    % this class stores a set of computed information about a KiloSortDataset,
    % particularly regarding template shapes, depths, and spiking drift
    % Note that this code borrows heavily from the math in github.com/cortex-lab/spikes
    
    properties(Transient) % not saved with metrics
        ks % Neuropixel.KiloSortDataset
    end
    
    properties(Dependent)
        nSpikes
        nChannelsSorted
        nTemplates
        nTemplateTimepoints
        nClusters
        nConcatenatedFiles
        maxTemplatesPerCluster
    end
    
    properties
        % copied over from ks
        fs
        channelMap
        channel_ids
        concatenatedSamples(:, 1) uint64
        concatenatedStarts(:, 1) uint64
        concatenatedNames(:, 1) string
        
        % per template properties
        
        % nTemplates x nTimePoints x nTemplateChannels
        template_unw single % unwhitened scaled templates
        template_scaled single % unwhitened, scaled templates
        
        % nTemplates x nSpatialCoord
        template_centerOfMass % single x,y,z, coords for each template's center of mass
        
        % nTemplates x 1
        template_is_localized(:, 1) logical
        
        % nTemplates x nTimepoints
        template_waveform single % single, maximum amplitude waveform
        
        % nTemplates x 1
        template_waveform_ch(:, 1) uint32 % channel where the maximum amplitude part of the template lives
        template_amplitude(:, 1) single % in uV
        template_ttp(:, 1) single % trough to peak
        
        % nTemplates x nChannelsSorted
        template_best_channels uint32 % nTemplates x nChannelsSorted matrix indicating the closest channels to the max (typically take first 20 cols)
        
        % per spike properties
        
        % nSpikes x 1
        spike_times(:,1) uint64
        
        % nSpikes x 1
        spike_amplitude(:, 1) single % in Uv,  product of template's largest amplitude over channels and scaling of each template (by ks.amplitudes)
        
        % nSpikes x nSpatialCoord
        spike_centerOfMass single % single x,y,z coords for each spike's center of mass based on private PCs
        
        % nSpikes x 1
        spike_templates(:, 1) uint32
        spike_clusters(:, 1) uint32
        
        % per template properties
        
        % unique cluster ids corresponding to each slot in the other cluster properties
        cluster_ids(:, 1) uint32
        
        % nClusters x 1
        cluster_template_mostUsed(:, 1) uint32 % nClusters x 1 which template is most used by each cluster
        
        % nClusters x maxTemplatesPerCluster
        cluster_template_list cell
        
        % nClusters x nTemplates
        cluster_template_useCount uint64 % nClusters x nTemplates number of spikes in cluster i using template j
        cluster_num_templates uint32 % nClusters x 1 number of templates
        
        % nClusters x nChannelsSorted channel ids
        cluster_best_channels uint32
        
        % nClusters x nSpatialCoord
        cluster_centerOfMass single
        
        % nClusters x 1 logical
        cluster_is_localized(:, 1) logical
        
        % nClusters x nTimePoint x maxTemplatesPerCluster
        cluster_waveform single
        
        % nClusters x maxTemplatesPerCluster
        cluster_waveform_ch uint32
        cluster_amplitude single
        cluster_ttp single
    end
    
    properties(Dependent)
        template_depth % column 2 (y) of center of mass
        spike_depth % column 2 (y) of center of mass
        spike_is_localized % logical
        cluster_depth
    end
    
    methods
        function m = KilosortMetrics(ks, varargin)
            p = inputParser();
            p.addParameter('ampThreshCenterOfMass', 0.3, @isscalar);
            p.addParameter('extThreshLocalizedTemplate', 0.5, @isscalar);
            p.addParameter('distThreshLocalizedTemplate', 200, @isscalar);
            p.parse(varargin{:});
            
            assert(isa(ks, 'Neuropixel.KiloSortDataset'));
            assert(ks.isLoaded, 'KiloSortDataset is not loaded');
            
            nSteps = 8;
            prog = Neuropixel.Utils.ProgressBar(nSteps, 'Computing metrics for KiloSort Dataset...\n');
            m.ks = ks;
            
            m.fs = ks.fsAP;
            m.spike_templates = ks.spike_templates;
            m.spike_clusters = ks.spike_clusters;
            m.cluster_ids = ks.cluster_ids;
            m.channelMap = ks.channelMap;
            m.channel_ids = ks.channel_ids;
            
            m.concatenatedSamples = ks.concatenatedSamples;
            m.concatenatedStarts = ks.concatenatedStarts;
            m.concatenatedNames = ks.concatenatedNames;
           
            m.spike_times = ks.spike_times;
            
            % A. compute unwhitened templates
            % ks.templates is potentially only specified on a subset of channels which may differ across templates.
            % although for Kilosort this is every channel
            ks = m.ks; ks.checkLoaded(); %#ok<*PROP>
            templates = ks.templates;
            template_unw = zeros(size(templates, 1), size(templates, 2), ks.nChannelsSorted, 'like', templates);
            wmi = single(ks.whitening_mat_inv);
            assert(size(wmi, 1) == size(template_unw, 3), 'dim 3 of templates must match whitening matrix inverse');
            for iT = 1:size(templates, 1)
                whichChannels = ks.templates_ind(iT, :);
                template_unw(iT, :, whichChannels) = templates(iT, :, :);
            end
            sz = size(template_unw);
            m.template_unw = reshape(reshape(template_unw, sz(1)*sz(2), sz(3)) * wmi, sz);
            
            % for each channel i, list other channels in order of spatial proximity
            % nChannelsSorted x nChannelsSorted, include each channel in its own closest list
            closest_lookup = [ks.channel_ids, ks.channelMap.getClosestChannels(ks.nChannelsSorted-1, ks.channel_ids, ks.channel_ids)];
            m.template_best_channels = nan(ks.nClusters, ks.nChannelsSorted);
            for iT = 1:ks.nTemplates
                [~, bestChannelInd] = max(range(m.template_unw(iT, :, :), 2));
                m.template_best_channels(iT, :) = closest_lookup(bestChannelInd, :)';
            end
             
            prog.increment('Template center of mass');
            % B. determine template center of mass
            %   1. compute template amp on each channel, zeroing out small (< 0.3 max) faraway channels
            templateUnscaledAmps = squeeze(max(m.template_unw,[],2) - min(m.template_unw,[],2)); % nTemplates x nTemplateChannels
            [templateUnscaledMaxAmp, templateMaxChInd] = max(templateUnscaledAmps, [], 2);
            threshAmp = templateUnscaledMaxAmp * p.Results.ampThreshCenterOfMass;
            templateUnscaledAmps(templateUnscaledAmps < threshAmp) = 0;
            
            %   2. compute template channel positions (nTemplates x nCh x nSpatialDim)
            templateChannelPos = reshape(ks.channel_positions(ks.templates_ind(:), :), [size(ks.templates_ind), size(ks.channel_positions, 2)]);
            
            %   3. compute template center of mass
            m.template_centerOfMass = Neuropixel.Utils.TensorUtils.squeezeDims(sum(templateUnscaledAmps .* templateChannelPos, 2) ./ ...
                sum(templateUnscaledAmps, 2), 2);
            
            % C. determine which templates are localized
            m.template_is_localized = false(ks.nTemplates, 1);
            for iT = 1:ks.nTemplates
                % find channels where the template has weight at least 0.5 of max
                threshExt = max(abs(m.template_unw(iT, :))) * p.Results.extThreshLocalizedTemplate;
                maskCh = squeeze(max(abs(m.template_unw(iT, :, :)), [], 2)) >= threshExt; % nChannelsSorted x 1
                
                if nnz(maskCh) < 2
                    m.template_is_localized(iT) = true;
                else
                    distMat = pdist(ks.channel_positions(maskCh, :));
                    m.template_is_localized(iT) = max(distMat(:)) < p.Results.distThreshLocalizedTemplate;
                end
            end
            
            % D. scale the spikes by both templates's largest intrinsic amplitude and then ks.amplitudes
            prog.increment('Spike / template amplitude scaling');
            m.spike_amplitude = ks.amplitudes .* templateUnscaledMaxAmp(ks.spike_templates) * ks.apScaleToUv;
            
            % E. determine template scale by averaging ks.amplitudes that use that template
            template_scale = accumarray(ks.spike_templates, ks.amplitudes, [ks.nTemplates, 1], @mean);
            m.template_scaled = m.template_unw .* template_scale * ks.apScaleToUv;
            
            % F. determine template largest scaled waveform
            prog.increment('Template waveforms');
            m.template_waveform_ch = nan(ks.nTemplates, 1);
            m.template_waveform = nan(ks.nTemplates, size(ks.templates, 2));
            for iT = 1:ks.nTemplates
                m.template_waveform_ch(iT) = ks.templates_ind(iT, templateMaxChInd(iT));
                m.template_waveform(iT, :) = m.template_scaled(iT, :, m.template_waveform_ch(iT));
            end
            m.template_amplitude = templateUnscaledMaxAmp .* template_scale * ks.apScaleToUv;
            
            % G. determine trough to peak time of each waveform
            [~, template_trough] = min(m.template_waveform, [], 2);
            [~, template_peak] = max(m.template_waveform, [], 2);
            m.template_ttp = (template_peak - template_trough) / ks.fsAP * 1000;
            m.template_ttp(m.template_ttp <= 0) = NaN;
            
            prog.increment('Individual spike center of mass');
            % H. use private pcs to determine spike center of mass a la driftmap
            %   1. weight channels based on squared projection onto 1st pc
            pc1proj = squeeze(ks.pc_features(:, 1, :)); % nSpikes x nPCFeatures
            pc1proj(pc1proj < 0) = 0; % only positive entries contribute
            pc1weight = pc1proj.^2;
            
            %   2. compute center of mass. (spikeFeatureChannelPos is nSpikes x nCh x nSpatialDim)
            spike_pcfeat_chind = ks.pc_feature_ind(ks.spike_templates, :);
            spikeFeatureChannelPos = reshape(ks.channel_positions(spike_pcfeat_chind(:), :), [size(spike_pcfeat_chind), size(ks.channel_positions, 2)]);
            m.spike_centerOfMass = Neuropixel.Utils.TensorUtils.squeezeDims(sum(pc1weight .* spikeFeatureChannelPos, 2) ./ ...
                sum(pc1weight, 2), 2);
            
            % I. compute cluster weighting over templates and list of templates used by each cluster, sorted by number of uses
            prog.increment('Cluster weighting over templates');
            [mask_spike_in_cluster, cluster_ind_each_spike] = ismember(ks.spike_clusters, ks.cluster_ids);
            m.cluster_template_useCount = accumarray([cluster_ind_each_spike(mask_spike_in_cluster), ks.spike_templates(mask_spike_in_cluster)], ...
                ones(ks.nSpikes, 1), [ks.nClusters, ks.nTemplates]);
            
            [~, m.cluster_template_mostUsed] = max(m.cluster_template_useCount, [], 2);
            m.cluster_best_channels = m.template_best_channels(m.cluster_template_mostUsed, :);
            
            m.cluster_num_templates = sum(m.cluster_template_useCount > 0, 2);
            maxTemplatesPerCluster = max(m.cluster_num_templates);
            m.cluster_template_list = cell(ks.nClusters, 1);
            for iC = 1:ks.nClusters
                inds = find(m.cluster_template_useCount(iC, :));
                [~, sortIdx] = sort(m.cluster_template_useCount(iC, inds), 'descend');
                m.cluster_template_list{iC} = inds(sortIdx)';
            end
            
            % J. cluster is localized if all templates used are localized
            m.cluster_is_localized = arrayfun(@(iC) all(m.template_is_localized(m.cluster_template_list{iC})), 1:ks.nClusters);
            
            prog.increment('Cluster center of mass');
            % K. find cluster center of mass
            cluster_template_weight = double(m.cluster_template_useCount) ./ sum(m.cluster_template_useCount, 2);
            m.cluster_centerOfMass = Neuropixel.Utils.TensorUtils.linearCombinationAlongDimension(m.template_centerOfMass, 1, cluster_template_weight);
            
            % L. cluster waveforms
            prog.increment('Cluster waveform');
            cluster_waveform = nan(ks.nClusters, size(m.template_waveform, 2), maxTemplatesPerCluster);
            cluster_waveform_ch = nan(ks.nClusters);
            for iC = 1:ks.nClusters
                cluster_waveform(iC, :, 1:m.cluster_num_templates(iC)) = permute(m.template_waveform(m.cluster_template_list{iC}, :), [3 2 1]);
                cluster_waveform_ch(iC, 1:m.cluster_num_templates(iC)) = m.template_waveform_ch(m.cluster_template_list{iC})';
            end
            m.cluster_waveform = cluster_waveform;
            m.cluster_waveform_ch = cluster_waveform_ch;
            % J. cluster_amplitudes and ttp - weighted mean of template amplitude / ttp
            prog.increment('Cluster amplitudes');
            m.cluster_amplitude = accumarray(cluster_ind_each_spike(mask_spike_in_cluster), m.spike_amplitude(mask_spike_in_cluster), [ks.nClusters, 1], @mean);
            template_mask = ~isnan(m.template_ttp);
            m.cluster_ttp = (double(m.cluster_template_useCount(:, template_mask)) * m.template_ttp(template_mask)) ./ sum(m.cluster_template_useCount(:, template_mask), 2);
            
            prog.finish();
        end
        
        function n = get.nSpikes(m)
            n = size(m.spike_times, 1);
        end
        
        function n = get.nChannelsSorted(m)
            n = size(m.template_unw, 3);
        end
        
        function n = get.nTemplates(m)
            n = size(m.template_unw, 1);
        end
        
        function n = get.nTemplateTimepoints(m)
            n = size(m.template_unw, 2);
        end 
        
        function n = get.nClusters(m)
            n = size(m.cluster_ids, 1);
        end
        
        function n = get.maxTemplatesPerCluster(m)
            n = size(m.cluster_waveform, 3);
        end
        
        function n = get.nConcatenatedFiles(m)
            n = numel(m.concatenatedSamples);
        end
        
        function d = get.spike_depth(m)
            d = m.spike_centerOfMass(:, 2);
        end
        
        function d = get.cluster_depth(m)
            d = m.cluster_centerOfMass(:, 2);
        end
        
        function tf = get.spike_is_localized(m)
            tf = m.template_is_localized(m.spike_templates);
        end
        
        function assertHasKs(m)
            assert(~isempty(m.ks), 'Must set .ks to KiloSortDataset');
        end
    end
    
    methods % Cluster existence over time
        function [fracBinsValid, fracBinsValidExcludingEdges] = computeClusterValidTimeSpan(m, varargin)
            % fracBinsValid is the total fraction of bins that had at least threshSpikes
            % Then we divide up the day into any zeros at the beginning and end of the file, and the middle region between these bands of zeros
            
            p = inputParser();
            p.addParameter('binsClose', 2, @isscalar);
            p.addParameter('threshSpikes', 5, @isscalar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            [counts, ~] = m.computeClusterBinnedCounts(p.Unmatched);
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
        
        function [counts, tvec] = computeClusterBinnedCounts(m, varargin)
            p = inputParser();
            p.addParameter('cluster_ids', [], @isvector);
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'Neuropixel.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('maskRegionsOutsideTrials', true, @islogical);
            p.addParameter('binWidth', 10, @isscalar); 
            p.parse(varargin{:});
            
            tsi = p.Results.tsi;
            if ~isempty(tsi) && p.Results.maskRegionsOutsideTrials
                mask = m.computeSpikeMaskWithinTrials(tsi);
            else
                mask = true(numel(m.spike_times), 1);
            end
            if ~isempty(p.Results.cluster_ids)
                cluster_ids = p.Results.cluster_ids;
                mask = mask & ismember(m.spike_clusters, p.Results.cluster_ids);
            else
                cluster_ids = m.cluster_ids;
            end
            
            spikeTimes = double(m.spike_times(mask)) / m.ks.fsAP; % convert to seconds
            spikeClu = m.spike_clusters(mask); %#ok<*PROPLC>
            
            nClu = numel(cluster_ids);
            
            edges = (0:p.Results.binWidth:max(spikeTimes))';
            tvec = edges(1:end-1);
            counts = zeros(nClu, numel(tvec));
            for iC = 1:nClu
                thisClu = spikeClu == cluster_ids(iC);
                
                times = spikeTimes(thisClu);
                
                nSp = histc(times, edges);
                nSp = nSp(1:end-1);
                
                counts(iC, :) = nSp;
            end
            
            % keep only bins within trial regions
            if ~isempty(tsi)
                mask = false(numel(tvec), 1);
                tvec_samples = tvec * m.fs;
                bin_width_samples = p.Results.binWidth * m.fs;
                [idxStart, idxStop, ~] = tsi.computeActiveRegions();
                for iR = 1:numel(idxStart)
                    mask(tvec_samples >= idxStart(iR) & tvec_samples+bin_width_samples <= idxStop(iR)) = true;
                end
                counts = counts(:, mask);
                tvec = tvec(mask);
            end
        end
        
        function meanFR = computeClusterMeanFR(m, varargin)
            p = inputParser();
            p.addParameter('binWidth', 10, @isscalar); 
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            counts = m.computeClusterBinnedCounts('binWidth', p.Results.binWidth, varargin{:});
            meanFR = mean(counts, 2) / p.Results.binWidth;
        end
    end
    
    methods % Plotting stability over time
        function mask = computeSpikeMaskWithinTrials(m, tsi)
            mask = false(numel(m.spike_times), 1);
            [idxStart, idxStop, ~] = tsi.computeActiveRegions();

            for iR = 1:numel(idxStart)
                mask(m.spike_times >= idxStart(iR) & m.spike_times <= idxStop(iR)) = true;
            end
        end
        
        function plotClusterDriftmap(m, varargin)
            p = inputParser();
            p.addParameter('cluster_ids', [], @isvector);
            p.addParameter('spike_mask', [], @(x) isempty(x) || isvector(x));
            p.addParameter('localizedOnly', true, @islogical);
            p.addParameter('minAmpQuantile', 0, @isscalar);
            p.addParameter('showIndividual', true, @islogical);
            p.addParameter('showSmooth', true, @islogical); % smooth amplitudes over spikes by this width
            p.addParameter('smoothWidthSeconds', 50, @isscalar); % smoothing kernel width in seconds
            p.addParameter('scaleWithAmp', false, @islogical);
            p.addParameter('colorByAmp', false, @islogical);
            p.addParameter('alpha', 1, @isscalar);
            p.addParameter('zShuffleClusters', true, @islogical);
            p.addParameter('maxClustersPlot', Inf, @isscalar);
            
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'Neuropixel.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('maskRegionsOutsideTrials', true, @islogical);
            p.addParameter('exciseRegionsOutsideTrials', false, @islogical);
            p.parse(varargin{:});
            
            mask = p.Results.spike_mask;
            tsi = p.Results.tsi;
            if ~isempty(tsi) && (p.Results.maskRegionsOutsideTrials || p.Results.exciseRegionsOutsideTrials)
                mask = m.computeSpikeMaskWithinTrials(tsi);
            else
                mask = true(numel(m.spike_times), 1);
            end
            if p.Results.localizedOnly
                mask = mask & m.spike_is_localized;
            end
            if ~isempty(p.Results.cluster_ids)
                mask = mask & ismember(m.spike_clusters, p.Results.cluster_ids);
            end
            if p.Results.minAmpQuantile > 0
                thresh = quantile(m.spike_amplitude, p.Results.minAmpQuantile);
                mask = mask & m.spike_amplitude >= thresh;
            end
            
            showSmooth = p.Results.showSmooth;
            showIndividual = p.Results.showIndividual;
            alpha = p.Results.alpha;
            if alpha == 1 && showSmooth
                alpha = 0.3;
            end
            colorByAmplitude = p.Results.colorByAmp;
            
            spikeTimes = m.spike_times(mask);
            if p.Results.exciseRegionsOutsideTrials
                timeShifts = tsi.computeShiftsExciseRegionsOutsideTrials();
                spikeTimes = timeShifts.shiftTimes(spikeTimes);
            end
            spikeTimes = double(spikeTimes) / m.ks.fsAP; % convert to seconds
            spikeAmps = m.spike_amplitude(mask);
            spikeYpos = m.spike_depth(mask);
            spikeClu = m.spike_clusters(mask); %#ok<*PROPLC>
            ampMax = quantile(spikeAmps, 0.95);
            
            uClu = unique(spikeClu);
            clusterInds = m.lookup_clusterIds(uClu);
            nClu = numel(uClu);
            clusterAmps = m.cluster_amplitude(clusterInds);
            [~, clusterAmpsSortOrder] = sort(clusterAmps, 'ascend');
            
            bgcolor = [0.92 0.92 0.95];
            
            if colorByAmplitude
                cluAmpMax = quantile(clusterAmps, 0.95);
                cluAmpNormalized = min(1, clusterAmps / cluAmpMax);  
%                 
%                 cmap = amplitudeCmap(clusterAmpsSortOrder / numel(clusterAmps));
                cmap = amplitudeCmap(cluAmpNormalized);
            else
                cmap = Neuropixel.Utils.distinguishable_colors(nClu, [1 1 1; 0 1 0]);
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
            spikeTimeBins = (0:smoothBinWidth:max(spikeTimes))';
            for iC_sorted = 1:nClu
                iC = clusterAmpsSortOrder(iC_sorted);
                thisClu = spikeClu == uClu(iC);
                
                x = spikeTimes(thisClu);
                y = spikeYpos(thisClu);
                
                if p.Results.scaleWithAmp
                    sz = spikeAmps(thisClu) / ampMax * 20;
                else
                    sz = 6;
                end
                
                ud = struct('cluster_id', uClu(iC), 'cluster_amplitude', sprintf('%.1f uV', clusterAmps(iC)), ...
                    'cluster_is_localized', m.cluster_is_localized(clusterInds(iC)), ...
                    'xname', 'Time', 'yname', 'Depth', 'xunits', 'sec', 'yunits', 'um');
                
                if showIndividual
                    h = plot(x, y, '.', 'Color', cmap(iC,:), 'MarkerSize', sz, 'UserData', ud);
                   TrialDataUtilities.Plotting.setMarkerOpacity(h, alpha);
 
                   hold on;
                end
                if showSmooth
%                     if showIndividual, width = 2; else, width = 1; end
                    width = 2;
                    [x, y] = binsmooth(x, y, spikeTimeBins, 5, smoothBinWidth); % min 10 spikes in 10 second increments
                    plot(x, y, '-', 'Color', cmap(iC,:), 'LineWidth', width, 'UserData', ud);
                    hold on;
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
            set(gcf, 'InvertHardcopy', 'off');
            axis tight;
            
            Neuropixel.Utils.configureDataTipsFromUserData(gcf);
                
            function [xb, yb] = binsmooth(x, y, edges, smoothBy, minSpikesPerBin)
                [nSp, whichBin] = histc(x, edges);
                xb = 0.5 * (edges(1:end-1) + edges(2:end));
                nSp = nSp(1:end-1);

                mask = whichBin > 0;
                yb = accumarray(whichBin(mask), y(mask), [numel(xb), 1], @mean, single(NaN));
                
                yb(nSp < minSpikesPerBin) = NaN;
                
                if smoothBy > 0
                    yb = smooth(yb, smoothBy, 'lowess', 2);
                end

                mask = nSp >= minSpikesPerBin;
                xb = xb(mask);
                yb = yb(mask);
            end  
            
            function cmap = amplitudeCmap(amp)
                N = numel(amp);
                cmap_hsl = cat(2, rand(N, 1), repmat(0.7, N, 1), repmat(0.5, N, 1));
                lerp = @(x, a, b) x*(b-a) + a; % map x from [0 1] to [a b]
                lum = lerp(amp, 0, 0.9);
%                 sat = lerp(amp, 0, 0.9);
%                 cmap_hsl(:, 2) = sat;
                cmap_hsl(:, 3) = lum;
                cmap = Neuropixel.Utils.hsl2rgb(cmap_hsl);
            end
        end
            
        % Inputs: spikeTimes, spikeAmps, spikeYpos - names self explanatory
        %         opt - optional, empty by default; 'mark' - will mark detected drifts, 'show' - will generate a different plot,
        %               where only large spikes are used, and the detection of drift locations is demonstrated
        function info = plotDriftmap(m, varargin)
            % Note: this is pretty much copied wholesale from cortex-lab/spikes, with some nice annotations added on
            
            p = inputParser();
            p.addParameter('mode', 'mark', @ischar); % or 'show' or 'mark'
            p.addParameter('spikeAmpQuantile', 0.966, @isscalar); % consider only spikes larger than quantile fo amplitude
            p.addParameter('segmentDepth', 800, @isscalar); % um to segment probe into
            p.addParameter('nAmpBins', 20, @isscalar);
            p.addParameter('localizedOnly', true, @islogical);
            p.addParameter('driftThreshold', 6, @isscalar); % um
            p.addParameter('driftTimeWindow', 10, @isscalar); % in seconds
            p.addParameter('minSpikesDrift', 50, @isscalar);
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'Neuropixel.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('maskRegionsOutsideTrials', false, @islogical);
            p.addParameter('exciseRegionsOutsideTrials', false, @islogical);
            p.addParameter('xOffset', 0, @isscalar);
            p.addParameter('ampRange', [], @(x) isvector(x) || isempty(x));
            
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
            
            % mask
            spike_times_samples = m.spike_times;
            if ~isempty(tsi) && (p.Results.maskRegionsOutsideTrials || p.Results.exciseRegionsOutsideTrials)
                mask = m.computeSpikeMaskWithinTrials(tsi);
            else
                mask = true(numel(spike_times_samples), 1);
            end
            if p.Results.localizedOnly
                mask = mask & m.spike_is_localized;
            end
            
            spikeTimes = m.spike_times(mask);
            if p.Results.exciseRegionsOutsideTrials
                timeShifts = tsi.computeShiftsExciseRegionsOutsideTrials();
                spikeTimes = timeShifts.shiftTimes(spikeTimes);
            else
                timeShifts = [];
            end
            spikeTimes = double(spikeTimes) / m.ks.fsAP; % convert to seconds
            spikeAmps = m.spike_amplitude(mask);
            spikeYpos = m.spike_depth(mask);
            
            if ~strcmpi(opt, 'show')
                m.plotSpikesByAmplitude('spike_mask', mask, 'time_shifts', timeShifts, 'nAmpBins', p.Results.nAmpBins, ...
                    'ampRange', p.Results.ampRange, ...
                    'localizedOnly', p.Results.localizedOnly, 'tsi', p.Results.tsi, 'xOffset', xOffset);
                
                ylim(m.channelMap.ylim);
                set(gca, 'XLimSpec', 'tight');
                box off;
                
                nD = floor(max(spikeYpos) / segDepth);
                driftEventsAll = cell(nD, 1);
                for iD = 1:nD % break the recording into 800 um segments
                    d = segDepth*(iD-1);
                    tmp = spikeAmps(spikeYpos >= d & spikeYpos < d+segDepth);
                    I = spikeAmps > mean(tmp) + spikeAmpStdThresh*std(tmp) & spikeYpos >= d & spikeYpos < d+segDepth; % large spikes in current segment
                    driftEvents = detectDriftEvents(spikeTimes(I), spikeYpos(I), strcmp(opt, 'show'));
                    driftEventsAll{iD} = driftEvents;
                    if strcmpi(opt, 'mark') && ~isempty(driftEvents)
                        plot(driftEvents(:,1) + xOffset, driftEvents(:,2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r')
                        % text(driftEvents(:,1)+1, driftEvents(:,2), num2str(round(driftEvents(:,3))), 'Color', 'r') % the magnitude of the drift
                    end
                end  
                
                driftEvents = cat(1, driftEventsAll{:});
                if ~isempty(driftEvents)
                    driftTimes = driftEvents(:, 1);
                    hold on;
                    Neuropixel.Utils.rugplot(driftTimes + xOffset, 'side', 'top', 'Color', [1 0.2 0.2]);
                else
                    driftTimes = [];
                end
            end
            
            if ~isempty(tsi)
                hold on;
                tsi.markTrialTicks('time_shifts', timeShifts, 'xOffset', xOffset, 'side', 'bottom', 'Color', [0.2 0.2 1]);
            end
            
            if p.Results.exciseRegionsOutsideTrials
                m.markExcisionBoundaries(timeShifts, 'xOffset', xOffset);
            end
            m.markConcatenatedFileBoundaries('time_shifts', timeShifts, 'xOffset', xOffset);
            
            hold off;
            
            info.driftTimes = driftTimes;
            info.xMax = max(spikeTimes) + xOffset;
            
            
            % driftEvents will contain a column of times, a column of depths, and a column of drift magnitudes
            function driftEvents = detectDriftEvents(spikeTimes, spikeDepths, doPlot)
                if nargin < 3
                    doPlot = false;
                end
                driftEvents = [];
                
                D = 2; % um
                bins = min(spikeDepths)-D:D:max(spikeDepths)+D;
                h = histc(spikeDepths, bins);
                
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
        
        function plotSpikesByAmplitude(m, varargin)
            % this plots each clusters' spikes in their own color, either as markers or as a line with smoothing
            p = inputParser();
            p.addParameter('spike_mask', [], @(x) isempty(x) || isvector(x));
            p.addParameter('localizedOnly', true, @islogical);
            p.addParameter('nAmpBins', 20, @isscalar);
            p.addParameter('ampRange', [], @(x) isvector(x) || isempty(x));
            
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec'));
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'Neuropixel.TrialSegmentationInfo')); % to mark trial boundaries
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            
            mask = p.Results.spike_mask;
            if isempty(mask)
                mask = true(numel(m.spike_times), 1);
            end
            if p.Results.localizedOnly
                mask = mask & m.spike_is_localized;
            end

            timeShifts =p.Results.time_shifts;
            spikeTimes = m.spike_times(mask);
            spikeTimesOrig = spikeTimes;
           
            if ~isempty(timeShifts)
                spikeTimes = timeShifts.shiftTimes(spikeTimes);
            end
            spikeTimes = double(spikeTimes) / m.ks.fsAP; % convert to seconds
            spikeClusters = m.spike_clusters(mask);
            spikeAmps = m.spike_amplitude(mask);
            spikeYpos = m.spike_depth(mask);
            
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

                h = plot(spikeTimes(theseSpikes) + xOffset, spikeYpos(theseSpikes), '.', 'Color', colors(b,:));
                
                if ~verLessThan('matlab', '9.6.0') % R2019a
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
                
                if numel(m.concatenatedStarts) > 1
                    % show original file name and ind
                    [fileInd, origInd] = m.lookup_sampleIndexInConcatenatedFile(spikeTimesOrig(theseSpikes));
                    
                    % grab associated file name
%                     fileNames = strings(numel(spikeTimesOrig), 1);
%                     mask = ~isnan(fileInd);
%                     fileNames(mask) = m.concatenatedNames(fileInd(mask));
                    
                    row = dataTipTextRow('Orig File Ind', fileInd, '%d');
                    h.DataTipTemplate.DataTipRows(end+1) = row;
                    row = dataTipTextRow('Orig Sample Ind', double(origInd), '%d');
                    h.DataTipTemplate.DataTipRows(end+1) = row;
                end
                
                hold on;
            end
            xlabel('time (sec)')
            ylabel('y position (um)')
            
            h = zoom;
            h.Motion = 'horizontal';
        end
        
        function h = markExcisionBoundaries(m, shifts, varargin)
            p = inputParser();
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            starts = double(shifts.idxShiftStart(2:end)) / m.fs;
            for iF = 1:numel(starts)
                h = xline(starts(iF) + xOffset, '-', sprintf('Exc %d', iF), 'Color', [0.9 0.3 0.9], 'LineWidth', 1, 'Interpreter', 'none');
%                 h.NodeChildren(1).NodeChildren(1).ColorData = uint8(255*[0.3 0.3 0.3 1]');
                h.NodeChildren(1).NodeChildren(1).BackgroundColor = uint8(255*[1 1 1 0.5]');
                hold on;
                
            end
        end
        
        function h = markConcatenatedFileBoundaries(m, varargin)
            p = inputParser();
            p.addParameter('time_shifts', [], @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec'));
            p.addParameter('xOffset', 0, @isscalar);
            p.parse(varargin{:});
            xOffset = p.Results.xOffset;
            
            
            timeShifts = p.Results.time_shifts;
            catStarts = m.concatenatedStarts;
            if ~isempty(timeShifts)
                catStarts = timeShifts.shiftTimes(catStarts);
            end
            starts = double(catStarts) / m.fs;
            for iF = 1:numel(starts)
                h = xline(starts(iF) + xOffset, '-', m.concatenatedNames{iF}, 'Color', [0.3 0.3 0.9], 'LineWidth', 1, 'Interpreter', 'none');
%                 h.NodeChildren(1).NodeChildren(1).ColorData = uint8(255*[0.3 0.3 0.3 1]');
                h.NodeChildren(1).NodeChildren(1).BackgroundColor = uint8(255*[1 1 1 0.5]');
                hold on;
            end
        end
        
        function [fileInds, origSampleInds] = lookup_sampleIndexInConcatenatedFile(m, inds)
           [fileInds, origSampleInds] = Neuropixel.Utils.lookup_sampleIndexInConcatenatedFile(m.concatenatedStarts, inds);
        end
    end
    
    methods % Plotting cluster waveforms
        
        function clusterInds = lookup_clusterIds(m, cluster_ids)
            [tf, clusterInds] = ismember(cluster_ids, m.cluster_ids);
            assert(all(tf), 'Some cluster ids were not found in ds.clusterids');
        end
        
        function clusterInds = lookup_channelIds(m, channel_ids)
            [tf, clusterInds] = ismember(channel_ids, m.channel_ids);
            assert(all(tf(:)), 'Some cluster ids were not found in ds.clusterids');
        end
        
        function plotClusterImage(m, cluster_ids, varargin)
            clusterInds = m.lookup_clusterIds(cluster_ids);
            templateLists = m.cluster_template_list(clusterInds);
            templateInds = cat(1, templateLists{:});
            cmap = Neuropixel.Utils.distinguishable_colors(numel(templateInds));
            
            m.plotTemplateImage(templateInds, 'colormap', cmap, varargin{:});
        end
        
        function plotTemplateImage(m, templateInds, varargin)
%             isholding = ishold;
            channel_ids_by_template = m.plotTemplateImageInternal(templateInds, varargin{:}); %#ok<NASGU>
%             channel_idx_all = unique(channel_ids_by_template(:)); 
%             hold on;
%             m.plotRecordingSites('channel_ids', channel_idx_all, 'showChannelLabels', false);
%             if ~isholding, hold off , end
        end
            
        function channel_ids_by_template = plotTemplateImageInternal(m, templateInds, varargin)
            p = inputParser();
            p.addParameter('xmag', 1.5, @isscalar);
            p.addParameter('ymag', 1.5, @isscalar);
            
            p.addParameter('cluster_ids', m.cluster_ids, @isvector);
            p.addParameter('colormap', @Neuropixel.Utils.distinguishable_colors, @(x) isa(x, 'function_handle') || ismatrix(x) || ischar(x));
            p.addParameter('templateLabels', {}, @iscell);
            
            % and ONE OR NONE of these to pick channels (or channels for each cluster)
            p.addParameter('channel_ids_by_template', [], @(x) isempty(x) || ismatrix(x));
            p.addParameter('best_n_channels', NaN, @isscalar); % or take the best n channels based on this clusters template when cluster_id is scalar
            
            p.parse(varargin{:});
            
            isholding = ishold;
            
            % figure out actual times requested
            if isfinite(p.Results.best_n_channels)
                channel_ids_by_template = m.template_best_channels(templateInds, 1:p.Results.best_n_channels);
            elseif ~isempty(p.Results.channel_ids_by_template)
                channel_ids_by_template = p.Results.channel_idx_by_cluster;
            else
                channel_ids_by_template = m.template_best_channels(templateInds, :);
            end
            channel_ind_by_template = m.lookup_channelIds(channel_ids_by_template);
            
            yspacing = m.channelMap.yspacing;
            xspacing = m.channelMap.xspacing;
            xmag = p.Results.xmag;
            ymag = p.Results.ymag;
            
            % plot relative time vector
            tvec = linspace(0, xspacing * xmag, size(m.template_scaled, 2));
            tvec_shift = tvec - mean(tvec);
            
            nTemp = numel(templateInds);
            nChannelsSorted = size(channel_ind_by_template, 2);
            data = nan(nTemp, size(m.template_scaled, 2), nChannelsSorted, 'single');
            for iT = 1:nTemp
                data(iT, :, :) = m.template_scaled(templateInds(iT), :, channel_ind_by_template(iT, :));
            end
            data = data - mean(data(:, :), 2); % center over time
            data = data ./ (max(data(:)) - min(data(:))) * yspacing * ymag; % normalize amplitudes
            
            cmap = p.Results.colormap;
            if isa(cmap, 'function_handle')
                cmap = cmap(nTemp);
            end
            templateLabels = p.Results.templateLabels;
            if isempty(templateLabels)
                templateLabels = arrayfun(@(ind) sprintf("template %d", ind), templateInds);
            end
            
            for iT = 1:nTemp
                this_channel_ind = channel_ind_by_template(iT, :);
%                 this_channel_ids = channel_ids_by_template(iT, :);
                xc = m.channelMap.xcoords(this_channel_ind);
                yc = m.channelMap.ycoords(this_channel_ind);
                cmapIdx = mod(iT-1, size(cmap, 1))+1;
                
                for iC = 1:nChannelsSorted
                    wave = Neuropixel.Utils.TensorUtils.squeezeDims(data(iT, :, iC), 1) + yc(iC);
                        
%                     ud = struct('template_ind', template_inds(iT), 'template_amplitude', sprintf('%.1f uV', m.template_amplitude(template_inds(iT)),
%                     'channel_id', this_channel_ids(iC), 'template_is_localized', m.template_is_localized(template_inds(iT)), ...
%                     'xname', 'Time', 'yname', 'Voltage', 'xoffset', xc(iC) - mean(tvec), 'yoffset', yc(iC) + dataOffset(iT, 1, iC), 'xscale', 1, 'yscale', waveScalingFactor_umtouV, 'xunits', 'ms', 'yunits', 'uV');
%                     
                    h = plot(tvec_shift + xc(iC), wave, 'Color', cmap(cmapIdx, :), 'LineWidth', 0.5);
                    hold on;
                    
                    if iC == 1 
                        Neuropixel.Utils.showFirstInLegend(h, templateLabels{iT});
                    else
                        Neuropixel.Utils.hideInLegend(h);
                    end
                end
                
                %drawnow;
            end
            
            axis off;
            axis tight;
            box off;
            if ~isholding, hold off, end
        end
        
        function plotRecordingSites(m, varargin)
            p = inputParser();
            p.addParameter('channel_ids', m.channel_ids, @isvector)
            p.addParameter('showChannelLabels', false, @islogical);
            p.addParameter('labelArgs', {}, @iscell);
            p.parse(varargin{:});
            
            channelInds = m.lookup_channelIds(p.Results.channel_ids);
            xc = m.channelMap.xcoords(channelInds);
            yc = m.channelMap.ycoords(channelInds);
            plot(xc, yc, '.', 'Marker', 's', 'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 5);
            if p.Results.showChannelLabels
                for iC = 1:numel(channelInds)
                    text(xc(iC), yc(iC), sprintf('ch %d', m.channel_ids(channelInds(iC))), ...
                            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'Background', 'none', ...
                            p.Results.labelArgs{:});
                end
            end
        end

        function plotClusterWaveformAtCenterOfMass(m, varargin)
            p = inputParser();
            p.addParameter('waveformScale', 10, @isscalar);
            p.addParameter('waveformWidth', m.channelMap.xspacing/20, @isscalar);
            p.addParameter('waveformHeight', m.channelMap.yspacing*3, @isscalar);
            p.addParameter('cluster_ids', m.cluster_ids, @isvector);
            p.addParameter('colormap', @Neuropixel.Utils.distinguishable_colors, @(x) isa(x, 'function_handle') || ismatrix(x) || ischar(x));
            
            p.addParameter('useAutoAxis', false, @islogical);
            p.parse(varargin{:});
            
            cluster_ids = p.Results.cluster_ids;
            clusterInds = m.lookup_clusterIds(cluster_ids);
            waves = m.cluster_waveform(clusterInds, :, 1); % nClusters x nTimepoints
            
            xvec = linspace(-p.Results.waveformWidth/2, p.Results.waveformWidth/2, size(waves, 2)) * p.Results.waveformScale;
            timeScaleFactor_umtoms = (size(waves, 2) / m.fs * 1000) / range(xvec);
            
            waveScalingFactor_umtouV = double((max(waves(:)) - min(waves(:))) / p.Results.waveformHeight / p.Results.waveformScale);
            waves = waves ./ waveScalingFactor_umtouV;
            
            colormap = p.Results.colormap;
            if isa(colormap, 'function_handle')
                colormap = colormap(size(waves, 1));
            end
            
            % nClusters x 2 or 3
            com = m.cluster_centerOfMass(clusterInds, :);
            
            holding = ishold;
            
            % plot recording sites
            m.plotRecordingSites();
            hold on
            
            xlim(double([min(com(:,1)) - m.channelMap.xspacing/2, max(com(:,1)) + m.channelMap.xspacing/2]));
            ylim(double([min(com(:,2)) - m.channelMap.yspacing*5, max(com(:,2)) + m.channelMap.yspacing*5]));
            
            for iC = 1:size(waves, 1)
                if ischar(colormap)
                    color = colormap;
                else
                    color = colormap(mod(iC-1, size(colormap, 1))+1, :);
                end
                
                ud = struct('cluster_id', cluster_ids(iC), 'cluster_amplitude', sprintf('%.1f uV', m.cluster_amplitude(clusterInds(iC))), ...
                    'cluster_is_localized', m.cluster_is_localized(clusterInds(iC)), ...
                    'xname', 'Time', 'yname', 'Voltage', 'xoffset', xvec(1) + com(iC, 1), 'yoffset', com(iC, 2), 'xscale', timeScaleFactor_umtoms, 'yscale', waveScalingFactor_umtouV, 'xunits', 'ms', 'yunits', 'uV');
                plot(xvec + com(iC, 1), waves(iC, :) + com(iC, 2), '-', 'Color', color, 'UserData', ud);
                hold on;    
            end
            
            if ~holding, hold off; end
            
            if p.Results.useAutoAxis
                AutoAxis.replaceScaleBars('xUnits', 'ms', 'xLength', 1, 'xScaleFactor', timeScaleFactor_umtoms, ...
                    'yUnits', 'uV', 'yLength', 200, 'yScaleFactor', waveScalingFactor_umtouV); 
            else
                xlabel('x (um)');
                ylabel('y (um)');
            end
%             
%             bgcolor = [0.9 0.9 0.9];
%             axh = gca;
%             axh.Color = bgcolor;
            box off;
            
            Neuropixel.Utils.configureDataTipsFromUserData(gcf);
        end
    end
    
    methods(Static) % multi plotting to simulate concatenated files
        function multiple_plotDriftmap(mSet, varargin)
            p = inputParser();
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'Neuropixel.TrialSegmentationInfo')); % to mark trial boundaries
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