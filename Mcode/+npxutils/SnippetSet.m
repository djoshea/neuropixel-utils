classdef SnippetSet < handle & matlab.mixin.Copyable
    % A collection of snippets of data collected from an IMEC data file
    % which includes the subset of channels which were sampled
    
    % data here should be thought of as a set of data for many snippets
    %
    
    properties (Transient)
        ks  % optional, handle to KilosortDataset for plotting waveforms and templates
    end
    
    properties
        data % (:, :, :, :) int16 % channels x time x snippets x layers - we deactivate the size checking here
        
        % optional, set if each snippet corresponds to a specific cluster
        cluster_ids (:,1) uint32 % nSnippets x 1, array indicating which cluster is extracted in each snippet, if this makes sense. otherwise will just be 1s
        
        group_ids (:,1)  % nSnippets x 1 metadata indicating which group a cluster belonged to
        
        channel_ids_by_snippet (:,:,:) uint32 % channels x snippets
        
        overlay_cluster_ids (:,1) uint32 % set of clusters whose waveforms or templates will be drawn on top of a given snippet, will be cluster_ids(i), or ks.cluster_ids
        
        % for manually specifying overlay matrices (can also be generated automatically with e.g. overlay_waveforms)
        overlay_labels (:,:,:,:) uint32 % channels x time x snippets x layers (label matrix corresponding to elements of data for the purpose of generating colored overlays)
        overlay_datatip_label (1,1) string
        overlay_datatip_values (:,1) % nLabels x 1
        
        sample_idx (:,1) uint64
        trial_idx (:,1) uint32
        window (:,2) int64 % in samples
        
        valid (:,1) logical
        
        channelMap % ChannelMap for coordinates
        scaleToUv (:,1) single % one value for each snippet (will auto inflate to nSnippets if set to be scalar
        fs
    end
    
    properties(Dependent)
        nChannels
        nTimepoints
        nSnippets
        nClusters
        nLayers
        
        data_valid
        time_samples (:,1) int64
        time_ms (:,1) double
        
        unique_cluster_ids
    end
    
    methods
        function ss = SnippetSet(ds, type)
            if nargin > 0
                if isa(ds, 'npxutils.KilosortTrialSegmentedDataset')
                    raw_dataset = ds.raw_dataset;
                    ks = ds.dataset;
                elseif isa(ds, 'npxutils.KilosortDataset')
                    raw_dataset = ds.raw_dataset;
                    ks = ds;
                elseif isa(ds, 'npxutils.ImecDataset')
                    raw_dataset = ds;
                    ks = [];
                else
                    error('Unknown arg type');
                end
                
                if nargin < 2
                    type = 'ap';
                end
                switch type
                    case 'ap'
                        ss.scaleToUv = raw_dataset.apScaleToUv;
                        ss.fs = raw_dataset.fsAP;
                    case 'lf'
                        ss.scaleToUv = raw_dataset.lfScaleToUv;
                        ss.fs = raw_dataset.fsLF;
                    otherwise
                        error('Unknown type %s', type);
                end
                
                ss.channelMap = raw_dataset.channelMap;
                ss.ks = ks;
            end
        end
        
        function n = get.nChannels(this)
            n = size(this.data, 1);
        end
        
        function n = get.nTimepoints(this)
            n = size(this.data, 2);
        end
        
        function n = get.nSnippets(this)
            n = size(this.data, 3);
        end
        
        function n = get.nLayers(this)
            n = size(this.data, 4);
        end
        
        function n = get.nClusters(this)
            n = numel(this.unique_cluster_ids);
        end
        
        function v = get.valid(this)
            if isempty(this.valid)
                v = true(size(this.data, 3), 1);
            else
                v = this.valid;
            end
        end
        
        function v = get.data_valid(this)
            v = this.data(:, :, this.valid);
        end
        
        function v = get.time_ms(this)
            v = double(this.window(1) : this.window(2)) / this.fs * 1000;
        end
        
        function v = get.time_samples(this)
            v = this.window(1) : this.window(2);
        end
        
        function v = get.unique_cluster_ids(this)
            v = unique(this.cluster_ids);
        end
        
        function v = get.overlay_cluster_ids(this)
            if isempty(this.overlay_cluster_ids)
                if isempty(this.unique_cluster_ids)
                    if isempty(this.ks)
                        v = [];
                    else
                        v = this.ks.cluster_ids;
                    end
                else
                    v = this.unique_cluster_ids;
                end
            else
                v = this.overlay_cluster_ids;
            end
        end
        
        function v = get.scaleToUv(this)
            v = this.scaleToUv;
            if isempty(v)
                v = NaN;
            end
            if numel(v) == 1
                v = repmat(v, this.nSnippets, 1);
            end
        end
        
        function colormap = generateOverlayColormap(this, overlay_labels, backgroundColor) %#ok<INUSL>
            if nargin < 3
                backgroundColor = [0 0 0; 1 1 1];
            end
            nLabels = max(overlay_labels(:));
            colormap = zeros(nLabels, 3);
            unique_labels = setdiff(unique(overlay_labels(:)), 0);
            label_found = ismember(1:nLabels, unique_labels);
            colormap(label_found, :) = npxutils.internal.graphics.distinguishable_colors(nnz(label_found), backgroundColor);
        end
        
        function this = selectClusters(this, cluster_ids)
            mask = ismember(this.cluster_ids, cluster_ids);
            this = this.selectData('maskSnippets', mask);
        end
        
        function this = selectData(this, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', this.valid, @isvector);
            p.addParameter('maskTime', true(this.nTimepoints, 1), @isvector);
            p.addParameter('maskChannels', true(this.nChannels, 1), @isvector);
            p.parse(varargin{:});
            
            maskSnippets = p.Results.maskSnippets;
            
            this = copy(this);
            this.valid = this.valid(maskSnippets); % first time auto generates based on data, so has to come first
            this.data = this.data(p.Results.maskChannels, p.Results.maskTime, maskSnippets);
            
            if ~isempty(this.overlay_labels)
                this.overlay_labels = this.overlay_labels(p.Results.maskChannels, p.Results.maskTime, maskSnippets);
            end
            
            this.cluster_ids = this.cluster_ids(maskSnippets);
            
            this.sample_idx = this.sample_idx(maskSnippets);
            if ~isempty(this.trial_idx)
                this.trial_idx = this.trial_idx(maskSnippets);
            end
            
            this.scaleToUv = this.scaleToUv(maskSnippets);
        end
        
        function [cluster_inds, cluster_ids] = lookup_clusterIds(this, cluster_ids)
            if islogical(cluster_ids)
                cluster_ids = this.unique_cluster_ids(cluster_ids);
            end
            [tf, cluster_inds] = ismember(cluster_ids, this.unique_cluster_ids);
            assert(all(tf), 'Some cluster_ids not found in ss.unique_cluster_ids');
        end
        
        function [channelInds, channelIds] = lookup_channelIdsForSnippet(ss, snippetInd, channelIds)
            assert(isscalar(snippetInd), 'snippetInd must be scalar');
            [tf, channelInds] = ismember(channelIds, ss.channel_ids_by_snippet(:, snippetInd));
            assert(all(tf), 'Not all channel_ids found in ss.channel_ids_by_snippet for this snippet''s cluster');
        end
    end
    
    methods % plotting
        function [hdata, settings] = plotAtProbeLocations(this, varargin)
            p = inputParser();
            % specify one of these
            p.addParameter('cluster_ids', [], @(x) isempty(x) || isvector(x));
            p.addParameter('maskSnippets', this.valid, @isvector);
            
            % and optionally these
            p.addParameter('maskTime', true(this.nTimepoints, 1), @isvector);
            p.addParameter('maskChannels', true(this.nChannels, 1), @isvector);
            p.addParameter('maxChannelsPerCluster', NaN, @isscalar);
            
            p.addParameter('maxPerCluster', Inf, @isscalar);
            
            p.addParameter('xmag', 0.8, @isscalar);
            p.addParameter('ymag', 0.8, @isscalar);
            p.addParameter('gain', NaN, @isscalar);
            
            p.addParameter('center', true, @islogical);
            p.addParameter('applyScaling', false, @islogical);
            
            p.addParameter('xoffset', 0, @isscalar);
            p.addParameter('yoffset', 0, @isscalar);
            
            p.addParameter('cluster_lags', [], @(x) isempty(x) || isvector(x));
            p.addParameter('xoffsetBetweenClusters', 0, @isscalar);
            p.addParameter('yoffsetBetweenClusters', 0, @isscalar);
            
            p.addParameter('showIndividual', true, @islogical);
            p.addParameter('showMean', false, @islogical);
            p.addParameter('showMedian', false, @islogical);
            p.addParameter('meanPlotArgs', {'LineWidth', 3}, @iscell);
            p.addParameter('medianPlotArgs', {'LineWidth', 3}, @iscell);
            
            p.addParameter('plotArgs', {'LineWidth', 0.5}, @iscell);
            
            p.addParameter('showChannelLabels', true, @islogical);
            p.addParameter('labelArgs', {}, @iscell);
            p.addParameter('alpha', 0.4, @isscalar);
            
            p.addParameter('xlim', [], @(x) true);
            p.addParameter('ylim', [], @(x) true);
            
            p.addParameter('colormap', get(groot, 'DefaultAxesColorOrder'), @(x) true);
            
            p.addParameter('suppressWarnings', false, @islogical);
            p.addParameter('axh', [], @(x) true);
            p.KeepUnmatched = false;
            p.parse(varargin{:});
            
            yspacing = this.channelMap.yspacing;
            xspacing = this.channelMap.xspacing;
            xmag = p.Results.xmag;
            ymag = p.Results.ymag;
            xoffset = p.Results.xoffset;
            yoffset = p.Results.yoffset;
            xoffsetBetweenClusters = p.Results.xoffsetBetweenClusters;
            yoffsetBetweenClusters = p.Results.yoffsetBetweenClusters;
            cluster_lags = double(p.Results.cluster_lags);
            if isempty(p.Results.cluster_ids) && ~isempty(cluster_lags)
                error('Cluster_ids must be specified when cluster_lags is specified');
            end
            specified_cluster_ids = p.Results.cluster_ids;
            
            axh = p.Results.axh;
            if isempty(axh), axh = gca; end
            
            % plot relative time vector
            tvec = linspace(0, xspacing * xmag, numel(this.window(1) : this.window(2)));
            dt = tvec(2) - tvec(1);
            
            if ~isempty(specified_cluster_ids)
                maskSnippets = ismember(this.cluster_ids, specified_cluster_ids) & this.valid;
            else
                maskSnippets = p.Results.maskSnippets;
            end
            
            hdata = struct();
            settings = struct();
            
            % center and scale each channel
            if ~any(maskSnippets)
                if ~p.Results.suppressWarnings
                    fprintf('Warning: no snippets for these clusters\n');
                end
                axis(axh, 'off');
                return;
            end
            
            if isinteger(this.data)
                data = single(this.data(p.Results.maskChannels, p.Results.maskTime, maskSnippets));
            else
                data = this.data(p.Results.maskChannels, p.Results.maskTime, maskSnippets);
            end
            
            if p.Results.applyScaling
                data = single(data) .* single(this.scaleToUv(maskSnippets));
            end
            
            if p.Results.center
                data = data - mean(data, 2); %#ok<*PROPLC>
            end
            
            gain = p.Results.gain;
            if isnan(gain)
                minmax = quantile(data, [0.025 0.975], 'all');
                gain = yspacing * ymag ./ (minmax(2) - minmax(1));
            end
            data = data .* gain;
            
            if ~isempty(specified_cluster_ids)
                % keep the list the same because colormap may be passed in
                unique_cluster_ids = specified_cluster_ids;
            else
                unique_cluster_ids = unique(this.cluster_ids(maskSnippets));
            end
            nUniqueClusters = numel(unique_cluster_ids);
            
            holding = ishold(axh);
            nChannels = size(data, 1);
            if isnan(p.Results.maxChannelsPerCluster)
                nChannelsMaxPerCluster = nChannels;
            else
                nChannelsMaxPerCluster = min(nChannels, p.Results.maxChannelsPerCluster);
            end
            handlesIndiv = cell(nChannelsMaxPerCluster, nUniqueClusters);
            handlesMean = gobjects(nChannelsMaxPerCluster, nUniqueClusters);
            handlesMedian = gobjects(nChannelsMaxPerCluster, nUniqueClusters);
            channels_plotted = false(this.channelMap.nChannels, 1);
            
            cmap = p.Results.colormap;
            if isa(cmap, 'function_handle')
                cmap = cmap(nUniqueClusters);
            end
            
            for iClu = 1:nUniqueClusters
                this_cluster_id = unique_cluster_ids(iClu);
                if isempty(cluster_lags)
                    this_cluster_lag = 0;
                else
                    this_cluster_lag = cluster_lags(find(specified_cluster_ids == this_cluster_id, 1));
                end
                
                % data is already sliced using maskSnippets, so we have a subset (size of nnz(maskSnippets)) and full mask (size of nSnippets)
                this_snippet_mask_subset = this.cluster_ids(maskSnippets) == this_cluster_id;
                if isfinite(p.Results.maxPerCluster)
                    idxKeep = find(this_snippet_mask_subset, p.Results.maxPerCluster, 'first');
                    if ~isempty(idxKeep)
                        this_snippet_mask_subset(idxKeep(end)+1:end) = false;
                    end
                end
                this_snippet_mask_full = maskSnippets;
                this_snippet_mask_full(~this_snippet_mask_subset) = false;
                
                % assume that all snippets for this cluster use the same channels
                if ~any(this_snippet_mask_full)
                    continue;
                end
                this_channel_ids = this.channel_ids_by_snippet(:, this_snippet_mask_full);
                this_channel_ids = this_channel_ids(:, 1);
                this_channel_ids = this_channel_ids(p.Results.maskChannels);
                this_channel_ids = this_channel_ids(1:nChannelsMaxPerCluster);
                
                this_channel_inds = this.channelMap.lookup_channelIds(this_channel_ids); % in channel map
                xc = this.channelMap.xcoords(this_channel_inds);
                yc = this.channelMap.ycoords(this_channel_inds);
                
                channels_plotted(this_channel_inds) = true;
                
                cmapIdx = mod(iClu-1, size(cmap, 1))+1;
                
                xc = xc + this_cluster_lag*dt + xoffset + xoffsetBetweenClusters*(iClu-1);
                yc = yc + yoffset + yoffsetBetweenClusters*(iClu-1);
                
                % arrange things carefully to plot all channels simultaneously spaced out at probe locations
                
                % ch x time x snippets
                data_clu = data(1:nChannelsMaxPerCluster, :, this_snippet_mask_subset) + yc(1:nChannelsMaxPerCluster);
                data_plot = permute(data_clu, [2 3 1]); % time x snippets x ch
                
                if p.Results.showIndividual
                    nSnippetsThis = size(data_plot, 2);
                    X = repmat(tvec', [1 nSnippetsThis nChannelsMaxPerCluster]) + shiftdim(xc, -2);
                    hthis = plot(axh, X(:, :), data_plot(:, :), 'Color', [cmap(cmapIdx, :), p.Results.alpha], ...
                        'AlignVertexCenters', true, p.Results.plotArgs{:}); % [snippets for ch1; snippets for ch2; ... ]
                    % split handles by channel
                    handlesIndiv(:, iClu) = mat2cell(hthis, repmat(nSnippetsThis, nChannelsMaxPerCluster, 1), 1);
                    npxutils.internal.graphics.hideInLegend(hthis);
                    hold(axh, 'on');
                end
                
                if p.Results.showMean
                    handlesMean(:, iClu) = plot(axh, tvec' + xc', squeeze(mean(data_plot, 2)), 'Color', cmap(cmapIdx, :), ...
                        'AlignVertexCenters', true, p.Results.meanPlotArgs{:});
                    hold(axh, 'on');
                end
                if p.Results.showMedian
                    handlesMedian(:, iClu) = plot(axh, tvec' + xc', squeeze(median(data_plot, 2)), 'Color', cmap(cmapIdx, :), ...
                        'AlignVertexCenters', true, p.Results.medianPlotArgs{:});
                    hold(axh, 'on');
                end
                
                if p.Results.showMean
                    npxutils.internal.graphics.showFirstInLegend(handlesMean(1, iClu), sprintf('cluster %d', this_cluster_id));
                elseif p.Results.showMedian
                    npxutils.internal.graphics.showFirstInLegend(handlesMedian(1, iClu), sprintf('cluster %d', this_cluster_id));
                else
                    npxutils.internal.graphics.showFirstInLegend(handlesIndiv{1, iClu}, sprintf('cluster %d', this_cluster_id));
                end
                %
                %                 for iC = 1:nChannelsMaxPerCluster
                %                     dmat = npxutils.internal.TensorUtils.squeezeDims(data(iC, :, this_snippet_mask_subset), 1) + yc(iC);
                %
                %                     if p.Results.showIndividual
                %                         handlesIndiv{iC, iClu} = plot(axh, tvec + xc(iC), dmat, 'Color', [cmap(cmapIdx, :), p.Results.alpha], ...
                %                             'AlignVertexCenters', true, p.Results.plotArgs{:});
                %                         npxutils.internal.graphics.hideInLegend(handlesIndiv{iC, iClu});
                %                         hold(axh, 'on');
                %                     end
                %                     if p.Results.showMean
                %                         handlesMean(iC, iClu) = plot(axh, tvec + xc(iC), mean(dmat, 2), 'Color', cmap(cmapIdx, :), ...
                %                             'AlignVertexCenters', true, p.Results.meanPlotArgs{:});
                %                         hold(axh, 'on');
                %                     end
                %                     if p.Results.showMedian
                %                         handlesMedian(iC, iClu) = plot(axh, tvec + xc(iC), median(dmat, 2), 'Color', cmap(cmapIdx, :), ...
                %                             'AlignVertexCenters', true, p.Results.medianPlotArgs{:});
                %                         hold(axh, 'on');
                %                     end
                %                     if iC == 1
                %                         if p.Results.showMean
                %                             npxutils.internal.graphics.showFirstInLegend(handlesMean(iC, iClu), sprintf('cluster %d', this_cluster_id));
                %                         elseif p.Results.showMedian
                %                             npxutils.internal.graphics.showFirstInLegend(handlesMedian(iC, iClu), sprintf('cluster %d', this_cluster_id));
                %                         else
                %                             npxutils.internal.graphics.showFirstInLegend(handlesIndiv{iC, iClu}, sprintf('cluster %d', this_cluster_id));
                %                         end
                %                     else
                %                         npxutils.internal.graphics.hideInLegend(handlesMean(iC, iClu));
                %                     end
                %                 end
                %
                %drawnow;
            end
            
            if p.Results.showChannelLabels
                xc = this.channelMap.xcoords + xoffset;
                yc = this.channelMap.ycoords + yoffset;
                
                for iC = 1:this.channelMap.nChannels
                    if channels_plotted(iC)
                        text(xc(iC), yc(iC), sprintf('ch %d', this.channelMap.channelIds(iC)), ...
                            'Parent', axh, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            p.Results.labelArgs{:});
                    end
                end
            end
            
            if ~holding
                hold(axh, 'off');
            end
            axis(axh, 'off');
            axis(axh, 'tight');
            if ~isempty(p.Results.xlim)
                xlim(axh, p.Results.xlim);
            end
            if ~isempty(p.Results.ylim)
                ylim(axh, p.Results.ylim);
            end
            box(axh, 'off');
            
            %             if p.Results.showIndividual
            %                 hall = cat(1, handlesIndiv{:});
            %                 for iH = 1:numel(hall)
            %                     hall(iH).Color(4) = p.Results.alpha;
            %                 end
            %             end
            
            hdata.waveforms = handlesIndiv;
            hdata.waveformMeans = handlesMean;
            
            settings.ymag = ymag;
            settings.xmag = xmag;
            settings.gain = gain;
            settings.xlim = xlim(axh);
            settings.ylim = ylim(axh);
            settings.xoffset = xoffset;
            settings.yoffset = yoffset;
            settings.xoffsetBetweenClusters = xoffsetBetweenClusters;
            settings.yoffsetBetweenClusters = yoffsetBetweenClusters;
            settings.cluster_lags = cluster_lags;
        end
        
        function [bg, traceColor] = setupFigureAxes(this, dark) %#ok<INUSL>
            if dark
                bg = [16 16 16]/255;
                traceColor = 1 - bg;
                set(gcf, 'Color', bg);
                set(gca, 'Color', 'none');
            else
                bg = [1 1 1];
                traceColor = [0 0 0];
                set(gcf, 'Color', get(groot, 'DefaultFigureColor'));
                set(gca, 'Color', get(groot, 'DefaultAxesColor'));
            end
        end
        
        function [overlay_labels, overlay_datatip_label, overlay_datatip_values] ...
                = buildWaveformOverlays(this, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', this.valid, @isvector);
            p.addParameter('cluster_ids', this.overlay_cluster_ids, @isvector);
            p.addParameter('sortChannelsVertically', false, @islogical);
            p.addParameter('best_n_channels', 24, @isscalar);
            p.parse(varargin{:});
            
            assert(~isempty(this.ks), 'KilosortDataset must be provided in .ks for waveform building');
            ks = this.ks;
            
            best_n_channels = p.Results.best_n_channels;
            
            m = ks.computeMetrics();
            [cluster_ind, cluster_ids] = ks.lookup_clusterIds(p.Results.cluster_ids);
            
            % loop through all the extracted snippets, and for each cluster, find that cluster's spikes, and highlight the surrounding
            % waveform window by creating an overlay that matches the data within this window and is NaN elsewhere
            
            cluster_best_channel_ids = m.cluster_best_channels(cluster_ind, 1:best_n_channels);
            if p.Results.sortChannelsVertically
                cluster_best_channel_ids = m.channelMap.sortChannelsVertically(cluster_best_channel_ids);
            end
            snippet_inds = find(p.Results.maskSnippets);
            nSnippets = numel(snippet_inds);
            times = this.sample_idx;
            labels = zeros(this.nChannels, this.nTimepoints, nSnippets, 1, 'uint32');
            
            for iiSn = 1:nSnippets
                iSn = snippet_inds(iiSn);
                channel_ids_this = this.channel_ids_by_snippet(:, iSn);
                
                % find all spike times in this window
                this_window = int64(times(iSn)) + this.window;
                mask_spikes = ks.spike_times >= this_window(1) & ks.spike_times <= this_window(2);
                
                % that belong to a cluster in cluster_ids
                mask_spikes(mask_spikes) = ismember(ks.spike_clusters(mask_spikes), cluster_ids);
                
                ind_spikes = find(mask_spikes);
                selected_cluster_ids = ks.spike_clusters(mask_spikes);
                [~, selected_cluster_inds] = ismember(selected_cluster_ids, cluster_ids);
                
                % mark each spike in the labels matrix
                for iSp = 1:numel(ind_spikes)
                    time_inds = ks.templateTimeRelative + int64(ks.spike_times(ind_spikes(iSp)))  - int64(times(iSn)) - int64(this.window(1)) + int64(1);
                    % chop off time inds before and after window
                    mask = time_inds >= 1 & time_inds <= this.nTimepoints;
                    time_inds = time_inds(mask);
                    
                    [ch_included, ch_inds] = ismember(cluster_best_channel_ids(selected_cluster_inds(iSp), :), channel_ids_this);
                    
                    labels(ch_inds(ch_included), time_inds, iSn) = selected_cluster_inds(iSp);
                end
            end
            
            overlay_labels = labels;
            overlay_datatip_label = 'cluster';
            overlay_datatip_values = cluster_ids;
        end
        
        function [templates, template_start_ind, spike_inds, cluster_ids] ...
                = buildTemplateOverlays(this, varargin)
            % temlates is ch x time x nTemplates reconstructed templates. the templates will be optionally nan'ed out at time
            % points that lie
            
            p = inputParser();
            p.addParameter('maskSnippets', this.valid, @isvector);
            p.addParameter('cluster_ids', this.overlay_cluster_ids, @isvector);
            p.addParameter('best_n_channels', 24, @isscalar);
            p.addParameter('sortChannelsVertically', false, @isscalar);
            p.addParameter('nanOutsideSnippet', true, @islogical);
            p.addParameter('nanOutsideBestChannels', true, @islogical);
            
            p.parse(varargin{:});
            
            assert(~isempty(this.ks), 'KilosortDataset must be provided in .ks for waveform building');
            ks = this.ks;
            
            nanOutsideSnippet = p.Results.nanOutsideSnippet;
            nanOutsideBestChannels = p.Results.nanOutsideBestChannels;
            
            if islogical(p.Results.maskSnippets)
                snippet_inds = find(p.Results.maskSnippets);
            else
                snippet_inds = p.Results.maskSnippets;
            end
            
            m = ks.computeMetrics();
            [cluster_ind, cluster_ids] = ks.lookup_clusterIds(p.Results.cluster_ids);
            best_n_channels = p.Results.best_n_channels;
            cluster_best_channel_ids = m.cluster_best_channels(cluster_ind, 1:best_n_channels);
            if p.Results.sortChannelsVertically
                cluster_best_channel_ids = m.channelMap.sortChannelsVertically(cluster_best_channel_ids);
            end
            [templates, template_start_ind, spike_inds, cluster_ids_by_spike] = deal(cell(numel(snippet_inds), 1));
            for iiSn = 1:numel(snippet_inds)
                iSn = snippet_inds(iiSn);
                this_window = int64(this.sample_idx(iSn)) + this.window;
                %                 cluster_ind_this = ss.lookup_clusterIds(ss.cluster_ids(iSn));
                ch_ids_this = this.channel_ids_by_snippet(:, iSn);
                
                [spike_inds_this, template_start_inds_this, templates_this] = ...
                    ks.findSpikesOverlappingWithWindow(this_window, 'channel_ids', ch_ids_this, 'cluster_ids', cluster_ids);
                
                % templates_this is scaled to uv, but our data is unscaled
                templates_this = templates_this ./ ks.apScaleToUv;
                
                cluster_ids_this = ks.spike_clusters(spike_inds_this);
                % NaN out regions of the templates that are outside the plotting window or not on the best channels
                % template_this is C x T x nSpikes
                for iSp = 1:numel(spike_inds_this)
                    [~, cluster_ind_this] = ismember(cluster_ids_this(iSp), cluster_ids);
                    
                    if nanOutsideBestChannels
                        ch_mask_this = ismember(ch_ids_this, cluster_best_channel_ids(cluster_ind_this, :));
                        templates_this(~ch_mask_this, :, iSp) = NaN;
                    end
                    
                    if nanOutsideSnippet
                        template_tvec_this = template_start_inds_this(iSp)+int64(1:ks.nTemplateTimepoints);
                        t_mask_this = template_tvec_this >= 1 & template_tvec_this <= this.nTimepoints;
                        templates_this(:, ~t_mask_this, iSp) = NaN;
                    end
                end
                
                mask_templates = squeeze(any(~isnan(templates_this), [1 2]));
                spike_inds{iiSn} = spike_inds_this(mask_templates);
                template_start_ind{iiSn} = template_start_inds_this(mask_templates);
                templates{iiSn} = templates_this(:, :, mask_templates);
                cluster_ids_by_spike{iiSn} = cluster_ids_this(mask_templates);
            end
            
            spike_inds = cat(1, spike_inds{:});
            template_start_ind = cat(1, template_start_ind{:});
            templates = cat(3, templates{:});
            cluster_ids = cat(1, cluster_ids_by_spike{:});
        end
        
        function [recon, cluster_ids] = buildTemplateReconstructions(this, varargin)
            [templates, template_start_ind, spike_inds, cluster_ids] = this.buildTemplateOverlays(varargin{:});
            
            recon = zeros(this.nChannels, this.nTimepoints, size(templates, 3), 'like', this.data);
            nTemplateTime = size(templates, 2);
            for iSp = 1:numel(spike_inds)
                insert_idx = int64(template_start_ind(iSp)) + int64(0:nTemplateTime-1);
                take_mask = insert_idx >= int64(1) & insert_idx <= size(recon, 2);
                insert_idx = insert_idx(take_mask);
                
                recon(:, insert_idx, iSp) = templates(:, take_mask, iSp);
            end
        end
        
        function plotStackedTraces(this, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', this.valid, @isvector);
            p.addParameter('maskChannels', true(this.nChannels, 1), @isvector);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 1.95, @isscalar);
            %             p.addParameter('car', false, @islogical);
            %             p.addParameter('downsample',1, @isscalar);
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.addParameter('dark', false, @islogical);
            p.addParameter('lineOpacity', 1, @isscalar);
            p.addParameter('showChannelDataTips', false, @islogical);
            p.parse(varargin{:});
            
            dark = p.Results.dark;
            [~, traceColor] = this.setupFigureAxes(dark);
            
            if p.Results.timeInMilliseconds
                time = this.time_ms;
            else
                time = this.time_samples;
            end
            mask_snippets = p.Results.maskSnippets;
            
            mask_channels = p.Results.maskChannels;
            data = this.data(mask_channels, :, mask_snippets, :);
            channel_ids = this.channel_ids_by_snippet(mask_channels, mask_snippets);
            
            % data is ch x time x snippets x layers
            % plot stacked traces wants time x ch x layers
            data_txcxl = permute(data(:, :, :), [2 1 3]);
            
            npxutils.internal.plotStackedTraces(time, data_txcxl, 'colors', traceColor, ...
                'lineWidth', 0.5, 'lineOpacity', p.Results.lineOpacity, ...
                'gain', p.Results.gain, 'invertChannels', this.channelMap.invertChannelsY, 'normalizeEach', false, ...
                'channel_ids', channel_ids, 'showChannelDataTips', p.Results.showChannelDataTips);
            
            axis off;
            hold off;
        end
        
        function plotStackedTracesWithOverlays(this, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', this.valid, @isvector);
            p.addParameter('maskChannels', true(this.nChannels, 1), @isvector);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 1.95, @isscalar);
            %             p.addParameter('car', false, @islogical);
            %             p.addParameter('downsample',1, @isscalar);
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.addParameter('dark', false, @islogical);
            p.addParameter('lineOpacity', 1, @isscalar);
            p.addParameter('showChannelDataTips', false, @islogical);
            p.addParameter('showOverlayDataTips', false, @islogical);
            p.addParameter('showTemplateDataTips', false, @islogical);
            
            p.addParameter('overlay_waveforms', false, @islogical);
            p.addParameter('overlay_templates', false, @islogical);
            p.addParameter('overlay_cluster_ids', this.overlay_cluster_ids, @isvector);
            p.addParameter('overlay_best_channels', 24, @isscalar);
            
            p.parse(varargin{:});
            
            dark = p.Results.dark;
            [bg, traceColor] = this.setupFigureAxes(dark);
            
            if p.Results.timeInMilliseconds
                time = this.time_ms;
            else
                time = this.time_samples;
            end
            mask_snippets = p.Results.maskSnippets;
            
            data = this.data(:, :, mask_snippets, :);
            
            overlay_cluster_ids = p.Results.overlay_cluster_ids;
            if isempty(overlay_cluster_ids) && ~isempty(this.ks)
                % overlay all clusters if none specified
                overlay_cluster_ids = this.ks.cluster_ids;
            end
            overlay_best_channels = p.Results.overlay_best_channels;
            overlay_cluster_colormap = this.generateOverlayColormap(1:numel(overlay_cluster_ids), [bg; traceColor]);
            if p.Results.overlay_waveforms
                [overlay_labels, overlay_datatip_label, overlay_datatip_values] = this.buildWaveformOverlays('maskSnippets', mask_snippets, ...
                    'cluster_ids', overlay_cluster_ids, ...
                    'best_n_channels', overlay_best_channels, 'sortChannelsVertically', true);
                overlay_colormap = overlay_cluster_colormap;
                
            elseif ~isempty(this.overlay_labels)
                overlay_labels = this.overlay_labels(:, :, mask_snippets, :);
                overlay_datatip_label = this.overlay_datatip_label;
                overlay_datatip_values = this.overlay_datatip_values;
                
                overlay_colormap = this.generateOverlayColormap(overlay_labels, [bg; traceColor]);
            else
                [overlay_labels, overlay_datatip_label, overlay_datatip_values] = deal([]);
                overlay_colormap = [];
            end
            
            % not supported by waveform or template overlays yet
            %             if p.Results.downsample > 1
            %                 data = data(:, 1:p.Results.downsample:end, :, :);
            %                 overlay_labels = overlay_labels(:, 1:p.Results.downsample:end, :, :);
            %                 time = time(1:p.Results.downsample:end);
            %             end
            
            % nCh x nSnippets
            channel_ids_ch_by_snippet = this.channel_ids_by_snippet(:, mask_snippets);
            
            % not supported by waveform or template overlays yet
            %             if p.Results.car
            %                 data = data - median(data, 1);
            %             end
            
            % data is ch x time x snippets x layers
            % plot stacked traces wants time x ch x layers
            data_txcxl = permute(data(:, :, :), [2 1 3]);
            overlay_labels_txcxl = permute(overlay_labels(:, :, :), [2 1 3]);
            
            [~, transform] = npxutils.internal.plotStackedTraces(time, data_txcxl, 'colors', traceColor, ...
                'lineWidth', 0.5, 'lineOpacity', p.Results.lineOpacity, ...
                'gain', p.Results.gain, 'invertChannels', this.channelMap.invertChannelsY, 'normalizeEach', false, ...
                'colorOverlayLabels', overlay_labels_txcxl, 'colorOverlayMap', overlay_colormap, ...
                'channel_ids', channel_ids_ch_by_snippet, 'showChannelDataTips', p.Results.showChannelDataTips, ...
                'showOverlayDataTips', p.Results.showOverlayDataTips, ...
                'colorOverlayDataTipLabel', overlay_datatip_label, 'colorOverlayDataTipValues', overlay_datatip_values);
            
            % plot templates if requested
            if p.Results.overlay_templates
                [templates, template_start_ind, spike_inds, cluster_ids] = this.buildTemplateOverlays('maskSnippets', mask_snippets, ...
                    'cluster_ids', overlay_cluster_ids, 'sortChannelsVertically', true, ...
                    'best_n_channels', p.Results.overlay_best_channels);
                
                [~, cluster_inds] = ismember(cluster_ids, overlay_cluster_ids);
                
                % generate expanded time vector we can insert templates into accounting for timepoints outside the snippet
                % where a template could start
                nTemplateTime = size(templates, 2);
                time_exp =  int64(this.window(1) + 1) - int64(nTemplateTime) : int64(this.window(2) - 1) + int64(nTemplateTime);
                if p.Results.timeInMilliseconds
                    time_exp = single(time_exp) / this.fs * 1000;
                end
                
                for iSp = 1:numel(spike_inds)
                    % template_start_ind corresponds to time, so we shift our start index to the right to account for the "pre-snippet"
                    % times we included in time_exp
                    time_this = time_exp( int64(template_start_ind(iSp)) + int64(nTemplateTime-1) + int64(0:nTemplateTime-1));
                    template_this = templates(:, :, iSp);
                    t_mask = any(~isnan(template_this), 1);
                    cluster_ind_this = cluster_inds(iSp);
                    if p.Results.showTemplateDataTips
                        dataTipArgs = {'dataTipLabel', 'cluster template', 'dataTipValues', cluster_ids(iSp), 'dataTipFormat', '%d'};
                    else
                        dataTipArgs = {};
                    end
                    hold on;
                    npxutils.internal.plotStackedTraces(time_this(t_mask), template_this(:, t_mask)', ...
                        'colors', overlay_cluster_colormap(cluster_ind_this, :), ...
                        'transformInfo', transform, ...
                        dataTipArgs{:});
                end
            end
            axis off;
            hold off;
        end
        
        function plotStackedTracesColorClusters(this, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', this.valid, @isvector);
            p.addParameter('maxPerCluster', 50, @isscalar);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 1.95, @isscalar);
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.addParameter('dark', false, @islogical);
            p.addParameter('lineOpacity', 1, @isscalar);
            p.parse(varargin{:});
            
            mask_snippets = p.Results.maskSnippets;
            unique_cluster_ids = this.unique_cluster_ids;
            
            maxPerCluster = p.Results.maxPerCluster;
            if isfinite(maxPerCluster)
                mask_inds = find(mask_snippets);
                cluster_ids = this.cluster_ids(mask_snippets);
                [~, cluster_inds] = ismember(cluster_ids, unique_cluster_ids);
                for iC = 1:numel(unique_cluster_ids)
                    this_cluster = find(cluster_inds == iC);
                    N = numel(this_cluster);
                    if N > maxPerCluster
                        drop = randsample(N, N - maxPerCluster, false);
                        mask_snippets(mask_inds(this_cluster(drop))) = false;
                    end
                end
            end
            
            dark = p.Results.dark;
            [bg, ~] = this.setupFigureAxes(dark);
            cluster_colormap = this.generateOverlayColormap(1:numel(unique_cluster_ids), bg);
            
            % nCh x nSnippets
            channel_ids_ch_by_snippet = this.channel_ids_by_snippet(:, mask_snippets);
            data = this.data(:, :, mask_snippets, :);
            data_txcxl = permute(data(:, :, :), [2 1 3]);
            
            cluster_ids = this.cluster_ids(mask_snippets);
            [~, cluster_inds] = ismember(cluster_ids, unique_cluster_ids);
            traceColors = cluster_colormap(cluster_inds, :);
            
            if p.Results.timeInMilliseconds
                time = this.time_ms;
            else
                time = this.time_samples;
            end
            
            npxutils.internal.plotStackedTraces(time, data_txcxl, 'layerColors', traceColors, ...
                'lineWidth', 0.5, 'lineOpacity', p.Results.lineOpacity, ...
                'gain', p.Results.gain, 'invertChannels', this.channelMap.invertChannelsY, 'normalizeEach', false, ...
                'channel_ids', channel_ids_ch_by_snippet, 'showChannelDataTips', false);
            
            axis off;
            hold off;
        end
        
        function plotHeatmap(this, snippet_ind, varargin)
            % plots one snippet as a heatmap with each of the clusters that occur during that period
            
            p = inputParser();
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.parse(varargin{:});
            
            if p.Results.timeInMilliseconds
                time = this.time_ms;
            else
                time = this.time_samples;
            end
            data = this.data(:, :, snippet_ind, 1);
            
            npxutils.internal.pmatbal(data, 'x', time);
            yline(0.5, '-', sprintf('snippet %d', snippet_ind), 'LabelVerticalAlignment', 'bottom');
        end
        
        function plotHeatmapWithTemplates(this, snippet_ind, varargin)
            % plots one snippet as a heatmap with each of the clusters that occur during that period
            
            p = inputParser();
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.addParameter('reconstruct_cluster_ids', this.overlay_cluster_ids, @isvector);
            p.addParameter('reconstruct_best_channels', 24, @isscalar);
            p.addParameter('successive_residuals', false, @islogical);
            p.parse(varargin{:});
            
            if p.Results.timeInMilliseconds
                time = this.time_ms;
            else
                time = this.time_samples;
            end
            data = this.data(:, :, snippet_ind, 1);
            
            reconstruct_cluster_ids = p.Results.reconstruct_cluster_ids;
            if isempty(reconstruct_cluster_ids) && ~isempty(this.ks)
                % overlay all clusters if none specified
                reconstruct_cluster_ids = this.ks.cluster_ids;
            end
            reconstruct_best_channels = p.Results.reconstruct_best_channels;
            
            % nCh x nSnippets
            %             cluster_inds_by_snippet = ss.lookup_clusterIds(ss.cluster_ids(snippet_ind));
            %             channel_ids_ch_by_snippet = ss.channel_ids_by_cluster(:, cluster_inds_by_snippet);
            
            [templates, cluster_ids] = this.buildTemplateReconstructions('maskSnippets', snippet_ind, ...
                'cluster_ids', reconstruct_cluster_ids, ...
                'best_n_channels', reconstruct_best_channels, 'sortChannelsVertically', true);
            
            residual = data;
            if p.Results.successive_residuals
                for iT = 1:size(templates, 3)
                    residual = residual - templates(:, :, iT);
                    templates(:, :, iT) = residual;
                end
            end
            templates_stacked = npxutils.internal.TensorUtils.reshapeByConcatenatingDims(templates, {[1 3], 2});
            
            
            data_stacked = cat(1, data, templates_stacked);
            
            npxutils.internal.pmatbal(data_stacked, 'x', time);
            
            yline(0.5, '-', sprintf('snippet %d', snippet_ind), 'LabelVerticalAlignment', 'bottom');
            for iC = 1:size(templates, 3)
                yline(this.nChannels * iC + 0.5, 'k-', sprintf('cluster %u', cluster_ids(iC)), 'LabelVerticalAlignment', 'bottom');
            end
            
            
        end
        
        function plotHeatmapWithReconstruction(this, snippet_ind, varargin)
            % plots one snippet as a heatmap with all of the overlapping clusters' reconstructed templates below
            
            default_cluster_reconstruct = this.ks.cluster_ids;
            hasClusterId = ~isempty(this.cluster_ids);
            if hasClusterId
                this_cluster_id = this.cluster_ids(snippet_ind);
                default_cluster_reconstruct = setdiff(default_cluster_reconstruct, this_cluster_id);
            end
            
            p = inputParser();
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.addParameter('reconstruct_cluster_ids', default_cluster_reconstruct, @isvector);
            p.addParameter('reconstruct_best_channels', this.ks.nChannelsSorted, @isscalar);
            p.parse(varargin{:});
            
            if p.Results.timeInMilliseconds
                time = this.time_ms;
            else
                time = this.time_samples;
            end
            data = this.data(:, :, snippet_ind, 1);
            
            reconstruct_cluster_ids = p.Results.reconstruct_cluster_ids;
            if isempty(reconstruct_cluster_ids) && ~isempty(this.ks)
                % overlay all clusters if none specified
                reconstruct_cluster_ids = this.ks.cluster_ids;
            end
            reconstruct_best_channels = p.Results.reconstruct_best_channels;
            
            if hasClusterId
                % build this cluster's tempalte
                this_template = this.buildTemplateReconstructions('maskSnippets', snippet_ind, ...
                    'cluster_ids', this_cluster_id, 'sortChannelsVertically', true, ...
                    'best_n_channels', reconstruct_best_channels);
            end
            
            % build other clusters' templates
            [templates, cluster_ids_recon] = this.buildTemplateReconstructions('maskSnippets', snippet_ind, ...
                'cluster_ids', reconstruct_cluster_ids, 'sortChannelsVertically', true, ...
                'best_n_channels', reconstruct_best_channels);
            reconstruction = sum(templates, 3);
            
            if hasClusterId
                data_stacked = cat(1, data, this_template, reconstruction);
            else
                data_stacked = cat(1, data, reconstruction);
            end
            npxutils.internal.pmatbal(data_stacked, 'x', time);
            
            if ~hasClusterId
                labelMain = sprintf('snippet %d', snippet_ind);
            else
                labelMain = sprintf('snippet %d (cluster %d)', snippet_ind, this.cluster_ids(snippet_ind));
            end
            yline(0.5, '-', labelMain, 'LabelVerticalAlignment', 'bottom');
            
            nOther = numel(unique(cluster_ids_recon));
            if nOther > 15
                other_clusters_str = sprintf('%d clusters', nOther);
            else
                other_clusters_str = strjoin(string(unique(cluster_ids_recon)), ',');
            end
            if hasClusterId
                yline(this.nChannels + 0.5, 'k-', sprintf('reconstruction cluster %d', this_cluster_id), 'LabelVerticalAlignment', 'bottom');
                yline(this.nChannels*2 + 0.5, 'k-', sprintf('reconstruction other clusters (%s)', other_clusters_str), 'LabelVerticalAlignment', 'bottom');
            else
                yline(this.nChannels + 0.5, 'k-', sprintf('reconstruction all clusters (%s)', other_clusters_str), 'LabelVerticalAlignment', 'bottom');
            end
        end
    end
end