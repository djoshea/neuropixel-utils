classdef SnippetSet < handle & matlab.mixin.Copyable
    % A collection of snippets of data collected from an IMEC data file
    % which includes the subset of channels which were sampled
    
    % data here should be thought of as a set of data for many snippets
    % 
    
    properties(Transient)
        ks  % optional, handle to KilosortDataset for plotting waveforms and templates
    end
    
    properties
        data % (:, :, :, :) int16 % channels x time x snippets x layers - we deactivate the size checking here 
        data_trust_mask logical % optional, same size as data indicating which data points are considered usable within data. If empty, all are valid
        
        % optional, set if each snippet corresponds to a specific cluster
        cluster_ids (:, 1) uint32 % nSnippets x 1, array indicating which cluster is extracted in each snippet, if this makes sense. otherwise will just be 1s
        
        group_ids (:, 1)  % nSnippets x 1 metadata indicating which group a cluster belonged to
        
        channel_ids_by_snippet (:, :, :) uint32 % channels x snippets
        
        overlay_cluster_ids (:, 1) uint32 % set of clusters whose waveforms or templates will be drawn on top of a given snippet, will be cluster_ids(i), or ks.cluster_ids
        
        % for manually specifying overlay matrices (can also be generated automatically with e.g. overlay_waveforms)
        overlay_labels(:, :, :, :) uint32 % channels x time x snippets x layers (label matrix corresponding to elements of data for the purpose of generating colored overlays)
        overlay_datatip_label (1, 1) string
        overlay_datatip_values(:, 1) % nLabels x 1
        
        sample_idx (:, 1) uint64
        trial_idx (:, 1) uint32
        window (:, 2) int64 % in samples
        
        valid (:, 1) logical
        
        channelMap % ChannelMap for coordinates
        scaleToUv (:, 1) single % one value for each snippet (will auto inflate to nSnippets if set to be scalar
        fs
    end
    
    properties(Dependent)
        nChannels
        nTimepoints
        nSnippets
        nClusters
        nLayers
        
        data_valid
        time_samples(:, 1) int64
        time_ms (:, 1) double
        
        unique_cluster_ids
    end
    
    methods
        function ss = SnippetSet(ds, type)
            if nargin > 0
                if isa(ds, 'Neuropixel.KilosortTrialSegmentedDataset')
                    raw_dataset = ds.raw_dataset;
                    ks = ds.dataset;
                elseif isa(ds, 'Neuropixel.KilosortDataset')
                    raw_dataset = ds.raw_dataset;
                    ks = ds;
                elseif isa(ds, 'Neuropixel.ImecDataset')
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
        
        function n = get.nChannels(ss)
            n = size(ss.data, 1);
        end
        
        function n = get.nTimepoints(ss)
            n = size(ss.data, 2);
        end
        
        function n = get.nSnippets(ss)
            n = size(ss.data, 3);
        end
        
        function n = get.nLayers(ss)
            n = size(ss.data, 4);
        end
        
        function n = get.nClusters(ss)
            n = numel(ss.unique_cluster_ids);
        end
        
        function v = get.valid(ss)
            if isempty(ss.valid)
                v = true(size(ss.data, 3), 1);
            else
                v = ss.valid;
            end
        end
        
        function v = get.data_valid(ss)
            v = ss.data(:, :, ss.valid);
        end
        
        function v = get.time_ms(ss)
            v = double(ss.window(1) : ss.window(2)) / ss.fs * 1000;
        end
        
        function v = get.time_samples(ss)
            v = ss.window(1) : ss.window(2);
        end
        
        function v = get.unique_cluster_ids(ss)
            v = unique(ss.cluster_ids);
        end
        
        function v = get.overlay_cluster_ids(ss)
            if isempty(ss.overlay_cluster_ids)
                if isempty(ss.unique_cluster_ids)
                    if isempty(ss.ks)
                        v = [];
                    else
                        v = ss.ks.cluster_ids;
                    end
                else
                    v = ss.unique_cluster_ids;
                end
            else
                v = ss.overlay_cluster_ids;
            end
        end
        
        function v = get.scaleToUv(ss)
            v = ss.scaleToUv;
            if isempty(v)
                v = NaN;
            end
            if numel(v) == 1
                v = repmat(v, ss.nSnippets, 1);
            end
        end
        
        function colormap = generateOverlayColormap(~, overlay_labels, backgroundColor)
            if nargin < 3
                backgroundColor = [0 0 0; 1 1 1];
            end
            nLabels = max(overlay_labels(:));
            colormap = zeros(nLabels, 3);
            unique_labels = setdiff(unique(overlay_labels(:)), 0);
            label_found = ismember(1:nLabels, unique_labels);
            colormap(label_found, :) = Neuropixel.Utils.distinguishable_colors(nnz(label_found), backgroundColor);
        end
        
        function ss = selectClusters(ss, cluster_ids)
            mask = ismember(ss.cluster_ids, cluster_ids);
            ss = ss.selectData('maskSnippets', mask);
        end 
        
        function ss = selectData(ss, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', ss.valid, @isvector);
            p.addParameter('maskTime', true(ss.nTimepoints, 1), @isvector);
            p.addParameter('maskChannels', true(ss.nChannels, 1), @isvector);
            p.parse(varargin{:});
            
            maskSnippets = p.Results.maskSnippets;
            
            ss = copy(ss);
            ss.valid = ss.valid(maskSnippets); % first time auto generates based on data, so has to come first
            ss.data = ss.data(p.Results.maskChannels, p.Results.maskTime, maskSnippets);
            
            if ~isempty(ss.overlay_labels)
                ss.overlay_labels = ss.overlay_labels(p.Results.maskChannels, p.Results.maskTime, maskSnippets);
            end
            
            if ~isempty(ss.cluster_ids)
                ss.cluster_ids = ss.cluster_ids(maskSnippets);
            end
        
            time_samples = ss.time_samples(p.Results.maskTime);
            ss.window = [time_samples(1) time_samples(end)];

            sample_offset = find(p.Results.maskTime, 1) - 1;
            ss.sample_idx = ss.sample_idx(maskSnippets) + sample_offset;
            
            if ~isempty(ss.trial_idx)
                ss.trial_idx = ss.trial_idx(maskSnippets);
            end
            
            ss.scaleToUv = ss.scaleToUv(maskSnippets);
        end
        
        function [cluster_inds, cluster_ids] = lookup_clusterIds(ss, cluster_ids)
            if islogical(cluster_ids)
                cluster_ids = ss.unique_cluster_ids(cluster_ids);
            end
            [tf, cluster_inds] = ismember(cluster_ids, ss.unique_cluster_ids);
            assert(all(tf), 'Some cluster_ids not found in ss.unique_cluster_ids');
        end
        
        function [channelInds, channelIds] = lookup_channelIdsForSnippet(ss, snippetInd, channelIds)
            assert(isscalar(snippetInd), 'snippetInd must be scalar');
            [tf, channelInds] = ismember(channelIds, ss.channel_ids_by_snippet(:, snippetInd));
            assert(all(tf), 'Not all channel_ids found in ss.channel_ids_by_snippet for this snippet''s cluster');
        end
    end
    
    methods % plotting
        function [hdata, settings, scale_info] = plotAtProbeLocations(ss, varargin)
            p = inputParser();
            % specify one of these 
            p.addParameter('cluster_ids', [], @(x) isempty(x) || isvector(x));
            p.addParameter('maskSnippets', ss.valid, @isvector);
            
            % and optionally these
            p.addParameter('maskTime', true(ss.nTimepoints, 1), @isvector);
            p.addParameter('maskChannels', true(ss.nChannels, 1), @isvector);
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
            p.addParameter('showSem', false, @islogical);
            p.addParameter('showStd', false, @islogical);

            p.addParameter('showMedian', false, @islogical);
            p.addParameter('showMad', false, @islogical);

            p.addParameter('meanPlotArgs', {'LineWidth', 1}, @iscell); 
            p.addParameter('medianPlotArgs', {'LineWidth', 1}, @iscell);
            
            p.addParameter('plotArgs', {'LineWidth', 0.5}, @iscell);
            
            p.addParameter('showChannelLabels', true, @islogical);
            p.addParameter('labelArgs', {}, @iscell);
            p.addParameter('alpha', 0.4, @isscalar);

            p.addParameter('markChannelIds', [], @isvector); % if specified, these channels may be marked with a round symbol 
            p.addParameter('markChannelColors', [0 0 0], @ismatrix); % nChannels x 3 colormap
            p.addParameter('markChannelSize', 5, @isscalar);
            p.addParameter('markChannelOffset', [0 0], @isvector);
            
            p.addParameter('xlim', [], @(x) true);
            p.addParameter('ylim', [], @(x) true);
            
            p.addParameter('colormap', get(groot, 'DefaultAxesColorOrder'), @(x) true);

            p.addParameter('suppressWarnings', false, @islogical);
            p.addParameter('axh', [], @(x) true);
            p.KeepUnmatched = false;
            p.parse(varargin{:});
            
            yspacing = ss.channelMap.yspacing;
            xspacing = ss.channelMap.xspacing;
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
            
            % plot relative time vector in ms
            nT = numel(ss.window(1) : ss.window(2));
            tvec = linspace(0, xspacing * xmag, nT);
            dt = tvec(2) - tvec(1); % dt is the actual delta x in plot

            scale_info.x_to_ms = (1000 / ss.fs) / dt; 
            
            if ~isempty(specified_cluster_ids)
                maskSnippets = ismember(ss.cluster_ids, specified_cluster_ids) & ss.valid;
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
            
            if isinteger(ss.data)
                data = single(ss.data(p.Results.maskChannels, p.Results.maskTime, maskSnippets));
            else
                data = ss.data(p.Results.maskChannels, p.Results.maskTime, maskSnippets);
            end
            
            if p.Results.applyScaling
                data = single(data) .* single(ss.scaleToUv(maskSnippets));
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
            scale_info.y_to_uv = 1 ./ gain;

            if ~isempty(ss.data_trust_mask)
                data_trust_mask = ss.data_trust_mask(p.Results.maskChannels, p.Results.maskTime, maskSnippets);
                data(~data_trust_mask) = NaN;
            end
            
            if ~isempty(specified_cluster_ids)
                % keep the list the same because colormap may be passed in
                unique_cluster_ids = specified_cluster_ids;
            else
                unique_cluster_ids = unique(ss.cluster_ids(maskSnippets));
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
            [handlesMean, handlesSem, handlesStd] = deal(gobjects(nChannelsMaxPerCluster, nUniqueClusters));
            [handlesMedian, handlesMad] = deal(gobjects(nChannelsMaxPerCluster, nUniqueClusters));
            channels_plotted = false(ss.channelMap.nChannels, 1);
            
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
                this_snippet_mask_subset = ss.cluster_ids(maskSnippets) == this_cluster_id;
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
                this_channel_ids = ss.channel_ids_by_snippet(:, this_snippet_mask_full);
                this_channel_ids = this_channel_ids(:, 1);
                this_channel_ids = this_channel_ids(p.Results.maskChannels);
                this_channel_ids = this_channel_ids(1:nChannelsMaxPerCluster);
                
                this_channel_inds = ss.channelMap.lookup_channelIds(this_channel_ids); % in channel map
                xc = ss.channelMap.xcoords(this_channel_inds);
                yc = ss.channelMap.ycoords(this_channel_inds);
                
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
                    Neuropixel.Utils.hideInLegend(hthis);
                    hold(axh, 'on');
                end
                
                if p.Results.showMean
                    data_plot_mean = squeeze(mean(data_plot, 2, 'omitnan'));

                    if p.Results.showSem
                        data_plot_sem = squeeze(Neuropixel.Utils.sem(data_plot, 2));
                        handlesSem(:, iClu) = Neuropixel.Utils.errorshadeInterval(tvec' + xc', data_plot_mean - data_plot_sem, data_plot_mean + data_plot_sem, ...
                            cmap(cmapIdx, :), 'axh', axh);
                        hold(axh, 'on');
                        Neuropixel.Utils.hideInLegend(handlesSem(:, iClu));
                    end

                    if p.Results.showStd
                        data_plot_std = squeeze(std(data_plot, 0, 2, 'omitnan'));
                        handlesStd(:, iClu) = Neuropixel.Utils.errorshadeInterval(tvec' + xc', data_plot_mean - data_plot_std, data_plot_mean + data_plot_std, ...
                            cmap(cmapIdx, :), 'axh', axh);
                        hold(axh, 'on');
                        Neuropixel.Utils.hideInLegend(handlesStd(:, iClu));
                    end
                    
                    handlesMean(:, iClu) = plot(axh, tvec' + xc', data_plot_mean, 'Color', cmap(cmapIdx, :), ...
                        'AlignVertexCenters', true, p.Results.meanPlotArgs{:});
                    hold(axh, 'on');
                    
                end
                if p.Results.showMedian
                    data_plot_median = squeeze(median(data_plot, 2, 'omitnan'));
                    if p.Results.showMad
                        data_plot_mad = squeeze(mad(data_plot, 1, 2, 'omitnan'));
                        handlesMad = Neuropixel.Utils.errorshadeInterval(tvec' + xc', data_plot_median - data_plot_mad, data_plot_median + data_plot_mad, ...
                            cmap(cmapIdx, :), 'axh', axh);
                        hold(axh, 'on');
                        Neuropixel.Utils.hideInLegend(handlesMad);
                    end
                    handlesMedian(:, iClu) = plot(axh, tvec' + xc', data_plot_median, 'Color', cmap(cmapIdx, :), ...
                        'AlignVertexCenters', true, p.Results.medianPlotArgs{:});
                    hold(axh, 'on');
                end
                
                if p.Results.showMean
                    Neuropixel.Utils.showFirstInLegend(handlesMean(1, iClu), sprintf('cluster %d', this_cluster_id));
                elseif p.Results.showMedian
                    Neuropixel.Utils.showFirstInLegend(handlesMedian(1, iClu), sprintf('cluster %d', this_cluster_id));
                else
                    Neuropixel.Utils.showFirstInLegend(handlesIndiv{1, iClu}, sprintf('cluster %d', this_cluster_id));
                end
% 
%                 for iC = 1:nChannelsMaxPerCluster
%                     dmat = Neuropixel.Utils.TensorUtils.squeezeDims(data(iC, :, this_snippet_mask_subset), 1) + yc(iC);
%                         
%                     if p.Results.showIndividual
%                         handlesIndiv{iC, iClu} = plot(axh, tvec + xc(iC), dmat, 'Color', [cmap(cmapIdx, :), p.Results.alpha], ...
%                             'AlignVertexCenters', true, p.Results.plotArgs{:});
%                         Neuropixel.Utils.hideInLegend(handlesIndiv{iC, iClu});
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
%                             Neuropixel.Utils.showFirstInLegend(handlesMean(iC, iClu), sprintf('cluster %d', this_cluster_id));
%                         elseif p.Results.showMedian
%                             Neuropixel.Utils.showFirstInLegend(handlesMedian(iC, iClu), sprintf('cluster %d', this_cluster_id));
%                         else
%                             Neuropixel.Utils.showFirstInLegend(handlesIndiv{iC, iClu}, sprintf('cluster %d', this_cluster_id));
%                         end
%                     else
%                         Neuropixel.Utils.hideInLegend(handlesMean(iC, iClu));
%                     end
%                 end
%                 
                %drawnow;
            end
            
            if p.Results.showChannelLabels
                xc = ss.channelMap.xcoords + xoffset;
                yc = ss.channelMap.ycoords + yoffset;
                
                for iC = 1:ss.channelMap.nChannels
                    if channels_plotted(iC)
                        text(xc(iC), yc(iC), sprintf('ch %d ', ss.channelMap.channelIds(iC)), ...
                            'Parent', axh, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Background', 'none',  ...
                            p.Results.labelArgs{:});
                    end
                end
            end

            if ~isempty(p.Results.markChannelIds)
                xc = ss.channelMap.xcoords + xoffset + p.Results.markChannelOffset(1);
                yc = ss.channelMap.ycoords + yoffset + yspacing/2 +  p.Results.markChannelOffset(2); 
                
                markChannels = p.Results.markChannelIds;
                markColors = p.Results.markChannelColors;
                nMark = numel(markChannels);
                [tf, markChannelInds] = ismember(markChannels, ss.channelMap.channelIds);
                assert(all(tf), 'Some markChannels could not be found in channelMap');
                if size(markColors, 1) == 1
                    markColors = repmat(markColors, nMark, 1);
                end

                handlesMarkChannels = gobjects(nMark, 1);
                for iM = 1:nMark
                    ind = markChannelInds(iM);
                    if channels_plotted(ind)
                        handlesMarkChannels(iM) = plot(xc(ind), yc(ind), 'o', MarkerFaceColor=markColors(iM, :), ...
                            MarkerEdgeColor="none", MarkerSize=p.Results.markChannelSize);
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
            if p.Results.showMean
                hdata.waveformMeans = handlesMean;
            end
            if p.Results.showSem
                hdata.waveformSem = handlesSem;
            end
            if p.Results.showStd
                hdata.waveformStd = handlesStd;
            end
            if p.Results.showMedian
                hdata.waveformMedian = handlesMedian;
            end
            if p.Results.showMad
                hdata.waveformMad = handlesMad;
            end
            
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
        
        function [bg, traceColor] = setupFigureAxes(~, dark)
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
        
        function [overlay_labels, overlay_datatip_label, overlay_datatip_values] = buildWaveformOverlays(ss, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', ss.valid, @isvector);
            p.addParameter('cluster_ids', ss.overlay_cluster_ids, @isvector);
            p.addParameter('sortChannelsVertically', false, @islogical);
            p.addParameter('best_n_channels', 24, @isscalar);
            p.parse(varargin{:});
            
            assert(~isempty(ss.ks), 'KilosortDataset must be provided in .ks for waveform building');
            ks = ss.ks;
            
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
             times = ss.sample_idx;
             labels = zeros(ss.nChannels, ss.nTimepoints, nSnippets, 1, 'uint32');
             
             for iiSn = 1:nSnippets
                 iSn = snippet_inds(iiSn);
                 channel_ids_this = ss.channel_ids_by_snippet(:, iSn);

                 % find all spike times in this window
                 this_window = int64(times(iSn)) + ss.window;
                 mask_spikes = ks.spike_times >= this_window(1) & ks.spike_times <= this_window(2);

                 % that belong to a cluster in cluster_ids
                 mask_spikes(mask_spikes) = ismember(ks.spike_clusters(mask_spikes), cluster_ids);

                 ind_spikes = find(mask_spikes);
                 selected_cluster_ids = ks.spike_clusters(mask_spikes);
                 [~, selected_cluster_inds] = ismember(selected_cluster_ids, cluster_ids);

                 % mark each spike in the labels matrix
                 for iSp = 1:numel(ind_spikes)
                    time_inds = ks.templateTimeRelative + int64(ks.spike_times(ind_spikes(iSp)))  - int64(times(iSn)) - int64(ss.window(1)) + int64(1);
                    % chop off time inds before and after window
                    mask = time_inds >= 1 & time_inds <= ss.nTimepoints;
                    time_inds = time_inds(mask);
                    
                    [ch_included, ch_inds] = ismember(cluster_best_channel_ids(selected_cluster_inds(iSp), :), channel_ids_this);

                    labels(ch_inds(ch_included), time_inds, iSn) = selected_cluster_inds(iSp);
                 end
             end

             overlay_labels = labels;
             overlay_datatip_label = 'cluster';
             overlay_datatip_values = cluster_ids;
        end
        
        function [templates, template_start_ind, spike_inds, cluster_ids] = buildTemplateOverlays(ss, varargin)
            % temlates is ch x time x nTemplates reconstructed templates. the templates will be optionally nan'ed out at time
            % points that lie 
            
            p = inputParser();
            p.addParameter('maskSnippets', ss.valid, @isvector);
            p.addParameter('cluster_ids', ss.overlay_cluster_ids, @isvector);
            p.addParameter('best_n_channels', 24, @isscalar);
            p.addParameter('sortChannelsVertically', false, @isscalar);
            p.addParameter('nanOutsideSnippet', true, @islogical);
            p.addParameter('nanOutsideBestChannels', true, @islogical);
             
            p.parse(varargin{:});
            
            assert(~isempty(ss.ks), 'KilosortDataset must be provided in .ks for waveform building');
            ks = ss.ks;
            
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
                this_window = int64(ss.sample_idx(iSn)) + ss.window;
%                 cluster_ind_this = ss.lookup_clusterIds(ss.cluster_ids(iSn));
                ch_ids_this = ss.channel_ids_by_snippet(:, iSn);
                 
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
                        t_mask_this = template_tvec_this >= 1 & template_tvec_this <= ss.nTimepoints;
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
        
        function [recon, cluster_ids] = buildTemplateReconstructions(ss, varargin)
            [templates, template_start_ind, spike_inds, cluster_ids] = ss.buildTemplateOverlays(varargin{:});
            
            recon = zeros(ss.nChannels, ss.nTimepoints, size(templates, 3), 'like', ss.data);
            nTemplateTime = size(templates, 2);
            for iSp = 1:numel(spike_inds)
                insert_idx = int64(template_start_ind(iSp)) + int64(0:nTemplateTime-1);
                take_mask = insert_idx >= int64(1) & insert_idx <= size(recon, 2);
                insert_idx = insert_idx(take_mask);
                
                recon(:, insert_idx, iSp) = templates(:, take_mask, iSp); 
            end
        end
        
        function plotStackedTraces(ss, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', ss.valid, @isvector);
            p.addParameter('maskChannels', true(ss.nChannels, 1), @isvector);
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
            [~, traceColor] = ss.setupFigureAxes(dark);

            if p.Results.timeInMilliseconds
                time = ss.time_ms;
            else
                time = ss.time_samples;
            end
            mask_snippets = p.Results.maskSnippets;

            mask_channels = p.Results.maskChannels;
            data = ss.data(mask_channels, :, mask_snippets, :);
            channel_ids = ss.channel_ids_by_snippet(mask_channels, mask_snippets);
            
            % data is ch x time x snippets x layers
            % plot stacked traces wants time x ch x layers
            data_txcxl = permute(data(:, :, :), [2 1 3]);
            
            Neuropixel.Utils.plotStackedTraces(time, data_txcxl, 'colors', traceColor, ...
                'lineWidth', 0.5, 'lineOpacity', p.Results.lineOpacity, ...
                'gain', p.Results.gain, 'invertChannels', ss.channelMap.invertChannelsY, 'normalizeEach', false, ...
                'channel_ids', channel_ids, 'showChannelDataTips', p.Results.showChannelDataTips);
             
            axis off;
            hold off;
        end
       
        function plotStackedTracesWithOverlays(ss, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', ss.valid, @isvector);
            p.addParameter('maskChannels', true(ss.nChannels, 1), @isvector);
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
            p.addParameter('overlay_cluster_ids', ss.overlay_cluster_ids, @isvector);
            p.addParameter('overlay_best_channels', 24, @isscalar);

            p.parse(varargin{:});
            
            dark = p.Results.dark;
            [bg, traceColor] = ss.setupFigureAxes(dark);

            if p.Results.timeInMilliseconds
                time = ss.time_ms;
            else
                time = ss.time_samples;
            end
            mask_snippets = p.Results.maskSnippets;

            data = ss.data(:, :, mask_snippets, :);
            
            overlay_cluster_ids = p.Results.overlay_cluster_ids;
            if isempty(overlay_cluster_ids) && ~isempty(ss.ks)
                % overlay all clusters if none specified
                overlay_cluster_ids = ss.ks.cluster_ids;
            end
            overlay_best_channels = p.Results.overlay_best_channels;
            overlay_cluster_colormap = ss.generateOverlayColormap(1:numel(overlay_cluster_ids), [bg; traceColor]);
            if p.Results.overlay_waveforms
                [overlay_labels, overlay_datatip_label, overlay_datatip_values] = ss.buildWaveformOverlays('maskSnippets', mask_snippets, ...
                    'cluster_ids', overlay_cluster_ids, ...
                    'best_n_channels', overlay_best_channels, 'sortChannelsVertically', true);
                overlay_colormap = overlay_cluster_colormap;
                
            elseif ~isempty(ss.overlay_labels)
                overlay_labels = ss.overlay_labels(:, :, mask_snippets, :);
                overlay_datatip_label = ss.overlay_datatip_label;
                overlay_datatip_values = ss.overlay_datatip_values;
                
                overlay_colormap = ss.generateOverlayColormap(overlay_labels, [bg; traceColor]);
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
            channel_ids_ch_by_snippet = ss.channel_ids_by_snippet(:, mask_snippets);
            
            % not supported by waveform or template overlays yet
%             if p.Results.car
%                 data = data - median(data, 1);
%             end
            
            % data is ch x time x snippets x layers
            % plot stacked traces wants time x ch x layers
            data_txcxl = permute(data(:, :, :), [2 1 3]);
            overlay_labels_txcxl = permute(overlay_labels(:, :, :), [2 1 3]);
            
            [~, transform] = Neuropixel.Utils.plotStackedTraces(time, data_txcxl, 'colors', traceColor, ...
                'lineWidth', 0.5, 'lineOpacity', p.Results.lineOpacity, ...
                'gain', p.Results.gain, 'invertChannels', ss.channelMap.invertChannelsY, 'normalizeEach', false, ...
                'colorOverlayLabels', overlay_labels_txcxl, 'colorOverlayMap', overlay_colormap, ...
                'channel_ids', channel_ids_ch_by_snippet, 'showChannelDataTips', p.Results.showChannelDataTips, ...
                'showOverlayDataTips', p.Results.showOverlayDataTips, ...
                'colorOverlayDataTipLabel', overlay_datatip_label, 'colorOverlayDataTipValues', overlay_datatip_values);
            
            % plot templates if requested
            if p.Results.overlay_templates
                [templates, template_start_ind, spike_inds, cluster_ids] = ss.buildTemplateOverlays('maskSnippets', mask_snippets, ...
                    'cluster_ids', overlay_cluster_ids, 'sortChannelsVertically', true, ...
                    'best_n_channels', p.Results.overlay_best_channels);
                
                [~, cluster_inds] = ismember(cluster_ids, overlay_cluster_ids);
                
                % generate expanded time vector we can insert templates into accounting for timepoints outside the snippet
                % where a template could start
                nTemplateTime = size(templates, 2);
                time_exp =  int64(ss.window(1) + 1) - int64(nTemplateTime) : int64(ss.window(2) - 1) + int64(nTemplateTime);
                if p.Results.timeInMilliseconds
                    time_exp = single(time_exp) / ss.fs * 1000;
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
                    Neuropixel.Utils.plotStackedTraces(time_this(t_mask), template_this(:, t_mask)', ...
                        'colors', overlay_cluster_colormap(cluster_ind_this, :), ...
                        'transformInfo', transform, ...
                        dataTipArgs{:});
                end
            end
            axis off;
            hold off;
        end
        
        function plotStackedTracesColorClusters(ss, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', ss.valid, @isvector);
            p.addParameter('maxPerCluster', 50, @isscalar);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 1.95, @isscalar);
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.addParameter('dark', false, @islogical);
            p.addParameter('lineOpacity', 1, @isscalar);
            p.parse(varargin{:});
            
            mask_snippets = p.Results.maskSnippets;
            unique_cluster_ids = ss.unique_cluster_ids;
            
            maxPerCluster = p.Results.maxPerCluster;
            if isfinite(maxPerCluster)
                mask_inds = find(mask_snippets);
                cluster_ids = ss.cluster_ids(mask_snippets);
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
            [bg, ~] = ss.setupFigureAxes(dark);
            cluster_colormap = ss.generateOverlayColormap(1:numel(unique_cluster_ids), bg);
            
            % nCh x nSnippets
            channel_ids_ch_by_snippet = ss.channel_ids_by_snippet(:, mask_snippets);
            data = ss.data(:, :, mask_snippets, :);
            data_txcxl = permute(data(:, :, :), [2 1 3]);
            
            cluster_ids = ss.cluster_ids(mask_snippets);
            [~, cluster_inds] = ismember(cluster_ids, unique_cluster_ids);
            traceColors = cluster_colormap(cluster_inds, :);
           
            if p.Results.timeInMilliseconds
                time = ss.time_ms;
            else
                time = ss.time_samples;
            end
            
            Neuropixel.Utils.plotStackedTraces(time, data_txcxl, 'layerColors', traceColors, ...
                'lineWidth', 0.5, 'lineOpacity', p.Results.lineOpacity, ...
                'gain', p.Results.gain, 'invertChannels', ss.channelMap.invertChannelsY, 'normalizeEach', false, ...
                'channel_ids', channel_ids_ch_by_snippet, 'showChannelDataTips', false);
            
            axis off;
            hold off;
        end
        
        function plotHeatmap(ss, snippet_ind, varargin)
            % plots one snippet as a heatmap with each of the clusters that occur during that period
            
            p = inputParser();
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.parse(varargin{:});
            
            if p.Results.timeInMilliseconds
                time = ss.time_ms;
            else
                time = ss.time_samples;
            end
            data = ss.data(:, :, snippet_ind, 1);
            
            Neuropixel.Utils.pmatbal(data, 'x', time);
            yline(0.5, '-', sprintf('snippet %d', snippet_ind), 'LabelVerticalAlignment', 'bottom');
        end
        
        function plotHeatmapWithTemplates(ss, snippet_ind, varargin)
            % plots one snippet as a heatmap with each of the clusters that occur during that period
            
            p = inputParser();
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.addParameter('reconstruct_cluster_ids', ss.overlay_cluster_ids, @isvector);
            p.addParameter('reconstruct_best_channels', 24, @isscalar);
            p.addParameter('successive_residuals', false, @islogical);
            p.parse(varargin{:});
            
            if p.Results.timeInMilliseconds
                time = ss.time_ms;
            else
                time = ss.time_samples;
            end
            data = ss.data(:, :, snippet_ind, 1);
            
            reconstruct_cluster_ids = p.Results.reconstruct_cluster_ids;
            if isempty(reconstruct_cluster_ids) && ~isempty(ss.ks)
                % overlay all clusters if none specified
                reconstruct_cluster_ids = ss.ks.cluster_ids;
            end
            reconstruct_best_channels = p.Results.reconstruct_best_channels;
            
            % nCh x nSnippets
%             cluster_inds_by_snippet = ss.lookup_clusterIds(ss.cluster_ids(snippet_ind));
%             channel_ids_ch_by_snippet = ss.channel_ids_by_cluster(:, cluster_inds_by_snippet);
            
             [templates, cluster_ids] = ss.buildTemplateReconstructions('maskSnippets', snippet_ind, ...
                'cluster_ids', reconstruct_cluster_ids, ...
                'best_n_channels', reconstruct_best_channels, 'sortChannelsVertically', true);
            
            residual = data;
            if p.Results.successive_residuals
                for iT = 1:size(templates, 3)
                    residual = residual - templates(:, :, iT);
                    templates(:, :, iT) = residual;
                end
            end
            templates_stacked = Neuropixel.Utils.TensorUtils.reshapeByConcatenatingDims(templates, {[1 3], 2});
            
            
            data_stacked = cat(1, data, templates_stacked);
            
            Neuropixel.Utils.pmatbal(data_stacked, 'x', time);
            
            yline(0.5, '-', sprintf('snippet %d', snippet_ind), 'LabelVerticalAlignment', 'bottom');
            for iC = 1:size(templates, 3)
                yline(ss.nChannels * iC + 0.5, 'k-', sprintf('cluster %u', cluster_ids(iC)), 'LabelVerticalAlignment', 'bottom');
            end
            
            
        end
        
        function plotHeatmapWithReconstruction(ss, snippet_ind, varargin)
            % plots one snippet as a heatmap with all of the overlapping clusters' reconstructed templates below
            
            default_cluster_reconstruct = ss.ks.cluster_ids;
            hasClusterId = ~isempty(ss.cluster_ids);
            if hasClusterId
                this_cluster_id = ss.cluster_ids(snippet_ind);
                default_cluster_reconstruct = setdiff(default_cluster_reconstruct, this_cluster_id);
            end
            
            p = inputParser();
            p.addParameter('timeInMilliseconds', false, @islogical);
            p.addParameter('reconstruct_cluster_ids', default_cluster_reconstruct, @isvector);
            p.addParameter('reconstruct_best_channels', ss.ks.nChannelsSorted, @isscalar);
            p.parse(varargin{:});
            
            if p.Results.timeInMilliseconds
                time = ss.time_ms;
            else
                time = ss.time_samples;
            end
            data = ss.data(:, :, snippet_ind, 1);
            
            reconstruct_cluster_ids = p.Results.reconstruct_cluster_ids;
            if isempty(reconstruct_cluster_ids) && ~isempty(ss.ks)
                % overlay all clusters if none specified
                reconstruct_cluster_ids = ss.ks.cluster_ids;
            end
            reconstruct_best_channels = p.Results.reconstruct_best_channels;
            
            if hasClusterId
                % build this cluster's tempalte
                this_template = ss.buildTemplateReconstructions('maskSnippets', snippet_ind, ...
                'cluster_ids', this_cluster_id, 'sortChannelsVertically', true, ...
                'best_n_channels', reconstruct_best_channels);
            end
                
            % build other clusters' templates
            [templates, cluster_ids_recon] = ss.buildTemplateReconstructions('maskSnippets', snippet_ind, ...
                'cluster_ids', reconstruct_cluster_ids, 'sortChannelsVertically', true, ...
                'best_n_channels', reconstruct_best_channels);
            reconstruction = sum(templates, 3);
            
            if hasClusterId
                data_stacked = cat(1, data, this_template, reconstruction);
            else
                data_stacked = cat(1, data, reconstruction);
            end
            Neuropixel.Utils.pmatbal(data_stacked, 'x', time);
            
            if ~hasClusterId
                labelMain = sprintf('snippet %d', snippet_ind);
            else
                labelMain = sprintf('snippet %d (cluster %d)', snippet_ind, ss.cluster_ids(snippet_ind));
            end
            yline(0.5, '-', labelMain, 'LabelVerticalAlignment', 'bottom');
            
            nOther = numel(unique(cluster_ids_recon));
            if nOther > 15
                other_clusters_str = sprintf('%d clusters', nOther);
            else
                other_clusters_str = strjoin(string(unique(cluster_ids_recon)), ',');
            end
            if hasClusterId
                yline(ss.nChannels + 0.5, 'k-', sprintf('reconstruction cluster %d', this_cluster_id), 'LabelVerticalAlignment', 'bottom');
                yline(ss.nChannels*2 + 0.5, 'k-', sprintf('reconstruction other clusters (%s)', other_clusters_str), 'LabelVerticalAlignment', 'bottom');
            else
                yline(ss.nChannels + 0.5, 'k-', sprintf('reconstruction all clusters (%s)', other_clusters_str), 'LabelVerticalAlignment', 'bottom');
            end
        end
    end
end