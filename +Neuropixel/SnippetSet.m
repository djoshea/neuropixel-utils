classdef SnippetSet < handle & matlab.mixin.Copyable
    % A collection of snippets of data collected from an IMEC data file
    % which includes the subset of channels which were sampled
    
    properties
        data (:, :, :) int16 % channels x time x snippets
        cluster_idx (:, 1) int32
        
        % one of these will be filled, with absolute channel inds
        channel_idx_by_cluster (:, :) uint32 % if each cluster draws from different channels
        unique_cluster_idx (:, 1) int32 % specifies the cluster_idx associated with each column here
        
        sample_idx (:, 1) uint64
        trial_idx (:, 1) uint32
        window (:, 2) int64 % in samples
        
        valid (:, 1) logical
        
        channelMap % ChannelMap for coordinates
        scaleToUv (1, 1) double
        fs
    end
    
    properties(Dependent)
        nChannels
        nTimepoints
        nSnippets
        nClusters
        
        data_valid
        time_ms (:, 1) double
    end
    
    methods
        function ss = SnippetSet(ds, type)
            if nargin > 0
                if isa(ds, 'Neuropixel.KiloSortTrialSegmentedDataset')
                    raw_dataset = ds.raw_dataset;
                elseif isa(ds, 'Neuropixel.KiloSortDataset')
                    raw_dataset = ds.raw_dataset;
                elseif isa(ds, 'Neuropixel.ImecDataset')
                    raw_dataset = ds;
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
        
        function n = get.nClusters(ss)
            n = numel(ss.unique_cluster_idx);
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
        
        function ss = selectClusters(ss, cluster_idx)
            mask = ismember(ss.cluster_idx, cluster_idx);
            ss = ss.selectData('maskSnippets', mask);
        end 
        
        function ss = selectData(ss, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', ss.valid, @isvector);
            p.addParameter('maskTime', true(ss.nTimepoints, 1), @isvector);
            p.addParameter('maskChannels', true(ss.nChannels, 1), @isvector);
            p.parse(varargin{:});
            
            mask = p.Results.maskSnippets;
            
            ss = copy(ss);
            ss.valid = ss.valid(mask); % first time auto generates based on data, so has to come first
            ss.data = ss.data(p.Results.maskChannels, p.Results.maskTime, mask);
            
            ss.cluster_idx = ss.cluster_idx(mask);
        
            ss.sample_idx = ss.sample_idx(mask);
            if ~isempty(ss.trial_idx)
                ss.trial_idx = ss.trial_idx(mask);
            end
        end
    end
    
    methods % plotting
        function hdata = plotAtProbeLocations(ss, varargin)
            p = inputParser();
            % specify one of these 
            p.addParameter('cluster_idx', [], @(x) isempty(x) || isvector(x));
            p.addParameter('maskSnippets', ss.valid, @isvector);
            
            % and optionally these
            p.addParameter('maskTime', true(ss.nTimepoints, 1), @isvector);
            p.addParameter('maskChannels', true(ss.nChannels, 1), @isvector);
            
            p.addParameter('maxPerCluster', Inf, @isscalar);
            
            p.addParameter('xmag', 1.5, @isscalar);
            p.addParameter('ymag', 1.5, @isscalar);
           
            p.addParameter('showIndividual', true, @islogical); 
            p.addParameter('showMean', false, @islogical);
            p.addParameter('meanPlotArgs', {'LineWidth', 3}, @iscell);
            p.addParameter('plotArgs', {'LineWidth', 0.5}, @iscell);
            
            p.addParameter('showChannelLabels', true, @islogical);
            p.addParameter('labelArgs', {}, @iscell);
            p.addParameter('alpha', 0.4, @isscalar);
            
            p.addParameter('colormap', get(groot, 'DefaultAxesColorOrder'), @(x) true);
            p.KeepUnmatched = false;
            p.parse(varargin{:});
            
            yspacing = ss.channelMap.yspacing;
            xspacing = ss.channelMap.xspacing;
            xmag = p.Results.xmag;
            ymag = p.Results.ymag;
            
            % plot relative time vector
            tvec = linspace(0, xspacing * xmag, numel(ss.window(1) : ss.window(2)));
            
            if ~isempty(p.Results.cluster_idx)
                maskSnippets = ismember(ss.cluster_idx, p.Results.cluster_idx) & ss.valid;
            else
                maskSnippets = p.Results.maskSnippets;
            end
           
            % center and scale each channel
            if ~any(maskSnippets)
                warning('No snippets for these clusters');
                return;
            end
            
            data = double(ss.data(p.Results.maskChannels, p.Results.maskTime, maskSnippets));
            data = data - mean(data(:, :), 2); %#ok<*PROPLC>
            data = data ./ (max(data(:)) - min(data(:))) * yspacing * ymag;
            
            cluster_idx = ss.cluster_idx(maskSnippets);
            unique_cluster_idx = unique(cluster_idx);
            nUniqueClusters = numel(unique_cluster_idx);
            
            holding = ishold;
            nChannels = size(data, 1);
            handlesIndiv = cell(nChannels, nUniqueClusters);
            handlesMean = gobjects(nChannels, nUniqueClusters);
            channels_plotted = false(ss.channelMap.nChannels, 1);
            
            cmap = p.Results.colormap;
            if isa(cmap, 'function_handle')
                cmap = cmap(nUniqueClusters);
            end
            
            for iClu = 1:nUniqueClusters
                this_cluster_idx = unique_cluster_idx(iClu);
                this_snippet_mask = cluster_idx == this_cluster_idx;
                [~, this_cluster_ind] = ismember(this_cluster_idx, ss.unique_cluster_idx);
                
                if isfinite(p.Results.maxPerCluster)
                    idxKeep = find(this_snippet_mask, p.Results.maxPerCluster, 'first');
                    this_snippet_mask(idxKeep(end)+1:end) = false;
                end
                
                this_channel_idx = ss.channel_idx_by_cluster(:, this_cluster_ind);
                this_channel_idx = this_channel_idx(p.Results.maskChannels);
                xc = ss.channelMap.xcoords(this_channel_idx);
                yc = ss.channelMap.ycoords(this_channel_idx);
                
                channels_plotted(this_channel_idx) = true;

                cmapIdx = mod(iClu-1, size(cmap, 1))+1;
                
                for iC = 1:nChannels
                    dmat = Neuropixel.Utils.TensorUtils.squeezeDims(data(iC, :, this_snippet_mask), 1) + yc(iC);
                        
                    if p.Results.showIndividual
                        handlesIndiv{iC, iClu} = plot(tvec + xc(iC), dmat, 'Color', [cmap(cmapIdx, :), p.Results.alpha], ...
                            p.Results.plotArgs{:});
                        Neuropixel.Utils.hideInLegend(handlesIndiv{iC, iClu});
                        hold on;
                    end
                    if p.Results.showMean
                        handlesMean(iC, iClu) = plot(tvec + xc(iC), mean(dmat, 2), 'Color', cmap(cmapIdx, :), ...
                            p.Results.meanPlotArgs{:});
                        hold on;
                    end
                    if iC == 1 
                        Neuropixel.Utils.showFirstInLegend(handlesMean(iC, iClu), sprintf('cluster %d', this_cluster_idx));
                    else
                        Neuropixel.Utils.hideInLegend(handlesMean(iC, iClu));
                    end
                end
                
                %drawnow;
            end
            
            if p.Results.showChannelLabels
                xc = ss.channelMap.xcoords;
                yc = ss.channelMap.ycoords;
                
                for iC = 1:ss.channelMap.nChannels
                    if channels_plotted(iC)
                        text(xc(iC), yc(iC), sprintf('ch %d', ss.channelMap.chanMap(iC)), ...
                            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            p.Results.labelArgs{:});
                    end
                end
            end
            
            if ~holding
                hold off;
            end
            axis off;
            axis tight;
            box off;
            
            hdata.waveforms = handlesIndiv;
            hdata.waveformMeans = handlesMean;
        end
        
        function plotStacked(ss, varargin)
            p = inputParser();
            p.addParameter('maskSnippets', ss.valid, @isvector);
            p.addParameter('maskTime', true(ss.nTimepoints, 1), @isvector);
            p.addParameter('maskChannels', true(ss.nChannels, 1), @isvector);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            channel_idx = TensorUtils.vectorMaskToIndices(p.Results.maskChannels);
            data = double(ss.data(p.Results.maskChannels, p.Results.maskTime, p.Results.maskSnippets));
            
            chNames = sprintfc('ch %d', channel_idx);
            ptstack(2, 1, ss.time_ms(p.Results.maskTime), data, ...
                'showLabels', true, 'labelsStacked', chNames, p.Unmatched);
            
        end
            
    end
end