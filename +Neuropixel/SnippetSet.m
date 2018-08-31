classdef SnippetSet < handle
    % A collection of snippets of data collected from an IMEC data file
    % which includes the subset of channels which were sampled
    
    properties
        data (:, :, :) int16 % channels x time x snippets
        cluster_idx (:, 1) uint32
        channel_idx (:, 1) uint32 % in absolute channel inds
        sample_idx (:, 1) uint64
        trial_idx (:, 1) uint32
        window (:, 2) int64 % in samples
        
        valid
        
        channelMap % ChannelMap for coordinates
        scaleToUv (1, 1) double
        fs
    end
    
    properties(Dependent)
        nChannels
        nTimepoints
        nSnippets
        
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
                elseif isa(ds, 'Neuropixel.ImecDataFile')
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
    end
    
    methods % plotting
        function plotAtProbeLocations(ss, varargin)
            p = inputParser();
            p.addParameter('xmag', 1.5, @isscalar);
            p.addParameter('ymag', 1.5, @isscalar);
            p.addParameter('showInLegendAs', 'snippets', @ischar);
            p.addParameter('showMean', false, @islogical);
            p.addParameter('meanPlotArgs', {'Color', 'r', 'LineWidth', 2}, @iscell);
            p.addParameter('plotArgs', {'Color', [0.3 0.3 0.3]}, @iscell);
            
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('labelArgs', {}, @iscell);
            p.addParameter('alpha', 0.5, @isscalar);
            p.KeepUnmatched = false;
            p.parse(varargin{:});
            
            yspacing = ss.channelMap.yspacing;
            xspacing = ss.channelMap.xspacing;
            xmag = p.Results.xmag;
            ymag = p.Results.ymag;
            xc = ss.channelMap.xcoords(ss.channel_idx);
            yc = ss.channelMap.ycoords(ss.channel_idx);
            
            % plot relative time vector
            tvec = linspace(0, xspacing * xmag, numel(ss.window(1) : ss.window(2)));
            
            % center each channel
            data = double(ss.data_valid);
            data = data - mean(ss.data(:, :), 2); %#ok<*PROPLC>
            data = data - min(data(:));
            data = data ./ max(data(:)) * yspacing * ymag;
            
            holding = ishold;
            handles = cell(ss.nChannels, 1);
            for iC = 1:ss.nChannels
                dmat = squeeze(data(iC, :, :)) + yc(iC);
                handles{iC} = plot(tvec + xc(iC), dmat, p.Results.plotArgs{:});
                if p.Results.alpha < 1
                    Neuropixel.Utils.setLineOpacity(handles{iC}, p.Results.alpha);
                end
                hold on;
                if p.Results.showMean
                    handles{iC}(end+1) = plot(tvec + xc(iC), mean(dmat, 2), p.Results.plotArgs{:}, p.Results.meanPlotArgs{:});
                end
                if p.Results.showLabels
                    text(xc(iC), yc(iC) + yspacing*0.5, sprintf('ch %d', ss.channel_idx(iC)), ...
                        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
                        p.Results.labelArgs{:});
                end
            end
            
            handles = cat(1, handles{:});
            Neuropixel.Utils.showFirstInLegend(handles, p.Results.showInLegendAs);
            
            if ~holding
                hold off;
            end
            axis off;
            box off;
        end
            
    end
end