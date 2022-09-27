classdef ChannelMap
% Author: Daniel J. O'Shea (2019)

    properties
        file  (1, 1) string
        name (1, 1) string
        
        channelIdsMapped (:, 1) uint32
        connected (:, 1) logical
        shankInd (:, 1)
        
        nSpatialDims = 2;
        xcoords (:, 1)
        ycoords (:, 1)
        zcoords (:, 1)
        
        syncChannelIndex (:, 1) uint32 % actual index in the AP & LF bin file
        syncChannelId (:, 1) uint32 % arbitrary channel id, typically the same as index

        % version 1 uses 32 ADC that sample 12 channels each
        channelsPerADC (1, 1) {mustBeInteger} = 32;
    end
    
    properties(Dependent)
        coords

        adcMap (:, 1) uint32 % which adc does each channel belong to (nChannelsMapped x 1)
        adcSampleShift (:, 1) double % fraction of the sample interval each channel is delayed by (nChannelsMapped x 1)
        
        syncInAPFile
        syncInLFFile
        channelIds
        nChannels
        nChannelsMapped
        connectedChannels
        referenceChannels
        
        invertChannelsY % plot first channel at bottom?
        
        yspacing
        xspacing
        zspacing
        
        xlim
        ylim
        zlim
    end
    
    methods
        function map = ChannelMap(spec)
            if nargin < 1 || isempty(spec)
                return;
            elseif isa(spec, 'Neuropixel.ChannelMap')
                % return the existing instance
                map = fname;
                return;
            elseif isstruct(spec)
                map = Neuropixel.ChannelMap.fromMeta(spec);
            elseif ischar(spec) || isstring(spec)
                if exist(spec, 'file') == 2
                    map = Neuropixel.ChannelMap.fromMatFile(spec);
                else
                    map = Neuropixel.ChannelMap.fromProbeName(spec);
                end
            else
                error('Unknown ChannelMap spec format');
            end
        end
    end
    
    methods(Static)
        function map = fromMatFile(fname)
            map = Neuropixel.ChannelMap();
            
            d = load(fname);
            map.file = fname;
            [~, map.name, ~] = fileparts(fname);
            map.channelIdsMapped = Neuropixel.Utils.makecol(d.chanMap);
            map.connected = Neuropixel.Utils.makecol(d.connected);
            if isfield(d, 'shankInd')
                map.shankInd = Neuropixel.Utils.makecol(d.shankInd);
            else
                map.shankInd = ones(size(map.connected));
            end
            map.xcoords = Neuropixel.Utils.makecol(d.xcoords);
            map.ycoords = Neuropixel.Utils.makecol(d.ycoords);
            if isfield(d, 'zcoords')
                map.zcoords = Neuropixel.Utils.makecol(d.zcoords);
                map.nSpatialDims = 3;
            else
                map.nSpatialDims = 2;
                map.zcoords = zeros(size(map.ycoords), 'like', map.ycoords);
            end
            
            if isfield(d, 'syncChannelIndex')
                map.syncChannelIndex = uint32(d.syncChannelId);
            else
                map.syncChannelIndex = uint32(d.chanMap(end)+ 1);
            end
            if isfield(d, 'syncChannelId')
                map.syncChannelId = uint32(d.syncChannelId);
            else
                map.syncChannelId = map.syncChannelIndex;
            end
        end
        
        function map = fromProbeName(key)
            key = lower(string(key));
            
            switch key
                case "phase3a"
                    fname = "neuropixPhase3A_kilosortChanMap.mat";
                case "phase3a4"
                    fname = "neuropixPhase3A_option4_kilosortChanMap.mat";
                case "phase3b1"
                    fname = "neuropixPhase3B1_kilosortChanMap.mat";
                case "phase3b2"
                    fname = "neuropixPhase3B2_kilosortChanMap.mat";
                otherwise 
                    error('Unknown channel map key or file %s. Valid options are phase3a, phase3a4, phase3b1, phase3b2.', key);
            end
            
            map = Neuropixel.ChannelMap.fromMatFile(fname);
        end
        
%         function map = fromMeta(meta)
%             % this function is from Nick Steinmetz: takes a .meta file from spikeglx and produces a kilosort-compatible channel map.
%             % In order to get the x and y coordinates correct in units of µm it needs a little additional info about the probe geometry
%             % which is hard-coded there for neuropixels 2.0, but that aspect of the function could be straightforwardly improved. 
%             
%             chanMap = [1:nCh(1)]'; 
%             chanMap0ind = chanMap-1;
%             connected = true(size(chanMap)); W
% 
%             shankSep = 250; 
%             rowSep = 15; 
%             colSep = 32;
% 
%             openParen = find(shankMap=='('); 
%             closeParen = find(shankMap==')'); 
%             for c = 1:nCh(1)
%                 thisMap = shankMap(openParen(c+1)+1: closeParen(c+1)-1); 
%                 thisMap(thisMap==':') = ',';
%                 n = str2num(thisMap); 
%                 xcoords(c) = (n(1)-1)*shankSep + (n(2)-1)*colSep; 
%                 ycoords(c) = (n(3)-1)*rowSep; 
%             end
% 
%             cm = struct();
%             cm.chanMap = chanMap; 
%             cm.chanMap0ind = chanMap0ind;
%             cm.xcoords = xcoords'; 
%             cm.ycoords = ycoords'; 
%             cm.connected = connected;
%             [~,name] = fileparts(m.imRoFile); 
%             cm.name = name;
%             
%             snsChanMap = char(snsChanMap);
%             map = Neuropixel.ChannelMap();
%             
%             
%         end
    end
    
    methods
        function adcMap = get.adcMap(map)
            % set adcBankProperty, based on https://github.com/int-brain-lab/ibl-neuropixel/blob/main/src/neuropixel.py 
            %     The sampling is serial within the same ADC, but it happens at the same time in all ADCs.
            %     The ADC to channel mapping is done per odd and even channels:
            %     ADC1: ch1, ch3, ch5, ch7...
            %     ADC2: ch2, ch4, ch6....
            %     ADC3: ch33, ch35, ch37...
            %     ADC4: ch34, ch36, ch38...
            %     Therefore, channels 1, 2, 33, 34 get sample at the same time

            %     version 1 uses 32 ADC that sample 12 channels each
            adc_channels = map.channelsPerADC;
            n_adc = map.nChannelsMapped / adc_channels;
            assert(n_adc == round(n_adc), 'Fix adc mapping for this version');

            cind = (0:map.nChannelsMapped-1)';
            adcMap = floor(cind / (adc_channels * 2)) * 2 + mod(cind, 2) + 1;
        end

        function adcSampleShift = get.adcSampleShift(map)
            adc_channels = map.channelsPerADC;

            adcSampleShift = zeros(map.nChannelsMapped, 1);
            adc = map.adcMap;
            n_adc = max(adc);
            for iA = 1:n_adc
                adcSampleShift(adc == iA) = (0:adc_channels-1)' / adc_channels;
            end
        end

        function tf = get.syncInAPFile(map)
            tf = ~isempty(map.syncChannelIndex);
        end
        
        function tf = get.syncInLFFile(map)
            tf = ~isempty(map.syncChannelIndex);
        end
        
        function zcoords = get.zcoords(map)
            if isempty(map.zcoords)
                zcoords = zeros(size(map.ycoords), 'like', map.ycoords);
            else
                zcoords = map.zcoords;
            end
        end
        
        function v = get.channelIds(map)
            v = map.channelIdsMapped;
            if ~isempty(map.syncChannelIndex)
                v(map.syncChannelIndex) = map.syncChannelId;
            end
        end
        
        function nChannels = get.nChannelsMapped(map)
            nChannels = numel(map.channelIdsMapped);
        end
        
        function nChannels = get.nChannels(map)
            nChannels = numel(map.channelIds);
        end
            
        function sites = get.connectedChannels(map)
            sites = map.channelIds(map.connected);
        end
        
        function sites = get.referenceChannels(map)
            sites = setdiff(map.channelIdsMapped, map.connectedChannels);
        end
        
        function zspacing = get.zspacing(map)
            zs = sort(map.zcoords(map.connected));
            dzs = diff(zs);
            dzs = dzs(dzs > 0);
            zspacing = min(dzs);
        end
        
        function yspacing = get.yspacing(map)
            ys = sort(map.ycoords(map.connected));
            dys = diff(ys);
            dys = dys(dys > 0);
            yspacing = min(dys);
        end
        
        function xspacing = get.xspacing(map)
            xs = sort(map.xcoords(map.connected));
            dxs = diff(xs);
            dxs = dxs(dxs > 0);
            xspacing = min(dxs);
        end
        
        function ylim = get.ylim(map)
            ys = map.yspacing;
            ylim = [min(map.ycoords) - ys, max(map.ycoords) + ys];
        end
        
        function xlim = get.xlim(map)
            xs = map.yspacing;
            xlim = [min(map.xcoords) - xs, max(map.xcoords) + xs];
        end
        
        function zlim = get.zlim(map)
            zs = map.yspacing;
            zlim = [min(map.zcoords) - zs, max(map.zcoords) + zs];
        end
        
        function coords = get.coords(map)
            if map.nSpatialDims == 3
                coords = cat(2, map.xcoords, map.ycoords, map.zcoords);
            else
                coords = cat(2, map.xcoords, map.ycoords);
            end
        end
        
        function tf = get.invertChannelsY(map)
            % bigger y coords are higher up on the probe. This is used when we want to plot stacked traces. If false, 
            % channel 1 belongs at the top (has largest y coord). If true, channel 1 belongs at the bottom (has smallest y coord)
            
            tf = map.ycoords(1) < map.ycoords(end);
        end
        
        function [channelInds, channelIds] = lookup_channelIds(map, channelIds)
             if islogical(channelIds)
                channelIds = map.channelIds(channelIds);
             end
            [tf, channelInds] = ismember(channelIds, map.channelIds);
            assert(all(tf), 'Some channel ids not found');
        end
        
        function closest_ids = getClosestChannels(map, nClosest, channel_ids, eligibleChannelIds)
            % channel_idx is is in 1:nChannels raw data indices
            % closest is numel(channel_idx) x nClosest
            
            if nargin < 3
                channel_ids = map.connectedChannels;
            end
            if nargin < 4
                eligibleChannelIds = map.connectedChannels;
            end
            
            eligibleChannelInds = map.lookup_channelIds(eligibleChannelIds);
                
            x = map.xcoords(eligibleChannelInds);
            y = map.ycoords(eligibleChannelInds);
            z = map.zcoords(eligibleChannelInds);
            N = numel(x);

            X = repmat(x, 1, N);
            Y = repmat(y, 1, N);
            Z = repmat(z, 1, N);

            % distance between all connected and non-connected channels
            distSq = (X - X').^2 + (Y - Y').^2 + (Z - Z').^2;
            distSq(logical(eye(N))) = Inf; % don't localize each channel to itself
            
            [tf, channel_inds] = ismember(channel_ids, eligibleChannelInds);
            assert(all(tf), 'Cannot localize non-connected channels or channels not in eligibleChannelIds');
            
            closest_ids = nan(numel(channel_inds), nClosest);
            for iC = 1:numel(channel_inds)
                [~, idxSort] = sort(distSq(channel_inds(iC), :), 'ascend');
                % idxSort will be in indices into eligible
                closest_ids(iC, :) = eligibleChannelIds(idxSort(1:nClosest));
            end
            
%             function idxFull = indicesIntoMaskToOriginalIndices(idxIntoMasked, mask)
%                 maskInds = find(mask);
%                 idxFull = maskInds(idxIntoMasked);
%             end
        end
        
        function channel_ids_sorted = sortChannelsVertically(map, channel_ids)
            channel_inds = map.lookup_channelIds(channel_ids);
            if size(map.coords, 2) == 3
                [~, sortIdx] = sortrows(map.coords(channel_inds, :), [-2 1 3]); % sort by y (high to low), x (low to high), then z (low to high)
            else
                [~, sortIdx] = sortrows(map.coords(channel_inds, :), [-2 1]); % sort by y (high to low), x (low to high)
            end
            channel_ids_sorted = channel_ids(sortIdx);
        end
        
        function plotRecordingSites(map, varargin)
            p = inputParser();
            p.addParameter('channel_ids', map.channelIdsMapped, @isvector);
            p.addParameter('goodChannels', [], @isvector);
            p.addParameter('badChannels', [], @isvector);
            
            p.addParameter('markerSize', 20, @isscalar);
            p.addParameter('showChannelLabels', false, @islogical);
            p.addParameter('labelArgs', {}, @iscell);
            
            p.addParameter('color_connected', [0.7 0.7 0.7], @(x) true);
            p.addParameter('color_bad', [1 0.7 0.7], @(x) true);
            p.addParameter('color_good', [0.7 1 0.7], @(x) true);
            p.addParameter('color_reference', [0.5 0.5 1], @(x) true);
            p.parse(varargin{:});
            
            [channelInds, channelIds] = map.lookup_channelIds(p.Results.channel_ids);
            xc = map.xcoords(channelInds);
            yc = map.ycoords(channelInds);
            
            is_connected = ismember(channelIds, map.connectedChannels);
            is_good = ismember(channelIds, p.Results.goodChannels);
            is_bad = ismember(channelIds, p.Results.badChannels);
            is_ref = ismember(channelIds, map.referenceChannels);
            cmap = zeros(numel(channelInds), 3);
            cmap(is_connected, :) = repmat(p.Results.color_connected, nnz(is_connected), 1);
            cmap(is_good, :) = repmat(p.Results.color_good, nnz(is_good), 1);
            cmap(is_bad, :) = repmat(p.Results.color_bad, nnz(is_bad), 1);
            cmap(is_ref, :) = repmat(p.Results.color_reference, nnz(is_ref), 1);
            
            channelTypes = strings(numel(channelInds), 1);
            channelTypes(is_connected) = "Connected";
            channelTypes(is_good) = "Good";
            channelTypes(is_bad) = "Bad";
            channelTypes(is_ref) = "Reference";
            
            h = scatter(xc, yc, p.Results.markerSize, cmap, 's', 'filled');
            h.DataTipTemplate.DataTipRows(1).Format = '%g um';
            h.DataTipTemplate.DataTipRows(2).Format = '%g um';
            h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('channel id', double(channelIds));
            h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('type', channelTypes);

            if p.Results.showChannelLabels
                for iC = 1:numel(channelInds)
                    text(xc(iC), yc(iC), sprintf('ch %d', m.channel_ids(channelInds(iC))), ...
                            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                            'Background', 'none', ...
                            p.Results.labelArgs{:});
                end
            end
        end
            
    end
end