classdef ChannelMap
% Author: Daniel J. O'Shea (2019)

    properties
        file 
        
        channelIdsMapped uint32
        connected logical
        shankInd
        xcoords
        ycoords
        zcoords
        
        syncChannelIndex uint32 % actual index in the AP bin file
        syncChannelId uint32 % arbitrary channel id, typically the same as index
    end
    
    properties(Dependent)
        syncInAPFile
        channelIds
        nChannels
        nChannelsMapped
        connectedChannels
        referenceChannels
        
        yspacing
        xspacing
        zspacing
    end
    
    methods
        function map = ChannelMap(fname)
            if isa(fname, 'Neuropixel.ChannelMap')
                map = fname;
                return;
            end
            
            d = load(fname);
            map.file = fname;
            map.channelIdsMapped = Neuropixel.Utils.makecol(d.chanMap);
            map.connected = Neuropixel.Utils.makecol(d.connected);
            map.shankInd = Neuropixel.Utils.makecol(d.shankInd);
            map.xcoords = Neuropixel.Utils.makecol(d.xcoords);
            map.ycoords = Neuropixel.Utils.makecol(d.ycoords);
            if isfield(d, 'zcoords')
                map.zcoords = Neuropixel.Utils.makecol(d.zcoords);
            else
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
        
        function tf = get.syncInAPFile(map)
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
            sites = map.channelIds(~map.connected);
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
            
    end
end