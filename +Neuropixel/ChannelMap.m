classdef ChannelMap
    properties
        file 
        
        chanMap
        chanMap0ind
        connected
        shankInd
        xcoords
        ycoords
        zcoords
    end
    
    properties(Dependent)
        nChannels
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
            map.chanMap = d.chanMap;
            map.chanMap0ind = d.chanMap0ind;
            map.connected = d.connected;
            map.shankInd = d.shankInd;
            map.xcoords = d.xcoords;
            map.ycoords = d.ycoords;
            if isfield(d, 'zcoords')
                map.zcoords = d.zcoords;
            else
                map.zcoords = zeros(size(map.ycoords), 'like', map.ycoords);
            end
        end
        
        function zcoords = get.zcoords(map)
            if isempty(map.zcoords)
                zcoords = zeros(size(map.ycoords), 'like', map.ycoords);
            else
                zcoords = map.zcoords;
            end
        end
        
        function nChannels = get.nChannels(map)
            nChannels = numel(map.chanMap);
        end
            
        function sites = get.connectedChannels(map)
            sites = map.chanMap(map.connected);
        end
        
        function sites = get.referenceChannels(map)
            sites = map.chanMap(~map.connected);
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
        
        function closest_idx = getClosestConnectedChannels(map, nClosest, channel_idx)
            % channel_idx is is in 1:nChannels raw data indices
            % closest is numel(channel_idx) x nClosest
            
            if nargin < 3
                channel_idx = map.connectedChannels;
            end
            
            x = map.xcoords(map.connected);
            y = map.ycoords(map.connected);
            z = map.zcoords(map.connected);
            N = numel(x);

            X = repmat(x, 1, N);
            Y = repmat(y, 1, N);
            Z = repmat(z, 1, N);

            % distance between all connected and non-connected channels
            distSq = (X - X').^2 + (Y - Y').^2 + (Z - Z').^2;
            distSq(logical(eye(N))) = Inf; % don't localize each channel to itself
            
            [tf, idx] = ismember(channel_idx, map.connectedChannels);
            assert(all(tf), 'Cannot localize non-connected channels');
            
            closest_idx = nan(numel(idx), nClosest);
            for iC = 1:numel(idx)
                [~, idxSort] = sort(distSq(iC, :), 'ascend');
                % idxSort will be in indices into connected
                closest_idx(iC, :) = indicesIntoMaskToOriginalIndices(idxSort(1:nClosest), map.connected);
            end
            
            function idxFull = indicesIntoMaskToOriginalIndices(idxIntoMasked, mask)
                maskInds = find(mask);
                idxFull = maskInds(idxIntoMasked);
            end
        end
            
    end
end