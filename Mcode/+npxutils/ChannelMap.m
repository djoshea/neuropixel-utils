classdef ChannelMap
    
    properties (Constant, Hidden)
        knownProbeNames = ["phase3a", "phase3a4", "phase3b1", "phase3b2"];        
    end
    
    properties
        file (1,1) string
        name (1,1) string
        
        channelIdsMapped (:,1) uint32
        connected (:,1) logical
        shankInd (:,1)
        
        nSpatialDims (1,1) double = 2;
        xcoords (:,1)
        ycoords (:,1)
        zcoords (:,1)
        
        % Actual index in the AP & LF bin file
        syncChannelIndex (:,1) uint32
        % Arbitrary channel id, typically the same as index
        syncChannelId (:,1) uint32
    end
    
    properties(Dependent)
        coords
        
        syncInAPFile
        syncInLFFile
        channelIds
        nChannels
        nChannelsMapped
        connectedChannels
        referenceChannels
        
        % Plot first channel at bottom?
        invertChannelsY
        
        yspacing
        xspacing
        zspacing
        
        xlim
        ylim
        zlim
    end
    
    methods
        function this = ChannelMap(spec)
            % Construct a new ChannelMap object
            %
            % map = npxutils.ChannelMap
            % map = npxutils.ChannelMap(spec)
            %
            % Spec may be:
            %   - A charvec or string
            %   - A struct
            %   - An npxutils.ChannelMap object
            %   - An empty of any type
            if nargin == 0 || isempty(spec)
                return;
            elseif isa(spec, 'npxutils.ChannelMap')
                % return the existing instance
                this = spec;
                return;
            elseif isstruct(spec)
                this = npxutils.ChannelMap.fromMeta(spec);
            elseif ischar(spec) || isstring(spec)
                if exist(spec, 'file') == 2
                    this = npxutils.ChannelMap.fromMatFile(spec);
                else
                    this = npxutils.ChannelMap.fromProbeName(spec);
                end
            else
                error('Unknown ChannelMap spec input format: %s', class(spec));
            end
        end
    end
    
    methods(Static)
        function map = fromMatFile(fname)
            % Read a ChannelMap from a MAT file.
            %
            % map = fromMatFile(fname)
            %
            % Reads a MAT file and constructs a ChannelMap from it.
            %
            % It looks like this assumes that you have added a directory with
            % probe map files to your Matlab path.
            %
            % Fname (string) is the path to a MAT file.
            %
            % The MAT file has to be in a particular format, which is not
            % documented here.
            %
            % Returns a scalar npxutils.ChannelMap object.
            map = npxutils.ChannelMap();
            
            d = load(fname);
            map.file = fname;
            [~, map.name, ~] = fileparts(fname);
            map.channelIdsMapped = npxutils.internal.makecol(d.chanMap);
            map.connected = npxutils.internal.makecol(d.connected);
            map.shankInd = npxutils.internal.makecol(d.shankInd);
            map.xcoords = npxutils.internal.makecol(d.xcoords);
            map.ycoords = npxutils.internal.makecol(d.ycoords);
            if isfield(d, 'zcoords')
                map.zcoords = npxutils.internal.makecol(d.zcoords);
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
            % Get the ChannelMap for a given probe name.
            %
            % map = npxutils.Channel.amp.fromProbeName(key)
            %
            % Key (string) is the name of a probe. Case insensitive.
            %
            % Looks like this needs you to have added a directory that contains
            % probe channel map files to your Matlab path before you call this.
            %
            % Returns a scalar ChannelMap.
            arguments
                key (1,1) string
            end
            
            switch lower(key)
                case "phase3a"
                    fname = "neuropixPhase3A_kilosortChanMap.mat";
                case "phase3a4"
                    fname = "neuropixPhase3A_option4_kilosortChanMap.mat";
                case "phase3b1"
                    fname = "neuropixPhase3B1_kilosortChanMap.mat";
                case "phase3b2"
                    fname = "neuropixPhase3B2_kilosortChanMap.mat";
                otherwise
                    error('Unknown probe name for channel map: %s. Valid options are: %s', ...
                        key, strjoin(npxutils.ChannelMap.knownProbeNames, ', '));
            end
            
            map = npxutils.ChannelMap.fromMatFile(fname);
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
        %             map = npxutils.ChannelMap();
        %
        %
        %         end
    end
    
    methods
        function tf = get.syncInAPFile(this)
            tf = ~isempty(this.syncChannelIndex);
        end
        
        function tf = get.syncInLFFile(this)
            tf = ~isempty(this.syncChannelIndex);
        end
        
        function zcoords = get.zcoords(this)
            if isempty(this.zcoords)
                zcoords = zeros(size(this.ycoords), 'like', this.ycoords);
            else
                zcoords = this.zcoords;
            end
        end
        
        function v = get.channelIds(this)
            v = this.channelIdsMapped;
            if ~isempty(this.syncChannelIndex)
                v(this.syncChannelIndex) = this.syncChannelId;
            end
        end
        
        function nChannels = get.nChannelsMapped(this)
            nChannels = numel(this.channelIdsMapped);
        end
        
        function nChannels = get.nChannels(this)
            nChannels = numel(this.channelIds);
        end
        
        function sites = get.connectedChannels(this)
            sites = this.channelIds(this.connected);
        end
        
        function sites = get.referenceChannels(this)
            sites = setdiff(this.channelIdsMapped, this.connectedChannels);
        end
        
        function zspacing = get.zspacing(this)
            zs = sort(this.zcoords(this.connected));
            dzs = diff(zs);
            dzs = dzs(dzs > 0);
            zspacing = min(dzs);
        end
        
        function yspacing = get.yspacing(this)
            ys = sort(this.ycoords(this.connected));
            dys = diff(ys);
            dys = dys(dys > 0);
            yspacing = min(dys);
        end
        
        function xspacing = get.xspacing(this)
            xs = sort(this.xcoords(this.connected));
            dxs = diff(xs);
            dxs = dxs(dxs > 0);
            xspacing = min(dxs);
        end
        
        function ylim = get.ylim(this)
            ys = this.yspacing;
            ylim = [min(this.ycoords) - ys, max(this.ycoords) + ys];
        end
        
        function xlim = get.xlim(this)
            xs = this.yspacing;
            xlim = [min(this.xcoords) - xs, max(this.xcoords) + xs];
        end
        
        function zlim = get.zlim(this)
            zs = this.yspacing;
            zlim = [min(this.zcoords) - zs, maz(this.zcoords) + zs];
        end
        
        function coords = get.coords(this)
            if this.nSpatialDims == 3
                coords = cat(2, this.xcoords, this.ycoords, this.zcoords);
            else
                coords = cat(2, this.xcoords, this.ycoords);
            end
        end
        
        function tf = get.invertChannelsY(this)
            % Get invertChannels
            
            % bigger y coords are higher up on the probe. This is used when we
            % want to plot stacked traces. If false, channel 1 belongs at the
            % top (has largest y coord). If true, channel 1 belongs at the
            % bottom (has smallest y coord)
            
            tf = this.ycoords(1) < this.ycoords(end);
        end
        
        function [channelInds, channelIds] = lookup_channelIds(this, channelIds)
            if islogical(channelIds)
                channelIds = this.channelIds(channelIds);
            end
            [tf, channelInds] = ismember(channelIds, this.channelIds);
            assert(all(tf), 'Some channel ids not found');
        end
        
        function closest_ids = getClosestChannels(this, nClosest, channelIds, eligibleChannelIds)
            % Get channels that are closest to specified channels.
            %
            % closest_ids = getClosestChannels(obj, nClosest, channelIds, eligibleChannelIds)
            %
            % ChannelIds is a list of channels. Defaults to
            % this.connectedChannels.
            %
            % EligibleChannelIds defaults to this.connectedChannels.
            
            % channel_idx is is in 1:nChannels raw data indices
            % closest is numel(channel_idx) x nClosest
            
            if nargin < 3
                channelIds = this.connectedChannels;
            end
            if nargin < 4
                eligibleChannelIds = this.connectedChannels;
            end
            
            eligibleChannelInds = this.lookup_channelIds(eligibleChannelIds);
            
            x = this.xcoords(eligibleChannelInds);
            y = this.ycoords(eligibleChannelInds);
            z = this.zcoords(eligibleChannelInds);
            N = numel(x);
            
            X = repmat(x, 1, N);
            Y = repmat(y, 1, N);
            Z = repmat(z, 1, N);
            
            % distance between all connected and non-connected channels
            distSq = (X - X').^2 + (Y - Y').^2 + (Z - Z').^2;
            distSq(logical(eye(N))) = Inf; % don't localize each channel to itself
            
            [tf, channel_inds] = ismember(channelIds, eligibleChannelInds);
            if ~all(tf)
                badChannelIds = unique(channelIds(~tf));
                error(['Cannot localize non-connected channels or channels not ' ...
                    'in eligibleChannelIds. Bad channel IDs: %s'], mat2str(badChannelIds));
            end
            
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
        
        function channel_ids_sorted = sortChannelsVertically(this, channel_ids)
            % Get this' channels, sorted vertically (by y, then x, then z
            %
            %
            channel_inds = this.lookup_channelIds(channel_ids);
            if size(this.coords, 2) == 3
                % sort by y (high to low), x (low to high), then z (low to high)
                [~, sortIdx] = sortrows(this.coords(channel_inds, :), [-2 1 3]);
            else
                % sort by y (high to low), x (low to high)
                [~, sortIdx] = sortrows(this.coords(channel_inds, :), [-2 1]);
            end
            channel_ids_sorted = channel_ids(sortIdx);
        end
        
        function plotRecordingSites(this, varargin)
            p = inputParser();
            p.addParameter('channel_ids', this.channelIdsMapped, @isvector);
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
            
            [channelInds, channelIds] = this.lookup_channelIds(p.Results.channel_ids);
            xc = this.xcoords(channelInds);
            yc = this.ycoords(channelInds);
            
            is_connected = ismember(channelIds, this.connectedChannels);
            is_good = ismember(channelIds, p.Results.goodChannels);
            is_bad = ismember(channelIds, p.Results.badChannels);
            is_ref = ismember(channelIds, this.referenceChannels);
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