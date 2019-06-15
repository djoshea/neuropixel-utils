classdef ImecDataset < handle
% Author: Daniel J. O'Shea (2019)

    properties(SetAccess = protected)
        pathRoot char = '';
        fileStem char = '';
        creationTime = NaN;
        nChannels = NaN;

        fileTypeAP = 'ap'; % typically ap or ap_CAR

        nSamplesAP = 0;
        nSamplesLF = 0;
        fsAP = NaN; % samples_per_second
        fsLF = NaN; % samples_per_second
        highPassFilterHz = NaN;
        apGain = NaN;
        apRange = [];
        lfGain = NaN;
        lfRange = []

        adcBits = 10;

        channelMap = []; % can be stored using set channelMap

        % see markBadChannels
        badChannels

        syncBitNames string;
        
        concatenationInfoAP
        concatenationInfoLF
    end

    properties
        % will be cached after loading, can also be cleared by user
        syncRaw int16 = [];
    end

    properties(Constant)
        bytesPerSample = 2;
    end

    properties(Dependent)
        hasAP
        hasLF

        channelMapFile
        mappedChannels
        mappedChannelInds
        nChannelsMapped % number of channels in the channel map (excludes sync)

        connectedChannels
        connectedChannelInds
        nChannelsConnected % excludes reference and sync channels

        goodChannels % connected channels sans badChannels
        goodChannelInds
        nGoodChannels
        
        channelIds % list of ids from ChannelMap
        channelNames % full list of channel names
        channelNamesPadded

        nSyncBits
        syncBitsNamed

        fileAP % .imec.ap.bin file without folder
        pathAP % .imec.ap.bin file with folder
        fileAPMeta
        pathAPMeta

        fileLF % without folder
        pathLF % .imec.lf.bin file with folder
        fileLFMeta
        pathLFMeta

        % if sync is stored in a separate file than AP
        fileSync
        pathSync % .imec.sync.bin file

        % after sync is cached to sync.mat file for faster reload
        fileSyncCached
        pathSyncCached % .sync.mat file (with cached sync)

        creationTimeStr

        apScaleToUv % multiply raw int16 by this to get uV
        lfScaleToUv

        % from channel map (although syncChannelIndex will be 1 if sync not in AP file)
        syncChannelIndex % if sync in AP file, at one index
        syncInAPFile % is the sync info in the ap file, or in a separate .sync file
    end

    methods
        function df = ImecDataset(fileOrFileStem, varargin)
            p = inputParser();
            p.addParameter('channelMap', [], @(x) true);
            p.addParameter('syncBitNames', [], @(x) isempty(x) || isstring(x) || iscellstr(x));
            p.parse(varargin{:})

            fileOrFileStem = char(fileOrFileStem);
            file = Neuropixel.ImecDataset.findImecFileInDir(fileOrFileStem, 'ap');
            if isempty(file)
                error('No AP Imec file found at or in %s', fileOrFileStem);
            end
            [df.pathRoot, df.fileStem, df.fileTypeAP] = Neuropixel.ImecDataset.parseImecFileName(file);
            if exist(df.pathAP, 'file')
                if ~exist(df.pathAPMeta, 'file')
                    error('Could not find AP meta file %s', df.pathAPMeta);
                end
                df.readInfo();
            else
                error('Could not find AP bin file %s', df.pathAP);
            end

            channelMapFile = p.Results.channelMap;
            if isempty(channelMapFile)
                channelMapFile = Neuropixel.Utils.getDefaultChannelMapFile(true);
            end
            df.channelMap = Neuropixel.ChannelMap(channelMapFile);
            assert(df.channelMap.nChannels <= df.nChannels, 'Channel count is less than number of channels in channel map');

            if ~isempty(p.Results.syncBitNames)
                df.setSyncBitNames(1:numel(p.Results.syncBitNames), p.Resuls.syncBitNames);
            end
        end

        function readInfo(df)
            meta = df.readAPMeta();
            df.nChannels = meta.nSavedChans;
            df.fsAP = meta.imSampRate;
            df.highPassFilterHz = meta.imHpFlt;
            df.creationTime = datenum(meta.fileCreateTime, 'yyyy-mm-ddTHH:MM:SS');

            if df.hasLF
                metaLF = df.readLFMeta();
                df.fsLF = metaLF.imSampRate;
            end

            % parse imroTable
            m = regexp(meta.imroTbl, '\(([\d, ]*)\)', 'tokens');
            gainVals = strsplit(m{2}{1}, ' ');
            df.apGain = str2double(gainVals{4});
            df.lfGain = str2double(gainVals{5});

            df.apRange = [meta.imAiRangeMin meta.imAiRangeMax];
            df.lfRange = [meta.imAiRangeMin meta.imAiRangeMax];

            % look at AP meta fields that might have been set by us
            if isfield(meta, 'badChannels')
                df.badChannels = union(df.badChannels, meta.badChannels);
            end
            if isfield(meta, 'syncBitNames')
                df.setSyncBitNames(1:numel(meta.syncBitNames), meta.syncBitNames);
            end

            if df.hasAP
                fid = df.openAPFile();
                fseek(fid, 0, 'eof');
                bytes = ftell(fid);
                fclose(fid);

                df.nSamplesAP = bytes / df.bytesPerSample / df.nChannels;
                assert(round(df.nSamplesAP) == df.nSamplesAP, 'AP bin file size is not an integral number of samples');
            end
            
            df.concatenationInfoAP = Neuropixel.ConcatenationInfo(df, 'ap', meta);

            if df.hasLF
                fid = df.openLFFile();
                fseek(fid, 0, 'eof');
                bytes = ftell(fid);
                fclose(fid);
                df.nSamplesLF = bytes / df.bytesPerSample / df.nChannels;
                assert(round(df.nSamplesAP) == df.nSamplesAP, 'LF bin file size is not an integral number of samples');
                
                df.concatenationInfoLF = Neuropixel.ConcatenationInfo(df, 'lf', meta);
            end
        end

        function setSyncBitNames(df, idx, names)
            % idx is the indices of which bits to set to the corresponding items from names
            assert(all(idx >= 1 & idx <= df.nSyncBits), 'Sync bit indices must be in [1 %d]', df.nSyncBits);
            if isscalar(idx) && ischar(names)
                df.syncBitNames{idx} = names;
            else
                names = string(names);
                df.syncBitNames(idx) = names;
            end
        end

        function idx = lookupSyncBitByName(df, names, ignoreNotFound)
            if nargin < 3
                ignoreNotFound = false;
            end
            if isnumeric(names)
                idx = names;
            else
                names = string(names);
                [tf, idx] = ismember(names, df.syncBitNames);
                if ignoreNotFound
                    idx(~tf) = NaN;
                elseif any(~tf)
                    error('Sync bit(s) %s not found', strjoin(names, ', '));
                end
            end
        end

        function newImec = copyToNewLocation(df, newRoot, newStem)
            if nargin < 3
                newStem = df.fileStem;
            end
            Neuropixel.Utils.mkdirRecursive(newRoot);

            f = @(suffix) fullfile(newRoot, [newStem suffix]);
            docopy(df.pathAP, f('.imec.ap.bin'));
            docopy(df.pathAPMeta, f('.imec.ap.meta'));
            docopy(df.pathLF, f('.imec.lf.bin'));
            docopy(df.pathLFMeta, f('.imec.lf.meta'));
            docopy(df.pathSync, f('.imec.sync.bin'));

            newImec = Neuropixel.ImecDataset(fullfile(newRoot, newStem), 'channelMap', df.channelMapFile);

            function docopy(from, to)
                if ~exist(from, 'file')
                    return;
                end
                fprintf('Copying to %s\n', to);
                [success, message, ~] = copyfile(from, to);
                if ~success
                    error('Error writing %s: %s', to, message);
                end
            end

        end
    end

    methods  % these functions read a contiguous block of samples over a contiguous band of channels
%         function data_ch_by_time = readAPChannelBand(df, chFirst, chLast, sampleFirst, sampleLast, msg)
%             if nargin < 4 || isempty(sampleFirst)
%                 sampleFirst = 1;
%             end
%             if nargin < 5 || isempty(sampleLast)
%                 sampleLast = df.nSamplesAP;
%             end
%             if nargin < 6 || isempty(msg)
%                 msg = 'Reading channels from neuropixel AP file';
%             end
% 
%             data_ch_by_time = df.readChannelBand('ap', chFirst, chLast, sampleFirst, sampleLast, msg);
%         end
% 
%         function data_ch_by_time = readLFChannelBand(df, chFirst, chLast, sampleFirst, sampleLast, msg)
%             if nargin < 4 || isempty(sampleFirst)
%                 sampleFirst = 1;
%             end
%             if nargin < 5 || isempty(sampleLast)
%                 sampleLast = df.nSamplesLF;
%             end
%             if nargin < 6 || isempty(msg)
%                 msg = 'Reading channels from neuropixel LF file';
%             end
% 
%             data_ch_by_time = df.readChannelBand('lf', chFirst, chLast, sampleFirst, sampleLast, msg);
%         end
% 
%         function data_by_time = readAPSingleChannel(df, ch, varargin)
%             data_by_time = df.readAPChannelBand(ch, ch, varargin{:})';
%         end
% 
%         function data_by_time = readLFSingleChannel(df, ch, varargin)
%             data_by_time = df.readLFChannelBand(ch, ch, varargin{:})';
%         end
    end
    
    methods % Sync channel read / cache
        function syncRaw = readSync(df, varargin)
            p = inputParser();
            p.addOptional('reload', false, @islogical); % if true, ignore cache in memory df.syncRaw 
            p.addParameter('ignoreCached', false, @islogical); % if true, ignore cache on disk in df.pathSyncCached
            p.parse(varargin{:});

            if isempty(df.syncRaw) || p.Results.reload
                if exist(df.pathSyncCached, 'file') && ~p.Results.ignoreCached
                    [~, f, e] = fileparts(df.pathSyncCached);
                    fprintf('Loading sync from cached %s%s\n', f, e);
                    ld = load(df.pathSyncCached);
                    df.syncRaw = ld.sync;
                else
                    % this will automatically redirect to a separate sync file
                    % or to the ap file depending on .syncInAPFile
                    fprintf('Loading sync channel (this will take some time)...\n');
                    mm = df.memmapSync_full();
                    df.syncRaw = mm.Data.x(df.syncChannelIndex, :)';

                    df.saveSyncCached();
                end
            end
            syncRaw = df.syncRaw;
        end

        function saveSyncCached(df)
            sync = df.readSync();
            save(df.pathSyncCached, 'sync');
        end

        function updateSyncCached(df, varargin)
            if exist(df.pathSyncCached, 'file')
                sync = df.readSync('reload', true, 'ignoreCached', true);
                save(df.pathSyncCached, 'sync');
            end
        end

        function tf = readSyncBit(df, bit)
            bit = df.lookupSyncBitByName(bit);
            tf = logical(bitget(df.readSync(), bit));
        end
        
        function vec = readSync_idx(df, idx)
            if ~isempty(df.syncRaw)
                vec = df.syncRaw(idx);
            else
                mm = df.memmapSync_full();
                vec = mm.Data.x(df.syncChannelIndex, idx)';
            end
        end
        
        function mat = readSyncBits_idx(df, bits, idx)
            % mat is nBits x nTime to match readAP_idx which is nChannels x nTime
            if isstring(bits) || ischar(bits)
                bits = df.lookupSyncBitByName(bits);
            end
            vec = df.readSync_idx(idx);
            mat = false(numel(bits), numel(vec));
            for iB = 1:numel(bits)
                mat(iB, :) = logical(bitget(vec, bits(iB)));
            end
        end
    end
    
    methods
        function sampleIdx = closestSampleAPForTime(df, timeSeconds)
            sampleIdx = round(timeSeconds * df.fsAP);
            sampleIdx(sampleIdx == 0) = 1;
            if any(sampleIdx < 0 | sampleIdx > df.nSamplesAP)
                error('Time seconds out of range');
            end 
        end
        
        function sampleIdx = closestSampleLFForTime(df, timeSeconds)
            sampleIdx = round(timeSeconds * df.fsLF);
            sampleIdx(sampleIdx == 0) = 1;
            if any(sampleIdx < 0 | sampleIdx > df.nSamplesLF)
                error('Time seconds out of range');
            end 
        end
        
        function data = readAP_idx(df, sampleIdx, varargin)
            p = inputParser();
            p.addParameter('applyScaling', true, @islogical); % convert to uV before processing
            p.parse(varargin{:});
            
            mm = df.memmapAP_full();
            
            if any(isnan(sampleIdx))
                mask = ~isnan(sampleIdx);
                data = single(mm.Data.x(:, sampleIdx)); % must be single to support NaNs
                data = Neuropixel.Utils.TensorUtils.inflateMaskedTensor(data, 2, mask, NaN);
            else
                data = mm.Data.x(:, sampleIdx);
            end

            if p.Results.applyScaling
                data = single(data);
                ch_conn_mask = df.lookup_channelIds(df.connectedChannels);
                data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(df.apScaleToUv);
            end
        end
       
        function [mat, sampleIdx] = readAP_timeWindow(df, timeWindowSec, varargin)
            idxWindow = df.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = df.readAP_idx(sampleIdx, varargin{:});
        end
        
        function [mat, sampleIdx] = readSyncBits_timeWindow(df, bits, timeWindowSec)
            idxWindow = df.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = readSyncBits_idx(bits, sampleIdx);
        end
        
        function [mat, sourceInds] = readAP_viaTimeShiftSpec(df, timeShifts, varargin)
            % read a matrix of AP data as requested by Neuropixel.TimeShiftSpec instance timeShifts
            sourceInds = timeShifts.computeSourceIndices(df.nSamplesAP);
            mat = df.readAP_idx(sourceInds, varargin{:});
        end
    end
    
    methods % Quick inspection
        function [channelInds, channelIds] = lookup_channelIds(df, channelIds)
             if islogical(channelIds)
                channelIds = df.channelIdx(channelIds);
             end
            [tf, channelInds] = ismember(channelIds, df.channelIds);
            assert(all(tf), 'Some channel ids not found');
        end
        
        function inspectAP_timeWindow(df, timeWindowSec, varargin)
            idxWindow = df.closestSampleAPForTime(timeWindowSec);
            df.inspectAP_idxWindow(idxWindow, 'timeInSeconds', true, varargin{:});
        end
        
        function inspectAP_idxWindow(df, idxWindow, varargin)
            p = inputParser();
            p.addParameter('channels', df.mappedChannels, @isvector);
            p.addParameter('invertChannels', false, @islogical);
            p.addParameter('showSync', true, @isvector);
            p.addParameter('syncBits', df.syncBitsNamed, @isvector);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 0.95, @isscalar);
            p.addParameter('car', false, @islogical);
            p.addParameter('downsample',1, @isscalar); 
            p.addParameter('timeInSeconds', false, @islogical);
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'Neuropixel.TrialSegmentationInfo')); % to mark trial boundaries
            p.parse(varargin{:});
            
            if numel(idxWindow) > 2
                idxWindow = [idxWindow(1), idxWindow(end)];
            end
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = df.readAP_idx(sampleIdx); % C x T
            labels = df.channelNamesPadded;
            
            [channelInds, channelIds] = df.lookup_channelIds(p.Results.channels); %#ok<*PROPLC>
            mat = mat(channelInds, :);
            labels = labels(channelInds);
            connected = ismember(channelIds, df.connectedChannels);
            bad = ismember(channelIds, df.badChannels);
            
            if p.Results.downsample > 1
                mat = mat(:, 1:p.Results.downsample:end);
                sampleIdx = sampleIdx(1:p.Results.downsample:end);
            end
            mat = double(mat);
            if p.Results.car
                mat = mat - median(mat, 1);
            end
            
            colors = zeros(size(mat, 1), 3);
            colors(~connected, 3) = 1; % not connected --> blue
            colors(bad & connected, 1) = 1; % bad --> red
            colors(bad & connected, 3) = 0; % bad --> red
            normalizeMask = true(size(mat, 1), 1);
            
            % append sync bit info to plot in purple
            syncBits = p.Results.syncBits;
            if ~isempty(syncBits) && p.Results.showSync
                syncBitMat = df.readSyncBits_idx(syncBits, sampleIdx);
                mat = cat(1, mat, syncBitMat);
                syncColor = [0.75 0 0.9];
                colors = cat(1, colors, repmat(syncColor, size(syncBitMat, 1), 1));
                labels = cat(1, labels, df.syncBitNames(syncBits));
                normalizeMask = cat(1, normalizeMask, false(size(syncBitMat, 1), 1));
            end

            if ~p.Results.showLabels
                labels = [];
            end
            
            if p.Results.timeInSeconds
                time = double(sampleIdx) / df.fsAP;
            else
                time = sampleIdx;
            end
            Neuropixel.Utils.plotStackedTraces(time, mat', 'colors', colors, 'labels', labels, ...
                'gain', p.Results.gain, 'invertChannels', p.Results.invertChannels, 'normalizeMask', normalizeMask, 'normalizeEach', false);
            
            tsi = p.Results.tsi;
            if ~isempty(tsi)
                tsi.markTrialTicks('sample_window', idxWindow, 'timeInSeconds', p.Results.timeInSeconds, ...
                    'side', 'bottom', 'Color', [0.2 0.2 1], 'expand_limits', true);
            end
        end
    end

    methods % Memory mapped read/write access to data
        function mm = memmapAP_by_sample(df)
            mm = memmapfile(df.pathAP, 'Format', {'int16', [df.nChannels 1], 'x'}, ...
               'Repeat', df.nSamplesAP);
        end

        function mm = memmapLF_by_sample(df)
            mm = memmapfile(df.pathLF, 'Format', {'int16', [df.nChannels 1], 'x'}, ...
               'Repeat', df.nSamplesLF);
        end

        function mm = memmapAP_by_chunk(df, nSamplesPerChunk)
            mm = memmapfile(df.pathAP, 'Format', {'int16', [df.nChannels nSamplesPerChunk], 'x'}, ...
               'Repeat', floor(df.nSamplesAP/nSamplesPerChunk));
        end

        function mm = memmapLF_by_chunk(df, nSamplesPerChunk)
            mm = memmapfile(df.pathLF, 'Format', {'int16', [df.nChannels nSamplesPerChunk], 'x'}, ...
               'Repeat', floor(df.nSamplesLF/nSamplesPerChunk));
        end

        function mm = memmapAP_full(df, varargin)
            p = inputParser();
            p.addParameter('Writable', false, @islogical);
            p.parse(varargin{:});

            mm = memmapfile(df.pathAP, 'Format', {'int16', [df.nChannels df.nSamplesAP], 'x'}, 'Writable', p.Results.Writable);
        end

        function mm = memmapLF_full(df, varargin)
            p = inputParser();
            p.addParameter('Writable', false, @islogical);
            p.parse(varargin{:});

            mm = memmapfile(df.pathLF, 'Format', {'int16', [df.nChannels df.nSamplesLF], 'x'}, 'Writable', p.Results.Writable);
        end

        function mm = memmapSync_full(df)
            if df.syncInAPFile
                % still has nChannels
                mm = memmapfile(df.pathSync, 'Format', {'int16', [df.nChannels df.nSamplesAP], 'x'});
            else
                % only sync channel
                mm = memmapfile(df.pathSync, 'Format', {'int16', [1 df.nSamplesAP], 'x'});
            end
        end
    end

    methods(Hidden) % Read data at specified times
        function [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet] = readSnippetsRaw(df, times, window, varargin)
            % for each sample index in times, read the window times + window(1):window(2)
            % of samples around this time from some channels

            p = inputParser();
            p.addParameter('source', 'ap', @ischar);
            
            % specify one of THESE (same channels every snippet or nChannels x numel(times))
            p.addParameter('channel_ids', [], @(x) isempty(x) || isvector(x));
            p.addParameter('channel_ids_by_snippet', [], @ismatrix);
            
            % or THESE (different channels for each group of snippets sharing the same cluster_ids
            p.addParameter('channel_ids_by_cluster', [], @ismatrix);
            p.addParameter('cluster_ids', [], @isvector); % same length as numel(times), corresponding to which cluster was pulled out (e.g. as a waveform)
            
            p.addParameter('car', false, @islogical); % subtract median over channels
            p.parse(varargin{:});

            channel_ids = p.Results.channel_ids;
            channel_ids_by_snippet = p.Results.channel_ids_by_snippet;
            channel_ids_by_cluster = p.Results.channel_ids_by_cluster;
            cluster_ids = p.Results.cluster_ids;
                
            if ~isempty(channel_ids_by_snippet)
                assert(size(channel_ids_by_snippet, 2) == numel(times), 'channel_ids must be nChannels x numel(times)');
                
            elseif ~isempty(channel_ids_by_cluster)
                assert(~isempty(cluster_ids), 'cluster_ids must be specified when channel_ids_by_cluster is used');
                unique_cluster_ids = unique(cluster_ids);
                [~, cluster_inds] = ismember(cluster_ids, unique_cluster_ids);
                channel_ids_by_snippet = channel_ids_by_cluster(:, cluster_inds);
                
            elseif ~isempty(channel_ids)
                channel_ids = Neuropixel.Utils.makecol(channel_ids);
                channel_ids_by_snippet = repmat(channel_ids, 1, numel(times));
                
            else
                error('Specify either channel_ids or channel_ids_by_cluster');
            end

            source = p.Results.source;
            switch source
                case 'ap'
                    mm = df.memmapAP_full();
                case 'lf'
                    mm = df.memmapLF_full();
                otherwise
                    error('Unknown source');
            end
            nC = size(channel_ids_by_snippet, 1);
            nS = numel(times);
            nT = numel(window(1):window(2));
            out = zeros(nC, nT, nS, 'int16');

            if numel(times) > 10
                prog = Neuropixel.Utils.ProgressBar(numel(times), 'Extracting %s snippets', upper(source));
            else
                prog = [];
            end
            for iS = 1:numel(times)
                if ~isempty(prog), prog.update(iS); end

                idx_start = times(iS)+window(1);
                idx_stop = idx_start + nT - 1;
                idx_request = idx_start:idx_stop;
                mask_idx_okay = idx_request >= uint64(1) & idx_request <= size(mm.Data.x, 2);
                idx_request = idx_request(mask_idx_okay);
                
                channel_idx = channel_ids_by_snippet(:, iS); % which channels for this spike
                if p.Results.car
                    extract = mm.Data.x(:, idx_request);
                    out(:, mask_idx_okay, iS) = extract(channel_idx, :) - median(extract(df.goodChannelInds, :), 1);
                else
                    out(:, mask_idx_okay, iS) = mm.Data.x(channel_idx, idx_request);
                end
            end
            data_ch_by_time_by_snippet = out;
            if ~isempty(prog), prog.finish(); end
        end
    end
    
    methods  % Read data at specified times
        function snippet_set = readAPSnippetSet(df, times, window, varargin)
            [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet] = ...
                df.readSnippetsRaw(times, window, 'source', 'ap', varargin{:});
            snippet_set = Neuropixel.SnippetSet(df, 'ap');
            snippet_set.data = data_ch_by_time_by_snippet;
            snippet_set.sample_idx = times;
            snippet_set.channel_ids_by_snippet = channel_ids_by_snippet;
            snippet_set.cluster_ids = cluster_ids;
            snippet_set.window = window;
        end

        function snippet_set = readLFSnippetSet(df, times, window, varargin)
            [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet] = ...
                df.readSnippetsRaw(times, window, 'source', 'lf', varargin{:});
            snippet_set = Neuropixel.SnippetSet(df, 'lf');
            snippet_set.data = data_ch_by_time_by_snippet;
            snippet_set.sample_idx = times;
            snippet_set.channel_ids_by_snippet = channel_ids_by_snippet;
            snippet_set.cluster_ids = cluster_ids;
            snippet_set.window = window;
        end

        function rms = computeRMSByChannel(df, varargin)
            p = inputParser();
            p.addParameter('sampleMaskFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % sampleMaskFn(data_ch_x_time, sample_idx_time) --> logical_time mask of time samples valid for use, useful if you have artifacts at known times
            p.addParameter('car', false, @islogical);
            p.addParameter('useChunks', 50, @isscalar);
            p.addParameter('chunkSize', 100000, @isscalar);
            p.parse(varargin{:});
            
            sampleMaskFn = p.Results.sampleMaskFn;
            
            % aim for the middle of the file
            chunkSize = min(df.fsAP, p.Results.chunkSize);
            mm = df.memmapAP_by_chunk(chunkSize);
            nChunks = numel(mm.Data);
            useChunks = min(nChunks, p.Results.useChunks);
            skipChunks = floor((nChunks-useChunks)/2);
            
            ch_mask = df.lookup_channelIds(df.mappedChannels); % for common average referencing

            sumByChunk = nan(df.nChannels, useChunks);
%             prog = Neuropixel.Utils.ProgressBar(useChunks, 'Computing RMS per channel');
            for iC =  1:useChunks
%                 prog.increment();
                data = mm.Data(iC+skipChunks).x;
                
                if p.Results.car
                    data(ch_mask, :) = data(ch_mask, :) - mean(data(ch_mask, :), 1, 'native');
                end
                
                if ~isempty(sampleMaskFn)
                    idx = (iC+skipChunks-1)*chunkSize + (1:chunkSize);
                    mask = sampleMaskFn(data, idx);
                    data = data(:, mask);
                end

                sumByChunk(:, iC) = sum((single(data) - mean(single(data), 2)).^2, 2);
            end
%             prog.finish();
            rms = sqrt(sum(sumByChunk, 2) ./ (useChunks * chunkSize));
            rms = rms * df.apScaleToUv;
        end
    end

    methods(Hidden)
        function fid = openAPFile(df)
            if ~exist(df.pathAP, 'file')
                error('RawDataFile: %s not found', df.pathAP);
            end
            fid = fopen(df.pathAP, 'r');

            if fid == -1
                 error('RawDataFile: Could not open %s', df.pathAP);
            end
        end

        function fid = openLFFile(df)
            if ~exist(df.pathLF, 'file')
                error('RawDataFile: %s not found', df.pathAP);
            end
            fid = fopen(df.pathLF, 'r');

            if fid == -1
                 error('RawDataFile: Could not open %s', df.pathAP);
            end
        end

        function fid = openSyncFile(df)
            if ~exist(df.pathSync, 'file')
                error('RawDataFile: %s not found', df.pathSync);
            end
            fid = fopen(df.pathSync, 'r');

            if fid == -1
                 error('RawDataFile: Could not open %s', df.pathSync);
            end
        end
    end

    methods % Dependent properties
        function tf = get.syncInAPFile(df)
            tf = df.channelMap.syncInAPFile;
        end
        
        function ind = get.syncChannelIndex(df)
            if df.syncInAPFile
                ind = df.channelMap.syncChannelIndex;
            else
                % if sync is in its own file, assume it's the first and only channel
                ind = uint32(1);
            end
        end
        
        function pathAP = get.pathAP(df)
            pathAP = fullfile(df.pathRoot, df.fileAP);
        end

        function fileAP = get.fileAP(df)
            fileAP = [df.fileStem '.imec.' df.fileTypeAP '.bin'];
        end

        function tf = get.hasAP(df)
            tf = exist(df.pathAP, 'file') == 2;
        end

        function fileAPMeta = get.fileAPMeta(df)
            fileAPMeta = [df.fileStem '.imec.ap.meta'];
        end

        function pathAPMeta = get.pathAPMeta(df)
            pathAPMeta = fullfile(df.pathRoot, df.fileAPMeta);
        end

        function pathLF = get.pathLF(df)
            pathLF = fullfile(df.pathRoot, df.fileLF);
        end

        function fileAP = get.fileLF(df)
            fileAP = [df.fileStem '.imec.lf.bin'];
        end

        function fileLFMeta = get.fileLFMeta(df)
            fileLFMeta = [df.fileStem '.imec.lf.meta'];
        end

        function pathLFMeta = get.pathLFMeta(df)
            pathLFMeta = fullfile(df.pathRoot, df.fileLFMeta);
        end

        function tf = get.hasLF(df)
            tf = exist(df.pathLF, 'file') == 2;
        end

        function fileSync = get.fileSync(df)
            if df.syncInAPFile
                fileSync = df.fileAP;
            else
                fileSync = [df.fileStem, '.imec.sync.bin'];
            end
        end

        function pathSync = get.pathSync(df)
            pathSync = fullfile(df.pathRoot, df.fileSync);
        end

        function fileSyncCached = get.fileSyncCached(df)
            fileSyncCached = [df.fileStem '.sync.mat'];
        end

        function pathSyncCached = get.pathSyncCached(df)
            pathSyncCached = fullfile(df.pathRoot, df.fileSyncCached);
        end

        function scale = get.apScaleToUv(df)
            scale = (df.apRange(2) - df.apRange(1)) / (2^df.adcBits) / df.apGain * 1e6;
        end

        function scale = get.lfScaleToUv(df)
            scale = (df.apRange(2) - df.apRange(1)) / (2^df.adcBits) / df.apGain * 1e6;
        end

        function file = get.channelMapFile(df)
            if isempty(df.channelMap)
                file = '';
            else
                file = df.channelMap.file;
            end
        end

        function list = get.mappedChannels(df)
            if isempty(df.channelMap)
                list = [];
            else
                list = df.channelMap.channelIds;
            end
        end
        
        function list = get.mappedChannelInds(df)
            list = df.lookup_channelIds(df.mappedChannels);
        end

        function list = get.connectedChannels(df)
            if isempty(df.channelMap)
                list = [];
            else
                list = df.channelMap.connectedChannels;
            end
        end
        
        function list = get.connectedChannelInds(df)
            list = df.lookup_channelIds(df.connectedChannels);
        end

        function n = get.nChannelsMapped(df)
            if isempty(df.channelMap)
                n = NaN;
            else
                n = df.channelMap.nChannels;
            end
        end

        function n = get.nChannelsConnected(df)
            if isempty(df.channelMap)
                n = NaN;
            else
                n = nnz(df.channelMap.connected);
            end
        end

        function ch = get.goodChannels(df)
            ch = setdiff(df.connectedChannels, df.badChannels);
        end
        
        function list = get.goodChannelInds(df)
            list = df.lookup_channelIds(df.goodChannels);
        end

        function n = get.nGoodChannels(df)
            n = numel(df.goodChannels);
        end
        
        function idx = get.channelIds(df)
            idx = df.channelMap.channelIds;
        end        
        
        function names = get.channelNames(df)
            names = strings(df.nChannels, 1);
            names(df.channelMap.channelIds) = string(sprintfc("ch %d", df.channelMap.channelIds));
            if ~isnan(df.syncChannelIndex)
                names(df.syncChannelIndex) = "sync";
            end
        end

        function names = get.channelNamesPadded(df)
            names = strings(df.nChannels, 1);
            names(df.channelMap.channelIds) = string(sprintfc("ch %03d", df.channelMap.channelIds));
            if ~isnan(df.syncChannelIndex)
                names(df.syncChannelIndex) = "sync";
            end
        end
        
        function n = get.nSyncBits(df)
            n = 8*df.bytesPerSample; % should be 16?
        end
        
        function bits = get.syncBitsNamed(df)
            names = df.syncBitNames;
            bits = find(names ~= "");
        end

        function names = get.syncBitNames(df)
            if isempty(df.syncBitNames)
                names = strings(df.nSyncBits, 1);
            else
                names = string(df.syncBitNames);
            end
        end

        function meta = readAPMeta(df)
            meta = Neuropixel.readINI(df.pathAPMeta);
        end

        function meta = generateModifiedAPMeta(df)
            meta = df.readAPMeta;

            meta.syncBitNames = df.syncBitNames;
            meta.badChannels = df.badChannels;
        end

        function writeModifiedAPMeta(df, varargin)
            p = inputParser();
            p.addParameter('extraMeta', struct(), @isstruct);
            p.parse(varargin{:});

            meta = df.generateModifiedAPMeta();

            % set extra user provided fields
            extraMeta = p.Results.extraMeta;
            extraMetaFields = fieldnames(extraMeta);
            for iFld = 1:numel(extraMetaFields)
                meta.(extraMetaFields{iFld}) = extraMeta.(extraMetaFields{iFld});
            end

            Neuropixel.writeINI([df.pathAPMeta], meta);
        end

        function meta = readLFMeta(df)
            meta = Neuropixel.readINI(df.pathLFMeta);
        end

        function str = get.creationTimeStr(df)
            str = datestr(df.creationTime);
        end
        
        function file = getAuxiliaryFileWithSuffix(df, suffix)
             suffix = char(suffix);
             file = fullfile(df.pathRoot, [df.fileStem, '.', suffix]);
        end
    end
    
    methods % Marking Channels as bad
            function [rmsBadChannels, rmsByChannel] = markBadChannelsByRMS(df, varargin)
            p = inputParser();
            p.addParameter('rmsRange', [3 100], @isvector);
            p.addParameter('sampleMaskFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % sampleMaskFn(data_ch_x_time, sample_idx_time) --> logical_time mask of time samples valid for use, useful if you have artifacts at known times
            p.parse(varargin{:});

            channelMask = true(df.nChannels, 1);

            channelMask(~df.channelMap.connected) = false;

            oldChannelMask = channelMask;

            rmsByChannel = df.computeRMSByChannel('sampleMaskFn', p.Results.sampleMaskFn);
            rmsMin = p.Results.rmsRange(1);
            rmsMax = p.Results.rmsRange(2);
            channelMask(rmsByChannel < rmsMin | rmsByChannel > rmsMax) = false;

            rmsBadChannels = find(~channelMask & oldChannelMask);
            df.markBadChannels(~channelMask);
        end

        function markBadChannels(df, list)
            % this adds to the set of bad channels, so multiple calls will
            % remove additional channels
            if islogical(list)
                list = find(list);
            end
            df.badChannels = union(df.badChannels, list);
        end

end

    methods % Modify bin data files in place
        function modifyAPInPlace(imec, varargin)
            imec.modifyInPlaceInternal('ap', varargin{:});
        end

        function modifyLFInPlace(imec, varargin)
            imec.modifyInPlaceInternal('lf', varargin{:});
        end

        function imecSym = symLinkAPIntoDirectory(imec, newFolder, varargin)
            p = inputParser();
            p.addParameter('relative', false, @islogical);
            p.parse(varargin{:});
            newFolder = char(newFolder);

            if ~exist(newFolder, 'dir')
                Neuropixel.Utils.mkdirRecursive(newFolder);
            end
            newAPPath = fullfile(newFolder, imec.fileAP);
            Neuropixel.Utils.makeSymLink(imec.pathAP, newAPPath, p.Results.relative);

            newAPMetaPath = fullfile(newFolder, imec.fileAPMeta);
            Neuropixel.Utils.makeSymLink(imec.pathAPMeta, newAPMetaPath, p.Results.relative);

            if ~imec.syncInAPFile && exist(imec.pathSync, 'file')
                newSyncPath = fullfile(newFolder, imec.fileSync);
                Neuropixel.Utils.makeSymLink(imec.pathSync, newSyncPath, p.Results.relative);
            end

            if exist(imec.pathSyncCached, 'file')
                newSyncCachedPath = fullfile(newFolder, imec.fileSyncCached);
                Neuropixel.Utils.makeSymLink(imec.pathSyncCached, newSyncCachedPath, p.Results.relative);
            end

            imecSym = Neuropixel.ImecDataset(newAPPath, 'channelMap', imec.channelMapFile);
        end
        
        function imecOut = saveTranformedDataset(df, outPath, varargin)
            p = inputParser();
            p.addParameter('transformAP', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (df, dataChunk) and return dataChunk someplace
            p.addParameter('transformLF', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (df, dataChunk) and return dataChunk someplace

            p.addParameter('gpuArray', false, @islogical);
            p.addParameter('applyScaling', false, @islogical); % convert to uV before processing

            p.addParameter('writeAP', true, @islogical);
            p.addParameter('writeSyncSeparate', false, @islogical); % true means ap will get only mapped channels, false will preserve channels as is
            p.addParameter('writeLF', false, @islogical);
            
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            p.addParameter('mappedChannelsOnly', false, @islogical);
            
            p.addParameter('chunkSize', 2^20, @isscalar);
            p.addParameter('chunkEdgeExtraSamples', [0 0], @isvector); 
            
            p.addParameter('timeShiftsAP', {}, @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            p.addParameter('timeShiftsLF', {}, @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            
            p.addParameter('extraMeta', struct(), @isstruct);
            p.parse(varargin{:});

            % this uses the same syntax as writeConcatenatedFileMatchGains
            imecOut = Neuropixel.ImecDataset.writeConcatenatedFileMatchGains({df}, outPath, p.Results);
         end
        
        function writeFolderForPhy(df, savePath, varargin)
            % writes the files needed for Phy template-gui to be able to inspect this file in Phy
            % essentially by generating a synthetic Kilosort output
            
            p = inputParser();
            p.addParameter('spikeTimes', [1 2], @(x) isempty(x) || isvector(x));
            p.addParameter('spikeClusters', [], @(x) isempty(x) || isvector(x));
            p.addParameter('clusterNames', [], @(x) isempty(x) || iscellstr(x) || isstring(x));
            p.parse(varargin{:});
            
            spikeTimes = Neuropixel.Utils.makecol(p.Results.spikeTimes);
            nSpikes = numel(spikeTimes);
            assert(nSpikes >= 2, 'At least 2 spikes are required for Phy');
                        
            if isempty(p.Results.spikeClusters)
                spikeClusters = zeros(nSpikes, 1, 'uint32');
            else
                spikeClusters = uint32(Neuropixel.Utils.makecol(p.Results.spikeClusters));
                assert(numel(spikeClusters) == nSpikes);
            end
            
            nTemplateTimepoints = 82;
            nCh = df.nGoodChannels;
            % these can't be 1 because of squeeze() inside Phy's code
            nTemplates = 2;
            nTemplateFeatures = 2;
            nFeaturesPerChannel = 2;
            nPCFeatures = 2;
            
            df.symLinkAPIntoDirectory(savePath);
            
            writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
            writeNPY(zeros(nSpikes, 1, 'uint32'), fullfile(savePath, 'spike_templates.npy'));
            writeNPY(spikeClusters, fullfile(savePath, 'spike_clusters.npy'));
    
            writeNPY(zeros(nSpikes, 1, 'double'), fullfile(savePath, 'amplitudes.npy'));
            
            templates = zeros(2, nTemplateTimepoints, nCh, 'single');
            writeNPY(templates, fullfile(savePath, 'templates.npy'));
    
            templatesInds = df.goodChannels';
            writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
            sortedInds = df.goodChannelInds;
            chanMap0ind = int32(df.channelMap.channelIdsMapped(sortedInds) - uint32(1));
            xcoords = df.channelMap.xcoords(sortedInds);
            ycoords = df.channelMap.ycoords(sortedInds);
            writeNPY(chanMap0ind, fullfile(savePath, 'channel_map.npy'));
            writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));
    
            templateFeatures = zeros([nTemplates nTemplateFeatures], 'single');
            writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
            
            templateFeatureInds = zeros(nTemplates, nTemplateFeatures, 'uint32');
            writeNPY(templateFeatureInds, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
            
            similarTemplates = zeros(nTemplates, nTemplates, 'single');
            writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
            
            pcFeatures = zeros([nSpikes, nFeaturesPerChannel, nPCFeatures], 'single');
            writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
            
            pcFeatureInds = zeros([nTemplates, nPCFeatures], 'uint32');
            writeNPY(pcFeatureInds, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
            whiteningMatrix = ones(nCh, nCh, 'double');
            writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
            whiteningMatrixInv = ones(nCh, nCh, 'double');
            writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
            
            % write params.py
            fid = fopen(fullfile(savePath,'params.py'), 'w');
            [~, fname, ext] = fileparts(df.fileAP);
            fprintf(fid,['dat_path = ''',fname ext '''\n']);
            fprintf(fid,'n_channels_dat = %i\n',df.nChannels);
            fprintf(fid,'dtype = ''int16''\n');
            fprintf(fid,'offset = 0\n');
            fprintf(fid,'sample_rate = %i\n', df.fsAP);
            fprintf(fid,'hp_filtered = False');
            fclose(fid);
            
            % write spike names column
            clusterIds = unique(spikeClusters);
            nClusters = numel(clusterIds);
            if isempty(p.Results.clusterNames)
                clusterNames = strings(nClusters, 1);
            else
                clusterNames = string(Neuropixel.Utils.makecol(p.Results.clusterNames));
                assert(numel(clusterNames) == nClusters, 'clusterNames must be nClusters long');
            end 
            fileID = fopen(fullfile(savePath, 'cluster_names.tsv'),'w');
            fprintf(fileID, 'cluster_id\tname\n');
            
            for iC = 1:nClusters
                fprintf(fileID, '%d\t%s\n', clusterIds(iC)-1, clusterNames{iC});
            end
            fclose(fileID);
    
        end
    end

    methods(Hidden)
        function [chInds, chIds] = build_channelSelectors_internal(imec, varargin)
            p = inputParser();
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            p.addParameter('mappedChannelsOnly', false, @islogical);
            p.parse(varargin{:});

            if p.Results.goodChannelsOnly
                [chInds, chIds] = imec.lookup_channelIds(imec.goodChannels);
                assert(~isempty(chInds), 'No channels marked good in dataset')

            elseif p.Results.connectedChannelsOnly
                [chInds, chIds] = imec.lookup_channelIds(imec.connectedChannels); % excludes sync channel
                assert(~isempty(chInds), 'No connected channels found in dataset');
                
            elseif p.Results.mappedChannelsOnly
                [chInds, chIds] = imec.lookup_channelIds(imec.mappedChannels); % excludes sync channel
                assert(~isempty(chInds), 'No mapped channels found in dataset');
            else
                chInds = 1:imec.nChannels;
                chIds = imec.channelIds(chInds);
            end
        end   
        
        function modifyInPlaceInternal(imec, mode, procFnList, varargin)
            p = inputParser();
            p.addParameter('chunkSize', 2^20, @isscalar);
            
            % [pre post] sextra samples to pass in with each chunk to the functions in procFnList, 
            % the output corresponding these times wont be saved and is simply used to avoid issues with processing at the edges of chunks
            p.addParameter('chunkEdgeExtraSamples', [0 0], @isvector); 
            p.addParameter('gpuArray', false, @islogical);
            p.addParameter('applyScaling', false, @islogical); % convert to uV before processing
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            p.addParameter('mappedChannelsOnly', false, @islogical);
            
            p.addParameter('extraMeta', struct(), @isstruct);
            p.addParameter('dryRun', false, @islogical); % for testing proc fn before modifying file
            p.addParameter('dryRunSampleInds', [], @(x) isempty(x) || isvector(x)); % if specified only include chunks with these sample inds 
            p.parse(varargin{:});

            chunkSize = p.Results.chunkSize;
            chunkExtra = p.Results.chunkEdgeExtraSamples;
            useGpuArray = p.Results.gpuArray;
            applyScaling = p.Results.applyScaling;
            dryRun = p.Results.dryRun;
            dryRunSampleInds = p.Results.dryRunSampleInds;

            if ~iscell(procFnList)
                procFnList = {procFnList};
            end
            if isempty(procFnList)
                error('No modification functions provided');
            end

            % open writable memmapfile
            switch mode
                case 'ap'
                    mm = imec.memmapAP_full('Writable', ~dryRun);
                case 'lf'
                    mm = imec.memmapLF_full('Writable', ~dryRun);
                otherwise
                    error('Unknown mode %s', mode);
            end

            % figure out which channels to keep
            [chInds, chIds] = imec.build_channelSelectors_internal('goodChannelsOnly', p.Results.goodChannelsOnly, ...
                'connectedChannelsOnly', p.Results.connectedChannelsOnly, ...
                'mappedChannelsOnly', p.Results.mappedChannelsOnly);

            dataSize = size(mm.Data.x, 2);
            nChunks = ceil(dataSize / chunkSize);
            prog = Neuropixel.Utils.ProgressBar(nChunks, 'Modifying %s file in place', mode);
            for iCh = 1:nChunks
                [idx, keepIdx] = Neuropixel.ImecDataset.determineChunkIdx(dataSize, iCh, nChunks, chunkSize, chunkExtra);
                
                if dryRun && ~isempty(dryRunSampleInds) && ~any(ismember(idx, dryRunSampleInds))
                    % skip this chunk unless some ind in idx should be processed
                    continue;
                end

                data = mm.Data.x(chInds, idx);

                % ch_connected_mask indicates which channels are
                % connected, which are the ones where scaling makes
                % sense. chIdx is all channels being modified by procFnList
                ch_conn_mask = ismember(chIds, imec.connectedChannels);

                % do additional processing here
                if applyScaling
                    % convert to uV and to single
                    switch mode
                        case 'ap'
                            data = single(data);
                            data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(imec.apScaleToUv);
                        case 'lf'
                            data = single(data);
                            data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(imec.lfScaleToUv);
                    end
                end

                if useGpuArray
                    data = gpuArray(data);
                end

                % apply each procFn sequentially
                for iFn = 1:numel(procFnList)
                    fn = procFnList{iFn};
                    data = fn(imec, data, chIds, idx);
                end

                if useGpuArray
                    data = gather(data);
                end

                if applyScaling
                    data(ch_conn_mask, :) = data(ch_conn_mask, :) ./ imec.scaleToUv;
                end

                data = int16(data);

                if ~dryRun
                    % slice off extra at edges
                    data = data(:, keepIdx);
                    mm.Data.x(chInds, idx) = data;
                end
                prog.increment();
            end
            prog.finish();

            if ~dryRun
                imec.writeModifiedAPMeta('extraMeta', p.Results.extraMeta);
            end
        end
    end

    methods(Static)
        function [selectIdx, keepIdx] = determineChunkIdx(dataSize, iCh, nChunks, chunkSize, chunkExtra)
            if iCh == nChunks
                % last chunk
                selectIdx = ((iCh-1)*(chunkSize)+1 - chunkExtra(1)) : dataSize;
                chunkExtraThis = [chunkExtra(1) 0];
            elseif iCh == 1
                % first chunk
                selectIdx = 1 : ((iCh-1)*(chunkSize) + chunkSize + chunkExtra(2));
                chunkExtraThis = [0 chunkExtra(2)];
            else
                % middle chunks
                selectIdx = ((iCh-1)*(chunkSize) + 1 - chunkExtra(1)) : ((iCh-1)*(chunkSize) + chunkSize + chunkExtra(2));
                chunkExtraThis = chunkExtra;
            end
            
            % trim invalid samples
            idxFirst = find(selectIdx >= 1, 1, 'first');
            if idxFirst > 1
                selectIdx = selectIdx(idxFirst:end);
                chunkExtraThis(1) = max(0, chunkExtraThis(1) - (idxFirst-1));
            end
            idxLast = find(selectIdx <= dataSize, 1, 'last');
            if idxLast <= numel(selectIdx)
                nDrop = numel(selectIdx) - idxLast;
                selectIdx = selectIdx(1:idxLast);
                chunkExtraThis(2) = max(0, chunkExtraThis(2) - nDrop);
            end
            
            nSelected = numel(selectIdx);
            keepIdx = 1+chunkExtraThis(1) : nSelected-chunkExtraThis(2);
        end   
        
        function [chIndsByFile, chIds] = multiFile_build_channelSelectors_internal(imecList, varargin)
            for iF = 1:numel(imecList)
                [~, chIdsThis] = imecList{iF}.build_channelSelectors_internal(varargin{:});
                if iF == 1
                    chIds = chIdsThis;
                else
                    chIds = intersect(chIds, chIdsThis);
                end
            end
            if isempty(chIds)
                error('No valid channels present across all datasets');
            end
            
            chIndsByFile = cellfun(@(imec) imec.lookup_channelIds(chIds), imecList, 'UniformOutput', false);
        end
        
        function [parent, leaf, ext] = filepartsMultiExt(file)
            % like fileparts, but a multiple extension file like file.test.meta
            % will end up with leaf = file and ext = .test.meta

            [parent, leaf, ext] = fileparts(file);
            if ~isempty(ext)
                [leaf, ext] = strtok([leaf, ext], '.');
            end
        end

        function tf = folderContainsDataset(fileOrFileStem)
            file = Neuropixel.ImecDataset.findImecFileInDir(fileOrFileStem, 'ap');
            if isempty(file)
                tf = false;
                return;
            end

            [pathRoot, fileStem, fileTypeAP] = Neuropixel.ImecDataset.parseImecFileName(file);
            pathAP = fullfile(pathRoot, [fileStem '.imec.' fileTypeAP '.bin']);
            pathAPMeta = fullfile(pathRoot, [fileStem '.imec.ap.meta']);

            tf = exist(pathAP, 'file') && exist(pathAPMeta, 'file');
            if exist(pathAP, 'file') && ~exist(pathAPMeta, 'file')
                warning('Found data file %s but not meta file %s', pathAP, pathAPMeta);
            end
        end

        function file = findImecFileInDir(fileOrFileStem, type)
            if nargin < 2
                type = 'ap';
            end
            fileOrFileStem = char(fileOrFileStem);

            if exist(fileOrFileStem, 'file') == 2
                file = fileOrFileStem;
                [~, ~, type] = Neuropixel.ImecDataset.parseImecFileName(file);
                switch type
                    case 'ap'
                        assert(ismember(type, {'ap', 'ap_CAR'}), 'Specify ap.bin or ap_CAR.bin file rather than %s file', type);
                    case 'lf'
                        assert(ismember(type, {'lf'}), 'Specify lf.bin file rather than %s file', type);
                end

            elseif exist(fileOrFileStem, 'dir')
                % it's a directory, assume only one imec file in directory
                path = fileOrFileStem;
                [~, leaf] = fileparts(path);

                switch type
                    case 'ap'
                        apFiles = Neuropixel.ImecDataset.listAPFilesInDir(path);
                        if ~isempty(apFiles)
                            if numel(apFiles) > 1
                                [tf, idx] = ismember([leaf '.imec.ap.bin'], apFiles);
                                if tf
                                    file = apFiles{idx};
                                    return
                                end
                                [tf, idx] = ismember([leaf '.imec.ap_CAR.bin'], apFiles);
                                if tf
                                    file = apFiles{idx};
                                    return
                                end
                                file = apFiles{1};
                                warning('Multiple AP files found in dir %s, choosing %s', path, file);
                            else
                                file = apFiles{1};
                            end
                        else
                            file = [];
                            return;
                        end

                    case 'lf'
                        lfFiles = Neuropixel.ImecDataset.listLFFilesInDir(path);
                        if ~isempty(lfFiles)
                            if numel(lfFiles) > 1
                                [tf, idx] = ismember([leaf '.imec.lf.bin'], lfFiles);
                                if tf
                                    file = lfFiles{idx};
                                    return
                                end
                                file = lfFiles{1};
                                warning('Multiple LF files found in dir, choosing %s', file);
                            else
                                file = lfFiles{1};
                            end
                        else
                            file = [];
                            return;
                        end
                    otherwise
                        error('Unknown type %s');
                end

                file = fullfile(path, file);
                
            else
                % not a folder or a file, but possibly pointing to the
                % stem of a file, e.g. '/path/data' pointing to
                % '/path/data.ap.imec.bin'
                [parent, leaf, ext] = fileparts(fileOrFileStem);
                if ~exist(parent, 'dir')
                    error('Folder %s does not exist', parent);
                end
                stem = [leaf, ext];
                
                % find possible matches
                switch type
                    case 'ap'
                        candidates = Neuropixel.ImecDataset.listAPFilesInDir(parent);
                        
                    case 'lf'
                        candidates = Neuropixel.ImecDataset.listLFFilesInDir(parent);
                        
                    otherwise
                        error('Unknown type %s');
                end
                mask = startsWith(candidates, stem);
                
                if ~any(mask)
                    error('No %s matches for %s* exist', type, fileOrFileStem);
                elseif nnz(mask) > 1
                    error('Multiple %s matches for %s* exist. Narrow down the prefix.', type, fileOrFileStem);
                end
                
                file = fullfile(parent, candidates{mask});
            end
        end

        function [pathRoot, fileStem, type] = parseImecFileName(file)
            if iscell(file)
                [pathRoot, fileStem, type] = cellfun(@Neuropixel.ImecDataset.parseImecFileName, file, 'UniformOutput', false);
                return;
            end
            file = char(file);

            [pathRoot, f, e] = fileparts(file);
            if isempty(e)
                error('No file extension specified on Imec file name');
            end
            file = [f, e];


            match = regexp(file, '(?<stem>\w+).imec.(?<type>\w+).bin', 'names', 'once');
            if ~isempty(match)
                type = match.type;
                fileStem = match.stem;
                return;
            end

            fileStem = file;
            type = '';
        end

        function apFiles = listAPFilesInDir(path)
            info = cat(1, dir(fullfile(path, '*.imec.ap.bin')), dir(fullfile(path, '*.imec.ap_CAR.bin')));
            apFiles = {info.name}';
        end

        function lfFiles = listLFFilesInDir(path)
            info = dir(fullfile(path, '*.imec.lf.bin'));
            lfFiles = {info.name}';
        end
        
        function clearDestinationStem(outPathStem)
            assert(~isempty(outPathStem));
            files = dir([outPathStem '*']);
            
            for iF = 1:numel(files)
                if ~files(iF).isdir
                    f = fullfile(files(iF).folder, files(iF).name);
                    delete(f);
                end
            end
        end

        function imecOut = writeConcatenatedFileMatchGains(imecList, outPath, varargin)
            p = inputParser();
            p.addParameter('writeAP', true, @islogical);
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('mappedChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            p.addParameter('writeSyncSeparate', false, @islogical); % true means ap will get only mapped channels, false will preserve channels as is
            p.addParameter('writeLF', false, @islogical);
            p.addParameter('chunkSize', 2^20, @isscalar);
            p.addParameter('chunkEdgeExtraSamples', [0 0], @isvector); 

            p.addParameter('gpuArray', false, @islogical);
            p.addParameter('applyScaling', false, @islogical); % convert to uV before processing

            p.addParameter('transformAP', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (df, dataChunk) and return dataChunk someplace
            p.addParameter('transformLF', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (df, dataChunk) and return dataChunk someplace

            p.addParameter('timeShiftsAP', {}, @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            p.addParameter('timeShiftsLF', {}, @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            
            p.addParameter('extraMeta', struct(), @isstruct);
            p.addParameter('dryRun', false, @islogical);
            p.parse(varargin{:});

            nFiles = numel(imecList);
            stemList = cellfun(@(imec) imec.fileStem, imecList, 'UniformOutput', false);
            dryRun = p.Results.dryRun;
            
            function s = lastFilePart(f)
                [~, f, e] = fileparts(f);
                s = [f, e];
            end

            [parent, leaf, ext] = Neuropixel.ImecDataset.filepartsMultiExt(outPath);
            if ~isempty(ext) && endsWith(ext, 'bin')
                % specified full file
                outPath = parent;
            else
                outPath = fullfile(parent, [leaf, ext]);
            end
            if ~exist(outPath, 'dir') && ~dryRun
                Neuropixel.Utils.mkdirRecursive(outPath);
            end

            % determine the gains that we will use
            function [multipliers, gain] = determineCommonGain(gains)
                uGains = unique(gains);

                if numel(uGains) == 1
                    gain = uGains;
                    multipliers = ones(nFiles, 1);
                    if numel(gains) > 1
                        fprintf('All files have common gain of %d\n', gain);
                    end
                else
                    % find largest gain that we can achieve by multiplying each
                    % file by an integer (GCD)
                    gain = uGains(1);
                    for g = uGains(2:end)
                        gain = lcm(gain, g);
                    end
                    multipliers = gain ./ gains;
                    assert(all(multipliers==round(multipliers)));
                    fprintf('Converting all files to gain of %d\n', gain);
                end

                multipliers = int16(multipliers);
            end

            % figure out which channels to keep
            [chIndsByFile, ~] = Neuropixel.ImecDataset.multiFile_build_channelSelectors_internal(imecList, 'goodChannelsOnly', p.Results.goodChannelsOnly, ...
                'connectedChannelsOnly', p.Results.connectedChannelsOnly, 'mappedChannelsOnly', p.Results.mappedChannelsOnly);

            chunkSize = p.Results.chunkSize;
            chunkEdgeExtraSamples = p.Results.chunkEdgeExtraSamples;

            useGpuArray = p.Results.gpuArray;
            applyScaling = p.Results.applyScaling;
            timeShiftsAP = p.Results.timeShiftsAP;
            timeShiftsLF = p.Results.timeShiftsLF;
            isConcatenation = numel(imecList) > 1;
            
            if ~dryRun
                Neuropixel.ImecDataset.clearDestinationStem(fullfile(outPath, leaf));
            end
            
            if p.Results.writeAP || ~isempty(p.Results.transformAP)
                gains = cellfun(@(imec) imec.apGain, imecList);
                [multipliers, gain] = determineCommonGain(gains);

                outFile = fullfile(outPath, [leaf '.imec.ap.bin']);
                metaOutFile = fullfile(outPath, [leaf '.imec.ap.meta']);

                % generate new meta file
                meta = imecList{1}.generateModifiedAPMeta();

                % adjust imroTabl to set gain correctly
                m = regexp(meta.imroTbl, '\(([\d, ]*)\)', 'tokens');
                pieces = cell(numel(m), 1);
                pieces{1} = m{1}{1};
                for iM = 2:numel(m)
                    gainVals = strsplit(m{iM}{1}, ' ');
                    gainVals{4} = sprintf('%d', gain);
                    pieces{iM} = strjoin(gainVals, ' ');
                end
                meta.imroTbl = ['(' strjoin(pieces, ')('), ')'];
                meta.fileName = [leaf '.imec.ap.bin'];
                
                % indicate concatenation time points in meta file
                if isConcatenation
                    meta.concatenated = strjoin(stemList, ':');
                    meta.concatenatedSamples = cellfun(@(imec) imec.nSamplesAP, imecList);
                    meta.concatenatedGains = gains;
                    meta.concatenatedMultipliers = multipliers;
                    meta.concatenatedAdcBits = cellfun(@(imec) imec.adcBits, imecList);
                    meta.concatenatedAiRangeMin = cellfun(@(imec) imec.apRange(1), imecList);
                    meta.concatenatedAiRangeMax = cellfun(@(imec) imec.apRange(2), imecList);
                end
                
                if ~isempty(timeShiftsAP) && isConcatenation
                    % log time shifts by file in meta
                    meta.concatenatedTimeShifts = strjoin(arrayfun(@(shift) shift.as_string(), timeShiftsAP), '; ');
                end

                % compute union of badChannels
                for iM = 2:numel(imecList)
                    meta.badChannels = union(meta.badChannels, imecList{iM}.badChannels);
                end

                % set extra user provided fields
                extraMeta = p.Results.extraMeta;
                extraMetaFields = fieldnames(extraMeta);
                for iFld = 1:numel(extraMetaFields)
                    meta.(extraMetaFields{iFld}) = extraMeta.(extraMetaFields{iFld});
                end

                fprintf('Writing AP meta file %s\n', lastFilePart(metaOutFile));
                if ~dryRun
                    Neuropixel.writeINI(metaOutFile, meta);
                end
                
                fprintf('Writing AP bin file %s\n', (outFile));
                writeCatFile(outFile, chIndsByFile, 'ap', multipliers, chunkSize, chunkEdgeExtraSamples, ...
                    p.Results.transformAP, timeShiftsAP, dryRun);
            end

            if p.Results.writeLF || ~isempty(p.Results.transformLF)
                gains = cellfun(@(imec) imec.lfGain, imecList);
                [multipliers, gain] = determineCommonGain(gains);

                outFile = fullfile(outPath, [leaf '.imec.lf.bin']);
                metaOutFile = fullfile(outPath, [leaf '.imec.lf.meta']);

                % generate new meta file
                meta = imecList{1}.readLFMeta();
                % adjust imroTabl to set gain correctly
                m = regexp(meta.imroTbl, '\(([\d, ]*)\)', 'tokens');
                pieces = cell(numel(m), 1);
                pieces{1} = m{1}{1};
                for iM = 2:numel(m)
                    gainVals = strsplit(m{iM}{1}, ' ');
                    gainVals{4} = sprintf('%d', gain);
                    pieces{iM} = strjoin(gainVals, ' ');
                end
                meta. imroTbl = ['(' strjoin(pieces, ')('), ')'];
                meta.fileName = [leaf '.imec.lf.bin'];
                
                % indicate concatenation time points in meta file
                if isConcatenation
                    meta.concatenated = strjoin(stemList, ':');
                    meta.concatenatedSamples = cellfun(@(imec) imec.nSamplesLF, imecList);
                    meta.concatenatedGains = gains;
                    meta.concatenatedMultipliers = multipliers;
                    meta.concatenatedAdcBits = cellfun(@(imec) imec.adcBits, imecList);
                    meta.concatenatedAiRangeMin = cellfun(@(imec) imec.lfRange(1), imecList);
                    meta.concatenatedAiRangeMax = cellfun(@(imec) imec.lfRange(2), imecList);
                end
                
                if ~isempty(timeShiftsLF) && isConcatenation
                    % log time shifts by file in meta
                    meta.concatenatedTimeShifts = strjoin(arrayfun(@(shift) shift.as_string(), timeShiftsLF, 'UniformOutput', false), '; ');
                end

                % compute union of badChannels
                for iM = 2:numel(imecList)
                    meta.badChannels = union(meta.badChannels, imecList{iM}.badChannels);
                end
                
                fprintf('Writing LF meta file %s\n', lastFilePart(metaOutFile));
                if ~dryRun
                    Neuropixel.writeINI(metaOutFile, meta);
                end

                fprintf('Writing LF bin file %s\n', lastFilePart(outFile));
                writeCatFile(outFile, chIndsByFile, 'lf', multipliers, chunkSize, chunkEdgeExtraSamples, ...
                    p.Results.transformLF, timeShiftsLF, dryRun);
            end

            if p.Results.writeSyncSeparate
                outFile = fullfile(outPath, [leaf '.imec.sync.bin']);
                fprintf('Writing separate sync bin file %s', lastFilePart(outFile));
                writeCatFile(outFile, imecList{1}.syncChannelIndex, 'sync', ones(nFiles, 1, 'int16'), chunkSize, chunkEdgeExtraSamples, ...
                    {}, timeShiftsAP, dryRun);
            end

            outFile = fullfile(outPath, [leaf '.imec.ap.bin']);
            imecOut = Neuropixel.ImecDataset(outFile, 'channelMap', imecList{1}.channelMapFile);

            function writeCatFile(outFile, chIndsByFile, mode, multipliers, chunkSize, chunkExtra, procFnList, timeShifts, dryRun)
                if ~iscell(procFnList)
                    procFnList = {procFnList};
                end
                multipliers = int16(multipliers);
                
                % generate new ap.bin file
                if ~dryRun
                    fidOut = fopen(outFile, 'w');
                    if fidOut == -1
                        error('Error opening output file %s', outFile);
                    end
                end

                for iF = 1:nFiles
                    fprintf("Writing contents of %s\n", imecList{iF}.fileStem);
                    
                    chInds = chIndsByFile{iF};
                    chIds = imecList{iF}.channelIds(chInds);

                    switch mode
                        case 'ap'
                            mm = imecList{iF}.memmapAP_full();
                            nSamplesSource = imecList{iF}.nSamplesAP;
                        case 'lf'
                            mm = imecList{iF}.memmapLF_full();
                            nSamplesSource = imecList{iF}.nSamplesLF;
                        case 'sync'
                            mm = imecList{iF}.memmapSync_full();
                            nSamplesSource = imecList{iF}.nSamplesAP;
                    end

                    % build idx vector
                    if isempty(timeShifts)
                        outSize = size(mm.Data.x, 2);
                        sourceIdxList = uint64(1):uint64(outSize);
                    else
                        sourceIdxList = timeShifts(iF).computeSourceIndices(nSamplesSource);
                        outSize = numel(sourceIdxList);
                    end
                    
                    nChunks = ceil(outSize / chunkSize);
                    prog = Neuropixel.Utils.ProgressBar(nChunks, 'Copying %s file %d / %d: %s', mode, iF, nFiles, imecList{iF}.fileStem);
                    
                    testChunkSize = true;
                    if testChunkSize
                        nOut = nan(nChunks, 1);
                        for iCh = 1:nChunks
                            [~, keepIdx] = Neuropixel.ImecDataset.determineChunkIdx(outSize, iCh, nChunks, chunkSize, chunkExtra);
                            nOut(iCh) = numel(keepIdx);
                        end
                        assert(sum(nOut) == outSize, 'Error with chunk selection');
                    end
                    
                    for iCh = 1:nChunks
                        [idx, keepIdx] = Neuropixel.ImecDataset.determineChunkIdx(outSize, iCh, nChunks, chunkSize, chunkExtra);

                        source_idx = sourceIdxList(idx);
                        data = mm.Data.x(chInds, source_idx);

                        % ch_connected_mask indicates which channels are
                        % connected, which are the ones where scaling makes
                        % sense. chIdx is all channels being written to
                        % output file
                        ch_conn_mask = ismember(chIds, imecList{iF}.connectedChannels);

                        if multipliers(iF) > 1
                            data(ch_conn_mask, :) = data(ch_conn_mask, :) * multipliers(iF);
                        end

                        % do additional processing here
                        if ~isempty(procFnList)
                            if applyScaling
                                % convert to uV and to single
                                switch mode
                                    case 'ap'
                                        data = single(data);
                                        data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(imecList{iF}.apScaleToUv);
                                    case 'lf'
                                        data = single(data);
                                        data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(imecList{iF}.lfScaleToUv);
                                end
                            end

                            if useGpuArray
                                data = gpuArray(data);
                            end

                            % apply each procFn sequentially
                            for iFn = 1:numel(procFnList)
                                fn = procFnList{iFn};
                                data = fn(imecList{iF}, data, chIds, source_idx);
                            end

                            if useGpuArray
                                data = gather(data);
                            end

                            if applyScaling
                                data(ch_conn_mask, :) = data(ch_conn_mask, :) ./ imecList{iF}.scaleToUv;
                            end

                            data = int16(data);
                        end

                        if ~dryRun
                            data = data(:, keepIdx);
                            fwrite(fidOut, data, 'int16');
                        end
                        prog.increment();
                    end
                    prog.finish();
                end
                
                if ~dryRun
                    fclose(fidOut);
                end
            end
        end
    end
end
