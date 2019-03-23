classdef ImecDataset < handle

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

        syncChannelIndex = NaN; % if in AP file, specify this
        syncInAPFile logical = true; % is the sync info in the ap file, or in a separate .sync file

        % see markBadChannels
        badChannels

        syncBitNames string;
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
        nChannelsMapped % number of channels in the channel map (excludes sync)

        connectedChannels
        nChannelsConnected % excludes reference and sync channels

        goodChannels % connected channels sans badChannels
        nGoodChannels
        
        channelIdx % 1 index channel list from channelMapFile
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
    end

    methods
        function df = ImecDataset(fileOrFileStem, varargin)
            p = inputParser();
            p.addParameter('channelMap', [], @(x) true);
            p.addParameter('syncInAPFile', true, @islogical);
            p.addParameter('syncChannelIndex', NaN, @isscalar);
            p.addParameter('syncBitNames', [], @(x) isempty(x) || isstring(x) || iscellstr(x));
            p.parse(varargin{:})

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

            if p.Results.syncInAPFile
                % check for cached sync file
                if isnan(p.Results.syncChannelIndex)
                    % assume last channel in ap file
                    df.syncChannelIndex = df.nChannels;
                else
                    df.syncChannelIndex = p.Results.syncChannelIndex;
                end
                df.syncInAPFile = true;
            else
                if isnan(p.Results.syncChannelIndex)
                    % assume first channel in sync file
                    df.syncChannelIndex = 1;
                else
                    df.syncChannelIndex = p.Results.syncChannelIndex;
                end
                df.syncInAPFile = false;
            end

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

            if df.hasLF
                fid = df.openLFFile();
                fseek(fid, 0, 'eof');
                bytes = ftell(fid);
                fclose(fid);
                df.nSamplesLF = bytes / df.bytesPerSample / df.nChannels;
                assert(round(df.nSamplesAP) == df.nSamplesAP, 'LF bin file size is not an integral number of samples');
            end
        end

        function setSyncBitNames(df, idx, names)
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
            mkdirRecursive(newRoot);

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
                debug('Copying to %s\n', to);
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
            p.addOptional('reload', false, @islogical);
            p.addParameter('ignoreCached', false, @islogical);
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

        function updateSyncCached(df)
            if exist(df.pathSyncCached, 'file')
                sync = df.readSync();
                save(df.pathSyncCached, 'sync');
            end
        end

        function tf = getSyncBit(df, bit)
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
            if isstring(bits)
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
            data = mm.Data.x(:, sampleIdx);

            if p.Results.applyScaling
                data = single(data);
                ch_conn_mask = df.lookup_channel_ids(df.connectedChannels);
                data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(df.apScaleToUv);
            end
        end
       
        function [mat, sampleIdx] = readAP_timeWindow(df, timeWindowSec, varargin)
            idxWindow = df.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = df.readAP_idx(sampleIdx, varargin{:});
        end
        
        function mat = readSyncBits_timeWindow(df, bits, timeWindowSec)
            idxWindow = df.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = readSyncBits_idx(bits, sampleIdx);
        end
    end
    
    methods % Quick inspection
        function [channelInds, channelIds] = lookup_channel_ids(df, channelIds)
             if islogical(channelIds)
                channelIds = df.channelIdx(channelIds);
             end
            [tf, channelInds] = ismember(channelIds, df.channelIdx);
            assert(all(tf), 'Some channel ids not found');
        end
        
        function inspectAP_timeWindow(df, timeWindowSec, varargin)
            idxWindow = df.closestSampleAPForTime(timeWindowSec);
            df.inspectAP_idxWindow(idxWindow, varargin{:});
        end
        
        function inspectAP_idxWindow(df, idxWindow, varargin)
            p = inputParser();
            p.addParameter('channels', df.mappedChannels, @isvector);
            p.addParameter('invertChannels', true, @islogical);
            p.addParameter('syncBits', df.syncBitsNamed, @isvector);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 0.95, @isscalar);
            p.addParameter('car', false, @islogical);
            p.addParameter('downsample',1, @isscalar); 
            p.parse(varargin{:});
            
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = df.readAP_idx(sampleIdx);
            labels = df.channelNamesPadded;
            
            [channelInds, channelIds] = df.lookup_channel_ids(p.Results.channels);
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
            
            if p.Results.invertChannels
                mat = flipud(mat);
                labels = flipud(labels);
                colors = flipud(colors);
            end
            
            % append sync bit info to plot in red
            syncBits = p.Results.syncBits;
            if ~isempty(syncBits)
                syncBitMat = df.readSyncBits_idx(syncBits, sampleIdx);
                mat = cat(1, mat, syncBitMat);
                colors = cat(1, colors, repmat([1 0 0], size(syncBitMat, 1), 1));
                labels = cat(1, labels, df.syncBitNames(syncBits));
                normalizeMask = cat(1, normalizeMask, false(size(syncBitMat, 1), 1));
            end

            if ~p.Results.showLabels
                labels = [];
            end
            Neuropixel.Utils.plotStackedTraces(sampleIdx, mat', 'colors', colors, 'labels', labels, ...
                'gain', p.Results.gain, 'normalizeMask', normalizeMask, 'normalizeEach', false);
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

    methods % Read data at specified times
        function [data_ch_by_time_by_snippet, cluster_idx, channel_idx_by_cluster, unique_cluster_idx] = readAPSnippetsRaw(df, times, window, varargin)
            % for each sample index in times, read the window times + window(1):window(2)
            % of samples around this time from some channels

            p = inputParser();
            p.addParameter('channel_idx_by_cluster', [], @ismatrix);
            p.addParameter('unique_cluster_idx', [], @isvector);
            p.addParameter('cluster_idx', ones(numel(times), 1), @isvector);
            p.addParameter('car', false, @islogical); % subtract median over channels
            p.parse(varargin{:});

            channel_idx_by_cluster = p.Results.channel_idx_by_cluster;
            unique_cluster_idx = p.Results.unique_cluster_idx;
            if isempty(unique_cluster_idx)
                unique_cluster_idx = unique(p.Results.cluster_idx);
            end
            if isempty(channel_idx_by_cluster)
                channel_idx_by_cluster = repmat(df.channelMap.chanMap, 1, numel(unique_cluster_idx));
            end
            if size(channel_idx_by_cluster, 1) == 1
                % same channels each cluster, repmat column to make matrix
                channel_idx_by_cluster = repmat(channel_idx_by_cluster, 1, numel(unique_cluster_idx));
            end
            assert(numel(unique_cluster_idx) == size(channel_idx_by_cluster, 2), ...
                'unique_cluster_idx must have same number of clusters as columns in channel_idx_by_cluster');

            cluster_idx = Neuropixel.Utils.makecol(p.Results.cluster_idx);
            if isscalar(cluster_idx)
                cluster_idx = repmat(cluster_idx, numel(times), 1);
            else
                assert(numel(cluster_idx) == numel(times), 'cluster_idx must have same length as requested times');
            end

            [tf, cluster_ind] = ismember(cluster_idx, unique_cluster_idx);
            assert(all(tf), 'Some cluster_idx were not found in unique_cluster_idx');

            mm = df.memmapAP_full();
            nC = size(channel_idx_by_cluster, 1);
            nS = numel(times);
            nT = numel(window(1):window(2));
            out = zeros(nC, nT, nS, 'int16');

            prog = ProgressBar(numel(times), 'Extracting AP snippets');
            for iS = 1:numel(times)
                prog.update(iS);

                idx_start = times(iS)+window(1);
                idx_stop = idx_start + nT - 1;

                channel_idx = channel_idx_by_cluster(:, cluster_ind(iS)); % which channels for this spike
                if p.Results.car
                    extract = mm.Data.x(:, idx_start:idx_stop);
                    out(:, :, iS) = extract(channel_idx, :) - median(extract, 1);
                else
                    out(:, :, iS) = mm.Data.x(channel_idx, idx_start:idx_stop);
                end
            end
            data_ch_by_time_by_snippet = out;
            prog.finish();
        end

        function snippet_set = readAPSnippetSet(df, times, window, varargin)
            [data_ch_by_time_by_snippet, cluster_idx, channel_idx_by_cluster, unique_cluster_idx] = ...
                df.readAPSnippetsRaw(times, window, varargin{:});
            snippet_set = Neuropixel.SnippetSet(df);
            snippet_set.data = data_ch_by_time_by_snippet;
            snippet_set.sample_idx = times;
            snippet_set.channel_idx_by_cluster = channel_idx_by_cluster;
            snippet_set.cluster_idx = cluster_idx;
            snippet_set.unique_cluster_idx = unique_cluster_idx;
            snippet_set.window = window;
        end

        function rms = computeRMSByChannel(df, varargin)
            p = inputParser();
            p.addParameter('sampleMask', [], @(x) isempty(x) || islogical(x));
            p.parse(varargin{:});
            % skip the first few chunks

            skipChunks = 5;
            useChunks = 5;
            chunkSize = 1000000;
            mm = df.memmapAP_by_chunk(chunkSize);

            sumByChunk = nan(df.nChannels, useChunks);
            prog = ProgressBar(useChunks, 'Computing RMS per channel');
            for iC =  1:useChunks
                prog.increment();
                data = mm.Data(iC+skipChunks).x;

                if ~isempty(p.Results.sampleMask)
                    idx = (iC+skipChunks-1)*chunkSize + (1:chunkSize);
                    mask = p.Results.sampleMask(idx);
                    data = data(:, mask);
                end

                sumByChunk(:, iC) = sum((single(data) - mean(single(data), 2)).^2, 2);
            end
            prog.finish();
            rms = sqrt(sum(sumByChunk, 2) ./ (useChunks * chunkSize));
            rms = rms * df.apScaleToUv;
        end

        function rmsBadChannels = markBadChannelsByRMS(df, varargin)
            p = inputParser();
            p.addParameter('rmsRange', [3 100], @isvector);
            p.addParameter('sampleMask', [], @(x) isempty(x) || islogical(x));
            p.parse(varargin{:});

            channelMask = true(df.nChannels, 1);

            channelMask(~df.channelMap.connected) = false;

            oldChannelMask = channelMask;

            rms = df.computeRMSByChannel('sampleMask', p.Results.sampleMask);
            rmsMin = p.Results.rmsRange(1);
            rmsMax = p.Results.rmsRange(2);
            channelMask(rms < rmsMin | rms > rmsMax) = false;

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

    methods(Hidden)
%         function data_ch_by_time = readChannelBand(df, type, chFirst, chLast, sampleFirst, sampleLast, msg)
%             if nargin < 5 || isempty(sampleFirst)
%                 sampleFirst = 1;
%             end
%             if nargin < 7 || isempty(msg)
%                 msg = 'Loading data from IMEC data file';
%             end
% 
%             switch type
%                 case {'ap', 'ap_CAR'}
%                     fid = df.openAPFile();
%                     nSamplesFull = df.nSamplesAP;
%                 case 'lf'
%                     fid = df.openLFFile();
%                     nSamplesFull = df.nSamplesLF;
%                 otherwise
%                     error('Unknown type %s', type);
%             end
%             if nargin < 6 || isempty(sampleLast)
%                 sampleLast = nSamplesFull;
%             end
% 
%             % skip to the sampleFirst channel
%             bytesOffsetSample = (sampleFirst-1)*df.nChannels*df.bytesPerSample;
% 
%             % skip to the chFirst channel
%             bytesOffsetChannel = (chFirst-1)*df.bytesPerSample;
% 
%             fseek(fid, bytesOffsetSample + bytesOffsetChannel,'bof');
% 
%             % read this channel only
%             nChRead = chLast - chFirst + 1;
%             skipBytes = (df.nChannels-nChRead)*df.bytesPerSample;
%             readStr = sprintf('%d*int16', nChRead);
%             nSamplesRead = sampleLast - sampleFirst + 1;
% 
%             % split into large reads
%             samplesPerSplit = 2^18;
%             nReadSplits = ceil(nSamplesRead / samplesPerSplit);
%             prog = ProgressBar(nReadSplits, msg);
%             data_ch_by_time = zeros(nChRead, nSamplesRead, 'int16');
%             sampleOffset = 0;
%             for iS = 1:nReadSplits
%                 prog.update(iS);
%                 nSamplesReadThis = min(samplesPerSplit, nSamplesRead - sampleOffset);
%                 data_ch_by_time(:, sampleOffset + (1:nSamplesReadThis)) = fread(fid, [nChRead, nSamplesReadThis], readStr, skipBytes);
%                 sampleOffset = sampleOffset + nSamplesReadThis;
%             end
%             fclose(fid);
%         end

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
                list = df.channelMap.chanMap;
            end
        end

        function list = get.connectedChannels(df)
            if isempty(df.channelMap)
                list = [];
            else
                list = df.channelMap.connectedChannels;
            end
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

        function n = get.nGoodChannels(df)
            n = numel(df.goodChannels);
        end
        
        function idx = get.channelIdx(df)
            idx = df.channelMap.chanMap;
        end
        
        function names = get.channelNames(df)
            names = strings(df.nChannels, 1);
            names(df.channelMap.chanMap) = string(sprintfc("ch %d", df.channelMap.chanMap));
            if ~isnan(df.syncChannelIndex)
                names(df.syncChannelIndex) = "sync";
            end
        end

        function names = get.channelNamesPadded(df)
            names = strings(df.nChannels, 1);
            names(df.channelMap.chanMap) = string(sprintfc("ch %03d", df.channelMap.chanMap));
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

            if ~exist(newFolder, 'dir')
                mkdirRecursive(newFolder);
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
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('writeSyncSeparate', false, @islogical); % true means ap will get only mapped channels, false will preserve channels as is
            p.addParameter('writeLF', false, @islogical);
            p.addParameter('chunkSize', 2^20, @isscalar);

            p.addParameter('extraMeta', struct(), @isstruct);
            p.parse(varargin{:});

            % this uses the same syntax as writeConcatenatedFileMatchGains
            imecOut = Neuropixel.ImecDataset.writeConcatenatedFileMatchGains({df}, outPath, p.Results);
        end
    end

    methods(Hidden)
        function modifyInPlaceInternal(imec, mode, procFnList, varargin)
            p = inputParser();
            p.addParameter('chunkSize', 2^20, @isscalar);
            p.addParameter('gpuArray', false, @islogical);
            p.addParameter('applyScaling', false, @islogical); % convert to uV before processing
            p.addParameter('goodChannelsOnly', true, @islogical);
            p.addParameter('mappedChannelsOnly', true, @islogical);
            p.addParameter('debug', false, @islogical); % for testing proc fn before modifying file

            p.addParameter('extraMeta', struct(), @isstruct);
            p.parse(varargin{:});

            chunkSize = p.Results.chunkSize;
            useGpuArray = p.Results.gpuArray;
            applyScaling = p.Results.applyScaling;
            debug = p.Results.debug;

            if ~iscell(procFnList)
                procFnList = {procFnList};
            end
            if isempty(procFnList)
                error('No modification functions provided');
            end

            % open writable memmapfile
            switch mode
                case 'ap'
                    mm = imec.memmapAP_full('Writable', ~debug);
                case 'lf'
                    mm = imec.memmapLF_full('Writable', ~debug);
                otherwise
                    error('Unknown mode %s', mode);
            end

            % figure out which channels to keep
            if p.Results.goodChannelsOnly
                chIdx = imec.goodChannels;
                assert(~isempty(chIdx), 'No goodChannels specified across all datasets')

            elseif p.Results.mappedChannelsOnly
                chIdx = imec.mappedChannels; % excludes sync channel
                assert(~isempty(chIdx), 'No mapped channels found in first dataset');
            else
                chIdx = 1:imec.nChannels;
            end

            nChunks = ceil(size(mm.Data.x, 2) / chunkSize);
            prog = ProgressBar(nChunks, 'Modifying %s file in place', mode);
            for iCh = 1:nChunks
                if iCh == nChunks
                    idx = (iCh-1)*(chunkSize)+1 : size(mm.Data.x, 2);
                else
                    idx = (iCh-1)*(chunkSize) + (1:chunkSize);
                end

                data = mm.Data.x(chIdx, idx);

                % ch_connected_mask indicates which channels are
                % connected, which are the ones where scaling makes
                % sense. chIdx is all channels being modified by procFnList
                ch_conn_mask = ismember(chIdx, imec.connectedChannels);

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
                    data = fn(imec, data, chIdx, idx);
                end

                if useGpuArray
                    data = gather(data);
                end

                if applyScaling
                    data(ch_conn_mask, :) = data(ch_conn_mask, :) ./ imec.scaleToUv;
                end

                data = int16(data);

                if ~debug
                    mm.Data.x(chIdx, idx) = data;
                end
                prog.increment();
            end
            prog.finish();

            imec.writeModifiedAPMeta('extraMeta', p.Results.extraMeta);
        end
    end

    methods(Static)
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
                                warning('Multiple AP files found in dir, choosing %s', file);
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
            p.addParameter('writeSyncSeparate', false, @islogical); % true means ap will get only mapped channels, false will preserve channels as is
            p.addParameter('writeLF', false, @islogical);
            p.addParameter('chunkSize', 2^20, @isscalar);

            p.addParameter('gpuArray', false, @islogical);
            p.addParameter('applyScaling', false, @islogical); % convert to uV before processing

            p.addParameter('transformAP', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (df, dataChunk) and return dataChunk someplace
            p.addParameter('transformLF', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (df, dataChunk) and return dataChunk someplace

            p.addParameter('timeShifts', {}, @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            
            p.addParameter('extraMeta', struct(), @isstruct);
            p.parse(varargin{:});

            nFiles = numel(imecList);
            stemList = cellfun(@(imec) imec.fileStem, imecList, 'UniformOutput', false);

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
            if ~exist(outPath, 'dir')
                mkdirRecursive(outPath);
            end

            % determine the gains that we will use
            function [multipliers, gain] = determineCommonGain(gains)
                uGains = unique(gains);

                if numel(uGains) == 1
                    gain = uGains;
                    multipliers = ones(nFiles, 1);
                    if numel(gains) > 1
                        debug('All files have common gain of %d\n', gain);
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
                    debug('Converting all files to gain of %d\n', gain);
                end

                multipliers = int16(multipliers);
            end

            % figure out which channels to keep
            imec1 = imecList{1};
            if p.Results.goodChannelsOnly
                goodMat = false(imec1.nChannels, nFiles);
                for i = 1:nFiles
                    goodMat(imecList{i}.goodChannels, i) = true;
                end
                chIdx = find(all(goodMat, 2));
                assert(~isempty(chIdx), 'No goodChannels specified across all datasets')

            elseif p.Results.writeSyncSeparate
                chIdx = imec1.mappedChannels; % excludes sync channel
                assert(~isempty(chIdx), 'No mapped channels found in first dataset');

            else
                chIdx = 1:imec1.nChannels;
            end

            chunkSize = p.Results.chunkSize;

            useGpuArray = p.Results.gpuArray;
            applyScaling = p.Results.applyScaling;
            timeShifts = p.Results.timeShifts;
            
            Neuropixel.ImecDataset.clearDestinationStem(fullfile(outPath, leaf));

            if p.Results.writeAP
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
                meta.concatenated = strjoin(stemList, ':');

                % indicate concatenation time points in meta file
                meta.concatenatedSamples = cellfun(@(imec) imec.nSamplesAP, imecList);
                meta.concatenatedGains = gains;
                meta.concatenatedMultipliers = multipliers;
                meta.concatenatedAdcBits = cellfun(@(imec) imec.adcBits, imecList);
                meta.concatenatedAiRangeMin = cellfun(@(imec) imec.apRange(1), imecList);
                meta.concatenatedAiRangeMax = cellfun(@(imec) imec.apRange(2), imecList);
                
                if ~isempty(timeShifts)
                    % log time shifts by file in meta
                    meta.concatenatedTimeShifts = strjoin(arrayfun(@(shift) shift.as_string(), timeShifts, 'UniformOutput', false), '; ');
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

                debug('Writing AP meta file %s\n', lastFilePart(metaOutFile));
                Neuropixel.writeINI(metaOutFile, meta);

                debug('Writing AP bin file %s\n', lastFilePart(outFile));
                writeCatFile(outFile, chIdx, 'ap', multipliers, chunkSize, p.Results.transformAP, timeShifts);
            end

            if p.Results.writeLF
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
                meta.imroTbl = ['(' strjoin(pieces, ')('), ')'];
                meta.fileName = [leaf '.imec.lf.bin'];
                meta.concatenated = strjoin(stemList, ':');

                % indicate concatenation time points in meta file
                meta.concatenatedSamples = cellfun(@(imec) imec.nSamplesLF, imecList);
                meta.concatenatedGains = gains;
                meta.concatenatedMultipliers = multipliers;
                meta.concatenatedAdcBits = cellfun(@(imec) imec.adcBits, imecList);
                meta.concatenatedAiRangeMin = cellfun(@(imec) imec.lfRange(1), imecList);
                meta.concatenatedAiRangeMax = cellfun(@(imec) imec.lfRange(2), imecList);

                debug('Writing LF meta file %s\n', lastFilePart(metaOutFile));
                Neuropixel.writeINI(metaOutFile, meta);

                debug('Writing LF bin file %s\n', lastFilePart(outFile));
                writeCatFile(outFile, chIdx, 'lf', multipliers, chunkSize);
            end

            if p.Results.writeSyncSeparate
                outFile = fullfile(outPath, [leaf '.imec.sync.bin']);
                debug('Writing separate sync bin file %s', lastFilePart(outFile));
                writeCatFile(outFile, imec1.syncChannelIndex, 'sync', ones(nFiles, 1, 'int16'), chunkSize, p.Results.transformAP, timeShifts);
            end

            outFile = fullfile(outPath, [leaf '.imec.ap.bin']);
            imecOut = Neuropixel.ImecDataset(outFile, 'channelMap', imec1.channelMapFile);

            function writeCatFile(outFile, chIdx, mode, multipliers, chunkSize, procFnList, timeShifts)
                if nargin < 6
                    procFnList = {};
                end
                if nargin < 7
                    timeShifts = {};
                end
                if ~iscell(procFnList)
                    procFnList = {procFnList};
                end
                multipliers = int16(multipliers);

                 % generate new ap.bin file
                fidOut = fopen(outFile, 'w');
                if fidOut == -1
                    error('Error opening output file %s', outFile);
                end

                for iF = 1:nFiles
                    debug("Writing contents of %s\n", imecList{iF}.fileStem);
                    switch mode
                        case 'ap'
                            mm = imecList{iF}.memmapAP_full();
                        case 'lf'
                            mm = imecList{iF}.memmapLF_full();
                        case 'sync'
                            mm = imecList{iF}.memmapSync_full();
                    end

                    % build idx vector
                    if isempty(timeShifts)
                        outSize = size(mm.Data.x, 2);
                        sourceIdxList = uint64(1):uint64(outSize);
                    else
                        sourceIdxList = timeShifts(iF).computeSourceIndices();
                        outSize = numel(sourceIdxList);
                    end
                    
                    nChunks = ceil(outSize / chunkSize);
                    prog = ProgressBar(nChunks, 'Copying %s file %d / %d: %s', mode, iF, nFiles, imecList{iF}.fileStem);
                    
                    for iCh = 1:nChunks
                        if iCh == nChunks
                            idx = (iCh-1)*(chunkSize)+1 : outSize;
                        else
                            idx = (iCh-1)*(chunkSize) + (1:chunkSize);
                        end

                        source_idx = sourceIdxList(idx);
                        data = mm.Data.x(chIdx, source_idx);

                        % ch_connected_mask indicates which channels are
                        % connected, which are the ones where scaling makes
                        % sense. chIdx is all channels being written to
                        % output file
                        ch_conn_mask = ismember(chIdx, imecList{iF}.connectedChannels);

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
                                data = fn(imecList{iF}, data, chIdx, source_idx);
                            end

                            if useGpuArray
                                data = gather(data);
                            end

                            if applyScaling
                                data(ch_conn_mask, :) = data(ch_conn_mask, :) ./ imecList{iF}.scaleToUv;
                            end

                            data = int16(data);
                        end

                        fwrite(fidOut, data, 'int16');
                        prog.increment();
                    end
                    prog.finish();
                end
            end
        end
    end
    
    methods
        function [chunkedOutputs, chunkSampleStartIdx] = testTransformInternal(imec, procFnList, mode, varargin)     
            p = inputParser();
            p.addParameter('timeIdxToInclude', [], @(x) isvector(x) || isempty(x));
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('writeSyncSeparate', false, @islogical); % true means ap will get only mapped channels, false will preserve channels as is
            p.addParameter('chunkSize', 2^20, @isscalar);
            p.addParameter('gpuArray', false, @islogical);
            p.addParameter('applyScaling', false, @islogical); % convert to uV before processing
            p.parse(varargin{:});

            if ~iscell(procFnList)
                procFnList = {procFnList};
            end
            chunkSize = p.Results.chunkSize;
            useGpuArray = p.Results.gpuArray;
            applyScaling = p.Results.applyScaling;

            % determine channels to pass along
            if p.Results.goodChannelsOnly
                goodMat = false(imec.nChannels, 1);
                goodMat(imec.goodChannels) = true;
                chIdx = find(all(goodMat, 2));
                assert(~isempty(chIdx), 'No goodChannels specified across all datasets')

            elseif p.Results.writeSyncSeparate
                chIdx = imec.mappedChannels; % excludes sync channel
                assert(~isempty(chIdx), 'No mapped channels found in first dataset');

            else
                chIdx = 1:imec.nChannels;
            end

            switch mode
                case 'ap'
                    mm = imec.memmapAP_full();
                case 'lf'
                    mm = imec.memmapLF_full();
                case 'sync'
                    mm = imec.memmapSync_full();
            end

            nChunks = ceil(size(mm.Data.x, 2) / chunkSize);

            chunkSampleStartIdx = 1 + chunkSize*(0:nChunks-1);
            chunkSampleEndIdx = min(chunkSampleStartIdx + chunkSize - 1, imec.nSamplesAP);

            timeIdxToInclude = p.Results.timeIdxToInclude;
            if isempty(timeIdxToInclude)
                maskRun = true(nChunks, 1);
            else
                maskRun = false(nChunks, 1);
                for iT = 1:numel(timeIdxToInclude)
                    maskRun(chunkSampleStartIdx <= timeIdxToInclude(iT) & timeIdxToInclude(iT) <= chunkSampleEndIdx) = true;
                end
            end
            
            if ~any(maskRun)
                error('timeIdxToInclude is out of range');
            end
            
            chunkedOutputs = cell(nChunks, 1);

            prog = ProgressBar(nChunks, 'Testing transform pipeline on %d / %d chunks', nnz(maskRun), nChunks);
            for iCh = 1:nChunks
                prog.update(iCh);
                if ~maskRun(iCh), continue, end
                
                if iCh == nChunks
                    idx = (iCh-1)*(chunkSize)+1 : size(mm.Data.x, 2);
                else
                    idx = (iCh-1)*(chunkSize) + (1:chunkSize);
                end

                data = mm.Data.x(chIdx, idx);

                % ch_connected_mask indicates which channels are
                % connected, which are the ones where scaling makes
                % sense. chIdx is all channels being written to
                % output file
                ch_conn_mask = ismember(chIdx, imec.connectedChannels);

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
                    data = fn(imec, data, chIdx, idx);
                end

                if useGpuArray
                    data = gather(data);
                end

                if applyScaling
                    data(ch_conn_mask, :) = data(ch_conn_mask, :) ./ imec.scaleToUv;
                end

                chunkedOutputs{iCh} = int16(data);
                
            end
            prog.finish();
            
            chunkedOutputs = chunkedOutputs(maskRun);
            chunkSampleStartIdx = chunkSampleStartIdx(maskRun); 
        end
    end
end
