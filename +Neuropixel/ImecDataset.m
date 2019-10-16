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
        badChannels (:, 1) uint32

        syncBitNames string;
        
        concatenationInfoAP
        concatenationInfoLF
        
        sourceDatasets % optional list of concatenated ImecDatasets that sourced the data for this file (via concatenationInfoAP)
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
        function imec = ImecDataset(fileOrFileStem, varargin)
            p = inputParser();
            p.addParameter('channelMap', [], @(x) true);
            p.addParameter('syncBitNames', [], @(x) isempty(x) || isstring(x) || iscellstr(x));
            p.addParameter('sourceDatasets', [], @(x) true);
            p.parse(varargin{:})

            fileOrFileStem = char(fileOrFileStem);
            file = Neuropixel.ImecDataset.findImecFileInDir(fileOrFileStem, 'ap');
            if isempty(file)
                error('No AP Imec file found at or in %s', fileOrFileStem);
            end
            [imec.pathRoot, imec.fileStem, imec.fileTypeAP] = Neuropixel.ImecDataset.parseImecFileName(file);
            if exist(imec.pathAP, 'file')
                if ~exist(imec.pathAPMeta, 'file')
                    error('Could not find AP meta file %s', imec.pathAPMeta);
                end
                imec.readInfo();
            else
                error('Could not find AP bin file %s', imec.pathAP);
            end

            channelMapFile = p.Results.channelMap;
            if isempty(channelMapFile)
                channelMapFile = Neuropixel.Utils.getDefaultChannelMapFile(true);
            end
            imec.channelMap = Neuropixel.ChannelMap(channelMapFile);
            assert(imec.channelMap.nChannels <= imec.nChannels, 'Channel count is less than number of channels in channel map');

            if ~isempty(p.Results.syncBitNames)
                imec.setSyncBitNames(1:numel(p.Results.syncBitNames), p.Resuls.syncBitNames);
            end
            
            if ~isempty(p.Results.sourceDatasets)
                assert(isa(p.Results.sourceDatasets, 'Neuropixel.ImecDataset'));
                imec.sourceDatasets = p.Results.sourceDatasets;
            end
        end

        function readInfo(imec)
            meta = imec.readAPMeta();
            imec.nChannels = meta.nSavedChans;
            imec.fsAP = meta.imSampRate;
            imec.highPassFilterHz = meta.imHpFlt;
            imec.creationTime = datenum(meta.fileCreateTime, 'yyyy-mm-ddTHH:MM:SS');

            if imec.hasLF
                metaLF = imec.readLFMeta();
                imec.fsLF = metaLF.imSampRate;
            end

            % parse imroTable
            m = regexp(meta.imroTbl, '\(([\d, ]*)\)', 'tokens');
            gainVals = strsplit(m{2}{1}, ' ');
            imec.apGain = str2double(gainVals{4});
            imec.lfGain = str2double(gainVals{5});

            imec.apRange = [meta.imAiRangeMin meta.imAiRangeMax];
            imec.lfRange = [meta.imAiRangeMin meta.imAiRangeMax];

            % look at AP meta fields that might have been set by us
            if isfield(meta, 'badChannels') && ~isempty(meta.badChannels)
                imec.badChannels = union(imec.badChannels, meta.badChannels);
            end
            if isfield(meta, 'syncBitNames')
                imec.setSyncBitNames(1:numel(meta.syncBitNames), meta.syncBitNames);
            end

            if imec.hasAP
                fid = imec.openAPFile();
                fseek(fid, 0, 'eof');
                bytes = ftell(fid);
                fclose(fid);

                imec.nSamplesAP = bytes / imec.bytesPerSample / imec.nChannels;
                assert(round(imec.nSamplesAP) == imec.nSamplesAP, 'AP bin file size is not an integral number of samples');
            end
            
            imec.concatenationInfoAP = Neuropixel.ConcatenationInfo(imec, 'ap', meta);

            if imec.hasLF
                fid = imec.openLFFile();
                fseek(fid, 0, 'eof');
                bytes = ftell(fid);
                fclose(fid);
                imec.nSamplesLF = bytes / imec.bytesPerSample / imec.nChannels;
                assert(round(imec.nSamplesAP) == imec.nSamplesAP, 'LF bin file size is not an integral number of samples');
              
                imec.concatenationInfoLF = Neuropixel.ConcatenationInfo(imec, 'lf', meta);
            end
        end

        function setSyncBitNames(imec, idx, names)
            % idx is the indices of which bits to set to the corresponding items from names
            assert(all(idx >= 1 & idx <= imec.nSyncBits), 'Sync bit indices must be in [1 %d]', imec.nSyncBits);
            if isscalar(idx) && ischar(names)
                imec.syncBitNames{idx} = names;
            else
                names = string(names);
                imec.syncBitNames(idx) = names;
            end
        end
        
        function setSourceDatasets(imec, imecList)
            assert(isempty(imecList) || isa(imecList, 'Neuropixel.ImecDataset'));
            imec.sourceDatasets = imecList;
        end

        function idx = lookupSyncBitByName(imec, names, ignoreNotFound)
            if nargin < 3
                ignoreNotFound = false;
            end
            if isnumeric(names)
                idx = names;
            else
                names = string(names);
                [tf, idx] = ismember(names, imec.syncBitNames);
                if ignoreNotFound
                    idx(~tf) = NaN;
                elseif any(~tf)
                    error('Sync bit(s) %s not found', strjoin(names, ', '));
                end
            end
        end

        function newImec = copyToNewLocation(imec, newRoot, newStem)
            if nargin < 3
                newStem = imec.fileStem;
            end
            Neuropixel.Utils.mkdirRecursive(newRoot);

            f = @(suffix) fullfile(newRoot, [newStem suffix]);
            docopy(imec.pathAP, f('.imec.ap.bin'));
            docopy(imec.pathAPMeta, f('.imec.ap.meta'));
            docopy(imec.pathLF, f('.imec.lf.bin'));
            docopy(imec.pathLFMeta, f('.imec.lf.meta'));
            docopy(imec.pathSync, f('.imec.sync.bin'));

            newImec = Neuropixel.ImecDataset(fullfile(newRoot, newStem), 'channelMap', imec.channelMapFile);

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
%         function data_ch_by_time = readAPChannelBand(imec, chFirst, chLast, sampleFirst, sampleLast, msg)
%             if nargin < 4 || isempty(sampleFirst)
%                 sampleFirst = 1;
%             end
%             if nargin < 5 || isempty(sampleLast)
%                 sampleLast = imec.nSamplesAP;
%             end
%             if nargin < 6 || isempty(msg)
%                 msg = 'Reading channels from neuropixel AP file';
%             end
% 
%             data_ch_by_time = imec.readChannelBand('ap', chFirst, chLast, sampleFirst, sampleLast, msg);
%         end
% 
%         function data_ch_by_time = readLFChannelBand(imec, chFirst, chLast, sampleFirst, sampleLast, msg)
%             if nargin < 4 || isempty(sampleFirst)
%                 sampleFirst = 1;
%             end
%             if nargin < 5 || isempty(sampleLast)
%                 sampleLast = imec.nSamplesLF;
%             end
%             if nargin < 6 || isempty(msg)
%                 msg = 'Reading channels from neuropixel LF file';
%             end
% 
%             data_ch_by_time = imec.readChannelBand('lf', chFirst, chLast, sampleFirst, sampleLast, msg);
%         end
% 
%         function data_by_time = readAPSingleChannel(imec, ch, varargin)
%             data_by_time = imec.readAPChannelBand(ch, ch, varargin{:})';
%         end
% 
%         function data_by_time = readLFSingleChannel(imec, ch, varargin)
%             data_by_time = imec.readLFChannelBand(ch, ch, varargin{:})';
%         end
    end
    
    methods % Sync channel read / cache
        function syncRaw = readSync(imec, varargin)
            p = inputParser();
            p.addOptional('reload', false, @islogical); % if true, ignore cache in memory imec.syncRaw 
            p.addParameter('ignoreCached', false, @islogical); % if true, ignore cache on disk in imec.pathSyncCached
            p.parse(varargin{:});

            if isempty(imec.syncRaw) || p.Results.reload
                if exist(imec.pathSyncCached, 'file') && ~p.Results.ignoreCached
                    [~, f, e] = fileparts(imec.pathSyncCached);
                    fprintf('Loading sync from cached %s%s\n', f, e);
                    ld = load(imec.pathSyncCached);
                    imec.syncRaw = ld.sync;
                else
                    % this will automatically redirect to a separate sync file
                    % or to the ap file depending on .syncInAPFile
                    fprintf('Loading sync channel (this will take some time)...\n');
                    mm = imec.memmapSync_full();
                    imec.syncRaw = mm.Data.x(imec.syncChannelIndex, :)';

                    imec.saveSyncCached();
                end
            end
            syncRaw = imec.syncRaw;
        end

        function saveSyncCached(imec)
            sync = imec.readSync();
            save(imec.pathSyncCached, 'sync');
        end

        function updateSyncCached(imec, varargin)
            if exist(imec.pathSyncCached, 'file')
                sync = imec.readSync('reload', true, 'ignoreCached', true);
                save(imec.pathSyncCached, 'sync');
            end
        end

        function tf = readSyncBit(imec, bit)
            bit = imec.lookupSyncBitByName(bit);
            tf = logical(bitget(imec.readSync(), bit));
        end
        
        function vec = readSync_idx(imec, idx, varargin)
            p = inputParser();
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.parse(varargin{:});
            
            fromSource = p.Results.fromSourceDatasets;
            
            if ~fromSource
                if ~isempty(imec.syncRaw)
                    vec = imec.syncRaw(idx);
                else
                    mm = imec.memmapSync_full();
                    vec = mm.Data.x(imec.syncChannelIndex, idx)';
                end
            else
                mmSet = imec.memmap_sourceAP_full();
                [sourceFileInds, sourceSampleInds] = imec.concatenationInfoAP.lookup_sampleIndexInSourceFiles(idx);
                vec = Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSet, sourceFileInds, sourceSampleInds, imec.syncChannelIndex);
            end
        end
        
        function mat = readSyncBits_idx(imec, bits, idx, varargin)
            % mat is nBits x nTime to match readAP_idx which is nChannels x nTime
            if isstring(bits) || ischar(bits)
                bits = imec.lookupSyncBitByName(bits);
            end
            vec = imec.readSync_idx(idx, varargin{:});
            mat = false(numel(bits), numel(vec));
            for iB = 1:numel(bits)
                mat(iB, :) = logical(bitget(vec, bits(iB)));
            end
        end
    end
    
    methods
        function sampleIdx = closestSampleAPForTime(imec, timeSeconds)
            sampleIdx = round(timeSeconds * imec.fsAP);
            sampleIdx(sampleIdx == 0) = 1;
            if any(sampleIdx < 0 | sampleIdx > imec.nSamplesAP)
                error('Time seconds out of range');
            end 
        end
        
        function sampleIdx = closestSampleLFForTime(imec, timeSeconds)
            sampleIdx = round(timeSeconds * imec.fsLF);
            sampleIdx(sampleIdx == 0) = 1;
            if any(sampleIdx < 0 | sampleIdx > imec.nSamplesLF)
                error('Time seconds out of range');
            end 
        end
        
        function data = readAP_idx(imec, sampleIdx, varargin)
            p = inputParser();
            p.addParameter('applyScaling', true, @islogical); % convert to uV before processing
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.parse(varargin{:});
            
            fromSource = p.Results.fromSourceDatasets;
            
            if ~fromSource
                mm = imec.memmapAP_full();
            
                if any(isnan(sampleIdx))
                    mask = ~isnan(sampleIdx);
                    data = single(mm.Data.x(:, sampleIdx)); % must be single to support NaNs
                    data = Neuropixel.Utils.TensorUtils.inflateMaskedTensor(data, 2, mask, NaN);
                else
                    data = mm.Data.x(:, sampleIdx);
                end
            else
                mmSet = imec.memmap_sourceAP_full();
                [sourceFileInds, sourceSampleInds] = imec.concatenationInfoAP.lookup_sampleIndexInSourceFiles(sampleIdx);
                data = Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSet, sourceFileInds, sourceSampleInds);
            end

            if p.Results.applyScaling
                data = single(data);
                ch_conn_mask = imec.lookup_channelIds(imec.connectedChannels);
                data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(imec.apScaleToUv);
            end
        end
       
        function [mat, sampleIdx] = readAP_timeWindow(imec, timeWindowSec, varargin)
            idxWindow = imec.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = imec.readAP_idx(sampleIdx, varargin{:});
        end
        
        function [mat, sampleIdx] = readSyncBits_timeWindow(imec, bits, timeWindowSec)
            idxWindow = imec.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = readSyncBits_idx(bits, sampleIdx);
        end
        
        function [mat, sourceInds] = readAP_viaTimeShiftSpec(imec, timeShifts, varargin)
            % read a matrix of AP data as requested by Neuropixel.TimeShiftSpec instance timeShifts
            sourceInds = timeShifts.computeSourceIndices(imec.nSamplesAP);
            mat = imec.readAP_idx(sourceInds, varargin{:});
        end
    end
    
    methods % Quick inspection
        function [channelInds, channelIds] = lookup_channelIds(imec, channelIds)
             if islogical(channelIds)
                channelIds = imec.channelIdx(channelIds);
             end
            [tf, channelInds] = ismember(channelIds, imec.channelIds);
            assert(all(tf, 'all'), 'Some channel ids not found');
        end
        
        function inspectAP_timeWindow(imec, timeWindowSec, varargin)
            idxWindow = imec.closestSampleAPForTime(timeWindowSec);
            imec.inspectAP_idxWindow(idxWindow, 'timeInSeconds', true, varargin{:});
        end
        
        function inspectAP_idxWindow(imec, idxWindow, varargin)
            p = inputParser();
            p.addParameter('channels', imec.mappedChannels, @isvector);
            p.addParameter('invertChannels', false, @islogical);
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            p.addParameter('showSync', true, @isvector);
            p.addParameter('syncBits', imec.syncBitsNamed, @isvector);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 0.95, @isscalar);
            p.addParameter('car', false, @islogical);
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.addParameter('downsample',1, @isscalar); 
            p.addParameter('timeInSeconds', false, @islogical);
            p.addParameter('timeRelativeTo', 0, @isscalar);
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'Neuropixel.TrialSegmentationInfo')); % to mark trial boundaries
            
            p.addParameter('markSampleIdx', [], @isvector);
            p.addParameter('markSampleMode', 'rug', @ischar);
            p.addParameter('markSampleColor', [0.5 0 0.5], @(x) true);
            p.parse(varargin{:});
            
            if numel(idxWindow) > 2
                idxWindow = [idxWindow(1), idxWindow(end)];
            end
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = imec.readAP_idx(sampleIdx, 'fromSourceDatasets', p.Results.fromSourceDatasets); % C x T
            labels = imec.channelNamesPadded;
            
            [channelInds, channelIds] = imec.lookup_channelIds(p.Results.channels); %#ok<*PROPLC>
            if p.Results.goodChannelsOnly
                mask = ismember(channelIds, imec.goodChannels);
                channelInds = channelInds(mask);
                channelIds = channelIds(mask);
            end
            if p.Results.connectedChannelsOnly
                mask = ismember(channelIds, imec.connectedChannels);
                channelInds = channelInds(mask);
                channelIds = channelIds(mask);
            end
            
            mat = mat(channelInds, :);
            labels = labels(channelInds);
            connected = ismember(channelIds, imec.connectedChannels);
            bad = ismember(channelIds, imec.badChannels);
            
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
                syncBitMat = imec.readSyncBits_idx(syncBits, sampleIdx, 'fromSourceDatasets', p.Results.fromSourceDatasets);
                mat = cat(1, mat, syncBitMat);
                syncColor = [0.75 0 0.9];
                colors = cat(1, colors, repmat(syncColor, size(syncBitMat, 1), 1));
                labels = cat(1, labels, imec.syncBitNames(syncBits));
                normalizeMask = cat(1, normalizeMask, false(size(syncBitMat, 1), 1));
            end

            if ~p.Results.showLabels
                labels = [];
            end
            
            timeRelativeTo = p.Results.timeRelativeTo;
            time = int64(sampleIdx) - int64(timeRelativeTo);
            if p.Results.timeInSeconds
                time = double(time) / imec.fsAP;
            end
            Neuropixel.Utils.plotStackedTraces(time, mat', 'colors', colors, 'labels', labels, ...
                'gain', p.Results.gain, 'invertChannels', p.Results.invertChannels, 'normalizeMask', normalizeMask, 'normalizeEach', false);
            
            tsi = p.Results.tsi;
            if ~isempty(tsi)
                tsi.markTrialTicks('sample_window', idxWindow, 'timeInSeconds', p.Results.timeInSeconds, ...
                    'side', 'bottom', 'Color', [0.2 0.2 1], 'expand_limits', true);
            end
            
            markSampleIdx = p.Results.markSampleIdx;
            if ~isempty(markSampleIdx)
                mask = markSampleIdx >= idxWindow(1) & markSampleIdx <= idxWindow(2);
                markSampleIdx = markSampleIdx(mask);
                assert(numel(markSampleIdx) < 100, 'Too many times to mark');
                
                markSampleIdx = int64(markSampleIdx) - int64(timeRelativeTo);
                if p.Results.timeInSeconds
                    markTimes = double(markSampleIdx) / imec.fsAP;
                else
                    markTimes = markSampleIdx;
                end
                
                if strcmp(p.Results.markSampleMode, 'rug')
                    Neuropixel.Utils.rugplot(markTimes, 'side', 'top', 'Color', p.Results.markSampleColor, 'expand_limits', true);
                else
                    for iM = 1:numel(markTimes)
                        xline(markTimes(iM), 'Color', p.Results.markSampleColor);
                    end
                end 
            end
            
            hold off;
        end
    end

    methods % Memory mapped read/write access to data
        function mm = memmapAP_by_sample(imec)
            mm = memmapfile(imec.pathAP, 'Format', {'int16', [imec.nChannels 1], 'x'}, ...
               'Repeat', imec.nSamplesAP);
        end

        function mm = memmapLF_by_sample(imec)
            mm = memmapfile(imec.pathLF, 'Format', {'int16', [imec.nChannels 1], 'x'}, ...
               'Repeat', imec.nSamplesLF);
        end

        function mm = memmapAP_by_chunk(imec, nSamplesPerChunk)
            mm = memmapfile(imec.pathAP, 'Format', {'int16', [imec.nChannels nSamplesPerChunk], 'x'}, ...
               'Repeat', floor(imec.nSamplesAP/nSamplesPerChunk));
        end

        function mm = memmapLF_by_chunk(imec, nSamplesPerChunk)
            mm = memmapfile(imec.pathLF, 'Format', {'int16', [imec.nChannels nSamplesPerChunk], 'x'}, ...
               'Repeat', floor(imec.nSamplesLF/nSamplesPerChunk));
        end

        function mm = memmapAP_full(imec, varargin)
            p = inputParser();
            p.addParameter('Writable', false, @islogical);
            p.parse(varargin{:});

            mm = memmapfile(imec.pathAP, 'Format', {'int16', [imec.nChannels imec.nSamplesAP], 'x'}, 'Writable', p.Results.Writable);
        end

        function mm = memmapLF_full(imec, varargin)
            p = inputParser();
            p.addParameter('Writable', false, @islogical);
            p.parse(varargin{:});

            mm = memmapfile(imec.pathLF, 'Format', {'int16', [imec.nChannels imec.nSamplesLF], 'x'}, 'Writable', p.Results.Writable);
        end

        function mm = memmapSync_full(imec)
            if imec.syncInAPFile
                % still has nChannels
                mm = memmapfile(imec.pathSync, 'Format', {'int16', [imec.nChannels imec.nSamplesAP], 'x'});
            else
                % only sync channel
                mm = memmapfile(imec.pathSync, 'Format', {'int16', [1 imec.nSamplesAP], 'x'});
            end
        end
        
        % refer back to the source datasets
        function mmSet = memmap_sourceAP_full(imec, varargin)
            assert(~isempty(imec.sourceDatasets));
            nSources = numel(imec.sourceDatasets);
            mmSet = cell(nSources, 1);
            for iF = 1:nSources
                mmSet{iF} = imec.sourceDatasets(iF).memmapAP_full(varargin{:});
            end
        end
        
        % refer back to the source datasets
        function mmSet = memmap_sourceLF_full(imec, varargin)
            assert(~isempty(imec.sourceDatasets));
            nSources = numel(imec.sourceDatasets);
            mmSet = cell(nSources, 1);
            for iF = 1:nSources
                mmSet{iF} = imec.sourceDatasets(iF).memmapLF_full(varargin{:});
            end
        end
    end
    
    methods(Static)
        function out = multi_mmap_extract_sample_idx(mmSet, fileInds, origSampleInds, chInds)
            % given [fileInds, origSampleInds] as returned by ConcatenationInfo/lookup_sampleIndexInSourceFiles
            % extract those samples from the set of memory mapped files in mmSet (returned by memmap**_all)
            
            
            assert(numel(fileInds) == numel(origSampleInds));
            nSamplesOut = numel(origSampleInds);
            cls = mmSet{1}.Format{1};
            nFiles = numel(mmSet);
            
            if nargin < 4
                nCh = mmSet{1}.Format{2}(1);
                chInds = 1:nCh;
            else
                if islogical(chInds)
                    nCh = nnz(chInds);
                else
                    nCh = numel(chInds);
                end
            end
            out = zeros(nCh, nSamplesOut, cls);
            
            for iF = 1:nFiles
                mask = fileInds == iF;
                out(:, mask) = mmSet{iF}.Data.x(chInds, origSampleInds(mask));
            end
        end
    end

    methods(Hidden) % Read data at specified times
        function [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet] = readSnippetsRaw(imec, times, window, varargin)
            % for each sample index in times, read the window times + window(1):window(2)
            % of samples around this time from some channels

            p = inputParser();
            p.addParameter('source', 'ap', @ischar);
            p.addParameter('fromSourceDatasets', false, @islogical);
            
            % specify one of THESE (same channels every snippet or nChannels x numel(times))
            p.addParameter('channel_ids', [], @(x) isempty(x) || isvector(x));
            p.addParameter('channel_ids_by_snippet', [], @ismatrix);
            
            % or THESE (different channels for each group of snippets sharing the same cluster_ids
            p.addParameter('channel_ids_by_cluster', [], @ismatrix);
            p.addParameter('unique_cluster_ids', [], @isvector);
            p.addParameter('cluster_ids_by_snippet', [], @isvector); % same length as numel(times), corresponding to which cluster was pulled out (e.g. as a waveform)
            
            p.addParameter('car', false, @islogical); % subtract median over channels
            p.parse(varargin{:});

            channel_ids = p.Results.channel_ids;
            channel_ids_by_snippet = p.Results.channel_ids_by_snippet;
            channel_ids_by_cluster = p.Results.channel_ids_by_cluster;
            cluster_ids = p.Results.cluster_ids_by_snippet;
            unique_cluster_ids = p.Results.unique_cluster_ids;
            if isempty(unique_cluster_ids), unique_cluster_ids = unique(cluster_ids); end
                
            if ~isempty(channel_ids_by_snippet)
                assert(size(channel_ids_by_snippet, 2) == numel(times), 'channel_ids must be nChannels x numel(times)');
                
            elseif ~isempty(channel_ids_by_cluster)
                assert(~isempty(cluster_ids), 'cluster_ids_by_snippet must be specified when channel_ids_by_cluster is used');
                [~, cluster_inds] = ismember(cluster_ids, unique_cluster_ids);
                channel_ids_by_snippet = channel_ids_by_cluster(:, cluster_inds);
                
            elseif ~isempty(channel_ids)
                channel_ids = Neuropixel.Utils.makecol(channel_ids);
                channel_ids_by_snippet = repmat(channel_ids, 1, numel(times));
                
            else
                error('Specify either channel_ids or channel_ids_by_cluster');
            end
            
            channel_inds_by_snippet = imec.lookup_channelIds(channel_ids_by_snippet);

            source = p.Results.source;
            fromSourceDatasets = p.Results.fromSourceDatasets;
            switch source
                case 'ap'
                    nSamples = imec.nSamplesAP;
                    if ~fromSourceDatasets
                        mm = imec.memmapAP_full();
                    else
                        mmSet = imec.memmap_sourceAP_full();
                        concatInfo = imec.concatenationInfoAP;
                    end
                case 'lf'
                    nSamples = imec.nSamplesLF;
                    if ~fromSourceDatasets
                        mm = imec.memmapLF_full();
                    else
                        mmSet = imec.memmap_sourceLF_full();
                        concatInfo = imec.concatenationInfoLF;
                    end 
                otherwise
                    error('Unknown source');
            end
            times = Neuropixel.Utils.makecol(uint64(times));
            nC = size(channel_ids_by_snippet, 1);
            nC_all = imec.nChannels;
            nS = numel(times);
            nT = numel(window(1):window(2));
            out = zeros(nC, nT, nS, 'int16');

            if numel(times) > 10
                if exist('ProgressBar', 'class') == 8
                    prog = ProgressBar(numel(times), 'Extracting %s snippets', upper(source));
                else
                    prog = Neuropixel.Utils.ProgressBar(numel(times), 'Extracting %s snippets', upper(source));
                end
            else
                prog = [];
            end
            
            good_ch_inds = imec.goodChannelInds;
            idx_request = int64(times') + int64(window(1):window(2))'; % nWindow x nTimes indices
            mask_idx_okay = idx_request >= int64(1) & idx_request <= nSamples;
            idx_request(~mask_idx_okay) = 1; % we'll clear out later
                
%             allAtOnce = true;
%             % two versions of the algorithm, one loops over snippets, one does all indexing at one time
%             if ~allAtOnce
%                 if ~fromSourceDatasets
%                     for iS = 1:numel(times)
%                         if mod(iS, 20) == 0 && ~isempty(prog), prog.update(iS); end
% 
%                         extract_all_ch = mm.Data.x(:, idx_request(:, iS));
%                         if p.Results.car
%                             ar = median(extract_all_ch(good_ch_inds, :), 1);
%                             out(:, :, iS) = extract_all_ch(channel_inds_by_snippet(:, iS), :) - ar; % which channels for this spike
%                         else
%                             out(:, :, iS) = extract_all_ch(channel_inds_by_snippet(:, iS), :);
%                         end
%                     end
% 
%                 else
%                     [sourceFileInds, sourceSampleInds] = concatInfo.lookup_sampleIndexInSourceFiles(idx_request);
%                     for iS = 1:numel(times)
%                         if mod(iS, 20) == 0 && ~isempty(prog), prog.update(iS); end
%                         extract_all_ch = Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSet, sourceFileInds(:, iS), sourceSampleInds(:, iS));
%                         if p.Results.car
%                             ar = median(extract_all_ch(good_ch_inds, :), 1);
%                             out(:, :, iS) = extract_all_ch(channel_inds_by_snippet(:, iS), :) - ar; % which channels for this spike
%                         else
%                             out(:, :, iS) = extract_all_ch(channel_inds_by_snippet(:, iS), :);
%                         end
%                     end
%                 end
%             else
%                 if ~fromSourceDatasets
%                     extract_all_ch = reshape(mm.Data.x(:, idx_request(:)), [nC_all nT, nS]);
%                 else
%                     [sourceFileInds, sourceSampleInds] = concatInfo.lookup_sampleIndexInSourceFiles(idx_request);
%                     extract_all_ch = reshape(Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSet, ...
%                         sourceFileInds(:), sourceSampleInds(:)), [nC_all nT nS]);
%                 end
%                 if p.Results.car
%                     ar = median(extract_all_ch(good_ch_inds, :, :), 1);
%                 else
%                     ar = zeros([1 nT nS], 'like', extract_all_ch);
%                 end
%                 for iS = 1:numel(times)
%                     if mod(iS, 20) == 0 && ~isempty(prog), prog.update(iS); end
%                     out(:, :, iS) = extract_all_ch(channel_inds_by_snippet(:, iS), :, iS) - ar(1, :, iS); % which channels for this spike
%                 end
%             end
            
            nPerGroup = 50; 
            nGroups = ceil(nS / nPerGroup);
            for iGroup = 1:nGroups
                idxS = nPerGroup*(iGroup-1) + 1 : min(nS,  nPerGroup*iGroup);
                nS_this = numel(idxS);
                idx_request_this = idx_request(:, idxS);
                if ~fromSourceDatasets
                    extract_all_ch = reshape(mm.Data.x(:, idx_request_this(:)), [nC_all nT, nS_this]);
                else
                    [sourceFileInds, sourceSampleInds] = concatInfo.lookup_sampleIndexInSourceFiles(idx_request_this);
                    extract_all_ch = reshape(Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSet, ...
                        sourceFileInds(:), sourceSampleInds(:)), [nC_all nT nS_this]);
                end
                if p.Results.car
                    ar = median(extract_all_ch(good_ch_inds, :, :), 1);
                else
                    ar = zeros([1 nT nS], 'like', extract_all_ch);
                end
                for iiS = 1:nS_this
                    out(:, :, idxS(iiS)) = extract_all_ch(channel_inds_by_snippet(:, idxS(iiS)), :, iiS) - ar(1, :, iiS); % which channels for this spike
                end
                if ~isempty(prog), prog.increment(nPerGroup); end
            end
            
            out(:, ~mask_idx_okay(:)) = 0;
            data_ch_by_time_by_snippet = out;
            if ~isempty(prog), prog.finish(); end
        end
    end
    
    methods  % Read data at specified times
        function snippet_set = readAPSnippetSet(imec, times, window, varargin)
            [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet] = ...
                imec.readSnippetsRaw(times, window, 'source', 'ap', varargin{:});
            snippet_set = Neuropixel.SnippetSet(imec, 'ap');
            snippet_set.data = data_ch_by_time_by_snippet;
            snippet_set.sample_idx = times;
            snippet_set.channel_ids_by_snippet = channel_ids_by_snippet;
            snippet_set.cluster_ids = cluster_ids;
            snippet_set.window = window;
        end

        function snippet_set = readLFSnippetSet(imec, times, window, varargin)
            [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet] = ...
                imec.readSnippetsRaw(times, window, 'source', 'lf', varargin{:});
            snippet_set = Neuropixel.SnippetSet(imec, 'lf');
            snippet_set.data = data_ch_by_time_by_snippet;
            snippet_set.sample_idx = times;
            snippet_set.channel_ids_by_snippet = channel_ids_by_snippet;
            snippet_set.cluster_ids = cluster_ids;
            snippet_set.window = window;
        end

        function rms = computeRMSByChannel(imec, varargin)
            % output will be nMappedChannels x 1 vector of rms
            p = inputParser();
            p.addParameter('sampleMaskFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % sampleMaskFn(data_ch_x_time, sample_idx_time) --> logical_time mask of time samples valid for use, useful if you have artifacts at known times
            p.addParameter('car', false, @islogical);
            p.addParameter('useChunks', 50, @isscalar);
            p.addParameter('chunkSize', 100000, @isscalar);
            p.parse(varargin{:});
            
            sampleMaskFn = p.Results.sampleMaskFn;
            
            % aim for the middle of the file
            chunkSize = min(imec.fsAP, p.Results.chunkSize);
            mm = imec.memmapAP_by_chunk(chunkSize);
            nChunks = numel(mm.Data);
            useChunks = min(nChunks, p.Results.useChunks);
            skipChunks = floor((nChunks-useChunks)/2);
            
            ch_mask = imec.lookup_channelIds(imec.mappedChannels); % for common average referencing

            sumByChunk = nan(imec.nChannels, useChunks);
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
            
            rms = rms(ch_mask); % only return mapped channels
            rms = rms * imec.apScaleToUv;
        end
    end

    methods(Hidden)
        function fid = openAPFile(imec)
            if ~exist(imec.pathAP, 'file')
                error('RawDataFile: %s not found', imec.pathAP);
            end
            fid = fopen(imec.pathAP, 'r');

            if fid == -1
                 error('RawDataFile: Could not open %s', imec.pathAP);
            end
        end

        function fid = openLFFile(imec)
            if ~exist(imec.pathLF, 'file')
                error('RawDataFile: %s not found', imec.pathAP);
            end
            fid = fopen(imec.pathLF, 'r');

            if fid == -1
                 error('RawDataFile: Could not open %s', imec.pathAP);
            end
        end

        function fid = openSyncFile(imec)
            if ~exist(imec.pathSync, 'file')
                error('RawDataFile: %s not found', imec.pathSync);
            end
            fid = fopen(imec.pathSync, 'r');

            if fid == -1
                 error('RawDataFile: Could not open %s', imec.pathSync);
            end
        end
    end

    methods % Dependent properties
        function tf = get.syncInAPFile(imec)
            tf = imec.channelMap.syncInAPFile;
        end
        
        function ind = get.syncChannelIndex(imec)
            if imec.syncInAPFile
                ind = imec.channelMap.syncChannelIndex;
            else
                % if sync is in its own file, assume it's the first and only channel
                ind = uint32(1);
            end
        end
        
        function pathAP = get.pathAP(imec)
            pathAP = fullfile(imec.pathRoot, imec.fileAP);
        end

        function fileAP = get.fileAP(imec)
            fileAP = [imec.fileStem '.imec.' imec.fileTypeAP '.bin'];
        end

        function tf = get.hasAP(imec)
            tf = exist(imec.pathAP, 'file') == 2;
        end

        function fileAPMeta = get.fileAPMeta(imec)
            fileAPMeta = [imec.fileStem '.imec.ap.meta'];
        end

        function pathAPMeta = get.pathAPMeta(imec)
            pathAPMeta = fullfile(imec.pathRoot, imec.fileAPMeta);
        end

        function pathLF = get.pathLF(imec)
            pathLF = fullfile(imec.pathRoot, imec.fileLF);
        end

        function fileAP = get.fileLF(imec)
            fileAP = [imec.fileStem '.imec.lf.bin'];
        end

        function fileLFMeta = get.fileLFMeta(imec)
            fileLFMeta = [imec.fileStem '.imec.lf.meta'];
        end

        function pathLFMeta = get.pathLFMeta(imec)
            pathLFMeta = fullfile(imec.pathRoot, imec.fileLFMeta);
        end

        function tf = get.hasLF(imec)
            tf = exist(imec.pathLF, 'file') == 2;
        end

        function fileSync = get.fileSync(imec)
            if imec.syncInAPFile
                fileSync = imec.fileAP;
            else
                fileSync = [imec.fileStem, '.imec.sync.bin'];
            end
        end

        function pathSync = get.pathSync(imec)
            pathSync = fullfile(imec.pathRoot, imec.fileSync);
        end

        function fileSyncCached = get.fileSyncCached(imec)
            fileSyncCached = [imec.fileStem '.sync.mat'];
        end

        function pathSyncCached = get.pathSyncCached(imec)
            pathSyncCached = fullfile(imec.pathRoot, imec.fileSyncCached);
        end

        function scale = get.apScaleToUv(imec)
            scale = (imec.apRange(2) - imec.apRange(1)) / (2^imec.adcBits) / imec.apGain * 1e6;
        end

        function scale = get.lfScaleToUv(imec)
            scale = (imec.apRange(2) - imec.apRange(1)) / (2^imec.adcBits) / imec.apGain * 1e6;
        end

        function file = get.channelMapFile(imec)
            if isempty(imec.channelMap)
                file = '';
            else
                file = imec.channelMap.file;
            end
        end

        function list = get.mappedChannels(imec)
            if isempty(imec.channelMap)
                list = [];
            else
                list = imec.channelMap.channelIdsMapped;
            end
        end
        
        function list = get.mappedChannelInds(imec)
            list = imec.lookup_channelIds(imec.mappedChannels);
        end

        function list = get.connectedChannels(imec)
            if isempty(imec.channelMap)
                list = [];
            else
                list = imec.channelMap.connectedChannels;
            end
        end
        
        function list = get.connectedChannelInds(imec)
            list = imec.lookup_channelIds(imec.connectedChannels);
        end

        function n = get.nChannelsMapped(imec)
            if isempty(imec.channelMap)
                n = NaN;
            else
                n = imec.channelMap.nChannelsMapped;
            end
        end

        function n = get.nChannelsConnected(imec)
            if isempty(imec.channelMap)
                n = NaN;
            else
                n = nnz(imec.channelMap.connected);
            end
        end

        function ch = get.goodChannels(imec)
            ch = setdiff(imec.connectedChannels, imec.badChannels);
        end
        
        function list = get.goodChannelInds(imec)
            list = imec.lookup_channelIds(imec.goodChannels);
        end

        function n = get.nGoodChannels(imec)
            n = numel(imec.goodChannels);
        end
        
        function idx = get.channelIds(imec)
            idx = imec.channelMap.channelIds;
        end        
        
        function names = get.channelNames(imec)
            names = strings(imec.nChannels, 1);
            names(imec.channelMap.channelIds) = string(sprintfc("ch %d", imec.channelMap.channelIds));
            if ~isnan(imec.syncChannelIndex)
                names(imec.syncChannelIndex) = "sync";
            end
        end

        function names = get.channelNamesPadded(imec)
            names = strings(imec.nChannels, 1);
            names(imec.channelMap.channelIds) = string(sprintfc("ch %03d", imec.channelMap.channelIds));
            if ~isnan(imec.syncChannelIndex)
                names(imec.syncChannelIndex) = "sync";
            end
        end
        
        function n = get.nSyncBits(imec)
            n = 8*imec.bytesPerSample; % should be 16?
        end
        
        function bits = get.syncBitsNamed(imec)
            names = imec.syncBitNames;
            bits = find(names ~= "");
        end

        function names = get.syncBitNames(imec)
            if isempty(imec.syncBitNames)
                names = strings(imec.nSyncBits, 1);
            else
                names = string(imec.syncBitNames);
            end
        end

        function meta = readAPMeta(imec)
            meta = Neuropixel.readINI(imec.pathAPMeta);
        end

        function meta = generateModifiedAPMeta(imec)
            meta = imec.readAPMeta;

            meta.syncBitNames = imec.syncBitNames;
            meta.badChannels = imec.badChannels;
        end

        function writeModifiedAPMeta(imec, varargin)
            p = inputParser();
            p.addParameter('extraMeta', struct(), @isstruct);
            p.parse(varargin{:});

            meta = imec.generateModifiedAPMeta();

            % set extra user provided fields
            extraMeta = p.Results.extraMeta;
            extraMetaFields = fieldnames(extraMeta);
            for iFld = 1:numel(extraMetaFields)
                meta.(extraMetaFields{iFld}) = extraMeta.(extraMetaFields{iFld});
            end

            Neuropixel.writeINI([imec.pathAPMeta], meta);
        end

        function meta = readLFMeta(imec)
            meta = Neuropixel.readINI(imec.pathLFMeta);
        end

        function str = get.creationTimeStr(imec)
            str = datestr(imec.creationTime);
        end
        
        function file = getAuxiliaryFileWithSuffix(imec, suffix)
             suffix = char(suffix);
             file = fullfile(imec.pathRoot, [imec.fileStem, '.', suffix]);
        end
    end
    
    methods % Marking Channels as bad
            function [rmsBadChannels, rmsByChannel] = markBadChannelsByRMS(imec, varargin)
            p = inputParser();
            p.addParameter('rmsRange', [3 100], @isvector);
            p.addParameter('sampleMaskFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % sampleMaskFn(data_ch_x_time, sample_idx_time) --> logical_time mask of time samples valid for use, useful if you have artifacts at known times
            p.parse(varargin{:});

            rmsByChannel = imec.computeRMSByChannel('sampleMaskFn', p.Results.sampleMaskFn);
            rmsMin = p.Results.rmsRange(1);
            rmsMax = p.Results.rmsRange(2);
            rmsBadMask = rmsByChannel < rmsMin | rmsByChannel > rmsMax;
            
            badMappedChannels = imec.mappedChannels(rmsBadMask);
            badConnectedChannels = badMappedChannels(ismember(badMappedChannels, imec.connectedChannels));
            imec.markBadChannels(badConnectedChannels);
            
            rmsBadChannels = badMappedChannels;
        end

        function markBadChannels(imec, list)
            % this adds to the set of bad channels, so multiple calls will
            % remove additional channels
            if islogical(list)
                list = find(list);
            end
            % filter for connected channels only
            badConnectedChannels = list(ismember(list, imec.connectedChannels));
            imec.badChannels = union(imec.badChannels, badConnectedChannels);
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
        
        function imecOut = saveTranformedDataset(imec, outPath, varargin)
            p = inputParser();
            p.addParameter('transformAP', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (imec, dataChunk) and return dataChunk someplace
            p.addParameter('transformLF', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (imec, dataChunk) and return dataChunk someplace

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
            imecOut = Neuropixel.ImecDataset.writeConcatenatedFileMatchGains({imec}, outPath, p.Results);
         end
        
        function writeFolderForPhy(imec, savePath, varargin)
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
            nCh = imec.nGoodChannels;
            % these can't be 1 because of squeeze() inside Phy's code
            nTemplates = 2;
            nTemplateFeatures = 2;
            nFeaturesPerChannel = 2;
            nPCFeatures = 2;
            
            imec.symLinkAPIntoDirectory(savePath);
            
            writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
            writeNPY(zeros(nSpikes, 1, 'uint32'), fullfile(savePath, 'spike_templates.npy'));
            writeNPY(spikeClusters, fullfile(savePath, 'spike_clusters.npy'));
    
            writeNPY(zeros(nSpikes, 1, 'double'), fullfile(savePath, 'amplitudes.npy'));
            
            templates = zeros(2, nTemplateTimepoints, nCh, 'single');
            writeNPY(templates, fullfile(savePath, 'templates.npy'));
    
            templatesInds = imec.goodChannels';
            writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
            sortedInds = imec.goodChannelInds;
            chanMap0ind = int32(imec.channelMap.channelIdsMapped(sortedInds) - uint32(1));
            xcoords = imec.channelMap.xcoords(sortedInds);
            ycoords = imec.channelMap.ycoords(sortedInds);
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
            [~, fname, ext] = fileparts(imec.fileAP);
            fprintf(fid,['dat_path = ''',fname ext '''\n']);
            fprintf(fid,'n_channels_dat = %i\n',imec.nChannels);
            fprintf(fid,'dtype = ''int16''\n');
            fprintf(fid,'offset = 0\n');
            fprintf(fid,'sample_rate = %i\n', imec.fsAP);
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
        
        function transformExtraArg = modifyInPlaceInternal(imec, mode, procFnList, varargin)
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
            
            p.addParameter('transformExtraArg', [], @(x) true);
            p.parse(varargin{:});

            chunkSize = p.Results.chunkSize;
            chunkExtra = p.Results.chunkEdgeExtraSamples;
            useGpuArray = p.Results.gpuArray;
            applyScaling = p.Results.applyScaling;
            dryRun = p.Results.dryRun;
            dryRunSampleInds = p.Results.dryRunSampleInds;
            transformExtraArg = p.Results.transformExtraArg;

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
                [source_idx, keepIdx] = Neuropixel.ImecDataset.determineChunkIdx(dataSize, iCh, nChunks, chunkSize, chunkExtra);
                
                if dryRun && ~isempty(dryRunSampleInds) && ~any(ismember(source_idx, dryRunSampleInds))
                    % skip this chunk unless some ind in idx should be processed
                    continue;
                end

                data = mm.Data.x(chInds, source_idx);

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
                    
                    % pass as many inputs to fn as it can handle, including user provided args at the end
                    extraArgs = {chIds, source_idx, transformExtraArg};
                    ninputs = nargin(fn);
                    if ninputs > 0 && ninputs < 2 + numel(extraArgs)
                        extraArgs = extraArgs(1:ninputs-2);
                    end

                    noutputs = nargout(fn);
                    if noutputs > 1
                        [data, transformExtraArg] = fn(imec, data, extraArgs{:});
                    else
                        data = fn(imec, data, extraArgs{:});
                    end
                                
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
                    mm.Data.x(chInds, source_idx) = data;
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
        
        function [tf, candidates] = pathPointsToSingleValidDataset(fileOrFileStem, type)
            if nargin < 2
                type = 'ap';
            end
            candidates = Neuropixel.ImecDataset.findImecFileInDir(fileOrFileStem, type, true, false);
            tf = numel(candidates) == 1;
        end

        function file = findImecFileInDir(fileOrFileStem, type, returnMultiple, errorIfNotFound)
            if nargin < 2
                type = 'ap';
            end
            if nargin < 3
                returnMultiple = false;
            end
            if nargin < 4
                errorIfNotFound = true;
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
                                if returnMultiple
                                    file = apFiles;
                                else
                                    file = apFiles{1};
                                    warning('Multiple AP files found in dir %s, choosing %s', path, file);
                                end
                            else
                                file = apFiles;
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
                                if returnMultiple
                                    file = lfFiles{1};
                                    warning('Multiple LF files found in dir, choosing %s', file);
                                else
                                    file = lfFiles;
                                end
                            else
                                file = lfFiles;
                            end
                        else
                            file = [];
                            return;
                        end
                    otherwise
                        error('Unknown type %s');
                end

                for iF = 1:numel(file)
                    file{iF} = fullfile(path, file{iF});
                end
                
            else
                % not a folder or a file, but possibly pointing to the
                % stem of a file, e.g. '/path/data' pointing to
                % '/path/data.ap.imec.bin'
                [parent, leaf, ext] = fileparts(fileOrFileStem);
                if ~exist(parent, 'dir')
                    if errorIfNotFound
                        error('Folder %s does not exist', parent);
                    else
                        file = strings(0, 1);
                        return;
                    end
                end
                stem = [leaf, ext];
                
                % find possible matches
                switch type
                    case 'ap'
                        candidates = Neuropixel.ImecDataset.listAPFilesInDir(parent);
                        bin_ext = '.imec.ap.bin';
                    case 'lf'
                        candidates = Neuropixel.ImecDataset.listLFFilesInDir(parent);
                        bin_ext = '.imec.lf.bin';
                    otherwise
                        error('Unknown type %s');
                end
                mask = startsWith(candidates, stem);
                
                if ~any(mask)
                    if errorIfNotFound
                        error('No %s matches for %s* exist', type, fileOrFileStem);
                    else
                        file = strings(0, 1);
                        return;
                    end
                        
                elseif nnz(mask) > 1
                    exactName = append(stem, bin_ext);
                    exactMatch = mask & strcmp(candidates, exactName);
                    if any(exactMatch)
                        mask = exactMatch;
                    elseif ~returnMultiple
                        error('Multiple %s matches for %s* exist, none exactly matches %s. Narrow down the prefix.', type, fileOrFileStem, exactName);
                    end
                end
                
                candidates = candidates(mask);
                file = strings(numel(candidates), 1);
                for iF = 1:numel(candidates)
                    file(iF) = fullfile(parent, candidates(iF));
                end
            end
            
            file = string(file);
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


            match = regexp(file, '(?<stem>[\w\.]+).imec.(?<type>\w+).bin', 'names', 'once');
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

        function [imecOut, transformExtraArg] = writeConcatenatedFileMatchGains(imecList, outPath, varargin)
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

            p.addParameter('transformAP', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (imec, dataChunk) and return dataChunk someplace
            p.addParameter('transformLF', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (imec, dataChunk) and return dataChunk someplace
            p.addParameter('transformExtraArg', [], @(x) true);
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
            transformExtraArg = p.Results.transformExtraArg;
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
%                 if isConcatenation
                    meta.concatenated = strjoin(stemList, ':');
                    meta.concatenatedSamples = cellfun(@(imec) imec.nSamplesAP, imecList);
                    meta.concatenatedGains = gains;
                    meta.concatenatedMultipliers = multipliers;
                    meta.concatenatedAdcBits = cellfun(@(imec) imec.adcBits, imecList);
                    meta.concatenatedAiRangeMin = cellfun(@(imec) imec.apRange(1), imecList);
                    meta.concatenatedAiRangeMax = cellfun(@(imec) imec.apRange(2), imecList);
%                 end
                
                if ~isempty(timeShiftsAP)
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
                transformExtraArg = writeCatFile(outFile, chIndsByFile, 'ap', multipliers, chunkSize, chunkEdgeExtraSamples, ...
                    p.Results.transformAP, timeShiftsAP, dryRun, transformExtraArg);
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
                transformExtraArg = writeCatFile(outFile, chIndsByFile, 'lf', multipliers, chunkSize, chunkEdgeExtraSamples, ...
                    p.Results.transformLF, timeShiftsLF, dryRun, transformExtraArg);
            end

            if p.Results.writeSyncSeparate
                outFile = fullfile(outPath, [leaf '.imec.sync.bin']);
                fprintf('Writing separate sync bin file %s', lastFilePart(outFile));
                transformExtraArg = writeCatFile(outFile, imecList{1}.syncChannelIndex, 'sync', ones(nFiles, 1, 'int16'), chunkSize, chunkEdgeExtraSamples, ...
                    {}, timeShiftsAP, dryRun, transformExtraArg);
            end

            outFile = fullfile(outPath, [leaf '.imec.ap.bin']);
            imecOut = Neuropixel.ImecDataset(outFile, 'channelMap', imecList{1}.channelMapFile);

            function transformExtraArg = writeCatFile(outFile, chIndsByFile, mode, multipliers, chunkSize, chunkExtra, procFnList, timeShifts, dryRun, transformExtraArg)
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
                                % pass as many inputs to fn as it can handle, including user provided args at the end
                                extraArgs = {chIds, source_idx};
                                if ~isempty(transformExtraArg)
                                    extraArgs{end+1} = transformExtraArg;
                                end
                                
                                ninputs = nargin(fn);
                                if ninputs > 0 && ninputs < 2 + numel(extraArgs)
                                    extraArgs = extraArgs(1:ninputs-2);
                                end
                                
                                noutputs = nargout(fn);
                                if noutputs > 1
                                    [data, transformExtraArg] = fn(imecList{iF}, data, extraArgs{:});
                                else
                                    data = fn(imecList{iF}, data, extraArgs{:});
                                end
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
