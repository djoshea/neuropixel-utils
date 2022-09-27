classdef ImecDataset < handle
% Author: Daniel J. O'Shea (2019)

    properties(SetAccess = protected)
        pathRoot char = '';
        fileStem char = '';
        fileImecNumber = NaN;
        creationTime = NaN;
        nChannels = NaN;

        fileTypeAP = 'ap'; % typically ap or ap_CAR
        fileTypeLF = 'lf'; % typically lf or lf_CAR
        
        nSamplesAP = 0;
        nSamplesLF = 0;
        fsAP = NaN; % samples_per_second
        fsLF = NaN; % samples_per_second
        fsSync = NaN;
        highPassFilterHz = NaN;
        apGain = NaN;
        apRange = [];
        lfGain = NaN;
        lfRange = []

        adcBits = 10;

        snsShankMap (1,1) string;
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
        syncRaw uint16 = [];
    end

    properties(Constant)
        bytesPerSample = 2;
    end

    properties(Dependent)
        hasAP
        hasLF
        hasSync
        
        hasSourceDatasets
        hasSourceAP
        hasSourceLF
        hasSourceSync

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

        % from channel map (although syncChannelIndex will be 1 if sync not in AP or LF file)
        syncChannelId
        syncChannelIndex % if sync in AP file, at one index
        
        % sync will come from either AP, LF, or separate file depending on these flags
        syncInAPFile % is the sync info in the ap file?
        syncInLFFile % is the sync info in the lf file?
    end

    methods
        function imec = ImecDataset(fileOrFileStem, varargin)
            p = inputParser();
            p.addParameter('channelMap', [], @(x) true);
            p.addParameter('syncBitNames', [], @(x) isempty(x) || isstring(x) || iscellstr(x));
            p.addParameter('sourceDatasets', [], @(x) true);
            p.parse(varargin{:})

            fileOrFileStem = char(fileOrFileStem);
            file = Neuropixel.ImecDataset.findImecFileInDir(fileOrFileStem, 'ap', true, false);
            if isempty(file)
                file = Neuropixel.ImecDataset.findImecFileInDir(fileOrFileStem, 'lf', false, false);
                if isempty(file)
                    error('No AP or LF Imec file found at or in %s', fileOrFileStem);
                else
                    isLFOnly = true;
                    [imec.pathRoot, imec.fileStem, imec.fileTypeLF, imec.fileImecNumber] = Neuropixel.ImecDataset.parseImecFileName(file);
                end
            elseif numel(file) == 1
                [imec.pathRoot, imec.fileStem, imec.fileTypeAP, imec.fileImecNumber] = Neuropixel.ImecDataset.parseImecFileName(file);
                isLFOnly= false;
            else
                for iF = 1:numel(file)
                    fprintf('Possible match: %s\n', file{iF});
                end
                error('Multiple imec datasets found in specified location, include file stem or a full path to refine the search');
            end
            
            if ~isLFOnly
                if exist(imec.pathAP, 'file')
                    if ~exist(imec.pathAPMeta, 'file')
                        error('Could not find AP meta file %s', imec.pathAPMeta);
                    end
                    imec.readInfo();
                else
                    error('Could not find AP bin file %s', imec.pathAP);
                end
            else
                if exist(imec.pathLF, 'file')
                    if ~exist(imec.pathLFMeta, 'file')
                        error('Could not find LF meta file %s', imec.pathLFMeta);
                    end
                    imec.readInfo();
                else
                    error('Could not find LF bin file %s', imec.pathLF);
                end
            end

            channelMap = p.Results.channelMap;
            if isempty(channelMap)
                channelMap = "";
            elseif ischar(channelMap)
                channelMap = string(channelMap);
            end
            
            if isa(channelMap, 'Neuropixel.ChannelMap')
                % manually specified channel map
                imec.channelMap = channelMap;
                
            elseif isstring(channelMap)
                if channelMap == ""
                    % try default paths with ks.path
                    channelMap = Neuropixel.Utils.searchForChannelMapInDirectory(imec.pathRoot);
                end
                if channelMap == ""
                    % use default channel map file
                    channelMap = Neuropixel.Utils.getDefaultChannelMapFile(true);
                end
                imec.channelMap = Neuropixel.ChannelMap(channelMap);
            end
            assert(imec.channelMap.nChannels <= imec.nChannels, 'Channel count is less than number of channels in channel map');

            if ~isempty(p.Results.syncBitNames)
                imec.setSyncBitNames(1:numel(p.Results.syncBitNames), p.Resuls.syncBitNames);
            end
            
            if ~isempty(p.Results.sourceDatasets)
                assert(isa(p.Results.sourceDatasets, 'Neuropixel.ImecDataset'));
                imec.setSourceDatasets(p.Results.sourceDatasets);
            end
        end

        function readInfo(imec)
            if imec.hasAP
                metaAP = imec.readAPMeta();
                if imec.hasLF
                    metaLF = imec.readLFMeta();
                else
                    metaLF = struct();
                end
                meta = metaAP;
                
            elseif imec.hasLF
                metaAP = struct();
                metaLF = imec.readLFMeta();
                meta = metaLF;
            else
                error('Must have either AP or LF data files');
            end
            
            imec.nChannels = meta.nSavedChans;
            imec.creationTime = datenum(meta.fileCreateTime, 'yyyy-mm-ddTHH:MM:SS');

            if imec.hasAP
                imec.fsAP = metaAP.imSampRate;
                if isfield(metaAP, 'imHpFlt')
                    imec.highPassFilterHz = metaAP.imHpFlt;
                end
            elseif imec.hasSourceDatasets
                imec.fsAP = imec.sourceDatasets(1).fsAP;
            end
            
            if imec.hasLF
                imec.fsLF = metaLF.imSampRate;
            elseif imec.hasSourceDatasets
                imec.fsLF = imec.sourceDatasets(1).fsLF;
            end

            % parse imroTable
            m = regexp(meta.imroTbl, '\(([\d, ]*)\)', 'tokens');
            gainVals = strsplit(m{2}{1}, ' ');
            imec.apGain = str2double(gainVals{4});
            imec.lfGain = str2double(gainVals{5});
            
            if imec.hasAP
                imec.apRange = [metaAP.imAiRangeMin metaAP.imAiRangeMax];
            end
            
            if imec.hasLF
                imec.lfRange = [metaLF.imAiRangeMin metaLF.imAiRangeMax];
            end
            
            % copy snsShankMap in case needed for building channel map
            if isfield(meta, 'snsShankMap')
                imec.snsShankMap = meta.snsShankMap;
            end

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
                if round(imec.nSamplesAP) ~= imec.nSamplesAP
                    warning('AP bin file size is not an integral number of samples, file data may not be fully copied, truncating nSamplesAP');
                    imec.nSamplesAP = floor(imec.nSamplesAP);
                end
                
                imec.concatenationInfoAP = Neuropixel.ConcatenationInfo(imec, 'ap', metaAP);
            end
            
            if imec.hasLF
                fid = imec.openLFFile();
                fseek(fid, 0, 'eof');
                bytes = ftell(fid);
                fclose(fid);
                imec.nSamplesLF = bytes / imec.bytesPerSample / imec.nChannels;
                if round(imec.nSamplesLF) ~= imec.nSamplesLF
                    warning('LF bin file size is not an integral number of samples, file data may not be fully copied, truncating nSamplesLF');
                    imec.nSamplesLF = floor(imec.nSamplesLF);
                end
              
                imec.concatenationInfoLF = Neuropixel.ConcatenationInfo(imec, 'lf', metaLF);
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
            
            % attempt to fill in missing metadata as well
            if isnan(imec.fsAP)
                imec.fsAP = imec.sourceDatasets(1).fsAP;
            end
            
            if isnan(imec.fsLF)
                imec.fsLF = imec.sourceDatasets(1).fsLF;
            end
            
            % attempt to cross infer concatenation info if missing
            if imec.hasAP && ~imec.hasLF && imec.hasSourceLF
                % infer LF concatenation info from AP info to enable source LF data access through concatInfo lookup
                if ~isempty(imecList)
                    imec.concatenationInfoLF = Neuropixel.ConcatenationInfo.inferFromOtherBand(imec, 'lf', 'ap');
                    imec.nSamplesLF = imec.concatenationInfoLF.nSamples;
                    imec.lfRange = imec.concatenationInfoLF.ranges(1, :); % we assume that the ranges are shared across all datasets
                else
                    % clear existing info
                    imec.concatenationInfoLF = [];
                    imec.nSamplesLF = 0;
                    imec.lfRange = [];
                end
            elseif imec.hasLF && ~imec.hasAP && imec.hasSourceAP
                if ~isempty(imecList)
                    % infer AP concatenation info from LF info to enable source AP data access through concatInfo lookup
                    imec.concatenationInfoAP = Neuropixel.ConcatenationInfo.inferFromOtherBand(imec, 'ap', 'lf');
                    imec.nSamplesAP = imec.concatenationInfoAP.nSamples;
                    imec.apRange = imec.concatenationInfoLF.ranges(1, :); % we assume that the ranges are shared across all datasets
                else
                    % clear existing info
                    imec.concatenationInfoAP = [];
                    imec.nSamplesAP = 0;
                    imec.apRange = [];
                end
            end
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

%     methods  % these functions read a contiguous block of samples over a contiguous band of channels
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
%     end
    
    methods % Sync channel read / cache
        function [syncRaw, fsSync] = readSync(imec, varargin)
            p = inputParser();
            p.addOptional('reload', false, @islogical); % if true, ignore cache in memory imec.syncRaw 
            p.addParameter('preferredBand', '', @Neuropixel.Utils.isstringlike); % prefer a band (ap or lf), though will only be obeyed if not already loaded / cached
            p.addParameter('ignoreCached', false, @islogical); % if true, ignore cache on disk in imec.pathSyncCached
            p.parse(varargin{:});
            preferredBand = string(p.Results.preferredBand);
            
            if isempty(imec.syncRaw) || p.Results.reload
                if exist(imec.pathSyncCached, 'file') && ~p.Results.ignoreCached
                    [~, f, e] = fileparts(imec.pathSyncCached);
                    fprintf('Loading sync from cached %s%s\n', f, e);
                    ld = load(imec.pathSyncCached);
                    imec.syncRaw = typecast(ld.sync, 'uint16');
                    if isfield(ld, 'fsSync')
                        imec.fsSync = ld.fsSync;
                    else
                        % assume from AP if not specified since I added sync from LF later
                        imec.fsSync = imec.fsAP;
                    end
                else
                    % this will automatically redirect to a separate sync file
                    % or to the ap file depending on .syncInAPFile
                    switch preferredBand 
                        case {"", "auto"}
                            fprintf('Loading sync channel auto (this will take some time)...\n');
                            [mm, imec.fsSync] = imec.memmapSync_full();
                            imec.syncRaw = typecast(mm.Data.x(imec.syncChannelIndex, :)', 'uint16');
                            imec.saveSyncCached();
                        case 'ap'
                            fprintf('Loading sync channel from AP band (this will take some time)...\n');
                            mm = imec.memmapAP_full();
                            imec.syncRaw = typecast(mm.Data.x(imec.syncChannelIndex, :)', 'uint16');
                            imec.fsSync = imec.fsAP;
                            imec.saveSyncCached();
                        case 'lf'
                            fprintf('Loading sync channel from LF band (this will take some time)...\n');
                            mm = imec.memmapLF_full();
                            imec.syncRaw = typecast(mm.Data.x(imec.syncChannelIndex, :)', 'uint16');
                            imec.fsSync = imec.fsLF;
                            imec.saveSyncCached();
                        otherwise
                            error('Unknown preferredBand %s', preferredBand);
                    end
                end
            end
            syncRaw = imec.syncRaw;
            fsSync = imec.fsSync;
        end

        function saveSyncCached(imec)
            sync = imec.readSync();
            fsSync = imec.fsSync; %#ok<PROP>
            save(imec.pathSyncCached, 'sync', 'fsSync');
        end
        
        function clearSyncCached(imec)
            if exist(imec.pathSyncCached, 'file') > 0
                fprintf('Deleting cached sync file %s\n', imec.pathSyncCached);
                delete(imec.pathSyncCached);
            end
            imec.syncRaw = [];
            imec.fsSync = NaN;
        end
        
        function updateSyncCached(imec, varargin)
            imec.syncRaw = [];
            imec.fsSync = [];
            if exist(imec.pathSyncCached, 'file')
                [sync, fsSync] = imec.readSync('reload', true, 'ignoreCached', true);
                save(imec.pathSyncCached, 'sync', 'fsSync');
                imec.syncRaw = sync;
                imec.fsSync = fsSync;
            end
        end

        function tf = readSyncBit(imec, bit)
            bit = imec.lookupSyncBitByName(bit);
            tf = logical(bitget(imec.readSync(), bit));
        end
        
        function vec = readSync_idx(imec, idx, varargin)
            p = inputParser();
            p.addParameter('band', '', @Neuropixel.Utils.isstringlike);
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.parse(varargin{:});
            
            band = string(p.Results.band);
            switch band
                case ""
                    fsRequested = NaN;
                case "auto"
                    fsRequested = NaN;
                case "lf"
                    fsRequested = imec.fsLF;
                case "ap"
                    fsRequested = imec.fsAP;
                otherwise
                    error('Unknown band %s', band);
            end
                    
            fromSource = p.Results.fromSourceDatasets;
            if ~fromSource
                % grab cached data if the sampling rate matches, otherwise use the memmap
                if ~isempty(imec.syncRaw) && imec.fsSync == fsRequested
                    vec = imec.syncRaw(idx);
                elseif ismember(band, ["", "auto"])
                    mm = imec.memmapSync_full();
                    vec = mm.Data.x(imec.syncChannelIndex, idx)';
                elseif band == "ap"
                    mm = imec.memmapAP_full();
                    vec = mm.Data.x(imec.syncChannelIndex, idx)';
                elseif band == "lf"
                    mm = imec.memmapLF_full();
                    vec = mm.Data.x(imec.syncChannelIndex, idx)';
                end
                
            elseif imec.syncInAPFile
                mmSet = imec.memmap_sourceAP_full();
                [sourceFileInds, sourceSampleInds] = imec.concatenationInfoAP.lookup_sampleIndexInSourceFiles(idx);
                vec = Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSet, sourceFileInds, sourceSampleInds, imec.syncChannelIndex);
            elseif imec.syncInLFFile
                mmSet = imec.memmap_sourceLF_full();
                [sourceFileInds, sourceSampleInds] = imec.concatenationInfoLF.lookup_sampleIndexInSourceFiles(idx);
                vec = Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSet, sourceFileInds, sourceSampleInds, imec.syncChannelIndex);
            else
                error('Cannot read source sync unless sync derives from AP or LF');
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
        
        function idxLF = closestSampleLFForAP(imec, idxAP)
            idxLF = floor(double(idxAP-1) * double(imec.fsLF) / double(imec.fsAP)) + 1;
        end
        
        function idxAP = closestSampleAPForLF(imec, idxLF)
            idxAP = floor(double(idxLF-1) * double(imec.fsAP) / double(imec.fsLF)) + 1;
        end
        
        function data = readAP_idx(imec, sampleIdx, varargin)
            data = imec.internal_read_idx(sampleIdx, 'band', 'ap', varargin{:});
        end
        
        function data = readAP_idx_alongside_source(imec, sampleIdx, varargin)
            data = imec.internal_read_idx_alongside_source(sampleIdx, 'band', 'ap', varargin{:});
        end
        
        function data = readLF_idx(imec, sampleIdx, varargin)
            data = imec.internal_read_idx(sampleIdx, 'band', 'lf', varargin{:});
        end
        
        function data = internal_read_idx_alongside_source(imec, sampleIdx, varargin)
            p = inputParser();
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.addParameter('scaleSourceToMatch', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            % returns data from this dataset as well as data fromSourceDatasets, concatenated along dim 3
            data_this = imec.internal_read_idx(sampleIdx, p.Unmatched);
            data_source = imec.internal_read_idx(sampleIdx, 'fromSourceDatasets', true, 'scaleSourceToMatch', p.Results.scaleSourceToMatch, p.Unmatched);
            
            data = cat(3, data_this, data_source);
        end
        
        function data = internal_read_idx(imec, sampleIdx, varargin)
            p = inputParser();
            p.addParameter('band', 'ap', @Neuropixel.Utils.isstringlike);
            p.addParameter('applyScaling', true, @islogical); % convert to uV before processing
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.addParameter('scaleSourceToMatch', false, @islogical); % when fromSourceDatasets is true, scale the raw data to match the scaling of the cleaned datsets
            p.parse(varargin{:});
            
            fromSource = p.Results.fromSourceDatasets;
            scaleSourceDataToMatch = p.Results.scaleSourceToMatch;
            band = string(p.Results.band);
            switch band
                case 'ap'
                    if ~fromSource
                        assert(imec.hasAP, 'ImecDataset does not have AP band');
                    else
                        assert(imec.hasSourceAP, 'ImecDataset does not have source datasets with AP band');
                    end
                case 'lf'
                    if ~fromSource
                        assert(imec.hasLF, 'ImecDataset does not have LF band');
                    else
                        assert(imec.hasSourceLF, 'ImecDataset does not have source datasets with LF band');
                    end
                otherwise
                    error('Unknown band %s', band);
            end
            
            ch_conn_mask = imec.lookup_channelIds(imec.connectedChannels);
            
            if ~fromSource
                if band == "ap"
                    mm = imec.memmapAP_full();
                    scaleToUv = imec.apScaleToUv;
                else
                    mm = imec.memmapLF_full();
                    scaleToUv = imec.lfScaleToUv;
                end
            
                if any(isnan(sampleIdx))
                    mask = ~isnan(sampleIdx);
                    data = single(mm.Data.x(:, sampleIdx)); % must be single to support NaNs
                    data = Neuropixel.Utils.TensorUtils.inflateMaskedTensor(data, 2, mask, NaN);
                else
                    data = mm.Data.x(:, sampleIdx);
                end
                
                if p.Results.applyScaling
                    data = single(data);
                    data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(scaleToUv);
                end 
            else
                if band == "ap"
                    mmSet = imec.memmap_sourceAP_full();
                    [sourceFileInds, sourceSampleInds] = imec.concatenationInfoAP.lookup_sampleIndexInSourceFiles(sampleIdx);
                    scalingByFile = imec.concatenationInfoAP.scaleToUvs;
                    scaleToUvThis = imec.apScaleToUv;
                else
                    mmSet = imec.memmap_sourceLF_full();
                    [sourceFileInds, sourceSampleInds] = imec.concatenationInfoLF.lookup_sampleIndexInSourceFiles(sampleIdx);
                    scalingByFile = imec.concatenationInfoLF.scaleToUvs;
                    scaleToUvThis = imec.lfScaleToUv;
                end
                
                % if applyScaling, use the per-file scalings to produce uv, scaling only the connected channels
                % otherwise we pass [] so that no scaling is done
                if ~p.Results.applyScaling
                    if scaleSourceDataToMatch
                        relativeScaling = scalingByFile ./ scaleToUvThis;
                        assert(all(round(relativeScaling) == relativeScaling), 'Non-integer scaling from ImecDataset source datasets');
                        scalingByFile = int16(relativeScaling);
                    else
                        scalingByFile = [];
                    end
                end
                data = Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSet, sourceFileInds, sourceSampleInds, imec.channelIds, scalingByFile, ch_conn_mask);
            end
        end
       
        function [mat, sampleIdx] = readAP_timeWindow(imec, timeWindowSec, varargin)
            idxWindow = imec.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = imec.readAP_idx(sampleIdx, varargin{:});
        end
        
        function [mat, sampleIdx] = readLF_timeWindow(imec, timeWindowSec, varargin)
            idxWindow = imec.closestSampleLFForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = imec.readLF_idx(sampleIdx, varargin{:});
        end
        
        function [mat, sampleIdx] = readSyncBits_timeWindow(imec, bits, timeWindowSec)
            idxWindow = imec.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = readSyncBits_idx(bits, sampleIdx);
        end
        
        function [mat, sampleIdx] = readSyncLFBits_timeWindow(imec, bits, timeWindowSec)
            idxWindow = imec.closestSampleLFForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = readSyncLFBits_idx(bits, sampleIdx);
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
            imec.internal_inspect_idxWindow(idxWindow, 'band', 'ap', varargin{:});
        end
        
        function inspectLF_idxWindow(imec, idxWindow, varargin)
            imec.internal_inspect_idxWindow(idxWindow, 'band', 'lf', varargin{:});
        end
        
        function inspectLF_timeWindow(imec, timeWindowSec, varargin)
            idxWindow = imec.closestSampleLFForTime(timeWindowSec);
            imec.inspectLF_idxWindow(idxWindow, 'timeInSeconds', true, varargin{:});
        end
        
        function internal_inspect_idxWindow(imec, idxWindow, varargin)
            p = inputParser();
            p.addParameter('band', 'ap', @ischar);
            p.addParameter('channels', imec.mappedChannels, @(x) isempty(x) || isvector(x));
            p.addParameter('invertChannels', imec.channelMap.invertChannelsY, @islogical);
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            p.addParameter('showSync', true, @isvector);
            p.addParameter('syncBits', imec.syncBitsNamed, @isvector);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 0.95, @isscalar);
            p.addParameter('car', false, @islogical); % subtract median of all channels at each time
            p.addParameter('center', false, @islogical); % subtract median of each channel over time
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.addParameter('compareSourceDatasets', false, @islogical); % if true, supercedes fromSourceDatasets and both will be plotted
            p.addParameter('syncFromSourceDatasets', [], @(x) isempty(x) || islogical(x));
            p.addParameter('downsample',1, @isscalar); 
            p.addParameter('timeInSeconds', false, @islogical);
            p.addParameter('timeRelativeTo', 0, @isscalar);
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'Neuropixel.TrialSegmentationInfo')); % to mark trial boundaries
            
            p.addParameter('markSampleIdx', [], @isvector);
            p.addParameter('markSampleMode', 'rug', @ischar);
            p.addParameter('markSampleColor', [0.5 0 0.5], @(x) true);
            
            p.addParameter('style', 'traces', @Neuropixel.Utils.isstringlike);
            
            p.parse(varargin{:});
            
            % by default, sync comes from same source v. processed data as the data being plotted, 
            % but this can be overrriden
            fromSource = p.Results.fromSourceDatasets;
            compareSource = p.Results.compareSourceDatasets;
            syncFromSource = p.Results.syncFromSourceDatasets;
            if isempty(syncFromSource)
                syncFromSource = fromSource;
            end
            
            band = string(p.Results.band);
            switch band
                case 'ap'
                    fsBand = imec.fsAP;
                case 'lf'
                    fsBand = imec.fsLF;
                otherwise
                    error('Unknown band %s', band);
            end
            
            if numel(idxWindow) > 2
                idxWindow = [idxWindow(1), idxWindow(end)];
            end
            sampleIdx = idxWindow(1):idxWindow(2);
            
            sampleIdx = floor(sampleIdx);
            
            if compareSource
                mat = imec.internal_read_idx_alongside_source(sampleIdx, 'band', band, 'applyScaling', true); % C x T x 2 (this, source)
                %labelsSuperimposed = ["processed", "source"];
            else
                mat = imec.internal_read_idx(sampleIdx, 'band', band, 'fromSourceDatasets', fromSource, 'applyScaling', true); % C x T
                if fromSource
                    %labelsSuperimposed = "source";
                else
                    %labelsSuperimpposed = "processed";
                end
            end
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
            
            mat = mat(channelInds, :, :);
            labels = labels(channelInds);
            connected = ismember(channelIds, imec.connectedChannels);
            bad = ismember(channelIds, imec.badChannels);
            
            if p.Results.downsample > 1
                mat = mat(:, 1:p.Results.downsample:end, :);
                sampleIdx = sampleIdx(1:p.Results.downsample:end);
            end
            mat = single(mat);
            if p.Results.center
                mat = mat - median(mat, 2);
            end
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
                normalizeMask = flipud(normalizeMask);
                colors = flipud(colors);
            end
            
            % append sync bit info to plot in purple
            showSync = p.Results.showSync;
            if (syncFromSource && ~imec.hasSourceSync) || (~syncFromSource && ~imec.hasSync)
                warning('Cannot show sync data since no sync data present');
                showSync = false;
            end
            syncBits = p.Results.syncBits;
            
            if showSync
                if isempty(syncBits)
                    syncBits = 1:imec.nSyncBits;
                end
                % make sure we grab sync from the matching band
                if compareSource
                    sync_this = imec.readSyncBits_idx(syncBits, sampleIdx, 'fromSourceDatasets', false, 'band', band);
                    sync_source = imec.readSyncBits_idx(syncBits, sampleIdx, 'fromSourceDatasets', true, 'band', band);
                    syncBitMat = cat(3, sync_this, sync_source);
                else
                    syncBitMat = imec.readSyncBits_idx(syncBits, sampleIdx, 'fromSourceDatasets', syncFromSource, 'band', band);
                end
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
                time = double(time) / fsBand;
            end
            
            style = string(p.Results.style);
            
            switch style
                case 'traces'
                    if compareSource
                        layer_colors = [0 0 0; 0.7 0.2 0.2];
                        Neuropixel.Utils.plotStackedTraces(time, permute(mat, [2 1 3]), 'layerColors', layer_colors, 'labels', labels, ...
                            'gain', p.Results.gain, 'invertChannels', false, ...
                            'normalizeMask', normalizeMask, 'normalizeEach', false); 
                    else
                        Neuropixel.Utils.plotStackedTraces(time, mat', 'colors', colors, 'labels', labels, ...
                            'gain', p.Results.gain, 'invertChannels', false, ...
                            'normalizeMask', normalizeMask, 'normalizeEach', false);
                    end
                case 'pmatbal'
                    if p.Results.invertChannels
                        mat = flipud(mat);
                        labelsF = flipud(labels);
                    else
                        labelsF = labels;
                    end
                    ytick = 1:size(mat, 1);
                    Neuropixel.Utils.pmatbal(mat, 'x', time, 'y', ytick);
                    set(gca, 'YTick', ytick, 'YTickLabel', labelsF);
                    
                otherwise
                    error('Unknown style %s', style);
            end
            
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
                    markTimes = double(markSampleIdx) / fsBand;
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
        
        function inspectSync_idxWindow(imec, idxWindow, varargin)
           imec.inspectAP_idxWindow(idxWindow, 'channels', [], 'showSync', true, varargin{:}); 
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

        function [mm, fsSync] = memmapSync_full(imec)
            if imec.syncInAPFile
                % still has nChannels
                mm = memmapfile(imec.pathSync, 'Format', {'int16', [imec.nChannels imec.nSamplesAP], 'x'});
                fsSync = imec.fsAP;
            elseif imec.syncInLFFile
                % still has nChannels
                mm = memmapfile(imec.pathSync, 'Format', {'int16', [imec.nChannels imec.nSamplesLF], 'x'});
                fsSync = imec.fsLF;
            else
                % only sync channel
                mm = memmapfile(imec.pathSync, 'Format', {'int16', [1 imec.nSamplesAP], 'x'});
                fsSync = imec.fsAP;
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
        
        function [mmSet, fsSync] = memmap_sourceSync_full(imec, varargin)
            assert(~isempty(imec.sourceDatasets));
            nSources = numel(imec.sourceDatasets);
            mmSet = cell(nSources, 1);
            fsSync = nan(nSources, 1);
            for iF = 1:nSources
                imecSrc = imec.sourceDatasets(iF);
                if imecSrc.syncInAPFile
                    % still has nChannels
                    mmSet{iF} = memmapfile(imecSrc.pathSync, 'Format', {'int16', [imec.nChannels imec.nSamplesAP], 'x'});
                    fsSync(iF) = imec.fsAP;
                elseif imecSrc.syncInLFFile
                    % still has nChannels
                    mmSet{iF} = memmapfile(imec.pathSync, 'Format', {'int16', [imec.nChannels imec.nSamplesLF], 'x'});
                    fsSync(iF) = imec.fsLF;
                else
                    % only sync channel
                    mmSet{iF} = memmapfile(imec.pathSync, 'Format', {'int16', [1 imec.nSamplesAP], 'x'});
                    fsSync(iF) = imec.fsAP;
                end
            end
            
            assert(numel(unique(fsSync)) == 1);
            fsSync = fsSync(1);
        end
    end
    
    methods(Static)
        function out = multi_mmap_extract_sample_idx(mmSet, fileInds, origSampleInds, chInds, scalingByFile, scalingChMask)
            % given [fileInds, origSampleInds] as returned by ConcatenationInfo/lookup_sampleIndexInSourceFiles
            % extract those samples from the set of memory mapped files in mmSet (returned by memmap**_all)
            
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
            
            % if scalingByFile is non-empty vector with length nFiles, samples will be scaled by these factors
            % on a per-file basis and converted to singles
            if nargin < 5
                scalingByFile = [];
            end
            if nargin < 6
                scalingChMask = true(numel(chInds), 1);
            end
            
            assert(numel(fileInds) == numel(origSampleInds));
            nSamplesOut = numel(origSampleInds);
            
            if isempty(scalingByFile) || isinteger(scalingByFile)
                cls = mmSet{1}.Format{1};
            else
                cls = 'single';
                assert(numel(scalingByFile) == numel(mmSet));
            end
            nFiles = numel(mmSet);
            
            out = zeros(nCh, nSamplesOut, cls);
            
            for iF = 1:nFiles
                mask = fileInds == iF;
                if isempty(scalingByFile)
                    out(:, mask) = mmSet{iF}.Data.x(chInds, origSampleInds(mask));
                elseif isinteger(scalingByFile)
                    % apply a scale factor but keep as int16
                    out(:, mask) = mmSet{iF}.Data.x(chInds, origSampleInds(mask));
                    out(scalingChMask, mask) = out(scalingChMask, mask) * cast(scalingByFile(iF), cls);
                else
                    % scale (typically to uV) while converting to single
                    out(:, mask) = single(mmSet{iF}.Data.x(chInds, origSampleInds(mask)));
                    out(scalingChMask, mask) = out(scalingChMask, mask) * single(scalingByFile(iF));
                end
            end
        end
    end

    methods(Hidden) % Read data at specified times
        function [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet, scaleToUv_by_snippet, group_ids, data_trust_mask] = readSnippetsRaw(imec, times, window, varargin)
            % for each sample index in times, read the window times + window(1):window(2)
            % of samples around this time from some channels

            p = inputParser();
            p.addParameter('band', 'ap', @Neuropixel.Utils.isstringlike);
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.addParameter('syncFromSourceDatasets', [], @(x) isempty(x) || islogical(x));
            
            % specify one of THESE (same channels every snippet or nChannels x numel(times))
            p.addParameter('channel_ids', [], @(x) isempty(x) || isvector(x));
            p.addParameter('channel_ids_by_snippet', [], @ismatrix);
            
            % or THESE (different channels for each group of snippets sharing the same cluster_ids
            p.addParameter('channel_ids_by_cluster', [], @ismatrix);
            p.addParameter('unique_cluster_ids', [], @isvector);
            p.addParameter('cluster_ids_by_snippet', [], @(x) isempty(x) || isvector(x)); % same length as numel(times), corresponding to which cluster was pulled out (e.g. as a waveform)
            p.addParameter('group_ids_by_snippet', [], @(x) isempty(x) || isvector(x)); % same length as numel(times), utility used for average_by_group_id
            
            p.addParameter('car', false, @islogical); % subtract median over channels
            p.addParameter('center', false, @islogical); % subtract median over time
            
            % two modes for averaging groups snippets together while extracting (saves memory vs doing this post hoc)
            p.addParameter('average_weight', [], @(x) isempty(x) || isvector(x)); % this will mutliply each snippet before averaging
            p.addParameter('average_by_cluster_id', false, @islogical); % average all snippsets from the same cluster_id together (cluster_ids_by_snippet must be passed)
            p.addParameter('average_by_group_id', false, @islogical); % average all snippets from the same group_id together (group_ids_by_snippet must be provided)
            % if both of these are true, the averaging will take 
            
            p.addParameter('applyScaling', false, @islogical); % scale to uV and return single
            p.addParameter('scaleSourceToMatch', false, @islogical); % when fromSourceDatasets is true, scale the raw data to match the scaling of the cleaned datsets

            p.addParameter('data_distrust_mask', [], @(x) isempty(x) || islogical(x)); % if provided, these samples will be flagged when included in snippets via ss.data_trust_mask
            p.parse(varargin{:});

            channel_ids = p.Results.channel_ids;
            channel_ids_by_snippet = p.Results.channel_ids_by_snippet;
            channel_ids_by_cluster = p.Results.channel_ids_by_cluster;
            cluster_ids = p.Results.cluster_ids_by_snippet;
            unique_cluster_ids = p.Results.unique_cluster_ids;
            if isempty(unique_cluster_ids), unique_cluster_ids = unique(cluster_ids); end
            
            average_by_cluster_id = p.Results.average_by_cluster_id;
            if isempty(cluster_ids) && average_by_cluster_id
                error('cluster_ids_by_snippet must be provided when average_by_cluster_id is true');
            end
            if ~isempty(cluster_ids)
                [~, cluster_inds] = ismember(cluster_ids, unique_cluster_ids);
            end
            
            group_ids = p.Results.group_ids_by_snippet;
            average_by_group_id = p.Results.average_by_group_id;
            if isempty(group_ids)
                if average_by_group_id
                    error('group_ids_by_snippet must be provided when average_by_group_id is true');
                end
            else
                [unique_group_ids, ~, group_inds] = unique(group_ids);
            end
            if average_by_group_id && ~average_by_cluster_id && ~isempty(cluster_ids)
                error('Cannot average_by_group_id and not average_by_cluster_id when cluster_ids is specified. Group-averaged snippets no longer correspond to clusters.');
            end
            
            if ~isempty(channel_ids_by_snippet)
                assert(size(channel_ids_by_snippet, 2) == numel(times), 'channel_ids must be nChannels x numel(times)');
                
            elseif ~isempty(channel_ids_by_cluster)
                assert(~isempty(cluster_ids), 'cluster_ids_by_snippet must be specified when channel_ids_by_cluster is used');
                channel_ids_by_snippet = channel_ids_by_cluster(:, cluster_inds);
                
            elseif ~isempty(channel_ids)
                channel_ids = Neuropixel.Utils.makecol(channel_ids);
                channel_ids_by_snippet = repmat(channel_ids, 1, numel(times));
                
            else
                error('Specify either channel_ids or channel_ids_by_cluster');
            end
            
            channel_inds_by_snippet = imec.lookup_channelIds(channel_ids_by_snippet);
            syncChannelIndex = imec.syncChannelIndex;
            syncInChannelInds = any(ismember(channel_inds_by_snippet, syncChannelIndex), 'all');

            band = string(p.Results.band);
            fromSourceDatasets = p.Results.fromSourceDatasets;
            syncFromSource = p.Results.syncFromSourceDatasets;
            if isempty(syncFromSource)
                syncFromSource = fromSourceDatasets;
            end

            data_distrust_mask = p.Results.data_distrust_mask;
            use_distrust_mask = ~isempty(data_distrust_mask);

            do_center = p.Results.center;
            do_car = p.Results.car;
            
            switch band
                case 'ap'
                    nSamples = imec.nSamplesAP;
                    if ~fromSourceDatasets
                        mm = imec.memmapAP_full();
                        scaleToUv = imec.apScaleToUv;
                    else
                        mmSet = imec.memmap_sourceAP_full();
                        concatInfo = imec.concatenationInfoAP;
                        scaleToUv_this = imec.apScaleToUv;
                        scaleToUv_by_source = cat(1, imec.sourceDatasets.apScaleToUv);
                    end     
                    
                    % do we need to specially handle the sync channel?
                    if syncFromSource ~= fromSourceDatasets && syncInChannelInds
                        handleSyncSeparately = true;
                        if ~syncFromSource
                            mmSync = imec.memmapAP_full();
                        else
                            mmSetSync = imec.memmap_sourceAP_full();
                        end
                    else
                        handleSyncSeparately = false;
                    end
                case 'lf'
                    nSamples = imec.nSamplesLF;
                    if ~fromSourceDatasets
                        mm = imec.memmapLF_full();
                        scaleToUv = imec.lfScaleToUv;
                    else
                        mmSet = imec.memmap_sourceLF_full();
                        concatInfo = imec.concatenationInfoLF;
                        scaleToUv_this = imec.lfScaleToUv;
                        scaleToUv_by_source = cat(1, imec.sourceDatasets.lfScaleToUv);
                    end
                    
                    % do we need to specially handle the sync channel?
                    if syncFromSource ~= fromSourceDatasets && syncInChannelInds
                        handleSyncSeparately = true;
                        if ~syncFromSource
                            mmSync = imec.memmapLF_full();
                        else
                            mmSetSync = imec.memmap_sourceLF_full();
                        end
                     else
                         handleSyncSeparately = false;
                     end
                    
                otherwise
                    error('Unknown source');
            end
            
            scaleSourceDataToMatch = p.Results.scaleSourceToMatch;
            if fromSourceDatasets
                if scaleSourceDataToMatch
                    sourceToDestScaling = scaleToUv_by_source ./ scaleToUv_this;
                end
            else
                scaleSourceDataToMatch = false;
            end
            
            assert(nSamples > 0, 'nSamples inferred for dataset is 0');
            
            times = Neuropixel.Utils.makecol(uint64(times));
            nC = size(channel_ids_by_snippet, 1);
            nC_all = imec.nChannels;
            nS = numel(times); % actual number of input snippets extracted

            if use_distrust_mask
                if size(data_distrust_mask, 1) == nC_all - 1
                    % doesn't include sync
                    nC_all_distrust = nC_all - 1;
                else
                    assert(size(data_distrust_mask, 1) == nC_all, 'data_distrust_mask channel count must be %d or %d', nC_all - 1, nC_all);
                    nC_all_distrust = nC;
                end
                assert(size(data_distrust_mask, 2) == size(mm.Data.x, 2), 'data_distrust_mask timepoints does not match data');
            end
            
            % nS out is the accumulated number of output snippets (equals nS or number of clusters if we're averaging)
            if average_by_cluster_id
                if average_by_group_id
                    % averaging within both cluster and group
                    % do unique conjunctions
                    combined_inds = [cluster_inds, group_inds];
                    [unique_inds, ~, averaging_dest_inds] = unique(combined_inds, 'rows'); 
                    nS_out = max(averaging_dest_inds);
                    
                    cluster_ids = unique_inds(:, 1); % associated with averaged output snippets
                    group_ids = unique_inds(:, 2);
                else
                    % averaging by cluster
                    nS_out = numel(unique_cluster_ids);
                    averaging_dest_inds = cluster_inds; % associated with averaged output snippets
                    cluster_ids = unique_cluster_ids;
                    group_ids = [];
                end
            elseif average_by_group_id
                % averaging within groups
                nS_out = numel(unique_group_ids);
                averaging_dest_inds = group_inds;
                cluster_ids = [];
                group_ids = unique_group_ids;
            else
                % no averaging
                nS_out = nS;
                averaging_dest_inds = (1:nS)';
            end
            nT = numel(window(1):window(2));
            
            outClass = 'int16';
            applyScaling  = p.Results.applyScaling;
            if applyScaling
                outClass = 'single';
            end
            if average_by_cluster_id || average_by_group_id
                outClass = 'single';
                assert(~fromSourceDatasets || applyScaling, 'applyScaling must be true for average_by_cluster_id or average_by_group_id with fromSourceDatasets since waveforms may have differing scales');
            end
            
            average_weight = p.Results.average_weight;
            if ~isempty(average_weight)
                average_weight = cast(average_weight, outClass);
            end
            
            out = zeros(nC, nT, nS_out, outClass);
            if use_distrust_mask
                accum_trust_counter = zeros(nC, nT, nS_out);
            else
                accum_counter = zeros(nS_out, 1);
            end
            
            if ~fromSourceDatasets
                scaleToUv_by_snippet = repmat(single(scaleToUv), nS_out, 1);
            else
                % populated below
                scaleToUv_by_snippet = nan(nS_out, 1, 'single');
            end
            if scaleSourceDataToMatch
                sourceToDestScaling_by_snippet = ones(nS_out, 1, 'int16');
            end

            if numel(times) > 10
                prog = Neuropixel.Utils.ProgressBar(numel(times), 'Extracting %s snippets', upper(band));
            else
                prog = [];
            end
            
            good_ch_inds = imec.goodChannelInds;
            idx_request = int64(times') + int64(window(1):window(2))'; % nWindow x nTimes indices
            mask_idx_okay = idx_request >= int64(1) & idx_request <= nSamples;
            idx_request(~mask_idx_okay) = 1; % we'll clear out later
            
            mask_request_okay = all(mask_idx_okay, 1)';
            
            nPerSegment = 50; 
            nSegments = ceil(nS / nPerSegment);
            
            for iSegment = 1:nSegments
                idxS = nPerSegment*(iSegment-1) + 1 : min(nS,  nPerSegment*iSegment);
                idxInsert = averaging_dest_inds(idxS); % this will match idxS if no averaging is being performed, otherwise it indicates which slot to insert into for each snippet

                nS_this = numel(idxS);
                idx_request_this = idx_request(:, idxS);
                mask_request_okay_this = mask_request_okay(idxS);
                if ~fromSourceDatasets
                    extract_all_ch = reshape(mm.Data.x(:, idx_request_this(:)), [nC_all nT, nS_this]);
                else
                    [sourceFileInds, sourceSampleInds] = concatInfo.lookup_sampleIndexInSourceFiles(idx_request_this);
                    extract_all_ch = reshape(Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSet, ...
                        sourceFileInds(:), sourceSampleInds(:)), [nC_all nT nS_this]);
                    scaleToUv_by_snippet(idxInsert) = scaleToUv_by_source(sourceFileInds(1, :)');
                    if scaleSourceDataToMatch
                        sourceToDestScaling_by_snippet(idxInsert) = sourceToDestScaling(sourceFileInds(1, :)');
                    end
                end
                
                if handleSyncSeparately
                    % sync is requested and coming from not the same source as the rest of the channels
                    if ~syncFromSource
                        extract_sync_ch = reshape(mmSync.Data.x(syncChannelIndex, idx_request_this(:)), [1 nT, nS_this]);
                    else
                        [sourceFileInds, sourceSampleInds] = concatInfo.lookup_sampleIndexInSourceFiles(idx_request_this);
                        extract_sync_ch = reshape(Neuropixel.ImecDataset.multi_mmap_extract_sample_idx(mmSetSync, ...
                            sourceFileInds(:), sourceSampleInds(:), syncChannelIndex), [1 nT nS_this]);
                    end
                end
                
                if do_car
                    ar = median(extract_all_ch(good_ch_inds, :, :), 1);
                else
                    ar = zeros([1 nT nS_this], 'like', extract_all_ch);
                end

                if use_distrust_mask
                    distrust_all_ch = reshape(full(data_distrust_mask(:, idx_request_this(:))), [nC_all_distrust, nT, nS_this]);
                end

                for iiS = 1:nS_this
                    this_extract =  extract_all_ch(channel_inds_by_snippet(:, idxS(iiS)), :, iiS);
                    if use_distrust_mask
                        this_distrust = distrust_all_ch(channel_inds_by_snippet(:, idxS(iiS)), :, iiS);
                    end
                    
                    % optionally insert sync from separate dataset, and track it so we don't apply numerical computations to sync line
                    mask_sync = channel_inds_by_snippet(:, idxS(iiS)) == syncChannelIndex;
                    if handleSyncSeparately
                        this_extract(mask_sync, :) = extract_sync_ch(:, :, iiS);
                    end
                    mask_scale = ~mask_sync;
                    this_extract(mask_scale, :) = this_extract(mask_scale, :) - ar(1, :, iiS);
                    
                    if applyScaling
                        % apply scaling to all but sync channels
                        this_extract = cast(this_extract, outClass);
                        this_extract(mask_scale, :) = this_extract(mask_scale, :) * scaleToUv_by_snippet(idxInsert(iiS));
                    end
                    if scaleSourceDataToMatch
                        % this would scale the source data to match the cleaned data while keeping everything int16
                        this_extract(mask_scale, :) = this_extract(mask_scale, :) * sourceToDestScaling_by_snippet(idxInsert(iiS));
                    end
                    if do_center
                        if use_distrust_mask
                            % exclude distrust from median calculation
                            temp = cast(this_extract(mask_scale, :), 'single');
                            temp(this_distrust(mask_scale, :)) = NaN;
                            this_extract(mask_scale, :) = this_extract(mask_scale, :) - cast(median(temp, 2, 'omitnan'), 'like', this_extract);
                        else
                            this_extract(mask_scale, :) = this_extract(mask_scale, :) - median(this_extract(mask_scale, :), 2);
                        end
                    end
                    % scale before accumulating
                    if ~isempty(average_weight)
                        this_extract(mask_scale, :) = this_extract(mask_scale, :) * average_weight(idxS(iiS));
                    end
                    if mask_request_okay_this(iiS)
                        if use_distrust_mask
                            % track a separate accumulator for each entry and zero the distrusted data before adding
                            this_extract(this_distrust) = 0;
                            out(:, :, idxInsert(iiS)) = out(:, :, idxInsert(iiS)) + cast(this_extract, outClass); % which channels for this spike
                            accum_trust_counter(:, :, idxInsert(iiS)) = accum_trust_counter(:, :, idxInsert(iiS)) + ~this_distrust;
                        else
                            out(:, :, idxInsert(iiS)) = out(:, :, idxInsert(iiS)) + cast(this_extract, outClass); % which channels for this spike
                            accum_counter(idxInsert(iiS)) = accum_counter(idxInsert(iiS)) + 1;
                        end
                    end
                end
               
                if ~isempty(prog), prog.increment(nPerSegment); end
            end
            
            if average_by_cluster_id || average_by_group_id
                 % divide the accumulator by the number of units, don't divide by the sum of average_weights, this is treated as normalized already
                 if use_distrust_mask
                     out = out ./ accum_trust_counter;
                 else
                    out = out ./ shiftdim(accum_counter, -2);
                 end
            else
                out(:, ~mask_idx_okay(:)) = 0;
            end
            if use_distrust_mask
                out_trust_mask = accum_trust_counter > 0;
                if isfloat(out)
                    out(~out_trust_mask) = NaN;
                else
                    out(~out_trust_mask) = 0;
                end
            else
                out_trust_mask = [];
            end
            
            if applyScaling
                % scaling already applied
                scaleToUv_by_snippet(:) = 1;
            end
            data_ch_by_time_by_snippet = out;
            data_trust_mask = out_trust_mask;
            if ~isempty(prog), prog.finish(); end
        end
    end
    
    methods  % Read data at specified times
        function snippet_set = readSnippetSet(imec, band, times, window, varargin)
            [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet, scaleToUv_by_snippet, group_ids, data_trust_mask] = ...
                imec.readSnippetsRaw(times, window, 'band', band, varargin{:});
            snippet_set = Neuropixel.SnippetSet(imec, band);
            snippet_set.data = data_ch_by_time_by_snippet;
            snippet_set.data_trust_mask = data_trust_mask;
            snippet_set.scaleToUv = scaleToUv_by_snippet;
            snippet_set.sample_idx = times;
            snippet_set.channel_ids_by_snippet = channel_ids_by_snippet;
            snippet_set.cluster_ids = cluster_ids;
            snippet_set.group_ids = group_ids;
            snippet_set.window = window;
        end
        
        function snippet_set = readAPSnippetSet(imec, times, window, varargin)
            snippet_set = imec.readSnippetSet('ap', times, window, varargin{:});
        end

        function snippet_set = readLFSnippetSet(imec, times, window, varargin)
            snippet_set = imec.readSnippetSet('lf', times, window, varargin{:});
        end

        function rms = computeRMSByChannel(imec, varargin)
            % output will be nMappedChannels x 1 vector of rms
            p = inputParser();
            p.addParameter('band', 'ap', @Neuropixel.Utils.isstringlike);
            p.addParameter('sampleMaskFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % sampleMaskFn(data_ch_x_time, sample_idx_time) --> logical_time mask of time samples valid for use, useful if you have artifacts at known times
            p.addParameter('car', false, @islogical);
            p.addParameter('useChunks', 50, @isscalar);
            p.addParameter('chunkSize', 100000, @isscalar);
            p.parse(varargin{:});
            
            sampleMaskFn = p.Results.sampleMaskFn;
            
            % aim for the middle of the file
            band = string(p.Results.band);

            switch band
                case "ap"
                    chunkSize = round(min(imec.fsAP, p.Results.chunkSize));
                    mm = imec.memmapAP_by_chunk(chunkSize);
                    scaleToUv = imec.apScaleToUv;
                case "lf"
                    chunkSize = round(min(imec.fsLF, p.Results.chunkSize));
                    mm = imec.memmapLF_by_chunk(chunkSize);
                    scaleToUv = imec.lfScaleToUv;
                otherwise
                    error("Unknown band");
            end
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
            rms = rms * scaleToUv;
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
            tf = imec.channelMap.syncInAPFile && imec.hasAP;
        end
        
        function tf = get.syncInLFFile(imec)
            tf = imec.channelMap.syncInLFFile && imec.hasLF;
        end
        
        function id = get.syncChannelId(imec)
            if imec.syncInAPFile || imec.syncInLFFile
                id = imec.channelMap.syncChannelId;
            else
                % if sync is in its own file, assume it's the first and only channel
                id = NaN;
            end
        end
        
        function ind = get.syncChannelIndex(imec)
            if imec.syncInAPFile || imec.syncInLFFile
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
            if isnan(imec.fileImecNumber)
                fileAP = [imec.fileStem '.imec.' imec.fileTypeAP '.bin'];
            else
                fileAP = [imec.fileStem, sprintf('.imec%d.', imec.fileImecNumber), imec.fileTypeAP, '.bin'];
            end
        end

        function tf = get.hasAP(imec)
            tf = exist(imec.pathAP, 'file') == 2;
        end

        function fileAPMeta = get.fileAPMeta(imec)
            if isnan(imec.fileImecNumber)
                fileAPMeta = [imec.fileStem '.imec.ap.meta'];
            else
                fileAPMeta = [imec.fileStem, sprintf('.imec%d.ap.meta', imec.fileImecNumber)];
            end
        end

        function pathAPMeta = get.pathAPMeta(imec)
            pathAPMeta = fullfile(imec.pathRoot, imec.fileAPMeta);
        end

        function pathLF = get.pathLF(imec)
            pathLF = fullfile(imec.pathRoot, imec.fileLF);
        end

        function fileLF = get.fileLF(imec)
            if isnan(imec.fileImecNumber)
                fileLF = [imec.fileStem '.imec.lf.bin'];
            else
                fileLF = [imec.fileStem, sprintf('.imec%d.', imec.fileImecNumber), 'lf.bin'];
            end
        end

        function fileLFMeta = get.fileLFMeta(imec)
            if isnan(imec.fileImecNumber)
                fileLFMeta = [imec.fileStem '.imec.lf.meta'];
            else
                fileLFMeta = [imec.fileStem, sprintf('.imec%d.lf.meta', imec.fileImecNumber)];
            end
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
            elseif imec.syncInLFFile
                fileSync = imec.fileLF;
            else
                fileSync = [imec.fileStem, '.imec.sync.bin'];
            end
        end

        function pathSync = get.pathSync(imec)
            pathSync = fullfile(imec.pathRoot, imec.fileSync);
        end
        
        function tf = get.hasSync(imec)
            tf = exist(imec.pathSync, 'file') == 2;
        end
        
        function fs = get.fsSync(imec)
            % defer to fsAP or fsLF, or to stored value if neither sources the sync signal
            if ~isempty(imec.fsSync)
                fs = imec.fsSync;
            elseif imec.syncInAPFile
                fs = imec.fsAP;
            elseif imec.syncInLFFile
                fs = imec.fsLF;
            else
                fs = NaN;
            end
        end         
        
        function tf = get.hasSourceSync(imec)
            tf = imec.hasSourceDatasets && all([imec.sourceDatasets.hasSync]);
        end

        function fileSyncCached = get.fileSyncCached(imec)
            fileSyncCached = [imec.fileStem '.sync.mat'];
        end

        function pathSyncCached = get.pathSyncCached(imec)
            pathSyncCached = fullfile(imec.pathRoot, imec.fileSyncCached);
        end

        function scale = get.apScaleToUv(imec)
            if isempty(imec.apRange)
                scale = NaN;
            else
                scale = (imec.apRange(2) - imec.apRange(1)) / (2^imec.adcBits) / imec.apGain * 1e6;
            end
        end

        function scale = get.lfScaleToUv(imec)
            if isempty(imec.lfRange)
                scale = NaN;
            else
                scale = (imec.lfRange(2) - imec.lfRange(1)) / (2^imec.adcBits) / imec.lfGain * 1e6;
            end
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

            meta.fileSizeBytes = imec.bytesPerSample * imec.nChannels * imec.nSamplesAP;
            meta.fileTimeSec = imec.nSamplesAP / imec.fsAP;
            meta.syncBitNames = imec.syncBitNames;
            meta.badChannels = imec.badChannels;
        end
        
        function tf = get.hasSourceDatasets(imec)
            tf = ~isempty(imec.sourceDatasets);
        end
        
        function tf = get.hasSourceAP(imec)
            tf = imec.hasSourceDatasets && all([imec.sourceDatasets.hasAP]);
        end
        
        function tf = get.hasSourceLF(imec)
            tf = imec.hasSourceDatasets && all([imec.sourceDatasets.hasLF]);
        end

        function writeModifiedAPMeta(imec, varargin)
            p = inputParser();
            p.addParameter('filename', imec.pathAPMeta, @isstringlike);
            p.addParameter('extraMeta', struct(), @isstruct);
            p.addParameter('prefixTildeFields', ["imroTbl", "snsChanMap", "snsShankMap"], @isstring); % these fields should have a tilde before them, ignoring this broke certain Python tools like Neo
            p.parse(varargin{:});

            meta = imec.generateModifiedAPMeta();

            % set extra user provided fields
            extraMeta = p.Results.extraMeta;
            extraMetaFields = string(fieldnames(extraMeta));
            for iFld = 1:numel(extraMetaFields)
                meta.(extraMetaFields{iFld}) = extraMeta.(extraMetaFields{iFld});
            end

            filename = string(p.Results.filename);
            Neuropixel.writeINI(filename, meta, 'prefixTildeFields', p.Results.prefixTildeFields);
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
            p.addParameter('band', 'ap', @Neuropixel.Utils.isstringlike);
            p.addParameter('rmsRange', [3 100], @isvector);
            p.addParameter('sampleMaskFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % sampleMaskFn(data_ch_x_time, sample_idx_time) --> logical_time mask of time samples valid for use, useful if you have artifacts at known times
            p.addParameter('printMessage', false, @islogical);
            p.parse(varargin{:});

            rmsByChannel = imec.computeRMSByChannel('band', p.Results.band, 'sampleMaskFn', p.Results.sampleMaskFn);
            rmsMin = p.Results.rmsRange(1);
            rmsMax = p.Results.rmsRange(2);
            rmsBadMask = rmsByChannel < rmsMin | rmsByChannel > rmsMax;

            if p.Results.printMessage
                fprintf('Marking %d / %d channels bad with RMS outside [%g %g] uV\n', nnz(rmsBadMask), numel(rmsByChannel), rmsMin, rmsMax);
            end
            
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
        
        function imecOut = saveTransformedDataset(imec, outPath, varargin)
            p = inputParser();
            p.addParameter('stem', "", @Neuropixel.Utils.isstringlike);
            
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
            p.addParameter('chunkEdgeExtraSamplesAP', [0 0], @isvector); 
            p.addParameter('chunkEdgeExtraSamplesLF', [0 0], @isvector); 
            
            
            p.addParameter('timeShiftsAP', {}, @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            p.addParameter('timeShiftsLF', {}, @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            
            p.addParameter('extraMeta', struct(), @isstruct);
            
            p.addParameter('dryRun', false, @islogical);
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
            
            Neuropixel.writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
            Neuropixel.writeNPY(zeros(nSpikes, 1, 'uint32'), fullfile(savePath, 'spike_templates.npy'));
            Neuropixel.writeNPY(spikeClusters, fullfile(savePath, 'spike_clusters.npy'));
    
            Neuropixel.writeNPY(zeros(nSpikes, 1, 'double'), fullfile(savePath, 'amplitudes.npy'));
            
            templates = zeros(2, nTemplateTimepoints, nCh, 'single');
            Neuropixel.writeNPY(templates, fullfile(savePath, 'templates.npy'));
    
            templatesInds = imec.goodChannels';
            Neuropixel.writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
            sortedInds = imec.goodChannelInds;
            chanMap0ind = int32(imec.channelMap.channelIdsMapped(sortedInds) - uint32(1));
            xcoords = imec.channelMap.xcoords(sortedInds);
            ycoords = imec.channelMap.ycoords(sortedInds);
            Neuropixel.writeNPY(chanMap0ind, fullfile(savePath, 'channel_map.npy'));
            Neuropixel.writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));
    
            templateFeatures = zeros([nTemplates nTemplateFeatures], 'single');
            Neuropixel.writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
            
            templateFeatureInds = zeros(nTemplates, nTemplateFeatures, 'uint32');
            Neuropixel.writeNPY(templateFeatureInds, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
            
            similarTemplates = zeros(nTemplates, nTemplates, 'single');
            Neuropixel.writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
            
            pcFeatures = zeros([nSpikes, nFeaturesPerChannel, nPCFeatures], 'single');
            Neuropixel.writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
            
            pcFeatureInds = zeros([nTemplates, nPCFeatures], 'uint32');
            Neuropixel.writeNPY(pcFeatureInds, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
            whiteningMatrix = ones(nCh, nCh, 'double');
            Neuropixel.writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
            whiteningMatrixInv = ones(nCh, nCh, 'double');
            Neuropixel.writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
            
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
                data_pre = data;
                
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

                changed = ~isequal(data, data_pre);
                if ~dryRun && changed
                    % slice off extra at edges
                    data = data(:, keepIdx);
                    mm.Data.x(chInds, source_idx) = data;
                end
                prog.increment();
            end
            prog.finish();

            if ~dryRun
                imec.writeModifiedAPMeta('extraMeta', p.Results.extraMeta);
                imec.clearSyncCached();
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
            if isempty(imecList)
                error('No imecList specified');
            end
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

            [parent, leaf, ext] = fileparts(char(file));
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

            [pathRoot, fileStem, fileTypeAP, imecNumber] = Neuropixel.ImecDataset.parseImecFileName(file);
            if isnan(imecNumber)
                pathAP = fullfile(pathRoot, [fileStem '.imec.' fileTypeAP '.bin']);
                pathAPMeta = fullfile(pathRoot, [fileStem '.imec.ap.meta']);
            else
                pathAP = fullfile(pathRoot, [fileStem sprintf('.imec%d.', imecNumber) fileTypeAP '.bin']);
                pathAPMeta = fullfile(pathRoot, [fileStem sprintf('.imec%d.ap.meta', imecNumber)]);
            end
            tf = exist(pathAP, 'file') && exist(pathAPMeta, 'file');
            if exist(pathAP, 'file') && ~exist(pathAPMeta, 'file')
                warning('Found data file %s but not meta file %s', pathAP, pathAPMeta);
            end
        end
        
        function [tf, candidates] = pathPointsToSingleValidDataset(fileOrFileStem, type)
            if nargin < 2
                type = "any";
            end
            type = string(type);
            
            if type == "any"
                candidates = Neuropixel.ImecDataset.findImecFileInDir(fileOrFileStem, 'ap', true, false);
                if numel(candidates) == 0
                    candidates = Neuropixel.ImecDataset.findImecFileInDir(fileOrFileStem, 'lf', true, false);
                end
            else
                candidates = Neuropixel.ImecDataset.findImecFileInDir(fileOrFileStem, type, true, false);
            end
            tf = numel(candidates) == 1;
        end

        function file = findImecFileInDir(fileOrFileStem, search_type, returnMultiple, errorIfNotFound)
            if nargin < 2
                search_type = 'ap';
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
                switch search_type
                    case 'ap'
                        if ismember(type, {'ap', 'ap_CAR'})
                            return;
                        else
                            % wrong type
                            file = [];
                        end
                        %assert(ismember(type, {'ap', 'ap_CAR'}), 'Specify ap.bin or ap_CAR.bin file rather than %s file', search_type);
                    case 'lf'
                        if ismember(type, {'lf'})
                            return;
                        else
                            % wrong type
                            file = [];
                        end
%                         assert(ismember(type, {'lf'}), 'Specify lf.bin file rather than %s file', search_type);
                end

            elseif exist(fileOrFileStem, 'dir')
                % it's a directory, assume only one imec file in directory
                path = fileOrFileStem;
%                 [~, leaf] = fileparts(path);

                switch search_type
                    case 'ap'
                        apFiles = Neuropixel.ImecDataset.listAPFilesInDir(path);
                        if ~isempty(apFiles)
                            if numel(apFiles) > 1
%                                 [tf, idx] = ismember([leaf '.ap.bin'], apFiles);
%                                 if tf
%                                     file = apFiles{idx};
%                                     return
%                                 end
%                                 [tf, idx] = ismember([leaf '.imec.ap_CAR.bin'], apFiles);
%                                 if tf
%                                     file = apFiles{idx};
%                                     return
%                                 end
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
%                                 [tf, idx] = ismember([leaf '.imec.lf.bin'], lfFiles);
%                                 if tf
%                                     file = lfFiles{idx};
%                                     return
%                                 end
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
                switch search_type
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
                        error('No %s matches for %s* exist', search_type, fileOrFileStem);
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
                        error('Multiple %s matches for %s* exist, none exactly matches %s. Narrow down the prefix.', search_type, fileOrFileStem, exactName);
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

        function [pathRoot, fileStem, type, imecNumber] = parseImecFileName(file)
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


            match = regexp(file, '(?<stem>[\w\-\.]+).imec(?<imecNumber>\d*).(?<type>\w+).bin', 'names', 'once');
            if ~isempty(match)
                type = match.type;
                fileStem = match.stem;
                if isempty(match.imecNumber)
                    imecNumber = NaN;
                else
                    imecNumber = str2double(match.imecNumber);
                end
                return;
            end

            fileStem = file;
            type = '';
            imecNumber = NaN;
        end

        function apFiles = listAPFilesInDir(path)
            info = cat(1, dir(fullfile(path, '*.ap.bin')), dir(fullfile(path, '*.ap_CAR.bin')));
            apFiles = {info.name}';
        end

        function lfFiles = listLFFilesInDir(path)
            info = dir(fullfile(path, '*.lf.bin'));
            lfFiles = {info.name}';
        end
        
        function checkClearDestinationStem(outPathStem)
            outPathStem = char(outPathStem);
            assert(~isempty(outPathStem));
            files = dir([outPathStem '*']);
            
            tf = false;
            for iF = 1:numel(files)
                if ~files(iF).isdir
                    f = fullfile(files(iF).folder, files(iF).name);
                    warning('This operation would overwrite file %s', f);
                    tf = true;
                end
            end
            if tf
                error('Refusing to proceed as one or more files would be overwritten by this operation');
            end
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

        function [imecOut, transformAPExtraArg, transformLFExtraArg] = writeConcatenatedFileMatchGains(imecList, outPath, varargin)
            p = inputParser();
            p.addParameter('stem', "", @Neuropixel.Utils.isstringlike);
            p.addParameter('writeAP', true, @islogical);
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('mappedChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            p.addParameter('writeSyncSeparate', false, @islogical); % true means ap will get only mapped channels, false will preserve channels as is
            p.addParameter('writeLF', false, @islogical);
            p.addParameter('chunkSize', 2^20, @isscalar);
            p.addParameter('chunkEdgeExtraSamplesAP', [0 0], @isvector); 
            p.addParameter('chunkEdgeExtraSamplesLF', [0 0], @isvector); 
            
            p.addParameter('gpuArray', false, @islogical);
            p.addParameter('applyScaling', false, @islogical); % convert to uV before processing

            p.addParameter('transformAP', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (imec, dataChunk) and return dataChunk someplace
            p.addParameter('transformLF', {}, @(x) iscell(x) || isa(x, 'function_handle')); % list of transformation functions that accept (imec, dataChunk) and return dataChunk someplace
            p.addParameter('transformAPExtraArg', struct(), @(x) true);
            p.addParameter('transformLFExtraArg', struct(), @(x) true);
            p.addParameter('timeShiftsAP', {}, @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            p.addParameter('timeShiftsLF', {}, @(x) isempty(x) || isa(x, 'Neuropixel.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            
            p.addParameter('extraMeta', struct(), @isstruct);
            p.addParameter('dryRun', false, @islogical);
            p.parse(varargin{:});

            nFiles = numel(imecList);
            assert(nFiles > 0);
            stemList = cellfun(@(imec) imec.fileStem, imecList, 'UniformOutput', false);
            dryRun = p.Results.dryRun;
            writeAP = p.Results.writeAP;
            writeLF = p.Results.writeLF;            
            
            function s = lastFilePart(f)
                [~, f, e] = fileparts(f);
                s = [f, e];
            end

            [parent, leaf, ext] = Neuropixel.ImecDataset.filepartsMultiExt(outPath);
            if strlength(ext) > 0 && endsWith(ext, 'bin')
                % specified full file
                outPath = parent;
            else
                outPath = fullfile(parent, [leaf, ext]);
            end
            if ~exist(outPath, 'dir') && ~dryRun
                Neuropixel.Utils.mkdirRecursive(outPath);
            end
            
            if p.Results.stem ~= ""
                leaf = char(p.Results.stem);
            else
                leaf = char(leaf);
            end

            % figure out which channels to keep
            [chIndsByFile, ~] = Neuropixel.ImecDataset.multiFile_build_channelSelectors_internal(imecList, 'goodChannelsOnly', p.Results.goodChannelsOnly, ...
                'connectedChannelsOnly', p.Results.connectedChannelsOnly, 'mappedChannelsOnly', p.Results.mappedChannelsOnly);

            chunkSize = p.Results.chunkSize;
            chunkEdgeExtraSamplesAP = p.Results.chunkEdgeExtraSamplesAP;
            chunkEdgeExtraSamplesLF = p.Results.chunkEdgeExtraSamplesLF;
            
            useGpuArray = p.Results.gpuArray;
            applyScaling = p.Results.applyScaling;
            
            % ensure all source files exist and all sampling rates match
            if writeAP
                assert(all(cellfun(@(imec) imec.hasAP, imecList)), 'All imecDatasets must have AP band to write concatenated AP');
                fsAP = cellfun(@(imec) imec.fsAP, imecList);
                assert(all(fsAP == fsAP(1)), 'All imecDatasets must have matching fsAP');
            end
            fsAP = imecList{1}.fsAP;
           
            if writeLF
                assert(all(cellfun(@(imec) imec.hasLF, imecList)), 'All imecDatasets must have LF band to write concatenated LF');
                fsLF = cellfun(@(imec) imec.fsLF, imecList);
                assert(all(fsLF == fsLF(1)), 'All imecDatasets must have matching fsLF');
            end
            fsLF = imecList{1}.fsLF;
            
            timeShiftsAP = p.Results.timeShiftsAP;
            timeShiftsLF = p.Results.timeShiftsLF;
            % exchange time shifts between LF and AP to ensure consistency
            if writeAP && isempty(timeShiftsAP) && ~isempty(timeShiftsLF)
                assert(~isnan(fsLF) && ~isnan(fsAP));
                nSamplesAPList = cellfun(@(imec) imec.nSamplesAP, imecList);
                timeShiftsAP = arrayfun(@(tsLF, maxSamples) tsLF.convertToDifferentSampleRate(fsLF, fsAP, maxSamples), timeShiftsLF, nSamplesAPList);
            elseif writeLF && isempty(timeShiftsLF) && ~isempty(timeShiftsAP)
                assert(~isnan(fsLF) && ~isnan(fsAP));
                nSamplesLFList = cellfun(@(imec) imec.nSamplesLF, imecList);
                timeShiftsLF = arrayfun(@(tsAP, maxSamples) tsAP.convertToDifferentSampleRate(fsAP, fsLF), timeShiftsAP, nSamplesLFList);
            elseif ~isempty(timeShiftsAP) && ~isempty(timeShiftsLF)
                warning('Both timeShiftsAP and timeShiftsLF are specified, which may lead to inconsistency in the concatenated output bands');
            end
        
            transformAPExtraArg = p.Results.transformAPExtraArg;
            transformLFExtraArg = p.Results.transformLFExtraArg;
            % throw an error if there are any collisions
            Neuropixel.ImecDataset.checkClearDestinationStem(fullfile(outPath, leaf));
            
            if ~dryRun
                % dont do this - dangerous
                %Neuropixel.ImecDataset.clearDestinationStem(fullfile(outPath, leaf));
            end
            
            if p.Results.writeAP || ~isempty(p.Results.transformAP)
                gains = cellfun(@(imec) imec.apGain, imecList);
                [multipliers, gain] = Neuropixel.ImecDataset.determineCommonGain(gains);

                outFile = fullfile(outPath, [leaf '.imec.ap.bin']);
                metaOutFile = fullfile(outPath, [leaf '.imec.ap.meta']);
                
                complainIfExtant(outFile);
                complainIfExtant(metaOutFile);

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
                if ~isfield(meta, 'badChannels')
                    meta.badChannels = imecList{1}.badChannels;
                end
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
                transformAPExtraArg = writeCatFile(outFile, chIndsByFile, 'ap', multipliers, chunkSize, chunkEdgeExtraSamplesAP, ...
                    p.Results.transformAP, timeShiftsAP, dryRun, transformAPExtraArg);
            end

            if p.Results.writeLF || ~isempty(p.Results.transformLF)
                gains = cellfun(@(imec) imec.lfGain, imecList);
                [multipliers, gain] = Neuropixel.ImecDataset.determineCommonGain(gains);

                outFile = fullfile(outPath, [leaf '.imec.lf.bin']);
                metaOutFile = fullfile(outPath, [leaf '.imec.lf.meta']);
                
                complainIfExtant(outFile);
                complainIfExtant(metaOutFile);

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
%                 if isConcatenation
                    meta.concatenated = strjoin(stemList, ':');
                    meta.concatenatedSamples = cellfun(@(imec) imec.nSamplesLF, imecList);
                    meta.concatenatedGains = gains;
                    meta.concatenatedMultipliers = multipliers;
                    meta.concatenatedAdcBits = cellfun(@(imec) imec.adcBits, imecList);
                    meta.concatenatedAiRangeMin = cellfun(@(imec) imec.lfRange(1), imecList);
                    meta.concatenatedAiRangeMax = cellfun(@(imec) imec.lfRange(2), imecList);
%                 end
                
                if ~isempty(timeShiftsLF)
                    % log time shifts by file in meta
                    meta.concatenatedTimeShifts = strjoin(arrayfun(@(shift) shift.as_string(), timeShiftsLF), '; ');
                end

                % compute union of badChannels
                if ~isfield(meta, 'badChannels')
                    meta.badChannels = imecList{1}.badChannels;
                end
                for iM = 2:numel(imecList)
                    meta.badChannels = union(meta.badChannels, imecList{iM}.badChannels);
                end
                
                fprintf('Writing LF meta file %s\n', lastFilePart(metaOutFile));
                if ~dryRun
                    Neuropixel.writeINI(metaOutFile, meta);
                end

                fprintf('Writing LF bin file %s\n', lastFilePart(outFile));
                transformLFExtraArg = writeCatFile(outFile, chIndsByFile, 'lf', multipliers, chunkSize, chunkEdgeExtraSamplesLF, ...
                    p.Results.transformLF, timeShiftsLF, dryRun, transformLFExtraArg);
            end

            if p.Results.writeSyncSeparate
                outFile = fullfile(outPath, [leaf '.imec.sync.bin']);
                fprintf('Writing separate sync bin file %s', lastFilePart(outFile));
                transformAPExtraArg = writeCatFile(outFile, imecList{1}.syncChannelIndex, 'sync', ones(nFiles, 1, 'int16'), chunkSize, chunkEdgeExtraSamplesAP, ...
                    {}, timeShiftsAP, dryRun, transformAPExtraArg);
            end

            outFile = fullfile(outPath, leaf);
            imecOut = Neuropixel.ImecDataset(outFile, 'channelMap', imecList{1}.channelMapFile, 'sourceDatasets', cat(1, imecList{:}));

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

                        % some elements of source_idx may be 0 (meaning they are not filled with source data)
                        source_idx = sourceIdxList(idx);
                        mask_source_idx = source_idx > 0;
                        
                        data = zeros(numel(chInds), numel(source_idx), 'int16');
                        data(:, mask_source_idx) = mm.Data.x(chInds, source_idx(mask_source_idx));

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
                                    extraArgs{end+1} = transformExtraArg; %#ok<AGROW>
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
            
            function complainIfExtant(dest_file)
                if exist(dest_file, 'file')
                    error('File %s would be overwritten by this operation', dest_file);
                end
            end
        end
        
        % determine the gains that will equalize scale across multiple datasets
        function [multipliers, gain] = determineCommonGain(gains, quiet)
            if nargin < 2
                quiet = false;
            end
            
            nFiles = numel(gains);
            uGains = unique(gains);

            if numel(uGains) == 1
                gain = uGains;
                multipliers = ones(nFiles, 1);
                if numel(gains) > 1 && ~quiet
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
                if ~quiet
                    fprintf('Converting all files to gain of %d\n', gain);
                end
            end

            multipliers = int16(multipliers);
        end
    end
end
