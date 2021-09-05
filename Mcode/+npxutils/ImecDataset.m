classdef ImecDataset < handle
    % An IMEC data set
    %
    % An IMEC data set is a set of data from a collection of related IMEC files.
    % This includes the raw data files and metadata files.
    %
    % The files in an IMEC data set must all be in the same directory on disk.
    %
    % This is a handle object, so it is pass-by-reference.
    
    properties (Constant)
        % Number of bytes per sample.
        bytesPerSample = 2;
    end
    
    properties(SetAccess = protected)
        % Full path to the parent directory of the files in this data set.
        pathRoot char = '';
        % Common base name of the files in this data set.
        fileStem char = '';
        fileImecNumber = NaN;
        % Data set creation time, as a datenum.
        creationTime = NaN;
        % Number of channels.
        nChannels = NaN;
        
        % Typically ap or ap_CAR
        fileTypeAP = 'ap';
        % Typically lf or lf_CAR
        fileTypeLF = 'lf';
        
        % Number of samples in the AP file.
        nSamplesAP = 0;
        % Number of samples in the LF file.
        nSamplesLF = 0;
        % ??? Maybe the frequency of samples in the AP file, in samples per second?
        fsAP = NaN;
        % ??? Maybe the frequency of samples in the LF file, in samples per second?
        fsLF = NaN;
        fsSync = NaN;
        highPassFilterHz = NaN;
        apGain = NaN;
        apRange = [];
        lfGain = NaN;
        lfRange = []
        
        adcBits = 10;
        
        snsShankMap (1,1) string;
        % Can be stored using set channelMap
        channelMap = [];
        
        % See markBadChannels
        badChannels (:, 1) uint32
        
        syncBitNames string;
        
        concatenationInfoAP
        concatenationInfoLF
        
        % Optional list of concatenated ImecDatasets that sourced the data for this file (via concatenationInfoAP)
        sourceDatasets
    end
    
    properties
        % Will be cached after loading, can also be cleared by user.
        syncRaw uint16 = [];
    end
    
    properties (Dependent)
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
        % Number of channels in the channel map (excludes sync)
        nChannelsMapped
        
        connectedChannels
        connectedChannelInds
        nChannelsConnected % Excludes reference and sync channels
        
        goodChannels % Connected channels sans badChannels
        goodChannelInds
        nGoodChannels
        
        channelIds % List of ids from ChannelMap
        channelNames % Full list of channel names
        channelNamesPadded
        
        nSyncBits
        syncBitsNamed
        
        fileAP % The .imec.ap.bin file, without folder
        pathAP % The .imec.ap.bin file, with folder
        fileAPMeta % The AP metadata file, without folder
        pathAPMeta % The AP Metadata file, with foler
        
        fileLF % The .imec.lf.bin file, without folder
        pathLF % The .imec.lf.bin file, with folder
        fileLFMeta % The .imec.lf.meta file, without folder
        pathLFMeta % The .imec.lf.meta file, with folder
        
        % The sync file, which may be a separate .imec.sync.bin file, without folder
        fileSync
        % The sync file, which may be a separate .imec.sync.bin file, with folder
        pathSync
        
        % After sync is cached to sync.mat file for faster reload, without folder
        fileSyncCached
        % After sync is cached to sync.mat file for faster reload, with folder
        pathSyncCached
        
        % Creation time, as a human-readable date string
        creationTimeStr
        
        % Multiply raw int16 by this to get uV
        apScaleToUv
        lfScaleToUv
        
        % From channel map (although syncChannelIndex will be 1 if sync not in AP or LF file)
        syncChannelId
        % If sync in AP file, at one index
        syncChannelIndex
        
        % Sync will come from either AP, LF, or separate file depending on these flags
        
        % Whether sync is in the ap file
        syncInAPFile
        % Whether sync is in the lf file
        syncInLFFile
    end
    
    methods
        function this = ImecDataset(fileOrFileStem, varargin)
            % Construct a new ImecDataset from a set of files.
            %
            % obj = npxutils.ImecDataset(fileOrFileStem, [Options ...])
            %
            % Constructing an ImecDataset object opens and reads an IMEC data
            % set from files on disk.
            %
            % FileOrFileStem is a scalar string specifying a file or file stem
            % that identifies the IMEC data set to load. It may be:
            %   - The .imec.ap.bin file
            %   - The base name of the .imec.* files (that is, the file name
            %     with no extension)
            %   - A common prefix of the base name of the file, as long as there
            %     is is no ambiguity with other files in that directory.
            %   - The parent directory of the data set, as long as there is only
            %     one .ap.bin file inside that directory.
            %
            % Options is a list of name/value options. Valid options:
            %
            %   channelMap - ???
            %   syncBitNames - (empty or string) ???
            %   sourceDatasets - ???
            
            p = inputParser();
            p.addParameter('channelMap', [], @(x) true);
            p.addParameter('syncBitNames', [], @(x) isempty(x) || isstring(x) || iscellstr(x));
            p.addParameter('sourceDatasets', [], @(x) true);
            p.parse(varargin{:})
            
            fileOrFileStem = char(fileOrFileStem);
            file = npxutils.ImecDataset.findImecFileInDir(fileOrFileStem, 'ap', true, false);
            if isempty(file)
                file = npxutils.ImecDataset.findImecFileInDir(fileOrFileStem, 'lf', false, false);
                if isempty(file)
                    error('No AP or LF Imec file found at or in %s', fileOrFileStem);
                else
                    isLFOnly = true;
                    [this.pathRoot, this.fileStem, this.fileTypeLF, this.fileImecNumber] = npxutils.ImecDataset.parseImecFileName(file);
                end
            elseif numel(file) == 1
                [this.pathRoot, this.fileStem, this.fileTypeAP, this.fileImecNumber] = npxutils.ImecDataset.parseImecFileName(file);
                isLFOnly = false;
            else
                for iF = 1:numel(file)
                    fprintf('Possible match: %s\n', file{iF});
                end
                error('Multiple imec datasets found in specified location, include file stem or a full path to refine the search');
            end
            
            if ~isLFOnly
                if exist(this.pathAP, 'file')
                    if ~exist(this.pathAPMeta, 'file')
                        error('Could not find AP meta file %s', this.pathAPMeta);
                    end
                    this.readInfo();
                else
                    error('Could not find AP bin file %s', this.pathAP);
                end
            else
                if exist(this.pathLF, 'file')
                    if ~exist(this.pathLFMeta, 'file')
                        error('Could not find LF meta file %s', this.pathLFMeta);
                    end
                    this.readInfo();
                else
                    error('Could not find LF bin file %s', this.pathLF);
                end
            end
            
            channelMap = p.Results.channelMap;
            if isempty(channelMap)
                channelMap = "";
            elseif ischar(channelMap)
                channelMap = string(channelMap);
            end
            
            if isa(channelMap, 'npxutils.ChannelMap')
                % manually specified channel map
                this.channelMap = channelMap;
            elseif isstring(channelMap)
                if channelMap == ""
                    % use default channel map file
                    channelMap = npxutils.util.getDefaultChannelMapFile(true);
                end
                this.channelMap = npxutils.ChannelMap(channelMap);
            end
            assert(this.channelMap.nChannels <= this.nChannels, ...
                'Channel count (%d) is less than number of channels in channel map (%d)', ...
                this.nChannels, this.channelMap.nChannels);
            
            if ~isempty(p.Results.syncBitNames)
                this.setSyncBitNames(1:numel(p.Results.syncBitNames), p.Resuls.syncBitNames);
            end
            
            if ~isempty(p.Results.sourceDatasets)
                assert(isa(p.Results.sourceDatasets, 'npxutils.ImecDataset'));
                this.setSourceDatasets(p.Results.sourceDatasets);
            end
        end
        
        function readInfo(this)
            % Read info from the files into this object.
            if this.hasAP
                metaAP = this.readAPMeta();
                if this.hasLF
                    metaLF = this.readLFMeta();
                else
                    metaLF = struct;
                end
                meta = metaAP;
            elseif this.hasLF
                metaAP = struct;
                metaLF = this.readLFMeta();
                meta = metaLF;
            else
                error('Must have either AP or LF data files');
            end
            
            this.nChannels = meta.nSavedChans;
            this.creationTime = datenum(meta.fileCreateTime, 'yyyy-mm-ddTHH:MM:SS');
            
            if this.hasAP
                this.fsAP = metaAP.imSampRate;
                if isfield(metaAP, 'imHpFlt')
                    this.highPassFilterHz = metaAP.imHpFlt;
                end
            elseif this.hasSourceDatasets
                this.fsAP = this.sourceDatasets(1).fsAP;
            end
            
            if this.hasLF
                this.fsLF = metaLF.imSampRate;
            elseif this.hasSourceDatasets
                this.fsLF = this.sourceDatasets(1).fsLF;
            end
            
            % parse imroTable
            m = regexp(meta.imroTbl, '\(([\d, ]*)\)', 'tokens');
            gainVals = strsplit(m{2}{1}, ' ');
            this.apGain = str2double(gainVals{4});
            this.lfGain = str2double(gainVals{5});
            
            if this.hasAP
                this.apRange = [metaAP.imAiRangeMin metaAP.imAiRangeMax];
            end
            
            if this.hasLF
                this.lfRange = [metaLF.imAiRangeMin metaLF.imAiRangeMax];
            end
            
            % copy snsShankMap in case needed for building channel map
            this.snsShankMap = meta.snsShankMap;
            
            % look at AP meta fields that might have been set by us
            if isfield(meta, 'badChannels') && ~isempty(meta.badChannels)
                this.badChannels = union(this.badChannels, meta.badChannels);
            end
            if isfield(meta, 'syncBitNames')
                this.setSyncBitNames(1:numel(meta.syncBitNames), meta.syncBitNames);
            end
            
            if this.hasAP
                fid = this.openAPFile();
                fseek(fid, 0, 'eof');
                bytes = ftell(fid);
                fclose(fid);
                
                this.nSamplesAP = bytes / this.bytesPerSample / this.nChannels;
                if round(this.nSamplesAP) ~= this.nSamplesAP
                    warning(['AP bin file size is not an integral number of samples, ' ...
                        'file data may not be fully copied, truncating nSamplesAP']);
                    this.nSamplesAP = floor(this.nSamplesAP);
                end
                
                this.concatenationInfoAP = npxutils.ConcatenationInfo(this, 'ap', metaAP);
            end
            
            if this.hasLF
                fid = this.openLFFile();
                fseek(fid, 0, 'eof');
                bytes = ftell(fid);
                fclose(fid);
                this.nSamplesLF = bytes / this.bytesPerSample / this.nChannels;
                if round(this.nSamplesLF) ~= this.nSamplesLF
                    warning('LF bin file size is not an integral number of samples, file data may not be fully copied, truncating nSamplesLF');
                    this.nSamplesLF = floor(this.nSamplesLF);
                end
                
                this.concatenationInfoLF = npxutils.ConcatenationInfo(this, 'lf', metaLF);
            end
        end
        
        function setSyncBitNames(this, idx, names)
            % Set the names of sync bits.
            %
            % setSyncBitNames(imec, idx, names)
            %
            % Idx is a vector if sync bit indices. Must be in range 1 ..
            % this.nSyncBits.
            %
            % Names is a string array or charvec of names corresponding to the
            % identified sync bits.
            
            % idx is the indices of which bits to set to the corresponding items from names
            assert(all(idx >= 1 & idx <= this.nSyncBits), ...
                'Sync bit indices must be in range [1..%d]', this.nSyncBits);
            if isscalar(idx) && ischar(names)
                this.syncBitNames{idx} = names;
            else
                names = string(names);
                this.syncBitNames(idx) = names;
            end
        end
        
        function setSourceDatasets(this, imecList)
            assert(isempty(imecList) || isa(imecList, 'npxutils.ImecDataset'));
            this.sourceDatasets = imecList;
            
            % attempt to fill in missing metadata as well
            if isnan(this.fsAP)
                this.fsAP = this.sourceDatasets(1).fsAP;
            end
            
            if isnan(this.fsLF)
                this.fsLF = this.sourceDatasets(1).fsLF;
            end
            
            % attempt to cross infer concatenation info if missing
            if this.hasAP && ~this.hasLF && this.hasSourceLF
                % infer LF concatenation info from AP info to enable source LF data access through concatInfo lookup
                if ~isempty(imecList)
                    this.concatenationInfoLF = npxutils.ConcatenationInfo.inferFromOtherBand(this, 'lf', 'ap');
                    this.nSamplesLF = this.concatenationInfoLF.nSamples;
                    this.lfRange = this.concatenationInfoLF.ranges(1, :); % we assume that the ranges are shared across all datasets
                else
                    % clear existing info
                    this.concatenationInfoLF = [];
                    this.nSamplesLF = 0;
                    this.lfRange = [];
                end
            elseif this.hasLF && ~this.hasAP && this.hasSourceAP
                if ~isempty(imecList)
                    % infer AP concatenation info from LF info to enable source AP data access through concatInfo lookup
                    this.concatenationInfoAP = npxutils.ConcatenationInfo.inferFromOtherBand(this, 'ap', 'lf');
                    this.nSamplesAP = this.concatenationInfoAP.nSamples;
                    this.apRange = this.concatenationInfoLF.ranges(1, :); % we assume that the ranges are shared across all datasets
                else
                    % clear existing info
                    this.concatenationInfoAP = [];
                    this.nSamplesAP = 0;
                    this.apRange = [];
                end
            end
        end
        
        function idx = lookupSyncBitByName(this, names, ignoreNotFound)
            if nargin < 3
                ignoreNotFound = false;
            end
            if isnumeric(names)
                idx = names;
            else
                names = string(names);
                [tf, idx] = ismember(names, this.syncBitNames);
                if ignoreNotFound
                    idx(~tf) = NaN;
                elseif any(~tf)
                    error('Sync bit(s) %s not found', strjoin(names, ', '));
                end
            end
        end
        
        function newImec = copyToNewLocation(this, newRoot, newStem)
            if nargin < 3
                newStem = this.fileStem;
            end
            npxutils.internal.mkdirRecursive(newRoot);
            
            f = @(suffix) fullfile(newRoot, [newStem suffix]);
            docopy(this.pathAP, f('.imec.ap.bin'));
            docopy(this.pathAPMeta, f('.imec.ap.meta'));
            docopy(this.pathLF, f('.imec.lf.bin'));
            docopy(this.pathLFMeta, f('.imec.lf.meta'));
            docopy(this.pathSync, f('.imec.sync.bin'));
            
            newImec = npxutils.ImecDataset(fullfile(newRoot, newStem), 'channelMap', this.channelMapFile);
            
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
    %         function data_ch_by_time = readAPChannelBand(this, chFirst, chLast, sampleFirst, sampleLast, msg)
    %             if nargin < 4 || isempty(sampleFirst)
    %                 sampleFirst = 1;
    %             end
    %             if nargin < 5 || isempty(sampleLast)
    %                 sampleLast = this.nSamplesAP;
    %             end
    %             if nargin < 6 || isempty(msg)
    %                 msg = 'Reading channels from neuropixel AP file';
    %             end
    %
    %             data_ch_by_time = this.readChannelBand('ap', chFirst, chLast, sampleFirst, sampleLast, msg);
    %         end
    %
    %         function data_ch_by_time = readLFChannelBand(this, chFirst, chLast, sampleFirst, sampleLast, msg)
    %             if nargin < 4 || isempty(sampleFirst)
    %                 sampleFirst = 1;
    %             end
    %             if nargin < 5 || isempty(sampleLast)
    %                 sampleLast = this.nSamplesLF;
    %             end
    %             if nargin < 6 || isempty(msg)
    %                 msg = 'Reading channels from neuropixel LF file';
    %             end
    %
    %             data_ch_by_time = this.readChannelBand('lf', chFirst, chLast, sampleFirst, sampleLast, msg);
    %         end
    %
    %         function data_by_time = readAPSingleChannel(this, ch, varargin)
    %             data_by_time = this.readAPChannelBand(ch, ch, varargin{:})';
    %         end
    %
    %         function data_by_time = readLFSingleChannel(imec, ch, varargin)
    %             data_by_time = this.readLFChannelBand(ch, ch, varargin{:})';
    %         end
    %     end
    
    methods % Sync channel read / cache
        function [syncRaw, fsSync] = readSync(this, varargin)
            p = inputParser();
            p.addOptional('reload', false, @islogical); % if true, ignore cache in memory imec.syncRaw
            p.addParameter('preferredBand', '', @npxutils.internal.isstringlike); % prefer a band (ap or lf), though will only be obeyed if not already loaded / cached
            p.addParameter('ignoreCached', false, @islogical); % if true, ignore cache on disk in imec.pathSyncCached
            p.parse(varargin{:});
            preferredBand = string(p.Results.preferredBand);
            
            if isempty(this.syncRaw) || p.Results.reload
                if exist(this.pathSyncCached, 'file') && ~p.Results.ignoreCached
                    [~, f, e] = fileparts(this.pathSyncCached);
                    fprintf('Loading sync from cached %s%s\n', f, e);
                    ld = load(this.pathSyncCached);
                    this.syncRaw = typecast(ld.sync, 'uint16');
                    if isfield(ld, 'fsSync')
                        this.fsSync = ld.fsSync;
                    else
                        % assume from AP if not specified since I added sync from LF later
                        this.fsSync = this.fsAP;
                    end
                else
                    % this will automatically redirect to a separate sync file
                    % or to the ap file depending on .syncInAPFile
                    switch preferredBand
                        case {"", "auto"}
                            fprintf('Loading sync channel auto (this will take some time)...\n');
                            [mm, this.fsSync] = this.memmapSync_full();
                            this.syncRaw = typecast(mm.Data.x(this.syncChannelIndex, :)', 'uint16');
                            this.saveSyncCached();
                        case 'ap'
                            fprintf('Loading sync channel from AP band (this will take some time)...\n');
                            mm = this.memmapAP_full();
                            this.syncRaw = typecast(mm.Data.x(this.syncChannelIndex, :)', 'uint16');
                            this.fsSync = this.fsAP;
                            this.saveSyncCached();
                        case 'lf'
                            fprintf('Loading sync channel from LF band (this will take some time)...\n');
                            mm = this.memmapLF_full();
                            this.syncRaw = typecast(mm.Data.x(this.syncChannelIndex, :)', 'uint16');
                            this.fsSync = this.fsLF;
                            this.saveSyncCached();
                        otherwise
                            error('Unknown preferredBand %s', preferredBand);
                    end
                end
            end
            syncRaw = this.syncRaw;
            fsSync = this.fsSync;
        end
        
        function saveSyncCached(this)
            sync = this.readSync();
            fsSync = this.fsSync; %#ok<PROP>
            save(this.pathSyncCached, 'sync', 'fsSync');
        end
        
        function clearSyncCached(this)
            if exist(this.pathSyncCached, 'file') > 0
                fprintf('Deleting cached sync file %s\n', this.pathSyncCached);
                delete(this.pathSyncCached);
            end
            this.syncRaw = [];
            this.fsSync = NaN;
        end
        
        function updateSyncCached(this, varargin)
            this.syncRaw = [];
            this.fsSync = [];
            if exist(this.pathSyncCached, 'file')
                [sync, fsSync] = this.readSync('reload', true, 'ignoreCached', true);
                save(this.pathSyncCached, 'sync', 'fsSync');
                this.syncRaw = sync;
                this.fsSync = fsSync;
            end
        end
        
        function tf = readSyncBit(this, bit)
            bit = this.lookupSyncBitByName(bit);
            tf = logical(bitget(this.readSync(), bit));
        end
        
        function vec = readSync_idx(this, idx, varargin)
            p = inputParser();
            p.addParameter('band', '', @npxutils.internal.isstringlike);
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.parse(varargin{:});
            
            band = string(p.Results.band);
            switch band
                case ""
                    fsRequested = NaN;
                case "auto"
                    fsRequested = NaN;
                case "lf"
                    fsRequested = this.fsLF;
                case "ap"
                    fsRequested = this.fsAP;
                otherwise
                    error('Unknown band %s', band);
            end
            
            fromSource = p.Results.fromSourceDatasets;
            if ~fromSource
                % grab cached data if the sampling rate matches, otherwise use the memmap
                if ~isempty(this.syncRaw) && this.fsSync == fsRequested
                    vec = this.syncRaw(idx);
                elseif ismember(band, ["", "auto"])
                    mm = this.memmapSync_full();
                    vec = mm.Data.x(this.syncChannelIndex, idx)';
                elseif band == "ap"
                    mm = this.memmapAP_full();
                    vec = mm.Data.x(this.syncChannelIndex, idx)';
                elseif band == "lf"
                    mm = this.memmapLF_full();
                    vec = mm.Data.x(this.syncChannelIndex, idx)';
                end
                
            elseif this.syncInAPFile
                mmSet = this.memmap_sourceAP_full();
                [sourceFileInds, sourceSampleInds] = this.concatenationInfoAP.lookup_sampleIndexInSourceFiles(idx);
                vec = npxutils.ImecDataset.multi_mmap_extract_sample_idx(mmSet, sourceFileInds, sourceSampleInds, this.syncChannelIndex);
            elseif this.syncInLFFile
                mmSet = this.memmap_sourceLF_full();
                [sourceFileInds, sourceSampleInds] = this.concatenationInfoLF.lookup_sampleIndexInSourceFiles(idx);
                vec = npxutils.ImecDataset.multi_mmap_extract_sample_idx(mmSet, sourceFileInds, sourceSampleInds, this.syncChannelIndex);
            else
                error('Cannot read source sync unless sync derives from AP or LF');
            end
        end
        
        function mat = readSyncBits_idx(this, bits, idx, varargin)
            % mat is nBits x nTime to match readAP_idx which is nChannels x nTime
            if isstring(bits) || ischar(bits)
                bits = this.lookupSyncBitByName(bits);
            end
            vec = this.readSync_idx(idx, varargin{:});
            mat = false(numel(bits), numel(vec));
            for iB = 1:numel(bits)
                mat(iB, :) = logical(bitget(vec, bits(iB)));
            end
        end
    end
    
    methods
        function sampleIdx = closestSampleAPForTime(this, timeSeconds)
            sampleIdx = round(timeSeconds * this.fsAP);
            sampleIdx(sampleIdx == 0) = 1;
            if any(sampleIdx < 0 | sampleIdx > this.nSamplesAP)
                error('Time seconds out of range');
            end
        end
        
        function sampleIdx = closestSampleLFForTime(this, timeSeconds)
            sampleIdx = round(timeSeconds * this.fsLF);
            sampleIdx(sampleIdx == 0) = 1;
            if any(sampleIdx < 0 | sampleIdx > this.nSamplesLF)
                error('Time seconds out of range');
            end
        end
        
        function idxLF = closestSampleLFForAP(this, idxAP)
            idxLF = floor(double(idxAP-1) * double(this.fsLF) / double(this.fsAP)) + 1;
        end
        
        function idxAP = closestSampleAPForLF(this, idxLF)
            idxAP = floor(double(idxLF-1) * double(this.fsAP) / double(this.fsLF)) + 1;
        end
        
        function data = readAP_idx(this, sampleIdx, varargin)
            data = this.internal_read_idx(sampleIdx, 'band', 'ap', varargin{:});
        end
        
        function data = readLF_idx(this, sampleIdx, varargin)
            data = this.internal_read_idx(sampleIdx, 'band', 'lf', varargin{:});
        end
        
        function data = internal_read_idx(this, sampleIdx, varargin)
            p = inputParser();
            p.addParameter('band', 'ap', @npxutils.internal.isstringlike);
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
                        assert(this.hasAP, 'ImecDataset does not have AP band');
                    else
                        assert(this.hasSourceAP, 'ImecDataset does not have source datasets with AP band');
                    end
                case 'lf'
                    if ~fromSource
                        assert(this.hasLF, 'ImecDataset does not have LF band');
                    else
                        assert(this.hasSourceLF, 'ImecDataset does not have source datasets with LF band');
                    end
                otherwise
                    error('Unknown band %s', band);
            end
            
            ch_conn_mask = this.lookup_channelIds(this.connectedChannels);
            
            if ~fromSource
                if band == "ap"
                    mm = this.memmapAP_full();
                    scaleToUv = this.apScaleToUv;
                else
                    mm = this.memmapLF_full();
                    scaleToUv = this.lfScaleToUv;
                end
                
                if any(isnan(sampleIdx))
                    mask = ~isnan(sampleIdx);
                    data = single(mm.Data.x(:, sampleIdx)); % must be single to support NaNs
                    data = npxutils.internal.TensorUtils.inflateMaskedTensor(data, 2, mask, NaN);
                else
                    data = mm.Data.x(:, sampleIdx);
                end
                
                if p.Results.applyScaling
                    data = single(data);
                    data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(scaleToUv);
                end
            else
                if band == "ap"
                    mmSet = this.memmap_sourceAP_full();
                    [sourceFileInds, sourceSampleInds] = this.concatenationInfoAP.lookup_sampleIndexInSourceFiles(sampleIdx);
                    scalingByFile = this.concatenationInfoAP.scaleToUvs;
                    scaleToUvThis = this.apScaleToUv;
                else
                    mmSet = this.memmap_sourceLF_full();
                    [sourceFileInds, sourceSampleInds] = this.concatenationInfoLF.lookup_sampleIndexInSourceFiles(sampleIdx);
                    scalingByFile = this.concatenationInfoLF.scaleToUvs;
                    scaleToUvThis = this.lfScaleToUv;
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
                data = npxutils.ImecDataset.multi_mmap_extract_sample_idx(mmSet, sourceFileInds, sourceSampleInds, this.channelIds, scalingByFile, ch_conn_mask);
            end
        end
        
        function [mat, sampleIdx] = readAP_timeWindow(this, timeWindowSec, varargin)
            idxWindow = this.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = this.readAP_idx(sampleIdx, varargin{:});
        end
        
        function [mat, sampleIdx] = readLF_timeWindow(this, timeWindowSec, varargin)
            idxWindow = this.closestSampleLFForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = this.readLF_idx(sampleIdx, varargin{:});
        end
        
        function [mat, sampleIdx] = readSyncBits_timeWindow(this, bits, timeWindowSec)
            idxWindow = this.closestSampleAPForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = readSyncBits_idx(bits, sampleIdx);
        end
        
        function [mat, sampleIdx] = readSyncLFBits_timeWindow(this, bits, timeWindowSec)
            idxWindow = this.closestSampleLFForTime(timeWindowSec);
            sampleIdx = idxWindow(1):idxWindow(2);
            mat = readSyncLFBits_idx(bits, sampleIdx);
        end
        
        function [mat, sourceInds] = readAP_viaTimeShiftSpec(this, timeShifts, varargin)
            % read a matrix of AP data as requested by npxutils.TimeShiftSpec instance timeShifts
            sourceInds = timeShifts.computeSourceIndices(this.nSamplesAP);
            mat = this.readAP_idx(sourceInds, varargin{:});
        end
    end
    
    methods % Quick inspection
        function [channelInds, channelIds] = lookup_channelIds(this, channelIds)
            if islogical(channelIds)
                channelIds = this.channelIdx(channelIds);
            end
            [tf, channelInds] = ismember(channelIds, this.channelIds);
            assert(all(tf, 'all'), 'Some channel ids not found');
        end
        
        function plotAP_timeWindow(this, timeWindowSec, varargin)
            idxWindow = this.closestSampleAPForTime(timeWindowSec);
            this.plotAP_idxWindow(idxWindow, 'timeInSeconds', true, varargin{:});
        end
        
        function plotAP_idxWindow(this, idxWindow, varargin)
            this.internal_plotIdxWindow(idxWindow, 'band', 'ap', varargin{:});
        end
        
        function plotLF_idxWindow(this, idxWindow, varargin)
            this.internal_plotIdxWindow(idxWindow, 'band', 'lf', varargin{:});
        end
        
        function plotLF_timeWindow(this, timeWindowSec, varargin)
            idxWindow = this.closestSampleLFForTime(timeWindowSec);
            this.plotLF_idxWindow(idxWindow, 'timeInSeconds', true, varargin{:});
        end
        
        function internal_plotIdxWindow(this, idxWindow, varargin)
            p = inputParser();
            p.addParameter('band', 'ap', @ischar);
            p.addParameter('channels', this.mappedChannels, @(x) isempty(x) || isvector(x));
            p.addParameter('invertChannels', this.channelMap.invertChannelsY, @islogical);
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            p.addParameter('showSync', true, @isvector);
            p.addParameter('syncBits', this.syncBitsNamed, @isvector);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('gain', 0.95, @isscalar);
            p.addParameter('car', false, @islogical); % subtract median of all channels at each time
            p.addParameter('center', false, @islogical); % subtract median of each channel over time
            p.addParameter('fromSourceDatasets', false, @islogical);
            p.addParameter('syncFromSourceDatasets', [], @(x) isempty(x) || islogical(x));
            p.addParameter('downsample',1, @isscalar);
            p.addParameter('timeInSeconds', false, @islogical);
            p.addParameter('timeRelativeTo', 0, @isscalar);
            p.addParameter('tsi', [], @(x) isempty(x) || isa(x, 'npxutils.TrialSegmentationInfo')); % to mark trial boundaries
            
            p.addParameter('markSampleIdx', [], @isvector);
            p.addParameter('markSampleMode', 'rug', @ischar);
            p.addParameter('markSampleColor', [0.5 0 0.5], @(x) true);
            
            p.addParameter('style', 'traces', @npxutils.internal.isstringlike);
            
            p.parse(varargin{:});
            
            % by default, sync comes from same source v. processed data as the data being plotted,
            % but this can be overrriden
            fromSource = p.Results.fromSourceDatasets;
            syncFromSource = p.Results.syncFromSourceDatasets;
            if isempty(syncFromSource)
                syncFromSource = fromSource;
            end
            
            band = string(p.Results.band);
            switch band
                case 'ap'
                    fsBand = this.fsAP;
                case 'lf'
                    fsBand = this.fsLF;
                otherwise
                    error('Unknown band %s', band);
            end
            
            if numel(idxWindow) > 2
                idxWindow = [idxWindow(1), idxWindow(end)];
            end
            sampleIdx = idxWindow(1):idxWindow(2);
            
            sampleIdx = floor(sampleIdx);
            
            mat = this.internal_read_idx(sampleIdx, 'band', band, 'fromSourceDatasets', fromSource, 'applyScaling', true); % C x T
            labels = this.channelNamesPadded;
            
            [channelInds, channelIds] = this.lookup_channelIds(p.Results.channels); %#ok<*PROPLC>
            if p.Results.goodChannelsOnly
                mask = ismember(channelIds, this.goodChannels);
                channelInds = channelInds(mask);
                channelIds = channelIds(mask);
            end
            if p.Results.connectedChannelsOnly
                mask = ismember(channelIds, this.connectedChannels);
                channelInds = channelInds(mask);
                channelIds = channelIds(mask);
            end
            
            mat = mat(channelInds, :);
            labels = labels(channelInds);
            connected = ismember(channelIds, this.connectedChannels);
            bad = ismember(channelIds, this.badChannels);
            
            if p.Results.downsample > 1
                mat = mat(:, 1:p.Results.downsample:end);
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
            if (syncFromSource && ~this.hasSourceSync) || (~syncFromSource && ~this.hasSync)
                warning('Cannot show sync data since no sync data present');
                showSync = false;
            end
            syncBits = p.Results.syncBits;
            if ~isempty(syncBits) && showSync
                % make sure we grab sync from the matching band
                syncBitMat = this.readSyncBits_idx(syncBits, sampleIdx, 'fromSourceDatasets', syncFromSource, 'band', band);
                mat = cat(1, mat, syncBitMat);
                syncColor = [0.75 0 0.9];
                colors = cat(1, colors, repmat(syncColor, size(syncBitMat, 1), 1));
                labels = cat(1, labels, this.syncBitNames(syncBits));
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
                    npxutils.internal.plotStackedTraces(time, mat', 'colors', colors, 'labels', labels, ...
                        'gain', p.Results.gain, 'invertChannels', false, ...
                        'normalizeMask', normalizeMask, 'normalizeEach', false);
                case 'pmatbal'
                    if p.Results.invertChannels
                        mat = flipud(mat);
                        labelsF = flipud(labels);
                    else
                        labelsF = labels;
                    end
                    ytick = 1:size(mat, 1);
                    npxutils.internal.pmatbal(mat, 'x', time, 'y', ytick);
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
                    npxutils.internal.graphics.rugplot(markTimes, 'side', 'top', 'Color', p.Results.markSampleColor, 'expand_limits', true);
                else
                    for iM = 1:numel(markTimes)
                        xline(markTimes(iM), 'Color', p.Results.markSampleColor);
                    end
                end
            end
            
            hold off;
        end
        
        function plotSync_idxWindow(this, idxWindow, varargin)
            this.plotAP_idxWindow(idxWindow, 'channels', [], 'showSync', true, varargin{:});
        end
    end
    
    methods % Memory mapped read/write access to data
        function mm = memmapAP_by_sample(this)
            mm = memmapfile(this.pathAP, 'Format', {'int16', [this.nChannels 1], 'x'}, ...
                'Repeat', this.nSamplesAP);
        end
        
        function mm = memmapLF_by_sample(this)
            mm = memmapfile(this.pathLF, 'Format', {'int16', [this.nChannels 1], 'x'}, ...
                'Repeat', this.nSamplesLF);
        end
        
        function mm = memmapAP_by_chunk(this, nSamplesPerChunk)
            mm = memmapfile(this.pathAP, 'Format', {'int16', [this.nChannels nSamplesPerChunk], 'x'}, ...
                'Repeat', floor(this.nSamplesAP/nSamplesPerChunk));
        end
        
        function mm = memmapLF_by_chunk(this, nSamplesPerChunk)
            mm = memmapfile(this.pathLF, 'Format', {'int16', [this.nChannels nSamplesPerChunk], 'x'}, ...
                'Repeat', floor(this.nSamplesLF/nSamplesPerChunk));
        end
        
        function mm = memmapAP_full(this, varargin)
            p = inputParser();
            p.addParameter('Writable', false, @islogical);
            p.parse(varargin{:});
            
            mm = memmapfile(this.pathAP, 'Format', {'int16', [this.nChannels this.nSamplesAP], 'x'}, 'Writable', p.Results.Writable);
        end
        
        function mm = memmapLF_full(this, varargin)
            p = inputParser();
            p.addParameter('Writable', false, @islogical);
            p.parse(varargin{:});
            
            mm = memmapfile(this.pathLF, 'Format', {'int16', [this.nChannels this.nSamplesLF], 'x'}, 'Writable', p.Results.Writable);
        end
        
        function [mm, fsSync] = memmapSync_full(this)
            if this.syncInAPFile
                % still has nChannels
                mm = memmapfile(this.pathSync, 'Format', {'int16', [this.nChannels this.nSamplesAP], 'x'});
                fsSync = this.fsAP;
            elseif this.syncInLFFile
                % still has nChannels
                mm = memmapfile(this.pathSync, 'Format', {'int16', [this.nChannels this.nSamplesLF], 'x'});
                fsSync = this.fsLF;
            else
                % only sync channel
                mm = memmapfile(this.pathSync, 'Format', {'int16', [1 this.nSamplesAP], 'x'});
                fsSync = this.fsAP;
            end
        end
        
        % refer back to the source datasets
        function mmSet = memmap_sourceAP_full(this, varargin)
            assert(~isempty(this.sourceDatasets));
            nSources = numel(this.sourceDatasets);
            mmSet = cell(nSources, 1);
            for iF = 1:nSources
                mmSet{iF} = this.sourceDatasets(iF).memmapAP_full(varargin{:});
            end
        end
        
        % refer back to the source datasets
        function mmSet = memmap_sourceLF_full(this, varargin)
            assert(~isempty(this.sourceDatasets));
            nSources = numel(this.sourceDatasets);
            mmSet = cell(nSources, 1);
            for iF = 1:nSources
                mmSet{iF} = this.sourceDatasets(iF).memmapLF_full(varargin{:});
            end
        end
        
        function [mmSet, fsSync] = memmap_sourceSync_full(this, varargin)
            assert(~isempty(this.sourceDatasets));
            nSources = numel(this.sourceDatasets);
            mmSet = cell(nSources, 1);
            fsSync = nan(nSources, 1);
            for iF = 1:nSources
                imecSrc = this.sourceDatasets(iF);
                if imecSrc.syncInAPFile
                    % still has nChannels
                    mmSet{iF} = memmapfile(imecSrc.pathSync, 'Format', {'int16', [this.nChannels this.nSamplesAP], 'x'});
                    fsSync(iF) = this.fsAP;
                elseif imecSrc.syncInLFFile
                    % still has nChannels
                    mmSet{iF} = memmapfile(this.pathSync, 'Format', {'int16', [this.nChannels this.nSamplesLF], 'x'});
                    fsSync(iF) = this.fsLF;
                else
                    % only sync channel
                    mmSet{iF} = memmapfile(this.pathSync, 'Format', {'int16', [1 this.nSamplesAP], 'x'});
                    fsSync(iF) = this.fsAP;
                end
            end
            
            assert(numel(unique(fsSync)) == 1);
            fsSync = fsSync(1);
        end
    end
    
    methods (Static)
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
    
    methods (Hidden) % Read data at specified times
        function [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet, scaleToUv_by_snippet, group_ids] ...
                = readSnippetsRaw(this, times, window, varargin)
            % for each sample index in times, read the window times + window(1):window(2)
            % of samples around this time from some channels
            
            p = inputParser();
            p.addParameter('band', 'ap', @npxutils.internal.isstringlike);
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
                channel_ids = npxutils.internal.makecol(channel_ids);
                channel_ids_by_snippet = repmat(channel_ids, 1, numel(times));
                
            else
                error('Specify either channel_ids or channel_ids_by_cluster');
            end
            
            channel_inds_by_snippet = this.lookup_channelIds(channel_ids_by_snippet);
            syncChannelIndex = this.syncChannelIndex;
            syncInChannelInds = any(ismember(channel_inds_by_snippet, syncChannelIndex), 'all');
            
            band = string(p.Results.band);
            fromSourceDatasets = p.Results.fromSourceDatasets;
            syncFromSource = p.Results.syncFromSourceDatasets;
            if isempty(syncFromSource)
                syncFromSource = fromSourceDatasets;
            end
            
            
            switch band
                case 'ap'
                    nSamples = this.nSamplesAP;
                    if ~fromSourceDatasets
                        mm = this.memmapAP_full();
                        scaleToUv = this.apScaleToUv;
                    else
                        mmSet = this.memmap_sourceAP_full();
                        concatInfo = this.concatenationInfoAP;
                        scaleToUv_this = this.apScaleToUv;
                        scaleToUv_by_source = cat(1, this.sourceDatasets.apScaleToUv);
                    end
                    
                    % do we need to specially handle the sync channel?
                    if syncFromSource ~= fromSourceDatasets && syncInChannelInds
                        handleSyncSeparately = true;
                        if ~syncFromSource
                            mmSync = this.memmapAP_full();
                        else
                            mmSetSync = this.memmap_sourceAP_full();
                        end
                    else
                        handleSyncSeparately = false;
                    end
                case 'lf'
                    nSamples = this.nSamplesLF;
                    if ~fromSourceDatasets
                        mm = this.memmapLF_full();
                        scaleToUv = this.lfScaleToUv;
                    else
                        mmSet = this.memmap_sourceLF_full();
                        concatInfo = this.concatenationInfoLF;
                        scaleToUv_this = this.lfScaleToUv;
                        scaleToUv_by_source = cat(1, this.sourceDatasets.lfScaleToUv);
                    end
                    
                    % do we need to specially handle the sync channel?
                    if syncFromSource ~= fromSourceDatasets  && syncInChannelInds
                        handleSyncSeparately = true;
                        if ~syncFromSource
                            mmSync = this.memmapLF_full();
                        else
                            mmSetSync = this.memmap_sourceLF_full();
                        end
                    else
                        handleSyncSeparately = true;
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
            
            times = npxutils.internal.makecol(uint64(times));
            nC = size(channel_ids_by_snippet, 1);
            nC_all = this.nChannels;
            nS = numel(times); % actual number of input snippets extracted
            
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
            accum_counter = zeros(nS_out, 1);
            
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
                prog = npxutils.internal.ProgressBar(numel(times), 'Extracting %s snippets', upper(band));
            else
                prog = [];
            end
            
            good_ch_inds = this.goodChannelInds;
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
                    extract_all_ch = reshape(npxutils.ImecDataset.multi_mmap_extract_sample_idx(mmSet, ...
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
                        extract_sync_ch = reshape(npxutils.ImecDataset.multi_mmap_extract_sample_idx(mmSetSync, ...
                            sourceFileInds(:), sourceSampleInds(:), syncChannelIndex), [1 nT nS_this]);
                    end
                end
                
                if p.Results.car
                    ar = median(extract_all_ch(good_ch_inds, :, :), 1);
                else
                    ar = zeros([1 nT nS_this], 'like', extract_all_ch);
                end
                
                for iiS = 1:nS_this
                    this_extract =  extract_all_ch(channel_inds_by_snippet(:, idxS(iiS)), :, iiS);
                    
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
                    if p.Results.center
                        this_extract(mask_scale, :) = this_extract(mask_scale, :) - median(this_extract(mask_scale, :), 2);
                    end
                    % scale before accumulating
                    if ~isempty(average_weight)
                        this_extract(mask_scale, :) = this_extract(mask_scale, :) * average_weight(idxS(iiS));
                    end
                    if mask_request_okay_this(iiS)
                        out(:, :, idxInsert(iiS)) = out(:, :, idxInsert(iiS)) + cast(this_extract, outClass); % which channels for this spike
                        accum_counter(idxInsert(iiS)) = accum_counter(idxInsert(iiS)) + 1;
                    end
                end
                
                if ~isempty(prog), prog.increment(nPerSegment); end
            end
            
            if average_by_cluster_id || average_by_group_id
                % divide the accumulator by the number of units, don't divide by the sum of average_weights, this is treated as normalized already
                out = out ./ shiftdim(accum_counter, -2);
            else
                out(:, ~mask_idx_okay(:)) = 0;
            end
            
            if applyScaling
                % scaling already applied
                scaleToUv_by_snippet(:) = 1;
            end
            data_ch_by_time_by_snippet = out;
            if ~isempty(prog), prog.finish(); end
        end
    end
    
    methods  % Read data at specified times
        function snippet_set = readSnippetSet(this, band, times, window, varargin)
            [data_ch_by_time_by_snippet, cluster_ids, channel_ids_by_snippet, scaleToUv_by_snippet, group_ids] = ...
                this.readSnippetsRaw(times, window, 'band', band, varargin{:});
            snippet_set = npxutils.SnippetSet(this, band);
            snippet_set.data = data_ch_by_time_by_snippet;
            snippet_set.scaleToUv = scaleToUv_by_snippet;
            snippet_set.sample_idx = times;
            snippet_set.channel_ids_by_snippet = channel_ids_by_snippet;
            snippet_set.cluster_ids = cluster_ids;
            snippet_set.group_ids = group_ids;
            snippet_set.window = window;
        end
        
        function snippet_set = readAPSnippetSet(this, times, window, varargin)
            snippet_set = this.readSnippetSet('ap', times, window, varargin{:});
        end
        
        function snippet_set = readLFSnippetSet(this, times, window, varargin)
            snippet_set = this.readSnippetSet('lf', times, window, varargin{:});
        end
        
        function rms = computeRMSByChannel(this, varargin)
            % output will be nMappedChannels x 1 vector of rms
            p = inputParser();
            p.addParameter('sampleMaskFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % sampleMaskFn(data_ch_x_time, sample_idx_time) --> logical_time mask of time samples valid for use, useful if you have artifacts at known times
            p.addParameter('car', false, @islogical);
            p.addParameter('useChunks', 50, @isscalar);
            p.addParameter('chunkSize', 100000, @isscalar);
            p.parse(varargin{:});
            
            sampleMaskFn = p.Results.sampleMaskFn;
            
            % aim for the middle of the file
            chunkSize = min(this.fsAP, p.Results.chunkSize);
            mm = this.memmapAP_by_chunk(chunkSize);
            nChunks = numel(mm.Data);
            useChunks = min(nChunks, p.Results.useChunks);
            skipChunks = floor((nChunks-useChunks)/2);
            
            ch_mask = this.lookup_channelIds(this.mappedChannels); % for common average referencing
            
            sumByChunk = nan(this.nChannels, useChunks);
            %             prog = npxutils.internal.ProgressBar(useChunks, 'Computing RMS per channel');
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
            rms = rms * this.apScaleToUv;
        end
    end
    
    methods (Hidden)
        function fid = openAPFile(this)
            if ~exist(this.pathAP, 'file')
                error('RawDataFile: %s not found', this.pathAP);
            end
            fid = fopen2(this.pathAP, 'r');
        end
        
        function fid = openLFFile(this)
            if ~exist(this.pathLF, 'file')
                error('RawDataFile: %s not found', this.pathAP);
            end
            fid = fopen2(this.pathLF, 'r');
        end
        
        function fid = openSyncFile(this)
            if ~exist(this.pathSync, 'file')
                error('RawDataFile: %s not found', this.pathSync);
            end
            fid = fopen2(this.pathSync, 'r');
        end
    
    end
    
        %% Dependent properties
    methods    
        function tf = get.syncInAPFile(this)
            tf = this.channelMap.syncInAPFile && this.hasAP;
        end
        
        function tf = get.syncInLFFile(this)
            tf = this.channelMap.syncInLFFile && this.hasLF;
        end
        
        function id = get.syncChannelId(this)
            if this.syncInAPFile || this.syncInLFFile
                id = this.channelMap.syncChannelId;
            else
                % if sync is in its own file, assume it's the first and only channel
                id = NaN;
            end
        end
        
        function ind = get.syncChannelIndex(this)
            if this.syncInAPFile || this.syncInLFFile
                ind = this.channelMap.syncChannelIndex;
            else
                % if sync is in its own file, assume it's the first and only channel
                ind = uint32(1);
            end
        end
        
        function pathAP = get.pathAP(this)
            pathAP = fullfile(this.pathRoot, this.fileAP);
        end
        
        function fileAP = get.fileAP(this)
            if isnan(this.fileImecNumber)
                fileAP = [this.fileStem '.imec.' this.fileTypeAP '.bin'];
            else
                fileAP = [this.fileStem, sprintf('.imec%d.', this.fileImecNumber), this.fileTypeAP, '.bin'];
            end
        end
        
        function tf = get.hasAP(this)
            tf = exist(this.pathAP, 'file') == 2;
        end
        
        function fileAPMeta = get.fileAPMeta(this)
            if isnan(this.fileImecNumber)
                fileAPMeta = [char(this.fileStem) '.imec.ap.meta'];
            else
                fileAPMeta = sprintf('%s.imec%d.ap.meta', this.fileStem, this.fileImecNumber);
            end
        end
        
        function pathAPMeta = get.pathAPMeta(this)
            pathAPMeta = fullfile(this.pathRoot, this.fileAPMeta);
        end
        
        function pathLF = get.pathLF(this)
            pathLF = fullfile(this.pathRoot, this.fileLF);
        end
        
        function fileLF = get.fileLF(this)
            if isnan(this.fileImecNumber)
                fileLF = [this.fileStem '.imec.lf.bin'];
            else
                fileLF = [this.fileStem, sprintf('.imec%d.', this.fileImecNumber), 'lf.bin'];
            end
        end
        
        function fileLFMeta = get.fileLFMeta(this)
            if isnan(this.fileImecNumber)
                fileLFMeta = [this.fileStem '.imec.lf.meta'];
            else
                fileLFMeta = [this.fileStem, sprintf('.imec%d.lf.meta', this.fileImecNumber)];
            end
        end
        
        function pathLFMeta = get.pathLFMeta(this)
            pathLFMeta = fullfile(this.pathRoot, this.fileLFMeta);
        end
        
        function tf = get.hasLF(this)
            tf = exist(this.pathLF, 'file') == 2;
        end
        
        function fileSync = get.fileSync(this)
            if this.syncInAPFile
                fileSync = this.fileAP;
            elseif this.syncInLFFile
                fileSync = this.fileLF;
            else
                fileSync = [this.fileStem, '.imec.sync.bin'];
            end
        end
        
        function pathSync = get.pathSync(this)
            pathSync = fullfile(this.pathRoot, this.fileSync);
        end
        
        function tf = get.hasSync(this)
            out = isfile(this.pathSync);
        end
        
        function fs = get.fsSync(this)
            % defer to fsAP or fsLF, or to stored value if neither sources the sync signal
            if ~isempty(this.fsSync)
                fs = this.fsSync;
            elseif this.syncInAPFile
                fs = this.fsAP;
            elseif this.syncInLFFile
                fs = this.fsLF;
            else
                fs = NaN;
            end
        end
        
        function tf = get.hasSourceSync(this)
            tf = this.hasSourceDatasets && all([this.sourceDatasets.hasSync]);
        end
        
        function out = get.fileSyncCached(this)
            out = [char(this.fileStem) '.sync.mat'];
        end
        
        function pathSyncCached = get.pathSyncCached(this)
            pathSyncCached = fullfile(this.pathRoot, this.fileSyncCached);
        end
        
        function scale = get.apScaleToUv(this)
            if isempty(this.apRange)
                scale = NaN;
            else
                scale = (this.apRange(2) - this.apRange(1)) / (2^this.adcBits) / this.apGain * 1e6;
            end
        end
        
        function scale = get.lfScaleToUv(this)
            if isempty(this.lfRange)
                scale = NaN;
            else
                scale = (this.lfRange(2) - this.lfRange(1)) / (2^this.adcBits) / this.lfGain * 1e6;
            end
        end
        
        function file = get.channelMapFile(this)
            if isempty(this.channelMap)
                file = '';
            else
                file = this.channelMap.file;
            end
        end
        
        function list = get.mappedChannels(this)
            if isempty(this.channelMap)
                list = [];
            else
                list = this.channelMap.channelIdsMapped;
            end
        end
        
        function list = get.mappedChannelInds(this)
            list = this.lookup_channelIds(this.mappedChannels);
        end
        
        function list = get.connectedChannels(this)
            if isempty(this.channelMap)
                list = [];
            else
                list = this.channelMap.connectedChannels;
            end
        end
        
        function list = get.connectedChannelInds(this)
            list = this.lookup_channelIds(this.connectedChannels);
        end
        
        function n = get.nChannelsMapped(this)
            if isempty(this.channelMap)
                n = NaN;
            else
                n = this.channelMap.nChannelsMapped;
            end
        end
        
        function n = get.nChannelsConnected(this)
            if isempty(this.channelMap)
                n = NaN;
            else
                n = nnz(this.channelMap.connected);
            end
        end
        
        function ch = get.goodChannels(this)
            ch = setdiff(this.connectedChannels, this.badChannels);
        end
        
        function list = get.goodChannelInds(this)
            list = this.lookup_channelIds(this.goodChannels);
        end
        
        function n = get.nGoodChannels(this)
            n = numel(this.goodChannels);
        end
        
        function idx = get.channelIds(this)
            idx = this.channelMap.channelIds;
        end
        
        function names = get.channelNames(this)
            names = strings(this.nChannels, 1);
            names(this.channelMap.channelIds) = string(sprintfc("ch %d", this.channelMap.channelIds));
            if ~isnan(this.syncChannelIndex)
                names(this.syncChannelIndex) = "sync";
            end
        end
        
        function names = get.channelNamesPadded(this)
            names = strings(this.nChannels, 1);
            names(this.channelMap.channelIds) = string(sprintfc("ch %03d", this.channelMap.channelIds));
            if ~isnan(this.syncChannelIndex)
                names(this.syncChannelIndex) = "sync";
            end
        end
        
        function n = get.nSyncBits(this)
            n = 8*this.bytesPerSample; % should be 16?
        end
        
        function bits = get.syncBitsNamed(this)
            names = this.syncBitNames;
            bits = find(names ~= "");
        end
        
        function names = get.syncBitNames(this)
            if isempty(this.syncBitNames)
                names = strings(this.nSyncBits, 1);
            else
                names = string(this.syncBitNames);
            end
        end
        
        function meta = readAPMeta(this)
            meta = npxutils.io.readINI(this.pathAPMeta);
        end
        
        function meta = generateModifiedAPMeta(this)
            meta = this.readAPMeta;
            
            meta.syncBitNames = this.syncBitNames;
            meta.badChannels = this.badChannels;
        end
        
        function tf = get.hasSourceDatasets(this)
            tf = ~isempty(this.sourceDatasets);
        end
        
        function tf = get.hasSourceAP(this)
            tf = this.hasSourceDatasets && all([this.sourceDatasets.hasAP]);
        end
        
        function tf = get.hasSourceLF(this)
            tf = this.hasSourceDatasets && all([this.sourceDatasets.hasLF]);
        end
        
        function writeModifiedAPMeta(this, varargin)
            p = inputParser();
            p.addParameter('extraMeta', struct(), @isstruct);
            p.parse(varargin{:});
            
            meta = this.generateModifiedAPMeta();
            
            % set extra user provided fields
            extraMeta = p.Results.extraMeta;
            extraMetaFields = fieldnames(extraMeta);
            for iFld = 1:numel(extraMetaFields)
                meta.(extraMetaFields{iFld}) = extraMeta.(extraMetaFields{iFld});
            end
            
            npxutils.io.writeINI([this.pathAPMeta], meta);
        end
        
        function meta = readLFMeta(this)
            meta = npxutils.io.readINI(this.pathLFMeta);
        end
        
        function str = get.creationTimeStr(this)
            str = datestr(this.creationTime);
        end
        
        function file = getAuxiliaryFileWithSuffix(this, suffix)
            suffix = char(suffix);
            file = fullfile(this.pathRoot, [this.fileStem, '.', suffix]);
        end
    
        %% Marking Channels as bad
        
        function [rmsBadChannels, rmsByChannel] = markBadChannelsByRMS(this, varargin)
            p = inputParser();
            p.addParameter('rmsRange', [3 100], @isvector);
            p.addParameter('sampleMaskFn', [], @(x) isempty(x) || isa(x, 'function_handle')); % sampleMaskFn(data_ch_x_time, sample_idx_time) --> logical_time mask of time samples valid for use, useful if you have artifacts at known times
            p.parse(varargin{:});
            
            rmsByChannel = this.computeRMSByChannel('sampleMaskFn', p.Results.sampleMaskFn);
            rmsMin = p.Results.rmsRange(1);
            rmsMax = p.Results.rmsRange(2);
            rmsBadMask = rmsByChannel < rmsMin | rmsByChannel > rmsMax;
            
            badMappedChannels = this.mappedChannels(rmsBadMask);
            badConnectedChannels = badMappedChannels(ismember(badMappedChannels, this.connectedChannels));
            this.markBadChannels(badConnectedChannels);
            
            rmsBadChannels = badMappedChannels;
        end
        
        function markBadChannels(this, list)
            % Mark specified channels as bad.
            %
            % markBadChannels(this, list)
            %
            % List is a logical or numeric index vector identifying the channels
            % to mark as bad. If you specify a non-connected or nonexistent
            % channel, it is silently ignored.
            %
            % This adds to the set of bad channels, so multiple calls will
            % remove additional channels.
            if islogical(list)
                list = find(list);
            end
            % filter for connected channels only
            badConnectedChannels = list(ismember(list, this.connectedChannels));
            this.badChannels = union(this.badChannels, badConnectedChannels);
        end
    
        %% Modify bin data files in place
        
        function modifyAPInPlace(this, varargin)
            this.modifyInPlaceInternal('ap', varargin{:});
        end
        
        function modifyLFInPlace(this, varargin)
            this.modifyInPlaceInternal('lf', varargin{:});
        end
        
        function imecSym = symLinkAPIntoDirectory(this, newFolder, varargin)
            p = inputParser();
            p.addParameter('relative', false, @islogical);
            p.parse(varargin{:});
            newFolder = char(newFolder);
            
            if ~exist(newFolder, 'dir')
                npxutils.internal.mkdirRecursive(newFolder);
            end
            newAPPath = fullfile(newFolder, this.fileAP);
            npxutils.internal.makeSymLink(this.pathAP, newAPPath, p.Results.relative);
            
            newAPMetaPath = fullfile(newFolder, this.fileAPMeta);
            npxutils.internal.makeSymLink(this.pathAPMeta, newAPMetaPath, p.Results.relative);
            
            if ~this.syncInAPFile && exist(this.pathSync, 'file')
                newSyncPath = fullfile(newFolder, this.fileSync);
                npxutils.internal.makeSymLink(this.pathSync, newSyncPath, p.Results.relative);
            end
            
            if exist(this.pathSyncCached, 'file')
                newSyncCachedPath = fullfile(newFolder, this.fileSyncCached);
                npxutils.internal.makeSymLink(this.pathSyncCached, newSyncCachedPath, p.Results.relative);
            end
            
            imecSym = npxutils.ImecDataset(newAPPath, 'channelMap', this.channelMapFile);
        end
        
        function imecOut = saveTransformedDataset(this, outPath, varargin)
            p = inputParser();
            p.addParameter('stem', "", @npxutils.internal.isstringlike);
            
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
            
            p.addParameter('timeShiftsAP', {}, @(x) isempty(x) || isa(x, 'npxutils.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            p.addParameter('timeShiftsLF', {}, @(x) isempty(x) || isa(x, 'npxutils.TimeShiftSpec')); % cell array of time shifts for each file, a time shift is a n x 3 matrix of idxStart, idxStop, newIdxStart. These are used to excise specific time windows from the file
            
            p.addParameter('extraMeta', struct(), @isstruct);
            
            p.addParameter('dryRun', false, @islogical);
            p.parse(varargin{:});
            
            % this uses the same syntax as writeConcatenatedFileMatchGains
            imecOut = npxutils.ImecDataset.writeConcatenatedFileMatchGains({this}, outPath, p.Results);
        end
        
        function writeFolderForPhy(this, savePath, varargin)
            % writes the files needed for Phy template-gui to be able to
            % inspect this file in Phy essentially by generating a
            % synthetic Kilosort output
            
            p = inputParser();
            p.addParameter('spikeTimes', [1 2], @(x) isempty(x) || isvector(x));
            p.addParameter('spikeClusters', [], @(x) isempty(x) || isvector(x));
            p.addParameter('clusterNames', [], @(x) isempty(x) || iscellstr(x) || isstring(x));
            p.parse(varargin{:});
            
            spikeTimes = npxutils.internal.makecol(p.Results.spikeTimes);
            nSpikes = numel(spikeTimes);
            assert(nSpikes >= 2, 'At least 2 spikes are required for Phy');
            
            if isempty(p.Results.spikeClusters)
                spikeClusters = zeros(nSpikes, 1, 'uint32');
            else
                spikeClusters = uint32(npxutils.internal.makecol(p.Results.spikeClusters));
                assert(numel(spikeClusters) == nSpikes);
            end
            
            nTemplateTimepoints = 82;
            nCh = this.nGoodChannels;
            % these can't be 1 because of squeeze() inside Phy's code
            nTemplates = 2;
            nTemplateFeatures = 2;
            nFeaturesPerChannel = 2;
            nPCFeatures = 2;
            
            this.symLinkAPIntoDirectory(savePath);
            
            npxutils.io.writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
            npxutils.io.writeNPY(zeros(nSpikes, 1, 'uint32'), fullfile(savePath, 'spike_templates.npy'));
            npxutils.io.writeNPY(spikeClusters, fullfile(savePath, 'spike_clusters.npy'));
            
            npxutils.io.writeNPY(zeros(nSpikes, 1, 'double'), fullfile(savePath, 'amplitudes.npy'));
            
            templates = zeros(2, nTemplateTimepoints, nCh, 'single');
            npxutils.io.writeNPY(templates, fullfile(savePath, 'templates.npy'));
            
            templatesInds = this.goodChannels';
            npxutils.io.writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
            
            sortedInds = this.goodChannelInds;
            chanMap0ind = int32(this.channelMap.channelIdsMapped(sortedInds) - uint32(1));
            xcoords = this.channelMap.xcoords(sortedInds);
            ycoords = this.channelMap.ycoords(sortedInds);
            npxutils.io.writeNPY(chanMap0ind, fullfile(savePath, 'channel_map.npy'));
            npxutils.io.writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));
            
            templateFeatures = zeros([nTemplates nTemplateFeatures], 'single');
            npxutils.io.writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
            
            templateFeatureInds = zeros(nTemplates, nTemplateFeatures, 'uint32');
            npxutils.io.writeNPY(templateFeatureInds, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
            
            similarTemplates = zeros(nTemplates, nTemplates, 'single');
            npxutils.io.writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
            
            pcFeatures = zeros([nSpikes, nFeaturesPerChannel, nPCFeatures], 'single');
            npxutils.io.writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
            
            pcFeatureInds = zeros([nTemplates, nPCFeatures], 'uint32');
            npxutils.io.writeNPY(pcFeatureInds, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
            
            whiteningMatrix = ones(nCh, nCh, 'double');
            npxutils.io.writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
            whiteningMatrixInv = ones(nCh, nCh, 'double');
            npxutils.io.writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
            
            % write params.py
            fid = fopen2(fullfile(savePath,'params.py'), 'w');
            [~, fname, ext] = fileparts(this.fileAP);
            fprintf(fid,['dat_path = ''',fname ext '''\n']);
            fprintf(fid,'n_channels_dat = %i\n',this.nChannels);
            fprintf(fid,'dtype = ''int16''\n');
            fprintf(fid,'offset = 0\n');
            fprintf(fid,'sample_rate = %i\n', this.fsAP);
            fprintf(fid,'hp_filtered = False');
            fclose2(fid);
            
            % write spike names column
            clusterIds = unique(spikeClusters);
            nClusters = numel(clusterIds);
            if isempty(p.Results.clusterNames)
                clusterNames = strings(nClusters, 1);
            else
                clusterNames = string(npxutils.internal.makecol(p.Results.clusterNames));
                assert(numel(clusterNames) == nClusters, 'clusterNames must be nClusters long');
            end
            fid2 = fopen2(fullfile(savePath, 'cluster_names.tsv'), 'w');
            fprintf(fid2, 'cluster_id\tname\n');
            
            for iC = 1:nClusters
                fprintf(fid2, '%d\t%s\n', clusterIds(iC)-1, clusterNames{iC});
            end
            fclose2(fid2);
            
        end
    end
    
    methods (Hidden)
        function [chInds, chIds] = build_channelSelectors_internal(this, varargin)
            p = inputParser();
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            p.addParameter('mappedChannelsOnly', false, @islogical);
            p.parse(varargin{:});
            
            if p.Results.goodChannelsOnly
                [chInds, chIds] = this.lookup_channelIds(this.goodChannels);
                assert(~isempty(chInds), 'No channels marked good in dataset')
                
            elseif p.Results.connectedChannelsOnly
                [chInds, chIds] = this.lookup_channelIds(this.connectedChannels); % excludes sync channel
                assert(~isempty(chInds), 'No connected channels found in dataset');
                
            elseif p.Results.mappedChannelsOnly
                [chInds, chIds] = this.lookup_channelIds(this.mappedChannels); % excludes sync channel
                assert(~isempty(chInds), 'No mapped channels found in dataset');
            else
                chInds = 1:this.nChannels;
                chIds = this.channelIds(chInds);
            end
        end
        
        function transformExtraArg = modifyInPlaceInternal(this, mode, procFnList, varargin)
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
                    mm = this.memmapAP_full('Writable', ~dryRun);
                case 'lf'
                    mm = this.memmapLF_full('Writable', ~dryRun);
                otherwise
                    error('Unknown mode %s', mode);
            end
            
            % figure out which channels to keep
            [chInds, chIds] = this.build_channelSelectors_internal('goodChannelsOnly', p.Results.goodChannelsOnly, ...
                'connectedChannelsOnly', p.Results.connectedChannelsOnly, ...
                'mappedChannelsOnly', p.Results.mappedChannelsOnly);
            
            dataSize = size(mm.Data.x, 2);
            nChunks = ceil(dataSize / chunkSize);
            prog = npxutils.internal.ProgressBar(nChunks, 'Modifying %s file in place', mode);
            for iCh = 1:nChunks
                [source_idx, keepIdx] = npxutils.ImecDataset.determineChunkIdx(dataSize, iCh, nChunks, chunkSize, chunkExtra);
                
                if dryRun && ~isempty(dryRunSampleInds) && ~any(ismember(source_idx, dryRunSampleInds))
                    % skip this chunk unless some ind in idx should be processed
                    continue;
                end
                
                data = mm.Data.x(chInds, source_idx);
                data_pre = data;
                
                % ch_connected_mask indicates which channels are
                % connected, which are the ones where scaling makes
                % sense. chIdx is all channels being modified by procFnList
                ch_conn_mask = ismember(chIds, this.connectedChannels);
                
                % do additional processing here
                if applyScaling
                    % convert to uV and to single
                    switch mode
                        case 'ap'
                            data = single(data);
                            data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(this.apScaleToUv);
                        case 'lf'
                            data = single(data);
                            data(ch_conn_mask, :) = data(ch_conn_mask, :) * single(this.lfScaleToUv);
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
                        [data, transformExtraArg] = fn(this, data, extraArgs{:});
                    else
                        data = fn(this, data, extraArgs{:});
                    end
                    
                end
                
                if useGpuArray
                    data = gather(data);
                end
                
                if applyScaling
                    data(ch_conn_mask, :) = data(ch_conn_mask, :) ./ this.scaleToUv;
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
                this.writeModifiedAPMeta('extraMeta', p.Results.extraMeta);
                this.clearSyncCached();
            end
        end
    end
    
    methods (Static)
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
            file = npxutils.ImecDataset.findImecFileInDir(fileOrFileStem, 'ap');
            if isempty(file)
                tf = false;
                return;
            end
            
            [pathRoot, fileStem, fileTypeAP] = npxutils.ImecDataset.parseImecFileName(file);
            pathAP = fullfile(pathRoot, [fileStem '.imec.' fileTypeAP '.bin']);
            pathAPMeta = fullfile(pathRoot, [fileStem '.imec.ap.meta']);
            
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
                candidates = npxutils.ImecDataset.findImecFileInDir(fileOrFileStem, 'ap', true, false);
                if numel(candidates) == 0
                    candidates = npxutils.ImecDataset.findImecFileInDir(fileOrFileStem, 'lf', true, false);
                end
            else
                candidates = npxutils.ImecDataset.findImecFileInDir(fileOrFileStem, type, true, false);
            end
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
                [~, ~, type] = npxutils.ImecDataset.parseImecFileName(file);
                switch type
                    case 'ap'
                        assert(ismember(type, {'ap', 'ap_CAR'}), 'Specify ap.bin or ap_CAR.bin file rather than %s file', type);
                    case 'lf'
                        assert(ismember(type, {'lf'}), 'Specify lf.bin file rather than %s file', type);
                end
                
            elseif exist(fileOrFileStem, 'dir')
                % it's a directory, assume only one imec file in directory
                path = fileOrFileStem;
                %                 [~, leaf] = fileparts(path);
                
                switch type
                    case 'ap'
                        apFiles = npxutils.ImecDataset.listAPFilesInDir(path);
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
                        lfFiles = npxutils.ImecDataset.listLFFilesInDir(path);
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
                switch type
                    case 'ap'
                        candidates = npxutils.ImecDataset.listAPFilesInDir(parent);
                        bin_ext = '.imec.ap.bin';
                    case 'lf'
                        candidates = npxutils.ImecDataset.listLFFilesInDir(parent);
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
        
        function [pathRoot, fileStem, type, imecNumber] = parseImecFileName(file)
            if iscell(file)
                [pathRoot, fileStem, type] = cellfun(@npxutils.ImecDataset.parseImecFileName, file, 'UniformOutput', false);
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
        
        function [imecOut, transformAPExtraArg, transformLFExtraArg] ...
                = writeConcatenatedFileMatchGains(imecList, outPath, varargin)
            p = inputParser();
            p.addParameter('stem', "", @npxutils.internal.isstringlike);
            p.addParameter('writeAP', true, @islogical);
            p.addParameter('goodChannelsOnly', false, @islogical);
            p.addParameter('mappedChannelsOnly', false, @islogical);
            p.addParameter('connectedChannelsOnly', false, @islogical);
            % true means ap will get only mapped channels, false will preserve channels as is
            p.addParameter('writeSyncSeparate', false, @islogical);
            p.addParameter('writeLF', false, @islogical);
            p.addParameter('chunkSize', 2^20, @isscalar);
            p.addParameter('chunkEdgeExtraSamplesAP', [0 0], @isvector);
            p.addParameter('chunkEdgeExtraSamplesLF', [0 0], @isvector);
            
            p.addParameter('gpuArray', false, @islogical);
            p.addParameter('applyScaling', false, @islogical); % convert to uV before processing
            
            % list of transformation functions that accept (imec, dataChunk) and return dataChunk someplace
            p.addParameter('transformAP', {}, @(x) iscell(x) || isa(x, 'function_handle'));
            % list of transformation functions that accept (imec, dataChunk) and return dataChunk someplace
            p.addParameter('transformLF', {}, @(x) iscell(x) || isa(x, 'function_handle'));
            p.addParameter('transformAPExtraArg', struct(), @(x) true);
            p.addParameter('transformLFExtraArg', struct(), @(x) true);
            % cell array of time shifts for each file, a time shift is a n x 3 matrix of
            % idxStart, idxStop, newIdxStart. These are used to excise specific
            % time windows from the file.
            p.addParameter('timeShiftsAP', {}, @(x) isempty(x) || isa(x, 'npxutils.TimeShiftSpec'));
            % cell array of time shifts for each file, a time shift is a n x 3 
            % matrix of idxStart, idxStop, newIdxStart. These are used to excise
            % specific time windows from the file
            p.addParameter('timeShiftsLF', {}, @(x) isempty(x) || isa(x, 'npxutils.TimeShiftSpec'));
            
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
            
            [parent, leaf, ext] = npxutils.ImecDataset.filepartsMultiExt(outPath);
            if strlength(ext) > 0 && endsWith(ext, 'bin')
                % specified full file
                outPath = parent;
            else
                outPath = fullfile(parent, [leaf, ext]);
            end
            if ~exist(outPath, 'dir') && ~dryRun
                npxutils.internal.mkdirRecursive(outPath);
            end
            
            if p.Results.stem ~= ""
                leaf = char(p.Results.stem);
            else
                leaf = char(leaf);
            end
            
            % figure out which channels to keep
            [chIndsByFile, ~] = npxutils.ImecDataset.multiFile_build_channelSelectors_internal(imecList, 'goodChannelsOnly', p.Results.goodChannelsOnly, ...
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
            npxutils.ImecDataset.checkClearDestinationStem(fullfile(outPath, leaf));
            
            if ~dryRun
                % dont do this - dangerous
                %npxutils.ImecDataset.clearDestinationStem(fullfile(outPath, leaf));
            end
            
            if p.Results.writeAP || ~isempty(p.Results.transformAP)
                gains = cellfun(@(imec) imec.apGain, imecList);
                [multipliers, gain] = npxutils.ImecDataset.determineCommonGain(gains);
                
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
                    npxutils.io.writeINI(metaOutFile, meta);
                end
                
                fprintf('Writing AP bin file %s\n', (outFile));
                transformAPExtraArg = writeCatFile(outFile, chIndsByFile, 'ap', multipliers, chunkSize, chunkEdgeExtraSamplesAP, ...
                    p.Results.transformAP, timeShiftsAP, dryRun, transformAPExtraArg);
            end
            
            if p.Results.writeLF || ~isempty(p.Results.transformLF)
                gains = cellfun(@(imec) imec.lfGain, imecList);
                [multipliers, gain] = npxutils.ImecDataset.determineCommonGain(gains);
                
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
                    npxutils.io.writeINI(metaOutFile, meta);
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
            imecOut = npxutils.ImecDataset(outFile, 'channelMap', imecList{1}.channelMapFile, 'sourceDatasets', cat(1, imecList{:}));
            
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
                    prog = npxutils.internal.ProgressBar(nChunks, 'Copying %s file %d / %d: %s', mode, iF, nFiles, imecList{iF}.fileStem);
                    
                    testChunkSize = true;
                    if testChunkSize
                        nOut = nan(nChunks, 1);
                        for iCh = 1:nChunks
                            [~, keepIdx] = npxutils.ImecDataset.determineChunkIdx(outSize, iCh, nChunks, chunkSize, chunkExtra);
                            nOut(iCh) = numel(keepIdx);
                        end
                        assert(sum(nOut) == outSize, 'Error with chunk selection');
                    end
                    
                    for iCh = 1:nChunks
                        [idx, keepIdx] = npxutils.ImecDataset.determineChunkIdx(...
                            outSize, iCh, nChunks, chunkSize, chunkExtra);
                        
                        % some elements of source_idx may be 0 (meaning they are
                        % not filled with source data)
                        source_idx = sourceIdxList(idx);
                        mask_source_idx = source_idx > 0;
                        
                        data = zeros(numel(chInds), numel(source_idx), 'int16');
                        data(:, mask_source_idx) = mm.Data.x(chInds, source_idx(mask_source_idx));
                        
                        % ch_connected_mask indicates which channels are
                        % connected, which are the ones where scaling makes
                        % sense. chIdx is all channels being written to output
                        % file.
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
                                        data(ch_conn_mask, :) = data(ch_conn_mask, :) ...
                                            * single(imecList{iF}.apScaleToUv);
                                    case 'lf'
                                        data = single(data);
                                        data(ch_conn_mask, :) = data(ch_conn_mask, :) ...
                                            * single(imecList{iF}.lfScaleToUv);
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
