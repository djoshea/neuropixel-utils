function thresholdNIDQAsExternalSyncFile(imec, niChannels, args)
    % performs a simple thresholding operation on the analog signals in the NIDQ signal to generate logical sync bits
    % more complex transformations can be performed by providing a thresholdFn directly
    % including creating fewer or additional sync bit signals
    % 
    % thresholdFn(data, fs=sampling rate, channelNames=string vector) --> (channels x time logical)
    % the output number of channels must match the length of assignSyncBits
    arguments
        imec Neuropixel.ImecDataset
        niChannels (:, 1) = (1:imec.nChannelsNI)';
        args.overwrite (1, 1) logical = false;
        args.errorIfExists (1, 1) logical = true;
        args.assignSyncBits (:, 1) uint16 = niChannels;
        args.thresholdSimple (:, 1) double = 2.5;
        args.thresholdFn (:, 1) function_handle; 
        args.chunkPadding (1, 2) uint32 = [0 0];
        args.chunkSize (1,1) uint64 = 2^20;
        args.dryRun (1,1) logical = true;
        args.plotChunks (:, 1) = [];
    end

    if isempty(imec.fileImecNumber)
        imec_str = "imec";
    else
        imec_str = sprintf("imec%d", imec.fileImecNumber);
    end
    sync_file = fullfile(imec.pathRoot, string(imec.fileStem) + "." + imec_str + ".sync.bin");
    if exist(sync_file, 'file') && ~args.overwrite
        if args.errorIfExists
            error("Sync file at %s already exists", sync_file);
        else
            return;
        end
    end
    sync_meta_file = fullfile(imec.pathRoot, string(imec.fileStem) + "." + imec_str + ".sync.meta");

    if isstring(niChannels)
        niChannels = imec.lookupNIChannelByName(niChannels, false);
    end
    nCh = numel(niChannels);
    niChannelNames = imec.niChannelNames(niChannels);

    thresholdSimple = args.thresholdSimple;
    if isscalar(thresholdSimple)
        thresholdSimple = repmat(thresholdSimple, nCh, 1);
    end
    if ~isempty(args.thresholdFn)
        assert(numel(thresholdSimple) == nCh);
    end
    
    nSamples = imec.nSamplesNI;
    nChunks = ceil(nSamples / args.chunkSize);
    
    mm = imec.memmapNI_full();

    assert(all(args.assignSyncBits >= 1 & args.assignSyncBits <= imec.nSyncBits));
    assert(numel(unique(args.assignSyncBits))==numel(args.assignSyncBits));

    % generate new ap.bin file
    if ~args.dryRun
        % save meta file
        meta.imSampRate = imec.fsNI;
        meta.nSavedChans = 1;
        meta.fileSizeBytes = imec.bytesPerSample * nSamples;
        Neuropixel.writeINI(sync_meta_file, meta);

        fidOut = fopen(sync_file, 'w');
        if fidOut == -1
            error('Error opening sync file file %s', sync_file);
        end
    end
    
    prog = ProgressBar(nChunks, "Thresholding %d NIDQ channels to external sync file", nCh);
    for iCh = 1:nChunks
        [idx, keepIdx] = Neuropixel.ImecDataset.determineChunkIdx(nSamples, iCh, nChunks, args.chunkSize, args.chunkPadding);

        niData = single(mm.Data.x(niChannels, idx)) * imec.niScaleToUv;
        
        if isempty(args.thresholdFn)
            niThresh = niData >= thresholdSimple;
        else
            niThresh = args.thresholdFn(niData, fs=imec.fsNI, channelNames=niChannelNames);
        end

        syncData = Neuropixel.Utils.bitpack(niThresh, args.assignSyncBits);

        if ismember(iCh, args.plotChunks)
            figure();
            ptstack(2, 1, cat(3, niData/1e6, niThresh));
            title(sprintf("chunk %d [ %d : %d ]", iCh, idx(1), idx(end)));
        end

        if ~args.dryRun
            syncData = syncData(:, keepIdx); % optional in case we do more sophisticated thresholding requiring padding 
            fwrite(fidOut, syncData, 'int16');
        end
        prog.increment();
    end
    prog.finish();
    
    if ~args.dryRun
        fclose(fidOut);

    end
end



