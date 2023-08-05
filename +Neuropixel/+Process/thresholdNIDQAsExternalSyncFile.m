function thresholdNIDQAsExternalSyncFile(imec, niChannels, args)
    % This function helps convert from a NIDQ-based sync / metadata setup to an internal sync bits based setup
    %
    % First, it performs a simple thresholding operation on the analog signals in 
    % the NIDQ signal to generate logical sync bits.
    % more complex transformations can be performed by providing a thresholdFn directly
    % including creating fewer or additional sync bit signals.
    % 
    % thresholdFn(data, fs=sampling rate, channelNames=string vector) --> (channels x time logical)
    % the output number of channels must match the length of assignSyncBits
    %
    % Second, it can optionally resample this sync file to the higher AP sampling rate. When doing this
    % it can optionally use a common syncronization signal present in both the internal AP sync 
    % and a NIDQ channel to do a better job of synchronizing. This process is akin to what catGT does
    % as explained at https://billkarsh.github.io/SpikeGLX/help/syncEdges/Sync_edges/

    arguments
        imec Neuropixel.ImecDataset
        niChannels (:, 1) = (1:imec.nChannelsNI)';
        args.overwrite (1, 1) logical = false;
        args.errorIfExists (1, 1) logical = true;
        args.assignSyncBits (:, 1) uint16 = niChannels;
        args.thresholdSimple (:, 1) double = 2.5;
        args.thresholdFn (:, 1) function_handle; % provide for more complex thresholding

        args.resampleToAP (1, 1) logical = false;
        args.commonSyncInternalBit (1, 1) uint16 = 0; % 1-index bit in AP sync channel
        args.commonSyncNIDQChannel (1, 1) uint16 = 0; % 1-index channel in NIDQ file
        args.commonSyncNIDQThresh (1, 1) double = 2.5 % threshold for the NIDQ channel to find edges

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

    assert(args.commonSyncInternalBit >= 0 && args.commonSyncInternalBit <= imec.nSyncBits, "commonSyncNIDQChannel out of range");
    assert(args.commonSyncNIDQChannel >= 0 && args.commonSyncNIDQChannel <= nCh, "commonSyncNIDQChannel out of range");

    thresholdSimple = args.thresholdSimple;
    if isscalar(thresholdSimple)
        thresholdSimple = repmat(thresholdSimple, nCh, 1);
    end
    if ~isempty(args.thresholdFn)
        assert(numel(thresholdSimple) == nCh);
    end

    assert(all(args.assignSyncBits >= 1 & args.assignSyncBits <= imec.nSyncBits));
    assert(numel(unique(args.assignSyncBits))==numel(args.assignSyncBits));

    mmNI = imec.memmapNI_full();

    if args.resampleToAP
        if args.commonSyncInternalBit > 0 && args.commonSyncNIDQChannel > 0
            mmAP = imec.memmapAP_full();
            niSyncThreshRaw = args.commonSyncNIDQThresh / imec.niScaleToUv * 1e6; % threshold in native scaling
            edgesAP = evalin('base', 'edgesAP');
            edgesNI = evalin('base', 'edgesNI');
            % [edgesAP, edgesNI] = extractAlignedEdges(mmAP, imec.syncChannelIndex, args.commonSyncInternalBit, mmNI, args.commonSyncNIDQChannel,niSyncThreshRaw, args.chunkSize);
        else
            warning('NIDQ is being resampled to AP, without using common sync signal. Provide commonSyncInternalBit and commonSyncNIDQChannel to specify common sync signal');
            [edgesAP, edgesNI] = deal([]);
        end

        % we'll run the loop over chunks in AP time
        nSamples = imec.nSamplesAP;
    else
        % we'll run the loop over chunks in NI time
        nSamples = imec.nSamplesNI;
    end
    nChunks = ceil(nSamples / args.chunkSize);

    % open new .bin file and write .meta
    if ~args.dryRun
        % save meta file
        if args.resampleToAP
            meta.imSampRate = imec.fsAP;
        else
            meta.imSampRate = imec.fsNI;
        end
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

        if args.resampleToAP
            idxNI = closestIdxNIforAP(idx, imec.fsAP, imec.fsNI, edgesAP, edgesNI);
        else
            idxNI = idx;
        end
        niData = single(mmNI.Data.x(niChannels, idxNI)) * imec.niScaleToUv;
        
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
   
function [edgesAP, edgesNI] = extractAlignedEdges(mmAP, apSyncChannel, apSyncBit, mmNI, nidqSyncChannel, nidqSyncThresh, chunkSize)

    % detect edges in AP by chunk
    nSamples = mmAP.Format{2}(2);
    nChunks = ceil(nSamples / chunkSize);
    [upAP, dnAP] = deal(cell(nChunks, 1));

    prog = ProgressBar(nChunks, "Extracting synchronization edges in AP ch %d bit %d", apSyncChannel, apSyncBit);
    for iCh = 1:nChunks
        idx = Neuropixel.ImecDataset.determineChunkIdx(nSamples, iCh, nChunks, chunkSize);
        this_sync = bitget(mmAP.Data.x(apSyncChannel, idx), apSyncBit) > 0;
        upAP{iCh} = uint64(find(diff(this_sync) == 1)) + idx(1);
        dnAP{iCh} = uint64(find(diff(this_sync) == -1)) + idx(1);
        prog.update(iCh);
    end
    prog.finish();

    % detect edges in NI by chunk
    nSamples = mmNI.Format{2}(2);
    nChunks = ceil(nSamples / chunkSize);
    [upNI, dnNI] = deal(cell(nChunks, 1));
    prog = ProgressBar(nChunks, "Extracting synchronization edges in NIDQ ch %d", nidqSyncChannel);
    for iCh = 1:nChunks
        idx = Neuropixel.ImecDataset.determineChunkIdx(nSamples, iCh, nChunks, chunkSize);
        this_sync = mmNI.Data.x(nidqSyncChannel, idx) >= nidqSyncThresh;
        upNI{iCh} = uint64(find(diff(this_sync) == 1)) + idx(1);
        dnNI{iCh} = uint64(find(diff(this_sync) == -1)) + idx(1);
        prog.update(iCh);
    end
    prog.finish();

    upAP = cat(2, upAP{:})';
    dnAP = cat(2, dnAP{:})';
    upNI = cat(2, upNI{:})';
    dnNI = cat(2, dnNI{:})';

    % check ordering of up/down is consistent
    if min(upAP) < min(dnAP)
        % up first 
        assert(min(upNI) < min(dnNI));
    else
        % down first
        assert(min(upNI) > min(dnNI));
    end
    
    edgesAP = sort(cat(1, upAP, dnAP));
    edgesNI = sort(cat(1, upNI, dnNI));

    % truncate so that lengths match
    if numel(edgesAP) > numel(edgesNI)
        edgesAP = edgesAP(1:numel(edgesNI));
    else
        edgesNI = edgesNI(1:numel(edgesAP));
    end
end

function idxNI = closestIdxNIforAP(idxAP, fsAP, fsNI, edgesAP, edgesNI)
    if isempty(edgesAP) || isempty(edgesNI)
        % simple scaling
        idxNI = floor(double(idxAP-1) * double(fsNI) / double(fsAP)) + 1;
    else
        % scaling relative to last common edge
        % here, edges and idxAP are sorted, so first we find the previous edge
        % for idxAP(1) and store its index in ie.
        nEdges = numel(edgesAP);
        ie = find(edgesAP > idxAP(1), 1, 'first');
        if isempty(ie)
            searchEdge = nEdges;
        else
            searchEdge = ie-1;
        end

        % then we loop over idxAP and store the previous edge for each entry
        % incrementing searchEdge as we go
        [refAP, refNI] = deal(zeros(size(idxAP), like=idxAP));
        for ii = 1:numel(idxAP)
            % increment searchEdge until the next edge occurs after idxAP(ii) or we run out of edges
            while searchEdge < nEdges && edgesAP(searchEdge + 1) <= idxAP(ii)
                searchEdge = searchEdge + 1;
            end

            if searchEdge == 0
                % this index occurs before the first edge, so reference off the start of the file
                refAP(ii) = 1;
                refNI(ii) = 1;
            elseif searchEdge >= nEdges
                % reached last edge, fill in the rest
                refAP(ii:end) = edgesAP(nEdges);
                refNI(ii:end) = edgesNI(nEdges);
                break;
            else
                % use this edge as this idx entry's reference.
                refAP(ii) = edgesAP(searchEdge);
                refNI(ii) = edgesNI(searchEdge);
            end
        end
        
        idxNI = cast(floor(double(idxAP - refAP) * double(fsNI) / double(fsAP)), like=idxAP) + refNI;
    end
end

