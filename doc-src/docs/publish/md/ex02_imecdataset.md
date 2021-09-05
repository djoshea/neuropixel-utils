# About ImecDataset


The `npxutils.ImecDataset` class wraps one individual recording session acquired with SpikeGLX. Currently, four files with extensions `.imec.ap.bin`, `.imec.ap.meta`, `.imec.lf.bin`, and `.imec.lf.meta` comprise one ImecDataset.


# Data setup

```matlab
homeDir = string(java.lang.System.getProperty('user.home'));
exDataDir = homeDir + 'work/npxutils/example-data';
```

# Constructing an ImecDataset


You construct a npxutils.ImecDataset by pointing it at the path to your dataset. How you specify the path is flexible. You can point directly at one of the files:



```matlab
imec = npxutils.ImecDataset(exDataDir + '/raw_datasets/neuropixel_01_g0_t0.imec.ap.bin');
```


```text
Error using npxutils.ImecDataset (line 173)
No AP or LF Imec file found at or in /Users/jankework/npxutils/example-data/raw_datasets/neuropixel_01_g0_t0.imec.ap.bin
```



Or specify the common prefix to the `.imec.`* files comprising the dataset:



```matlab
imecB = npxutils.ImecDataset(exDataDir + '/raw_datasets/neuropixel_01_g0_t0');
```



The common prefix can be shorter as long as there is no ambiguity:



```matlab
imecC = npxutils.ImecDataset(exDataDir + '/raw_datasets/neuropixel_01');
```



You can also point at the parent directory as long as only one `.ap.bin` file is contained within:



```matlab
imecD = npxutils.ImecDataset(exDataDir + '/data/raw_datasets/');
```



These are all equivalent, in that the resulting `imec` instance will wrap both AP and LF bin and meta files. (Though it’s okay if the LF files are missing).


## Specifying a Channel Map


When constructing the ImecDataset, you can specify a channel map. If you don’t specify one, the default will be returned by `npxutils.util.getDefaultChannelMapFile()`, which in turn will look for the file pointed to by the environment variables `'NEUROPIXEL_MAP_FILE'` or `'NPIX_MAP_FILE'`. You can set these using Matlab’s `setenv` function.




These .mat files are expected to be in the same format as found on the [neuropixels repo](https://github.com/cortex-lab/neuropixels). For the phase 3A probe with 384 channels, the file `neuropixPhase3A_kilosortChanMap.mat` contains:



```matlab
ld = load('neuropixPhase3A_kilosortChanMap.mat')
```



You can construct a ChannelMap directly by pointing at the `.mat` file, although every function within `neuropixel-utils` will also accept the filename and construct the `ChannelMap` for you:



```matlab
probeName = 'neuropixPhase3B1';
probeMatFilePath = npxutils.globals.distroot + '/map_files/' + probeName + '_kilosortChanMap.mat';
chanMapObject = npxutils.ChannelMap(probeMatFilePath)
```



Then you can use this map for an ImecDataset using either of the following:



```matlab
imecE = npxutils.ImecDataset('channelMap', chanMapObject);
imecF = npxutils.ImecDataset('channelMap', probeMatFilePath);
```

# Exploring metadata


When the ImecDataset is created, the metadata are loaded from the `.ap.meta` file. You can request the full metadata using:



```matlab
meta = imec.readAPMeta()
```



The most commonly accessed metadata is stored in properties of the ImecDataset:



```matlab
disp(imec)
```

## Setting sync bit names


For convenience, if you use the sync bits for individual TTL signals during your recordings, you can set their names for the recording so that subsequent processing can refer to the bits by name rather than index. You may also find it useful to store additional data in some sync bits during subsequent processing, such as marking regions of time where an artifact was detected, such that specifying the bit by name is useful. Bits are ordered like `bitget`, i.e. 1 is the least significant bit.



```matlab
imec.setSyncBitNames(1, "trialStart");
```



The full set of sync bits is found in `imec.syncBitNames`, or you can lookup the bit corresponding to a given name using:



```matlab
idx = imec.lookupSyncBitByName("trialStart");
```

## Marking bad channels


You can manually mark specific channels as bad using the `markBadChannels` function and passing a list of IDs to it:



```matlab
channelIds = [17 23 52];
imec.markBadChannels(channelIds);
```

### WARNING: Use channel IDs, not indices


Note that like all channel lists, `channelIds` is specified using the actual unique ids of each channel (as specified in the ChannelMap), which is not necessarily their index into the channel map if the channels are not contiguously numbered.




You can use `lookup_channelIds` to find the channel indices for a given set of channel ids if needed:



```matlab
[channelInds, channelIds] = imec.lookup_channelIds(channelIds)
```



One common task is marking channels as bad if their RMS voltage lies outside a reasonable range:



```matlab
goodRmsRange = [3 100]; % [low high] range of RMS in uV
imec.markBadChannelsByRMS('rmsRange', goodRmsRange);
```



The remaining “good” channel ids can always be accessed using `imec.goodChannels` which will be the set of connected channels excluding the bad channels.


## Writing modified metadata back to disk


After making any modifications to the metadata, such as setting the sync bit names or marking bad channels, you can write it back to disk (in the `.imec.ap.meta` file) such that it will be reloaded automatically the next time you create an ImecDataset instance for that file.



```matlab
imec.writeModifiedAPMeta();
```



You can also append whatever fields you want to the meta file using the `extraMeta` parameter:



```matlab
extraMeta.cleaned = true;
extraMeta.cleaningAlgorithm = 'v1';
imec.writeModifiedAPMeta('extraMeta', extraMeta);
```

# Accessing data
## Raw memory maps


The most convenient way to access the data is to request a memory map:



```matlab
mmap = imec.memmapAP_full();
value = mmap.Data.x(ch_index, sample_index) % access a specific sample
```



Equivalent functionality is available for LF files using `imec.memmapLF_full()`. If you wish to modify the underlying data file directly, you can also request a Writeable version of the memory map:



```matlab
mmap = imec.memmapAP_full('Writiable', true);
mmap.Data.x(ch_index, sample_index) = new_value; % overwrite a specific sample
```

## Reading specific time window


You can also request a specific sample or time window directly:



```matlab
idxWindow = [idxFirst idxLast];
[partialData, sampleIdx] = imec.readAP_idx(idxWindow);

timeWindow = [secStart, secStop]; % in seconds
[partialDataB, sampleIdxB] = imec.readAP_timeWindow(timeWindow);
```

## Plotting specific time windows


You can also quickly generate a stacked traces plot of a specific time window, optionally selecting which channels to plot. Take care to select a reasonable time window to avoid overwhelming your system. All channels are individually centered and then collectively normalized by the maximum value before plotting. You can change the global scaling factor by specifying `gain` > 1.




`imec``.``inspectAP_idxWindow``(``idxWindow``,`` ``...``)`




`imec``.``inspectAP_timeWindow``(``timeWindow``,`` ``...``)`




There are additional optional parameters you can specify:



   -  `channels`: channel indices to plot, defaults to `imec.mappedChannels` 
   -  `syncBits` : which sync bits to plot individually, defaults to `imec.syncBitsNamed` 
   -  `gain`: global scaling factor, values larger then 1 will magnify the individual channels, defaults to 0.95 
   -  `car`: perform common average referencing before plotting, defaults to false 
   -  `downsample` : take every nth sample to speed up plotting, defaults to 1 



Good channels are plotted in black, non-connected channels in blue, and bad channels in red. Sync bits are also shown in red and are not affected by the normalization gain.



```matlab
% TODO: Some code that actually generates a plot here!
```

## Reading sync channel


You can access the sync channel all at once using:



```matlab
sync = imec.readSync();
trialStartBit = imec.readSyncBit("trialStart");
```



Or a specific sample or time window using:



```matlab
[partialSync, sampleIdx] = imec.readSync_idx(idxWindow);
[partialSyncB, sampleIdxB] = imec.readSync_timeWindow(timeWindow);
```



You can also access the logical values of specific bits either by bit index or name using:



```matlab
partial = imec.readSyncBits_idx(bits_or_names, idxWindow); %  nTime x nBits
trialStart_partial = imec.readSyncBits_idx("trialStart", idxWindow); % nTime x 1
```

# Building a preprocessing pipeline


If the raw `.imec.ap.bin` file must be processed in some way before running Kilosort, e.g. to remove artifacts, you can implement this efficiently by writing a transformation function that will act on chunks of the data. One example is found in `Neuropixel.DataProcessFn.commonAverageReference`:



```matlab
function [data, extra] = commonAverageReference(imec, data, chIds, sampleIdx) %#ok<INUSD>
    % chIds are the channel ids and imec.goodChannels is a list of channel ids
    % so this will lookup the channel indices of those channels both in chIdx that
    % are also marked as good. This prevents us from computing the reference
    % from the sync, reference, or bad channels.
    chanMask = imec.lookup_channelIds(intersect(chIds, imec.goodChannels));

    % subtract median of each channel over time
    data(chanMask, :) = bsxfun(@minus, data(chanMask, :), median(data(chanMask, :), 2));

    % subtract median across good channels
    data(chanMask, :) = bsxfun(@minus, data(chanMask, :), median(data(chanMask, :), 1));
end
```



Essentially, your transform function can perform any modifications to the data matrix and return the resulting data matrix. Here, `imec` will be the ImecDataset being transformed, `data` will be the `nChannels x nTime` chunk of data being processed. `chIds` will be the channel ids (not necessarily their indices into data but the ids assigned by the channel map), and will typically be the full list of channel ids present in the data file. `sampleIdx` is the sample indices in the current chunk. `extra` is an optional output argument that allows you to store any additional metadata. After the transformation function has been run on every chunk of the dataset, these individual `extra` outputs will be accumulated in a cell array by chunk.


## Modifying a dataset during copy to new location


Once you’ve written your transform function (or functions), you can run them on the dataset using:



```matlab
imecOut = imec.saveTransformedDataset(outPath, 'transformAP', ...\
  {cell of function handles}, 'transformLF', {cell of function handles});
```



Here, `outPath` should include the folder where the new datasets should be written. By default, the file stem (preceeding `.ap.bin`) will match the leaf directory in `outPath`, but this can be specified manually by passing a `stem` parameter:



```matlab
imecOut = imec.saveTransformedDataset('/path/to/datasets', 'stem', 'modifiedDataset', ...)
% creates /path/to/datasets/modifiedDataset.ap.bin, .ap.meta, etc.
```



You can provide one or more function handles (e.g. `@Neuropixel.DataProcessFn.commonAverageReference`) that will be applied sequentially. Other optional parameters include:



   -  `dryRun`: if true, no actual data will be modified on disk, facilitating testing or step by step debugging of the transform functions before writing data. (default false) 
   -  `gpuArray`: if true, the data chunks will be copied to the GPU and the transformation functions will receive and return GPU arrays 
   -  `applyScaling`: if true, the data will be converted to floating point values with the correct analog scale. if false (default) the data will remain in the original, unscaled `int16` quantization. 
   -  `writeAP`: if true, the AP file will be transformed and copied 
   -  `writeLF`: if true, the LF file will be transformed and copied (default false but will be set true automatically if any transformLF is non-empty) 
   -  `goodChannelsOnly`: send only the channels marked good to the transform function 
   -  `connectedChannelsOnly`: send only the connected channels to the transform function 
   -  `mappedChannelsOnly`: send only the mapped channels to the transform function 
   -  `chunkSize`: specify the number of time samples sent to transform functions at once 
   -  `extraMeta`: a struct with extra meta fields to include or overwrite with the output file 
   -  `timeShifts`: a `Neuropixel.TimeShiftSpec` instance used to excise time windows, see [excising time windows](https://djoshea.github.io/neuropixel-utils/imec_dataset/#excising-time-windows) 

### WARNING: Save transformed data to a new folder!


Ensure that `outPath` refers to separate directory so that you make a copy of the dataset rather than writing over the same location. An error will be thrown if any existing files would be overwritten by this call.


## Modifying a dataset in place


Rather than generate a copy, you can also modify a file in place if you’re short on disk space, but be careful, as there’s no undo if anything goes wrong. The same additional parameters are accepted, and you may wish to test your code first by passing `'dryRun', true`.



```matlab
imec.modifyAPInPlace(outPath, {cell of function handles}, ...);
imec.modifyLFInPlace(outPath, {cell of function handles}, ...);
```

## Concatenating multiple files together


If you have multiple separate recording files that you wish to process and sort together, you can concatenate them together during the copy. The code will also scale the datasets up to a common gain factor if the file gains differ from each other. This will not of course increase the resolution of low-gain files, but it will ensure that the signal amplitudes match across files with different gains.



```matlab
imecList = {imec1, imec2, ...};
imecOut = npxutils.ImecDataset.writeConcatenatedFileMatchGains(imecList, outPath, ...
        'transformAP', {cell of function handles}, ...
        'transformLF', {cell of function handles}, ...);
```

# Making copies and symbolic links


You can generate a copy of a dataset using:



```matlab
[imecOut] = imec.saveTransformedDataset(outPath, 'writeAP', true, 'writeLF', true);
```



Alternatively, you can create symbolic links to the AP bin and meta files in a new location using:



```matlab
imecLinked = imec.symLinkAPIntoDirectory(outPath);
```



This is useful for running Kilosort with varying parameters, since each run would ideally live in its own directory but there’s no need for a real copy of the raw data.


## Excising time windows


Occasionally it can be beneficial to remove certain time windows from a file, or to omit them while plotting data. This may be accomplished using `npxutils.TimeShiftSpec` instances to indicate which windows of time to keep and how to shift them so as to remove gaps. A `TimeShiftSpec` specifies a list of sample intervals bounded by a start and stop index in properties `idxStart` and `idxStop`. The start index in `idxStart` will be shifted to lie at sample index `idxShiftStart`. You can calculate these shifts directly, but it is typically easier to specify only the sample intervals you wish to keep and then construct the `TimeShiftSpec` using:



```matlab
timeShifts = npxutils.TimeShiftSpec.buildToExciseGaps(idxStartList, idxStopList);
```



If you have known trial boundaries in your file (see Segmenting a Kilosort dataset into trials for more information), you can also excise the regions of time far from trial boundaries using the `TrialSegmentationInfo` instance. I’ve found this to be useful to exclude time windows where the subject was asleep from further analysis.



```matlab
timeShifts = tsi.computeShiftsExciseRegionsOutsideTrials('maxPauseSec', 20);
```



You can then pass along this `npxutils.TimeShiftSpec` to any of the data transform functions. Depending on whether the `timeShifts` object was created in indices of `AP` band sample rate or `LF` band sample rate, you should pass it along as `timeShiftsAP` or `timeShiftsLF`. The conversion to the other sampling rate will be handled automatically so that the excision is performed on both datasets appropriately.



```matlab
imecOut = imec.saveTransformedDataset(outPath, 'timeShiftsAP', timeShifts, ...);
```



A cell array of `npxutils.TimeShiftSpec` instances can be provided when concatenating multiple files:



```matlab
 imecOut = npxutils.ImecDataset.writeConcatenatedFileMatchGains(outPath, imecList, ...
        'timeShiftsAP', {timeShift1, timeShift2, ... }, ...);
```

## Referring to Source Datasets


If helpful, when loading a derived `ImecDataset`, you can specify the `sourceDatasets` parameter to provide an array of `ImecDataset` instances corresponding to the original, pre-processed source datasets. Here, this would be the set of raw datasets provided as the `imecList` argument above. For certain methods, you can then pass parameter `fromSourceDatasets`, true and the corresponding window of time from the source datasets will be plotted instead of the processed data. This will automatically handle any time shifts and excisions performed; consequently it is helpful for debugging processing pipelines to see the before and after data.



```matlab
imecProcessed = npxutils.ImecDataset.writeConcatenatedFileMatchGains(outPath, imecList);
% ... or:
imecProcessedB = npxutils.ImecDataset(outPath, 'sourceDatasets', {imecRaw1, imecRaw2});

imecProcessed.inspectAP_timeWindow([1 2], 'fromSourceDatasets', true);
```

