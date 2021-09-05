# Example 1: Walkthrough
# Data setup

```matlab
homeDir = string(java.lang.System.getProperty('user.home'));
exDataDir = homeDir + 'work/npxutils/example-data';
ksDataDir = homeDir + 'work/npxutils/kilosort';
```

# Create ImecDataset


Create an ImecDataset from a specified raw data file and probe name:



```matlab
channelMapFile = 'neuropixPhase3A_kilosortChanMap.mat';
apBinFile = exDataDir + '/Vinnie/raw_datasets/2018-08-17/Vinnie_20180817_All.imec.ap.bin';
imec = Neuropixel.ImecDataset(apBinFile, 'channelMap', channelMapFile);
disp(imec)
```



Preprocess it a bit:



```matlab
% Mark individual channels as bad based on RMS voltage
rmsBadChannels = imec.markBadChannelsByRMS('rmsRange', [3 100]);

% Specify names for the individual bits in the sync channel
imec.setSyncBitNames([1 2 3], {'trialInfo', 'trialStart', 'stim'});

% Save the bad channels and Sync bit names to the .imec.ap.meta file so they are loaded next time
imec.writeModifiedAPMeta();
```



Common average referencing:



```matlab
% Perform common average referencing on the file and save the results
% to a separate "cleaned" directory
cleanedDir = exDataDir + '/cleaned_datasets';
cleanedApBinFile = cleanedDir + 'Vinnie_20180817_All_cleaned.imec.ap.bin';
extraMeta = struct();
extraMeta.commonAverageReferenced = true;
fnList = {@npxutils.dataprocess.commonAverageReference};
imec = imec.saveTransformedDataset(cleanedApBinFile, 'transformAP', fnList, 'extraMeta', extraMeta);
```



Get it ready for sending to Kilosort:



```matlab
% Symlink the cleaned dataset into a separate directory for Kilosort2
ksPath = ksDataDir + '/Vinnie_20180817_All_cleaned.imec.ap.bin';
imec = imec.symLinkAPIntoDirectory(ksPath);
```

# Visual inspection

```matlab
% Inspect the raw IMEC traces
imec.inspectAP_timeWindow([200 201]); % 200-201 seconds into the recording
```

# Kilosort time!

```matlab
% Run Kilosort2
npxutils.runKilosort2(imec);

% Load the Kilosort2 results
ks = npxutils.KilosortDataset();
ks.load()

% Define how the data are segmented into trials
% (You must implement this function yourself)
% tsi = computeTrialSegmentation(imec)

% TODO: Figure out how to implement a trial segmentation function for these
% examples so the Live Script can run through.

% Segment the KilosortDataset into trials based on start and stop idx.
% trial_ids are specified according to the data structure being merged into.
% If a trial is included in trial_ids but not found in tsi,
% its contents would be blank and trial_has_data(i) would be set false.
%
% Each of seg's properties are now nTrials x ... cells containing the data
% corresponding to that trial
trial_ids = min(tsi.trialId):max(tsi.trialId);
seg = npxutils.KilosortTrialSegmentedDataset(ks, tsi, trial_ids);
disp(seg)

% Compute useful stats about each template and cluster
metrics = ks.getMetrics()

% Plot a drift map, annotated with trial start markers
metrics.plotDriftmap('tsi', tsi);
```



Note that the cluster structure looks distinct during a few time windows near the beginning and end of the recording, corresponding to regions when no trials were being performed (blue ticks near the bottom).



```matlab
% Extract raw waveforms for a specific cluster id at the 24 largest amplitude channels
% Clean these waveforms by subtracting the contribution of other clusters spiking within the same time window
ss = ks.getWaveformsFromRawData('cluster_id', 255, 'num_waveforms', 100, ...
  'best_n_channels', 24,  'subtractOtherClusters', true)

% Plot these waveforms at their physical coordinates on the neuropixel
ss.plotAtProbeLocations()
```

