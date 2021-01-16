# Running Kilosort 1/2 and loading the results

`KilosortDataset` is a wrapper around the output of Kilosort or Kilosort2, which will load the output files back into Matlab for further analysis. Most of these fields are explained in detail in the [Phy documentation](https://phy-contrib.readthedocs.io/en/latest/template-gui/) but we document them here for convenience.

## Running Kilosort

To run Kilosort or Kilosort2 on an ImecDataset:

```matlab
Neuropixel.runKilosort1(imec, ...);
```

Or for Kilosort 2:

```matlab
Neuropixel.runKilosort2(imec, ...);
```

By default, the standard configuration settings will be used. For Kilosort1, these are hardcoded based on `configFiles/StandardConfig_MOVEME.m`. For Kilosort2, the script `configFiles/configFile384.m` will be run to produce the `ops` struct, unless a different configuration file is set in the environment variable `KILOSORT_CONFIG_FILE`, which must be on the path. Default configuration settings can be overridden by passing in extra parameters, e.g.

```matlab
Neuropixel.runKilosort1(imec, 'Th', [4 10], 'GPU', false);
Neuropixel.runKilosort2(imec, 'minfr_goodchannels', 0.1);
```

## Loading Kilosort results

You can create a KilosortDataset instance by pointing at the folder containing the Kilosort output:

```matlab
ks = Neuropixel.KilosortDataset(pathToKilosortOutput();
ks.load();
```

The constructor will optionally take an 'imecDataset' parameter providing the `Neuropixel.ImecDataset` instance if there is no `.imec.ap.bin` file in the Kilosort directory, and a 'channelMap' parameter in case the default is not correct. The results can then be loaded using `ks.load()`.

The descriptions of each property can be found in the `+Neuropixel/KilosortDataset.m` code, copied here for convenience, originally described in the [Phy documentation](https://phy-contrib.readthedocs.io/en/latest/template-gui/):

```matlab
>> ks

KilosortDataset with properties:

                  path: '/data/kilosort/neuropixel_01'
           raw_dataset: [1×1 Neuropixel.ImecDataset]
            channelMap: [1×1 Neuropixel.ChannelMap]
                  fsAP: 30000
           apScaleToUv: 2.3438
                  meta: [1×1 struct]
              pathLeaf: 'neuropixel_01'
              isLoaded: 1
         hasRawDataset: 1
             nChannels: 371
               nSpikes: 8181228
             nClusters: 592
            nTemplates: 653
           nPCFeatures: 32
   nFeaturesPerChannel: 3
          syncBitNames: [16×1 string]
              dat_path: 'neuropixel_01.imec.ap.bin'
        n_channels_dat: 385
                 dtype: 'int16'
                offset: 0
           sample_rate: 30000
           hp_filtered: 0
            amplitudes: [nSpikes × 1 double]
           channel_ids: [nChannels × 1 uint32]
     channel_positions: [nChannels × 2 double]
           pc_features: [nSpikes × nFeaturesPerChannel × nPCFeatures single]
        pc_feature_ind: [nTemplates × nPCFeatures uint32]
     similar_templates: [nTemplates × nTemplates single]
       spike_templates: [nSpikes × 1 uint32]
           spike_times: [nSpikes × 1 uint64]
     template_features: [nSpikes × nPCFeatures single]
  template_feature_ind: [nTemplates × nPCFeatures uint32]
             templates: [nTemplates × nTimePoints × nChannels single]
         templates_ind: [nTemplates × nChannels double]
         whitening_mat: [nChannels × nChannels double]
     whitening_mat_inv: [nChannels × nChannels double]
        spike_clusters: [nSpikes × 1 uint32]
        cluster_groups: [nClusters × 1 categorical]
           cluster_ids: [nClusters × 1 uint32]
         clusters_good: [nClustersGood × 1 uint32]
          clusters_mua: [nClustersMUA × 1 uint32]
        clusters_noise: [nClustersNoise × 1 uint32]
     clusters_unsorted: [nClustersUnsorted × 1 uint32]
```

* `nChannels` : number of channels used by Kilosort

* `nSpikes` : number of spikes extracted

* `nClusters` : number of unique clusters

* `nTemplates` : number of spike templates

* `nPCFeatures` number of spatiotemporal PC features used for templates

* `nFeaturesPerChannel` : number of PC features used for each channel

* `amplitudes` - `[nSpikes]` double vector with the amplitude scaling factor that was applied to the template when extracting that spike

* `channel_ids` - `[nChannels]` uint32 vector with the channel ids used for sorting

* `channel_positions` - `[nChannels, 2]` double matrix with each row giving the x and y coordinates of that channel. Together with the channel map, this determines how waveforms will be plotted in WaveformView (see below).

* `pc_features` - `[nSpikes, nFeaturesPerChannel, nPCFeatures]` single matrix giving the PC values for each spike. The channels that those features came from are specified in pc_features_ind. E.g. the value at pc_features[123, 1, 5] is the projection of the 123rd spike onto the 1st PC on the channel given by pc_feature_ind[5].

* `pc_feature_ind` - `[nTemplates, nPCFeatures]` uint32 matrix specifying which channels contribute to each entry in dim 3 of the pc_features matrix

* `similar_templates` - `[nTemplates, nTemplates]` single matrix giving the similarity score (larger is more similar) between each pair of templates
similar_templates(:, :) single

* `spike_templates` - `[nSpikes]` uint32 vector specifying the identity of the template that was used to extract each spike

* `spike_times` - `[nSpikes]` uint64 vector giving the spike time of each spike in samples. To convert to seconds, divide by sample_rate from params.py.

* `template_features` - `[nSpikes, nTempFeatures]` single matrix giving the magnitude of the projection of each spike onto nTempFeatures other features. Which other features is specified in template_feature_ind.

* `template_feature_ind` - `[nTemplates, nTempFeatures]` uint32 matrix specifying which templateFeatures are included in the template_features matrix.

* `templates` - `[nTemplates, nTimePoints, nTemplateChannels]` single matrix giving the template shapes on the channels given in templates_ind

* `templates_ind` - `[nTemplates, nTempChannels]` double matrix specifying the channels on which each template is defined. In the case of Kilosort templates_ind is just the integers from 0 to nChannels-1, since templates are defined on all channels.

* `whitening_mat` - `[nChannels, nChannels]` double whitening matrix applied to the data during automatic spike sorting

* `whitening_mat_inv` - `[nChannels, nChannels] `double, the inverse of the whitening matrix.

* `spike_clusters` - `[nSpikes]` uint32 vector giving the cluster identity of each spike.

* `cluster_groups` - `[nClusters]` categorical vector giving the "cluster group" of each cluster (noise, mua, good, unsorted)

* `cluster_ids` - `[nClusters]` unique clusters in spike_clusters

## Segmenting a Kilosort dataset into trials

`ks.spike_times` contains the times for each spike in samples from the beginning of the file, but there is a more useful representation for data collected with a trial structure: split the spikes into separate groups based on which trial they occurred in, and convert the times to milliseconds since the start of the trial.

### TrialSegmentationInfo
In order to do this, you need to figure out where trials start and stop. You'll need to write this code, since this will differ for each experimental setup. Essentially, you need to create a `Neuropixel.TrialSegmentationInfo` instance and populate its fields with the correct values:

```matlab
tsi = Neuropixel.TrialSegmentationInfo(nTrials, fsAP);
tsi.idxStart = [list of start sample indices]
tsi.idxStop = [list of stop sample indices];
tsi.trialId = [list of trial ids];
```

Here is an example script that uses the sync channel to determine where trials begin and end. It expects one bit (named `'trialStart'`) to contain TTL pulses each time a trial starts, and another bit (named `'trialInfo'`) to contain ASCII-serialized bits of text occurring at the start of each trial. For example, the string `id=1;c=2` would correspond to `trialId=1` and `conditionId=2`. It also assumes that a trial ends when the next trial begins (or at the end of the file). Long trials can be subsequently truncated using `tsi.truncateTrialsLongerThan(maxDurationSeconds)`.

```matlab
function tsi = parseTrialInfoFromSync(syncRaw, fs, syncBitNames)
    % fs is in samples per second
    % parses the sync line of an neuropixel .imec.ap.bin data file
    % and produces a scalar TrialSegmentationInfo

    % parse the sync
    if isempty(syncBitNames)
        trialInfoBitNum = 1;
        trialStartBitNum = 2;
    else
        [tf, trialInfoBitNum] = ismember('trialInfo', syncBitNames);
        if ~tf, trialInfoBitNum = 1; end
        [tf, trialStartBitNum] = ismember('trialStart', syncBitNames);
        if ~tf, trialStartBitNum = 1; end
    end

    serialBit = bitget(syncRaw, trialInfoBitNum);
    trialStart = bitget(syncRaw, trialStartBitNum);

    % trials start when going high
    idxStart = find(diff(trialStart) == 1) + 1;
    nTrials = numel(idxStart);

    tsi = Neuropixel.TrialSegmentationInfo(nTrials, fs);

    samplesEachBit = round(fs / 1000); % each bit delivered per ms
    for iR = 1:nTrials
        if iR < nTrials
            idxNext = idxStart(iR+1) - 1;
        else
            idxNext = numel(serialBit);
        end
        bitsByTrial = uint8(serialBit(floor(samplesEachBit/2) + idxStart(iR) : samplesEachBit : idxNext));

        lastHigh = find(bitsByTrial, 1, 'last');
        lastHigh = ceil(lastHigh / 8) * 8;
        bitsByTrial = bitsByTrial(1:lastHigh);

        infoThis = parseInfoString(bitsToString(bitsByTrial));

        if isfield(infoThis, 'id')
            tsi.trialId(iR) = str2double(infoThis.id); %#ok<*AGROW>
        else
            tsi.trialId(iR) = NaN;
        end
        if isfield(infoThis, 'c')
            tsi.conditionId(iR) = str2double(infoThis.c);
        else
            tsi.conditionId(iR) = NaN;
        end

        tsi.idxStart(iR) = idxStart(iR);
        tsi.idxStop(iR) = idxNext;
    end

    function out = bitsToString(bits)
        nChar = numel(bits) / 8;
        assert(nChar == round(nChar), 'Bit length must be multiple of 8');
        out = blanks(nChar);
        for iC = 1:nChar
            idx = (1:8) + (iC-1)*8;
            out(iC)  = char(bin2dec(sprintf('%u', bits(idx))));
        end
    end

    function out = parseInfoString(str)
        keyval = regexp(str, '(?<key>\w+)=(?<value>[\d\.]+)', 'names');
        if isempty(keyval)
            warning('Could not parse info string "%s"', str);
            out = struct();
        else
            for i = 1:numel(keyval)
                out.(keyval(i).key) = keyval(i).value;
            end
        end
    end
end
```

### KilosortTrialSegmentedDataset

Once you have the trial boundaries stored in your `TrialSegmentationInfo` instance, you can split the properties of the `KilosortDataset` into each trial, resulting in a `Neuropixel.KilosortTrialSegmentedDataset` instance. To facilitate merging this into another data structure later, you will need to specify the ultimate `trialId` order you want the `KilosortTrialSegmentedDataset` to have. For example, if you have a behavioral data structure, you can extract the list of trial ids from that so that your `KilosortTrialSegmentedDataset` will have a matching trial sequence.

```matlab
trialIds = cat(1, behaviorStruct.trialId);
```

Any trials not found in the `TrialSegmentationInfo` will simply be empty in the `KilosortTrialSegmentedDataset`. If you simply want to preserve the trials in the order they are in `tsi`, you can simply use:

```matlab
trialIds = tsi.trialIds;
```

You can then segment the KilosortDataset using:

```matlab
>> seg = Neuropixel.KilosortTrialSegmentedDataset(ks, tsi, trial_ids)

KilosortTrialSegmentedDataset with properties:

                   dataset: [1×1 Neuropixel.KilosortDataset] % unsegmented KilosortDataset
                 trial_ids: [nTrials × 1 uint32] % trial ids
            trial_has_data: [nTrials × 1 logical] % indicator if trial found in tsi
               trial_start: [nTrials × 1 uint64] % nTrials start sample idx copied from tsi
                trial_stop: [nTrials × 1 uint64] % nTrials start sample idx copied from tsi
                 spike_idx: {nTrials × nClusters cell} % nTrials x nClusters lists of indices into ks.spike_times array for each spike
               cluster_ids: [nClusters × 1 uint32] % copied from ks
            cluster_groups: [nClusters × 1 categorical] % copied from ks
                      sync: {1×1 cell} %  segmented contents of sync channel
              syncBitNames: [16×1 string] % copied from ImecDataset
               raw_dataset: [1×1 Neuropixel.ImecDataset] % original ImecDataset copied from ks
                   nTrials: 1092 % number of trials as numel(trialIds)
           nTrialsHaveData: 1092 % number of trials with matching trialIds in tsi
                 nClusters: 592 % number of clusters (sorted units)
           nChannelsSorted: 385 % number of channels
               channel_ids: [nChannelsSorted × 1 unit32] % list of channel ids used for sorting
         trial_duration_ms: [nTrials × 1 double]
                      fsAP: 30000
                amplitudes: {nTrials × nClusters cell} % each of these will contain a vector with one entry for each spike from that cluster on that trial
               pc_features: {nTrials × nClusters cell}
               spike_times: {nTrials × nClusters cell} % raw sample times from ks
  spike_times_ms_rel_start: {nTrials × nClusters cell} % times from trial start in milliseconds
         template_features: {nTrials × nClusters cell}
           spike_templates: {nTrials × nClusters cell}
```
