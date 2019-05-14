% dsPath = '/Volumes/pixin/will_sort/Pierre_20180609_08';
kilosortDatasetPath = '/path/to/kilosort_output';

% build path to channel map file
channelMapFile '/path/to/neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap.mat');

% load the kilosort dataset
ds = Kilosort.Dataset(kilosortDatasetPath, channelMapFile);
ds.load();

% Figure out where trials start and stop
% tsi is TrialSegmentationInfo that stores this info
% this will load trialSegmentationInfo.mat file if present, but the first time will parse the sync from raw imec.ap.bin file
imecAPBinFile = '/path/to/data.imec.ap.bin';
tsi = NeuropixelExpt.DataLoad.readTrialInfo(imecAPBinFile);

% trial_ids is the list of trial ids for the behavioral data you want to merge into
% e.g. from TrialData object
trial_ids = td.getParam('trialId');

% segment the kilosort dataset into trials, with the final arrays matching
% trial_ids (and missing trials being empty)
% so seg will have numel(trial_ids) trials in it regardless of how many trials are found in the ImecDataFile
seg = Kilosort.TrialSegmentedDataset(ds08, tsi, trial_ids);

% merge into trial data
td = td.dropSpikeChannels();
td = seg.addSpikesToTrialData(td, 'npix');
