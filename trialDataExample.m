% dsPath = '/Volumes/pixin/will_sort/Pierre_20180609_08';
dsPath = '/path/to/kilosort_output';

ds = KiloSort.Dataset(dsPath);
ds.load();

% Figure out where trials start and stop
% tsi is TrialSegmentationInfo
% load hopefully cached trialSegmentationInfo.mat file, else parse sync from raw imec.ap.bin file
imecAPBinFile = '/path/to/data.imec.ap.bin';
tsi = NeuropixelExpt.DataLoad.readTrialInfo(imecAPBinFile);

% trial_ids is the list of trial ids for the behavioral data you want to merge into
% e.g. from TrialData object
trial_ids = td.getParam('trialId');
seg = KiloSort.TrialSegmentedDataset(ds08, tsi, trial_ids);

td = td.dropSpikeChannels();
td = seg.addSpikesToTrialData(td, 'npix');
