function rezFull = runKilosort2(imec, varargin)
% Run Kilosort2 on an IMEC dataset.
%
% rezFull = runKilosort2(imec, ...)
%
% Imec (npxutils.ImecDataset) is the IMEC data set.
%
% This is probably broken.
%
% Options:
%   saveDir
%   workingDir - path to working dir for intermediate files. Should be on a 
%       fast SSD; defaults to tempdir().
%   
% Requires npy-matlab to be installed on your system.


p = inputParser();
p.addParameter('saveDir', imec.pathRoot, @ischar);
p.addParameter('workingDir', tempdir, @ischar); % should be on fast SSD
p.KeepUnmatched = true;
p.parse(varargin{:});

if exist('writeNPY', 'file') ~= 2
    error('npy-matlab was not found on the Matlab path');
end

ops = defaultConfig();
ops.fig = false; % default avoid plotting in main loop, can be overriden as parameter to runKilosort2
ops.trange = [0 Inf];

% custom params added locally
ops.markSplitsOnly = false; % custom parameter for local working version of Kilosort2
ops.spikeThreshBothDirs = false; % apply spkTh threshold from above and below

ops.fproc = fullfile(p.Results.workingDir, sprintf('temp_wh_%s.dat', imec.fileStem));

flds = fieldnames(p.Unmatched);
for iF = 1:numel(flds)
    fld = flds{iF};
    if isfield(ops, fld)
        ops.(fld) = p.Unmatched.(fld);
    else
        error('Unknown option %s', fld);
    end
end
assert(ops.spkTh < 0, 'Option spkTh should be negative');

ops.root = imec.pathRoot;
ops.fs = imec.fsAP;        % sampling rate		(omit if already in chanMap file)
ops.fbinary = imec.pathAP;
ops.scaleToUv = imec.apScaleToUv;
if ~exist(ops.fbinary, 'file')
    error('Imec data file not found.');
end

ops.saveDir = p.Results.saveDir;
if ~exist(ops.saveDir, 'dir')
    mkdir(ops.saveDir);
end

% build channel map for kilosort, providing coordinates only for good channels
goodChannels = imec.goodChannels;
assert(~isempty(goodChannels), 'Must mark good channels');

map = imec.channelMap;
chanMap = map.channelIdsMapped;
xcoords = map.xcoords;
ycoords = map.ycoords;
% this is a mask over mapped channels that is true if channel is good
connected = ismember(1:imec.nChannelsMapped, goodChannels);
kcoords = map.shankInd;
ops.chanMap = fullfile(ops.root,'chanMap.mat');
ops.Nchan = nnz(connected);
ops.NchanTOT  = imec.nChannels;           % total number of channels (omit if already in chanMap file)

fs = imec.fsAP;
save(ops.chanMap, 'chanMap', 'xcoords', 'ycoords', 'connected', 'kcoords', 'fs');

fprintf('Kilosort2: preprocessing\n');
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
fprintf('Kilosort2: Time reordering\n');
rez = clusterSingleBatches(rez);

% main optimization
fprintf('Kilosort2: Main optimization\n');
rez = learnAndSolve8b(rez);

% final merges
fprintf('Kilosort2: Merging clusters\n');
rez = find_merges(rez, 1);

% final splits by SVD
fprintf('Kilosort2: Final splits by SVD\n');
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
fprintf('Kilosort2: Final splits by amplitudes');
rez = splitAllClusters(rez, 0);

% decide on cutoff
fprintf('Kilosort2: Deciding on cutoff');
rez = set_cutoff(rez);

fprintf('Kilosort2: Found %d / %d good units \n', nnz(rez.good>0), numel(rez.good));

% write to Phy
fprintf('Kilosort2: Saving results for Phy\n')
rezToPhy(rez, ops.saveDir);

fprintf('Kilosort2: Saving rez to rez.mat\n')
npxutils.io.exportRezToMat(rez, ops.saveDir);

%     % mark templates that maybe should be split
%     fprintf('Marking split candidates')
%     rez = markSplitCandidates(rez);
%
%     % merge templates that should be merged into new clusters
%     rez = mergeTemplatesIntoClusters(rez);
%
%     % mark split candidates and merged clusters as "mua"
%     rez = assignClusterGroups(rez);

% remove temporary file
delete(ops.fproc);
end

function ops = defaultConfig() %#ok<STOUT>
configFile = getenv('KILOSORT_CONFIG_FILE');
if isempty(configFile)
    configFile = 'configFile384';
end
if ~exist(configFile, 'file')
    error('Could not find Kilosort2 config file %s', configFile);
end

% treat either as full path or as m file on the path
run(configFile);

assert(exist('ops', 'var') > 0, 'Config file %s did not produce ops variable', configFile);
end
