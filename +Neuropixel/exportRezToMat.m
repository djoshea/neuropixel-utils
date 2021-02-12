function exportRezToMat(rez, savePath)

if nargin < 2
    savePath = rez.ops.saveDir;
end

% write cutoff spikes to separate files to keep rez.mat sufficiently small
if isfield(rez, 'cProj_cutoff_invalid') && ~isempty(rez.cProj_cutoff_invalid)
    templateFeatures = rez.cProj_cutoff_invalid;
    writeNPY(templateFeatures, fullfile(savePath, 'cutoff_template_features.npy'));
end

if isfield(rez, 'cProjPC_cutoff_invalid') && ~isempty(rez.cProjPC_cutoff_invalid)
    pcFeatures = rez.cProjPC_cutoff_invalid;
    writeNPY(pcFeatures, fullfile(savePath, 'cutoff_pc_features.npy'));
end

% drop features from rez, too large
% deciding to keep  so that we can reconstruct the spikes here
rez = clearFields(rez, {'temp', 'cProj_cutoff_invalid', 'cProjPC_cutoff_invalid', 'cProj', 'cProjPC', 'dWUA', 'DATA', 'distrust_batched'});
if isfield(rez, 'ops')
    rez.ops = clearFields(rez.ops, {'gui', 'distrust_data_mask'});
end

% sort spike times
[~, isort]   = sort(rez.st3(:,1), 'ascend');
rez.st3      = rez.st3(isort, :);

% gather all gpuArrays
flds = fieldnames(rez);
for iF = 1:numel(flds)
    val = rez.(flds{iF});
    if isa(val, 'gpuArray')
        rez.(flds{iF}) = gather(val);
    end
end

% save final results as rez
fname = fullfile(savePath, 'rez.mat');
save(fname, 'rez', '-v7.3');

end

function s = clearFields(s, flds)
    for iF = 1:numel(flds)
        fld = flds{iF};
        if isfield(s, fld)
            s = rmfield(s, fld);
        end
    end
end  