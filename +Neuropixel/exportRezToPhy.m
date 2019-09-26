function exportRezToPhy(rez, savePath, varargin)
% based on a rez dataset, save it's contents to disk as a set of .npy files
% which are compatible both with Phy/Phy2 sorting and subsequent loading as a Neuropixel.KilosortDataset
    if nargin < 2 || isempty(savePath)
        savePath = rez.ops.saveDir;
    end

    p = inputParser();
    p.addParameter('export_cutoff_as_offset_clusters', false, @islogical); % export the cutoff spikes with cluster numbers 10000 + original, only for debugging, neuropixel-utils directly handles the cutoff spikes directly
    p.addParameter('export_cutoff_hidden', true, @islogical); % export the cutoff spikes with cluster numbers 10000 + original, only for debugging, neuropixel-utils directly handles the cutoff spikes directly
    p.parse(varargin{:});
    export_cutoff_as_offset_clusters = p.Results.export_cutoff_as_offset_clusters;
    export_cutoff_hidden = p.Results.export_cutoff_hidden;
    if export_cutoff_hidden
        assert(isfield(rez, 'st3_cutoff_invalid'), 'Missing rez.st3_cutoff_invalid');
    end

    cluster_col = 7; % use post split cluster assignments
    markSplitsOnly = getOr(rez.ops, 'markSplitsOnly', false);

    % spikeTimes will be in samples, not seconds
    Wphy = gather(single(rez.Wphy));
    U = gather(single(rez.U));
    rez.mu = gather(single(rez.mu));

    CUTOFF_CLUSTER_OFFSET = 10000;
    if export_cutoff_as_offset_clusters
        nNonCutoff = size(rez.st3, 1);
        rez.st3 = cat(1, rez.st3, rez.st3_cutoff_invalid);
        rez.st3(nNonCutoff+1:end, cluster_col) = rez.st3(nNonCutoff+1:end, cluster_col) + CUTOFF_CLUSTER_OFFSET; % offset the cluster numbers and export as usual
        rez.cProj = cat(1, rez.cProj, rez.cProj_cutoff_invalid);
        rez.cProjPC = cat(1, rez.cProjPC, rez.cProjPC_cutoff_invalid);
    end

    [~, isort]   = sort(rez.st3(:,1), 'ascend');
    rez.st3      = rez.st3(isort, :);
    rez.cProj    = rez.cProj(isort, :);
    rez.cProjPC  = rez.cProjPC(isort, :, :);

    % ix = rez.st3(:,4)>12;
    % rez.st3 = rez.st3(ix, :);
    % rez.cProj = rez.cProj(ix, :);
    % rez.cProjPC = rez.cProjPC(ix, :,:);

    fs = dir(fullfile(savePath, '*.npy'));
    for i = 1:length(fs)
       delete(fullfile(savePath, fs(i).name)); 
    end
    if exist(fullfile(savePath, '.phy'), 'dir')
        rmdir(fullfile(savePath, '.phy'), 's');
    end

    spikeTimes = uint64(rez.st3(:,1));
    % [spikeTimes, ii] = sort(spikeTimes);
    spikeTemplatesPreSplit = uint32(rez.st3(:,2));
    spikeTemplates = uint32(rez.st3(:,6));
    spikeClusters = uint32(rez.st3(:,cluster_col));
    amplitudes = rez.st3(:,3);

    if export_cutoff_hidden
        cutoff_spikeTimes = uint64(rez.st3_cutoff_invalid(:,1));
        cutoff_spikeTemplatesPreSplit = uint32(rez.st3_cutoff_invalid(:,2));
        cutoff_spikeTemplates = uint32(rez.st3_cutoff_invalid(:,6));
        cutoff_spikeClusters = uint32(rez.st3_cutoff_invalid(:,cluster_col));
        cutoff_amplitudes = rez.st3_cutoff_invalid(:,3);
    end

    Nchan = rez.ops.Nchan;

    xcoords     = rez.xcoords(:);
    ycoords     = rez.ycoords(:);
    chanMap     = rez.ops.chanMap(:);
    chanMap0ind = chanMap - 1;

    nt0_phy = size(rez.Wphy,1);

    Nfilt = size(Wphy,2);

    templates = zeros(Nchan, nt0_phy, Nfilt, 'single');
    for iNN = 1:size(templates,3)
       templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(Wphy(:,iNN,:))'; 
    end
    templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
    templatesInds = repmat(0:size(templates,3)-1, size(templates,1), 1); % we include all channels so this is trivial

    templateFeatures = rez.cProj;
    templateFeatureInds = uint32(rez.iNeigh);
    pcFeatures = rez.cProjPC;
    pcFeatureInds = uint32(rez.iNeighPC);

    if export_cutoff_hidden
        cutoff_templateFeatures = rez.cProj_cutoff_invalid;
        cutoff_pcFeatures = rez.cProjPC_cutoff_invalid;
    end

    whiteningMatrix = rez.Wrot/rez.ops.scaleproc;
    whiteningMatrixInv = whiteningMatrix^-1;

    % here we compute the amplitude of every template...

    % unwhiten all the templates
    tempsUnW = zeros(size(templates));
    for t = 1:size(templates,1)
        tempsUnW(t,:,:) = squeeze(templates(t,:,:))*whiteningMatrixInv;
    end

    % The amplitude on each channel is the positive peak minus the negative
    tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

    % The template amplitude is the amplitude of its largest channel 
    tempAmpsUnscaled = max(tempChanAmps,[],2);

    % assign all spikes the amplitude of their template multiplied by their
    % scaling amplitudes
    spikeAmps = tempAmpsUnscaled(spikeTemplates).*amplitudes;

    % take the average of all spike amps to get actual template amps (since
    % tempScalingAmps are equal mean for all templates)

    % important that we ignore the +10000 clusters added when export_cutoff_as_offset_clusters is true
    cluster_offset = 1;
    nClusters = length(rez.good);
    cluster_ids = (1:nClusters) - cluster_offset;
    cluster_subs = spikeClusters;
    if export_cutoff_as_offset_clusters
        cluster_subs(cluster_subs > CUTOFF_CLUSTER_OFFSET) = cluster_subs(cluster_subs > CUTOFF_CLUSTER_OFFSET) - CUTOFF_CLUSTER_OFFSET;
    end
    clusterAmps = accumarray(cluster_subs, spikeAmps, [nClusters, 1], @mean);
    gain = getOr(rez.ops, 'gain', 1);
    clusterAmps = gain * clusterAmps;

    % ta = clusterAverage(spikeTemplates, spikeAmps);
    % tids = unique(spikeTemplates);
    % tempAmps(tids) = ta; % because ta only has entries for templaters that had at least one spike
    % gain = getOr(rez.ops, 'gain', 1);
    % tempAmps = gain*tempAmps'; % for consistency, make first dimension template number

    %% write primary files needed by phy 
    
    writeNPY_local(spikeTimes, 'spike_times.npy');
    writeNPY_local(uint32(spikeTemplates-cluster_offset), 'spike_templates.npy'); % -1 for zero indexing
    writeNPY_local(uint32(spikeTemplatesPreSplit-cluster_offset), 'spike_templates_preSplit.npy'); % -1 for zero indexing
    writeNPY_local(uint32(spikeClusters-cluster_offset), 'spike_clusters.npy'); % -1 for zero indexing
    writeNPY_local(uint32(spikeClusters-cluster_offset), 'spike_clusters_ks2orig.npy'); % -1 for zero indexing (this is a copy made that won't be touched by Phy / Unit merge tool)

    writeNPY_local(amplitudes, 'amplitudes.npy');
    writeNPY_local(templates, 'templates.npy');
    writeNPY_local(templatesInds, 'templates_ind.npy');

    chanMap0ind = int32(chanMap0ind);

    writeNPY_local(chanMap0ind, 'channel_map.npy');
    writeNPY_local([xcoords ycoords], 'channel_positions.npy');

    writeNPY_local(templateFeatures, 'template_features.npy');
    writeNPY_local(templateFeatureInds'-cluster_offset, 'template_feature_ind.npy');% -1 for zero indexing
    writeNPY_local(pcFeatures, 'pc_features.npy');
    writeNPY_local(pcFeatureInds'-cluster_offset, 'pc_feature_ind.npy');% -1 for zero indexing

    if export_cutoff_hidden
        writeNPY_local(cutoff_spikeTimes, 'cutoff_spike_times.npy');
        writeNPY_local(uint32(cutoff_spikeTemplates-cluster_offset), 'cutoff_spike_templates.npy'); % -1 for zero indexing
        writeNPY_local(uint32(cutoff_spikeTemplatesPreSplit-cluster_offset), 'cutoff_spike_templates_preSplit.npy'); % -1 for zero indexing
        writeNPY_local(uint32(cutoff_spikeClusters-cluster_offset), 'cutoff_spike_clusters.npy'); % -1 for zero indexing
        writeNPY_local(cutoff_amplitudes, 'cutoff_amplitudes.npy');
        writeNPY_local(cutoff_templateFeatures, 'cutoff_template_features.npy');
        writeNPY_local(cutoff_pcFeatures, 'cutoff_pc_features.npy');
    end

    writeNPY_local(whiteningMatrix, 'whitening_mat.npy');
    writeNPY_local(whiteningMatrixInv, 'whitening_mat_inv.npy');

    if isfield(rez, 'simScore')
        similarTemplates = rez.simScore;
        writeNPY_local(similarTemplates, 'similar_templates.npy');
    end

    %% write cluster metadata for cluster table in phy
    ks_labels = strings(nClusters, 1);
    ks_labels(rez.good == 1) = "good";
    ks_labels(rez.good == 0) = "mua";
    rez.est_contam_rate(isnan(rez.est_contam_rate)) = 1;
    est_contam_rate = rez.est_contam_rate*100;
    splitsrc = rez.splitsrc - cluster_offset;
    splitdst = rez.splitdst - cluster_offset;
    splitauc = rez.splitauc;
    mergecount = rez.mergecount;
    
    % write MergeSplit description
    merge_split_desc = strings(nClusters, 1);
    for j = 1:nClusters
        % construct note indicating merge or split
        str = '';
        if rez.mergecount(j) > 1
            str = [str, sprintf('merged %d; ', rez.mergecount(j))]; %#ok<AGROW>
        end
        if markSplitsOnly
            if rez.split_candidate(j)
                str = [str, 'split candidate; ']; %#ok<AGROW>
            end
        else
            if ~isnan(rez.splitsrc(j)) && rez.splitsrc(j) ~= j
                str = [str, sprintf('split from %d; ', rez.splitsrc(j) - cluster_offset)]; %#ok<AGROW> % templates are written for phy as 0 indexed, so subtract 1
            end
            if ~isnan(rez.splitdst(j)) && rez.splitdst(j) ~= j
                str = [str, sprintf('split to %d; ', rez.splitdst(j) - cluster_offset)]; %#ok<AGROW> % templates are written for phy as 0 indexed, so subtract 1
            end
        end
        if length(str) > 2
            str = str(1:end-2); % strip trailing '; '
        end

        merge_split_desc(j) = str;
    end

    writeTSV('cluster_KSLabel.tsv', 'KSLabel', ks_labels, '%s');
    writeTSV('cluster_ContamPct.tsv', 'ContamPct', est_contam_rate, '%.1f');
    writeTSV('cluster_Amplitude.tsv', 'Amplitude', clusterAmps, '%.1f');
    writeTSV('cluster_MergeSplit.tsv', 'MergeSplit', merge_split_desc, '%s');
    
    %% export other things not needed by phy but needed by KilosortDataset
    ops = rez.ops;
    save(fullfile(savePath, 'ks_ops.mat'), 'ops');
    
    nBatches  = ops.Nbatch;
    ntPerBatch  	= uint64(ops.NT - ops.ntbuff);
    batch_starts = (uint64(1) : ntPerBatch : (ntPerBatch*uint64(nBatches) + uint64(1)))';
    writeNPY_local(uint64(batch_starts), 'batch_starts.npy');
    
    writeNPY_local(uint32(splitsrc), 'cluster_splitsrc.npy');
    writeNPY_local(uint32(splitdst), 'cluster_splitdst.npy');
    writeNPY_local(splitauc, 'cluster_splitauc.npy');
    writeNPY_local(mergecount, 'cluster_mergecount.npy');
    writeNPY_local(rez.est_contam_rate, 'cluster_est_contam_rate.npy');
    writeNPY_local(rez.split_candidate, 'cluster_split_candidate.npy');
    writeNPY_local(rez.split_orig_template, 'cluster_split_orig_template.npy');
    
    writeNPY_local(rez.W, 'template_W.npy');
    writeNPY_local(rez.U, 'template_U.npy');
    writeNPY_local(rez.ccb, 'batchwise_ccb.npy');
    writeNPY_local(rez.iorig, 'batch_sort_order.npy');

    nTemplateRank = size(rez.WA, 3);
    nTemplatePCs = size(rez.W_a, 2);
    nTemplates = size(templates, 1);
    nTemplateTimepoints = size(rez.WA, 1);
    nChannelsSorted = rez.ops.Nchan;

    writeNPY_local(single(rez.WA), 'template_W_batch.npy');
    writeNPY_local(reshape(single(rez.W_a), [nTemplateTimepoints, nTemplateRank, nTemplatePCs, nTemplates]), 'template_W_batch_US.npy');
    writeNPY_local(single(rez.W_b), 'template_W_batch_V.npy');
    
    writeNPY_local(single(rez.UA), 'template_U_batch.npy');
    writeNPY_local(reshape(single(rez.U_a), [nChannelsSorted, nTemplateRank, nTemplatePCs, nTemplates]), 'template_U_batch_US.npy');
    writeNPY_local(single(rez.U_b), 'template_U_batch_V.npy');
    
    %% write params.py file
    fid = fopen(fullfile(savePath,'params.py'), 'w');
    [~, fname, ext] = fileparts(rez.ops.fbinary);
    fprintf(fid,['dat_path = ''',fname ext '''\n']);
    fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
    fprintf(fid,'dtype = ''int16''\n');
    fprintf(fid,'offset = 0\n');
    if mod(rez.ops.fs,1)
        fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
    else
        fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
    end
    fprintf(fid,'hp_filtered = False');
    fclose(fid);
    
    function writeNPY_local(vals, fname)
        writeNPY(gather(vals), fullfile(savePath, fname));
    end

    function writeTSV(fname, col_name, vals, fmat)
        tsvfid = fopen(fullfile(savePath, fname),'w');
        if tsvfid == -1
            error('Error opening %s for writing', fname);
        end
        fprintf(tsvfid, 'cluster_id\t%s\r\n', col_name);

        fmatStr = append('%d\t', fmat, '\r\n');
        for ii = 1:numel(cluster_ids)
            fprintf(tsvfid, fmatStr, cluster_ids(ii), vals(ii));
        end

        fclose(tsvfid); 
    end    

end