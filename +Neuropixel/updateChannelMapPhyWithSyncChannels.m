function updateChannelMapPhyWithSyncChannels(imec)
% update the channel map file so that phy shows the sync channel as well

savePath = imec.pathRoot;
map = imec.channelMap;
syncChannelIndex = imec.syncChannelIndex;
assert(isscalar(syncChannelIndex));

outputs = {'channel_map.npy', 'channel_positions.npy'};

% fs = dir(fullfile(savePath, '*.npy'));
% for i = 1:length(fs)
%     fname = fs(i).name;
%     % don't delete .npy files which have nothing to do with us
%     if find(strcmp(fname, outputs))
%         delete(fullfile(savePath, fname));
%     end
% end

connected   = map.connected(:);
xcoords     = map.xcoords(connected);
ycoords     = map.ycoords(connected);

chanMap     = map.chanMap(connected);

if ~ismember(syncChannelIndex, chanMap)
    chanMap(end+1) = syncChannelIndex;
    xcoords(end+1) = min(xcoords);
    ycoords(end+1) = min(ycoords) - max(diff(sort(map.ycoords)));
end

chanMap0ind = int32(chanMap - 1);

whitening = readNPY(fullfile(savePath, 'whitening_mat.npy'));

writeNPY(chanMap0ind, fullfile(savePath, 'channel_map.npy'));
writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));

nCh = nnz(connected);

if size(whitening, 1) == nCh
    % expand whitening and whitening inv
    % need to expand by 1 for sync
    whitening(end+1, :) = 0;
    whitening(:, end+1) = 0;
    whitening(end, end) = 1;

    whitening_inv = whitening^(-1);
    
    writeNPY(whitening, fullfile(savePath, 'whitening_mat.npy'));
    writeNPY(whitening_inv, fullfile(savePath, 'whitening_mat_inv.npy'));
end



end
