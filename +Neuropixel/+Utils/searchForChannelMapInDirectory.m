function file = searchForChannelMapInDirectory(path)
    % try searching the directory for a file called channel_map.mat
    file = "";
    candidates = ["channel_map.mat", "chanMap.mat"];
    for iC = 1:numel(candidates)
        candidateChannelMapFile = fullfile(path, candidates(iC));
        if exist(candidateChannelMapFile, 'file')
            file = candidateChannelMapFile;
            break;
        end
    end
end