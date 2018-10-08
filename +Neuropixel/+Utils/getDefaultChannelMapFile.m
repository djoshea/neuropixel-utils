function file = getDefaultChannelMapFile()
    file = getenv('NEUROPIXEL_MAP_FILE');
    if isempty(file) || ~exist(file, 'file')
        file = getenv('NPIX_MAP_FILE');
    end
end