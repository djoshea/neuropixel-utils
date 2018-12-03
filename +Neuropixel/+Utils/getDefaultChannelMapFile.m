function file = getDefaultChannelMapFile(complainIfEmpty)
    file = getenv('NEUROPIXEL_MAP_FILE');
    if isempty(file) || ~exist(file, 'file')
        file = getenv('NPIX_MAP_FILE');
    end
    
    if nargin < 1
        complainIfEmpty = false;
    end
    
    if complainIfEmpty
        assert(~isempty(file), 'No default channel map file is set. Call setenv to set either NEUROPIXEL_MAP_FILE or NPIX_MAP_FILE');
    end
end