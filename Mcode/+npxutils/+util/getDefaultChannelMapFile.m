function file = getDefaultChannelMapFile(complainIfEmpty)
% Get the path to the default channel map file for this npxutils session.
%
% file = getDefaultChannelMapFile(complainIfEmpty)
%
% ComplainIfEmpty (logical,false*) is whether to raise an error if the default
% map file is not set.
%
% This is drawn from the NEUROPIXEL_MAP_FILE or NPIX_MAP_FILE environment
% variables, in that order. If NEUROPIXEL_MAP_FILE is defined, but the file does
% not exist on disk, this falls back to NPIX_MAP_FILE. If the file in
% NPIX_MAP_FILE does not exist on disk, this uses that value anyway.
arguments
    complainIfEmpty (1,1) logical = false
end

file = getenv('NEUROPIXEL_MAP_FILE');
if isempty(file) || ~exist(file, 'file')
    file = getenv('NPIX_MAP_FILE');
end

if isempty(file) && complainIfEmpty
    error(['No default channel map file is set. Call setenv to set either ' ...
        'NEUROPIXEL_MAP_FILE or NPIX_MAP_FILE.']);
end

end