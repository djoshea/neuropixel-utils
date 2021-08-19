function path = generatePath(subject, type, date, fileStem, varargin)
% Get the path to a file for given subject/type/etc.
%
% path = npxutils.generatePath(subject, type, date, fileStem, ...)
%
% Options:
%
%   dataRoot (string) - data root path; defaults to using the
%       NEUROPIXEL_DATAROOT environment variable.
%
% Returns the full, absolute path to the specified data file.

p = inputParser();
p.addParameter('dataRoot', '', @ischar);
p.parse(varargin{:});
dataRoot = p.Results.dataRoot;
if isempty(dataRoot)
    dataRoot = getenv('NEUROPIXEL_DATAROOT');
    if isempty(dataRoot)
        error('Must specify either dataRoot or setenv NEUROPIXEL_DATAROOT');
    end
end

if nargin < 2
    type = '';
end
if nargin < 3
    date = '';
end
if nargin < 4
    fileStem = '';
end

if ischar(date)
    dateStr = date;
else
    dateStr = datestr(date, 'YYYY-mm-dd');
end

path = fullfile(dataRoot, subject, type, dateStr, fileStem);

end