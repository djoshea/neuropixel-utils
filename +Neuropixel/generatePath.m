function path = generatePath(subject, type, date, fileStem, varargin)

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