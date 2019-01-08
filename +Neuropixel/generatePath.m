function path = generatePath(subject, type, date, fileStem)

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

path = fullfile(getenv('NEUROPIXEL_DATAROOT'), subject, type, dateStr, fileStem);

end