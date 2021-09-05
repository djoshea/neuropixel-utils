function [names, details] = dir2(dirPath)
% Like DIR, but better
%
% [names, details] = dir2(path)
%
% Like DIR, but better: Doesn't include . and .. in the results. Returns the
% details as a table, and time values as datetimes.
%
% Path is the path to the directory, file, or fileglob pattern to list.
%
% Returns:
%    names - a list of the child files and directories under path, as string
%        array
%    details - entry details, as a table array, with variables:
%        name - base name of file (string)
%        path - full path to file (string)
%        mtime - last modification time of file (datetime with TimeZone)
%        bytes - size of file in bytes (double)
%        isdir - whether file is a directory (logical)

d = dir(dirPath);
tfDots = ismember({d.name}, {'.', '..'});
d(tfDots) = [];

names = string({d.name});

if nargout > 1
  name = names';
  folder = string({d.folder}');
  path = fullfile(folder, name);
  mtime = datetime([d.datenum]', 'ConvertFrom','datenum');
  mtime.TimeZone = datetime.SystemTimeZone;
  bytes = [d.bytes]';
  isdir = [d.isdir]';
  details = table(name, mtime, bytes, isdir, path, folder);
  % TODO: add ishidden var to details?
end

end