function rmrf(files)
% Recursively delete files and directories
%
% npxutils.io.rmrf(files)
%
% Files (string) is a list of paths to files or directories. All are deleted.
%
% If any deletion fails, an error is raised.
arguments
  files string
end
npxutils.internal.util.rmrf(files);
end