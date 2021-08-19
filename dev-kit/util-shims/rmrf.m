function rmrf(files)
% Recursively delete files and directories
%
% rmrf(files)
arguments
  files string
end
npxutils.internal.util.rmrf(files);
end