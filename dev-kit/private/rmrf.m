function rmrf(files)
% Recursively delete files and directories
%
% rmrf(files)
arguments
  files string
end

%#ok<*STRNU>

RAII.warn = withwarnoff('MATLAB:DELETE:FileNotFound');
% Sigh. We have to glob out the files ourselves because Matlab doesn't have a
% delete operation that will work on both files and directories.
for iGlob = 1:numel(files)
  [~,details] = dir2(files(iGlob));
  paths = details.path;
  for file = paths'
    if isfile(file)
      delete(file);
    elseif isfolder(file)
      rmdir2(file, 's');
    else
      warning('Dir result "%s" is neither a file nor folder. Huh???', file);
    end
  end
end

end