function out = npxutils_repo_root
% Get the path to the root of this repo

thisDir = fileparts(mfilename('fullpath'));
out = fileparts(fileparts(thisDir));
out = string(out);

end
