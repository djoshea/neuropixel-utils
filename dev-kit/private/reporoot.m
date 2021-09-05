function out = reporoot
% The root dir of the neuropixel-utils repo
out = string(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end