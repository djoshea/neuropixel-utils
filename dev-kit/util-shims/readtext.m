function out = readtext(file, encoding)
% Read the contents of a text file as a string
%
% This is analagous to Matlab's readcsv and readtable, and exists because Matlab
% doesn't provide a basic file-slurping mechanism.
arguments
  file (1,1) string
  encoding (1,1) string = 'UTF-8' % TODO: auto-detect file encoding via sniffing
end
out = npxutils.internal.util.readtext(file, encoding);
end
