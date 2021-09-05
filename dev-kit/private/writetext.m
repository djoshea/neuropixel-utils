function writetext(text, file, encoding)
% Write text to a file
%
% writetext(text, file, encoding)
%
% Encoding deaults to UTF-8.
%
% Replaces the original file contents.
arguments
  text (1,1) string
  file (1,1) string
  encoding (1,1) string = 'UTF-8'
end
[fid,msg] = fopen(file, 'w', 'n', encoding);
if fid < 1
  error('Failed opening file %s: %s', file, msg);
end
RAII.fh = onCleanup(@() fclose(fid));
fprintf(fid, '%s', text);
end
