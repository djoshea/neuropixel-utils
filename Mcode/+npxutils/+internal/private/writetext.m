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
npxutils.internal.util.writetext(text, file, encoding);
end
