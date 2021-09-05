function copyfile2(src, dest)
% A version of copyfile that raises an error on failure
arguments
  src (1,:) string
  dest (1,1) string
end
for file = src
  [ok,msg] = copyfile(file, dest);
  if ~ok
    error('Failed copying file "%s" to "%s": %s', file, dest, msg);
  end
end
end
