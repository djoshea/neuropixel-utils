function copyfile2(src, dest)
% A version of copyfile that raises an error on failure
arguments
  src (1,:) string
  dest (1,1) string
end
npxutils.internal.util.copyfile2(src, dest);
end
