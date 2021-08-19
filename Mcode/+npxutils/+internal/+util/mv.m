function mv(source, dest)
% A version of movefile that raises an error on failure
[ok,msg] = movefile(source, dest);
if ~ok
  error('Failed moving "%s" to "%s": %s', source, dest, msg);
end
end