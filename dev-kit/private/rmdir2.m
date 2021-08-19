function rmdir2(dir, varargin)
% A version of rmdir that raises errors on failure
[ok,msg] = rmdir(dir, varargin{:});
if ~ok
  error('rmdir of "%s" failed: %s', dir, msg);
end
end