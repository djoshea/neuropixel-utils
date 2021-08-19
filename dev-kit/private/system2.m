function out = system2(cmd)
% A version of system that raises an error on failure

if nargout == 0
  status = system(cmd);
else
  [status,out] = system(cmd);
end

if status ~= 0
  error('Command failed (exit status %d). Command: %s', status, cmd);
end

end