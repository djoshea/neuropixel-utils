function fseek2(fid, offset, origin)
% A version of fseek that errors on failure
status = fseek(fid, offset, origin);
if status ~= 0
  file = fopen(fid);
  error('Failed doing fseek on file %s. Error message unavailable.', file);
end
end