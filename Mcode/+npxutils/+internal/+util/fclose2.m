function fclose2(fid)
% A version of fclose that errors on failure
status = fclose(fid);
if status ~= 0
  if isequal(fid, 'all')
    file = '<all filehandles>';
  else
    file = fopen(fid);
  end
  error('Failed doing fclose on file %s. Error message unavailable.', file);
end
end
