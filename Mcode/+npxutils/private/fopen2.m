function out = fopen2(filename, permission, machinefmt, encodingIn)
% A version of fopen that raises an error on failure and defaults to Unicode
arguments
  filename (1,1) string
  permission (1,1) string = 'r'
  machinefmt (1,1) string = 'n'
  encodingIn (1,1) string = 'UTF-8'
end
out = npxutils.internal.util.fopen2(filename, permission, machinefmt, encodingIn);
end