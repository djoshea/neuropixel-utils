function out = fopen2(filename, permission, machinefmt, encodingIn)
% A version of fopen that raises an error on failure and defaults to Unicode
arguments
  filename (1,1) string
  permission (1,1) string = 'r'
  machinefmt (1,1) string = 'n'
  encodingIn (1,1) string = 'UTF-8'
end
[out,msg] = fopen(filename, permission, machinefmt, encodingIn);
if out < 0
  error('Failed opening file "%s": %s', filename, msg);
end
end