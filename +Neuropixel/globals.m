classdef globals
  % Global library properties and settings for Neuropixel Utils
  
  properties (Constant)
    % Path to the root directory of this Neuropixel Utils library
    distroot = string(fileparts(fileparts(mfilename('fullpath'))));
  end
  
  methods (Static)
    
    function out = version
      % The version of this Neuropixel Utils library
      %
      % Returns a string.
      persistent val
      if isempty(val)
        versionFile = fullfile(Neuropixel.globals.distroot, 'VERSION');
        val = strtrim(readtext(versionFile));
      end
      out = val;
    end
    
  end
  
end

function out = readtext(file, encoding)
% Read the contents of a text file as a string
%
% This is analagous to Matlab's readcsv and readtable, and exists because Matlab
% doesn't provide a basic file-slurping mechanism.
  
arguments
  file (1,1) string
  encoding (1,1) string = 'UTF-8' % TODO: auto-detect file encoding via sniffing
end
[fid,msg] = fopen(file, 'r', 'n', encoding);
if fid < 1
  error('Failed opening file %s: %s', file, msg);
end
RAII.fh = onCleanup(@() fclose(fid));
c = fread(fid, Inf, 'uint8=>char');
out = string(c');

end
