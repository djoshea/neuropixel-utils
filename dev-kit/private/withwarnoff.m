function out = withwarnoff(warningId)
% Temporarily disable warnings
arguments
  warningId string
end
origWarnState = warning;
out.RAII = onCleanup(@() warning(origWarnState));
for i = 1:numel(warningId)
  warning('off', warningId(i));
end
end
