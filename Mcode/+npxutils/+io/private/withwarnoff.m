function out = withwarnoff(warningId)
% Temporarily disable warnings
arguments
  warningId string
end
out = npxutils.internal.util.withwarnoff(warningId);
end
