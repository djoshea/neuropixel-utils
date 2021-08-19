function npxutils_batch_package_toolbox
% An entry point for building the toolbox from a "matlab -batch" call
%
% This has error handling to set matlab's exit status appropriately.
%
% You don't need to call this directly in an interactive Matlab session; for
% that, call `npxutils_make test` instead.

try
  load_npxutils;
  npxutils_package_toolbox;
catch err
  fprintf('Error occurred:\n');
  fprintf('%s\n', getReport(err));
  exit(1);
end
