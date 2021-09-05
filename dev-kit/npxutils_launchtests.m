function npxutils_launchtests
% Entry point for running full test suite from command line or automation
%
% This function is intended for being called by Continuous Integration tools or
% other automated contexts, via a "matlab -batch" call. It loads the 
% neuropixel-utils library and runs the full test suite.
%
% You don't need to call this directly in an interactive Matlab session; for
% that, call `npxutils_make test` instead.

load_npxutils;

results = npxutils.test.runtests %#ok<NOPRT,NASGU>

end
