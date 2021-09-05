function out = npxutils_toolbox_info
% Info about your toolbox, for use during the build
%
% This is the file you edit to have your toolbox-specific information.
%
% It is expected to be run while your toolbox is loaded, so you can use your
% toolbox functions inside it.
%
% It must return a struct with the fields:
%   name    - a charvec 
%   version - a charvec
%
% The version may have "-[\w+]" suffixes in addition to the standard
% "0.0[.0[.0]]" numeric version format.

out.name = 'neuropixel-utils';
out.version = npxutils.globals.version;
