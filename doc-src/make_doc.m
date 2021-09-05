function make_doc(varargin)
% Build these doc sources into the final doc/ directory
%
% make_doc
% make_doc --preview
% make_doc --build-only
%
% This will require special configuration on Windows to get it working. And this
% is probably broken on Mac, because you need to use the Anaconda version of
% mkdocs there.
%
% Requires MkDocs and its "material" theme to be installed. See https://www.mkdocs.org/.

action = "install";
args = string(varargin);
if ismember("--preview", args)
  action = "preview";
elseif ismember("--build-only", args)
  action = "build";
end

%#ok<*STRNU>

RAII.cd = withcd(fileparts(mfilename('fullpath')));

if action == "preview"
  % Use plain system() and quash because we expect this to error out when user Ctrl-C's it
  system('mkdocs serve');
else
  system2('mkdocs build');
  if action == "install"
    rmdir2('../docs', 's');
    copyfile2('_site/*.*', '../docs');
    delete('../docs/make_doc.m');
  end
end

end

function out = withcd(dir)
% Temporarily change to a new directory
origDir = pwd;
cd(dir);
out = onCleanup(@() cd(origDir));
end

function out = system2(cmd)
% A version of system that raises an error on failure

if nargout == 0
  status = system(cmd);
else
  [status,out] = system(cmd);
end

if status ~= 0
  error('Command failed (exit status %d). Command: %s', status, cmd);
end

end
