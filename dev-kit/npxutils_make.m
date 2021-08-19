function npxutils_make(target, varargin)
% Build tool for the neuropixel-utils library
%
% This is the main build tool for doing all the build and packaging operations
% for the Neuropixel Utils library. It's intended to be called as a command in
% most cases. This is what you will use to build & package the distribution
% files for a release of the library.
%
% Operations:
%   npxutils_make build        - "build" the source code
%   npxutils_make test         - run the tests
%   npxutils_make dist         - build the distribution files
%   
%   npxutils_make toolbox      - build the Matlab Toolbox .mltbx installer file
%   npxutils_make clean        - delete all the derived artifacts
%   npxutils_make doc          - build the project doco (probably broken)
%   npxutils_make doc-preview  - live-preview the project doco (probably broken)
%
% The "doc" build targets are probably broken for you, because they require 
% MkDocs with its "material" theme, and the only way (apparently) to get that
% installed and working (at least on Mac) is with Anaconda. The "doc" targets
% here don't pull in Anaconda; they assume you have a working `mkdocs` command
% on the path. Instead, open a shell, cd in to `doc-src`, and run the `make_doc`
% command there.
%
% The "dist" target does not automatically build the docs for you. You need to
% do that first, if you want the docs in the distribution files to be up to date
% (which you probably do).

%#ok<*STRNU>

arguments
  target (1,1) string
end
arguments (Repeating)
  varargin
end

if target == "build"
  npxutils_build;
elseif target == "doc-src"
  make_package_docs --src
elseif target == "doc"
  make_package_docs;
elseif target == "doc-preview"
  preview_docs;
elseif target == "m-doc"
  make_mdoc;
elseif target == "toolbox"
  npxutils_make m-doc;
  npxutils_package_toolbox;
elseif target == "clean"
  make_clean
elseif target == "test"
  npxutils_launchtests
elseif target == "dist"
  npxutils_make build
  npxutils_make m-doc
  make_dist
elseif target == "util-shim"
  pkg = varargin{1};
  make_util_shim(pkg);
else
  error("Undefined target: %s", target);
end

end

function make_mdoc
rmrf('build/M-doc')
mkdir2('build/M-doc')
copyfile2('docs/*', 'build/M-doc')
if isfile('build/M-doc/feed.xml')
  delete('build/M-doc/feed.xml')
end
end

function preview_docs
import npxutils.internal.util.*;
RAII.cd = withcd('doc-src');
make_doc --preview
end

function make_dist
program = "neuropixel-utils";
distName = program + "-" + npxutils_version;
verDistDir = fullfile("dist", distName);
distfiles = ["build/Mcode" "docs" "examples" "README.md" "LICENSE" "LICENSE-ThirdParty.md" "CHANGES.md"];
rmrf([verDistDir, distName+".tar.gz", distName+".zip"])
if ~isfolder('dist')
  mkdir2('dist')
end
mkdir2(verDistDir)
% Ugh, copyfile() doesn't work right with directories. Guess we'll just call
% a system utility, and require devs to run this on Unix.
quotedFiles = strcat("""", distfiles, """");
copyCmd = sprintf("cp -R %s %s", strjoin(quotedFiles, ' '), verDistDir);
system2(copyCmd);
RAII.cd = withcd('dist');
tar(distName+".tar.gz", distName)
zip(distName+".zip", distName)
end

function make_clean
rmrf(strsplit("dist/* build doc-src/_site M-doc test-output", " "));
end

function make_package_docs(varargin)
doOnlySrc = ismember('--src', varargin);
build_docs;
if ~doOnlySrc
  build_doc;
end
end

function build_doc_src
% Build the generated parts of the doc-src sources
%
% This includes the live scripts from the examples directory.
RAII.cd = withcd(reporoot);
docSrcDir = fullfile(reporoot, 'doc-src');
% TODO: Build stuff from live scripts
end

function build_doc
% Build the final doc files
RAII.cd = withcd(fullfile(reporoot, 'doc-src'));
make_doc;
end

function make_util_shim(pkg)
shimsDir = fullfile(reporoot, 'dev-kit', 'util-shims');
relpkgpath = strjoin(strcat("+", strsplit(pkg, ".")));
pkgdir = fullfile(fullfile(reporoot, 'Mcode'), relpkgpath);
if ~isfolder
  error('Package folder does not exist: %s', pkgdir);
end
privateDir = fullfile(pkgdir, 'private');
if ~isfolder(privateDir)
  mkdir(privateDir);
end
copyfile2(fullfile(shimsDir, '*.m'), privateDir);
fprintf('Util-shimmed package %s', pkg);
end
