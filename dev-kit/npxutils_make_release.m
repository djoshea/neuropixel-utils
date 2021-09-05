function npxutils_make_release(newVersion)
% Cut a release of neuropixel-utils
%
% npxutils_make_release <newVersion>
%
% This does most of the steps of performing a release for you. It will:
%    - run the tests
%    - build the project
%    - build the distribution artifacts
%    - update the VERSION file
%    - tag a new release in git
%
% NewVersion (string) is the version number of the release you want to make, like
% "1.2.3".

%#ok<*STRNU>

arguments
  newVersion (1,1) string
end

% Enforce SemVer version formatting
% TODO: Fix this to prohibit leading zeros, but still allow "0" as a version
% number component.
validVersionPat = '^[0-9]+\.[0-9]+\.[0-9]+[+a-zA-Z-]*$';
if isempty(regexp(newVersion, validVersionPat, 'once'))
  error('Invalid version number: "%s". Please use a valid SemVer version number.', newVersion);
end
echo('Cutting release %s', newVersion)

RAII.cd = withcd(reporoot);

echo("Checking prerequisites...")
rslt = system2('git status --porcelain');
if ~isempty(rslt)
  error("Error: Your repo is not clean! There are local changes. No release for you!");
end

% TODO: Figure out how to check for pending remote changes on the tracked branch
% in git. Diffing origin/main is incorrect, because we may be releasing a patch
% release from a different branch.
%must git remote update
%if [[ ! -z "$(git diff origin/main)" ]]; then
%  echo >&2 "Error: There are pending remote changes in git. No release for you!"
%  exit 1
%fi

% First run the tests before we make any changes

echo("Running tests...")
% TODO

% Okay, let's do the release

markVersion(newVersion);

echo("Regenerating doco...");
npxutils_make doc

echo("Building dist...")
% TODO: Implement "make dist"
echo("Building toolbox...")
npxutils_package_toolbox

echo("Tagging release...")
system2(sprintf('git commit -a -m "Version %s"', newVersion));
system2(sprintf('git tag "v%s"', newVersion));
markVersion(newVersion+"+")
system2(sprintf('git commit -a -m "Open development for next release"'))
system2(sprintf('git push'))
system2(sprintf('git push --tags'))

echo("Release is pushed! Now go to GitHub and draft the actual Release:")
echo("   https://github.com/djoshea/neuropixel-utils/releases")
if ismac
  % TODO: Figure out how to do this on Windows and Linux
  system('open "https://github.com/djoshea/neuropixel-utils/releases"')
end

end

function markVersion(version)
writetext(sprintf('%s\n'), version), 'VERSION');
end

function echo(fmt, varargin)
if nargin == 0
  fprintf('\n');
  return
end
fprintf([char(fmt) '\n'], varargin{:});
end

