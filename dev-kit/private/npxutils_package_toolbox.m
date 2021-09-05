function npxutils_package_toolbox
% Packages this toolbox as a Matlab Toolbox .mltbx file
%
% npxutils_package_toolbox
%
% The package must be loaded on to the Matlab path in order for this to work.
%
% This must be run with your cwd set to the root of the repo.

tbxInfo = npxutils_toolbox_info;
tbxName = tbxInfo.name;

if ~isfolder('dist')
  mkdir('dist');
end

% Munge the project file
% Toolboxes don't support "-<pre>" or "+" suffixes in versions
tbxVer = tbxInfo.version;
baseTbxVer = strrep(regexprep(tbxVer, '-.*', ''), '+', '');
fprintf('Packaging %s %s (as %s)\n', tbxName, tbxVer, baseTbxVer);

prjInFile = sprintf('%s.prj.in', tbxName);
prjFile = sprintf('%s.prj', tbxName);
prjTxt = fileread(prjInFile);
prjTxt = strrep(prjTxt, '${PROJECT_VERSION}', baseTbxVer);
spew(prjFile, prjTxt);

% I can't control the output file name from the project file, so we have to move
% it ourselves
builtFile = [tbxName '.mltbx'];
targetFile = sprintf('dist/%s-%s.mltbx', tbxName, tbxVer);
if isfile(targetFile)
  delete(targetFile);
end

matlab.addons.toolbox.packageToolbox(prjFile);
movefile(builtFile, targetFile);
delete(prjFile);

fprintf('%s %s packaged to %s (as %s)\n', tbxName, tbxVer, targetFile, baseTbxVer);

end

function spew(file, txt)
[fid,msg] = fopen(file, 'w');
RAII.fid = onCleanup(@() fclose(fid));
if fid < 1
  error('Failed opening file %s for writing: %s', file, msg);
end
fprintf(fid, '%s', txt);
end