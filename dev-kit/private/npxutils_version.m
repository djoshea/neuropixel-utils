function out = npxutils_version
versionFile = fullfile(reporoot, 'VERSION');
txt = readtext(versionFile);
txt = regexprep(txt, '\r?\n.*', '');
out = string(txt);
end