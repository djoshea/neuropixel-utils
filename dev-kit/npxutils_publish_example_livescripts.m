function npxutils_publish_example_livescripts
% Publish all the Live Scripts in examples to doc-src
%
% Does not actually rebuild the final doco in the docs folder! You must
% manually do that by opening a terminal window and running MkDocs with the
% doc-src/make_doc script.

% See: https://www.mathworks.com/matlabcentral/answers/282820-programmatically-run-and-export-live-script

%#ok<*NASGU>
%#ok<*ASGLU>

rootDir = npxutils_repo_root;
exDir = rootDir + '/examples';
docSrcDir = rootDir + '/doc-src';
publishDir = docSrcDir + '/docs/publish';
htmlDir = publishDir + '/html';
mdDir = publishDir + '/md';

[fnames, fdeets] = dir2(exDir + '/*.mlx');
n = numel(fnames);

fprintf('Exporting %d example Live Scripts to HTML and Markdown\n', n);

for i = 1:n
  mlxFileBase = fnames(i);
  mlxFilePath = exDir + '/' + mlxFileBase;
  [pDir, fileStem, mlxExt] = fileparts(mlxFilePath);
  fileStemPath = fullfile(pDir, fileStem);
  fprintf('Exporting: %s to %s.*\n', mlxFileBase, fileStemPath);
  
  % Export to HTML
  %
  % Yes, we're using an undocumented Matlab internal here. Gotta have it.
  htmlFile = htmlDir + '/' + fileStem + '.html';
  matlab.internal.liveeditor.openAndConvert(char(mlxFilePath), char(htmlFile));
  fprintf('Exported:  -> %s\n', htmlFile);
  
  % Export to Markdown/LaTeX
  latexFile = mdDir + '/' + fileStem + '.tex';
  matlab.internal.liveeditor.openAndConvert(char(mlxFilePath), char(latexFile));
  mdFile = mdDir + '/' + fileStem + '.md';
  latexFileStemPath = mdDir + '/' + fileStem;
  npxutils.internal.livescript2markdown.latex2markdown(latexFileStemPath);
  fprintf('Exported:  -> %s\n', mdFile);
end

fprintf('All example Live Scripts published\n');

end